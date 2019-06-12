// -*- c-basic-offset: 4; related-file-name: "devicebalancer.hh" -*-
/*
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, subject to the conditions
 * listed in the Click LICENSE file. These conditions include: you must
 * preserve this copyright notice, and you cannot mention the copyright
 * holders in advertising related to the Software without their permission.
 * The Software is provided WITHOUT ANY WARRANTY, EXPRESS OR IMPLIED. This
 * notice is a summary of the Click LICENSE file; the license in that file is
 * legally binding.
 */

#include <click/config.h>
#include <click/args.hh>
#include <click/error.hh>
#include <click/dpdkdevice.hh>
#include <click/multithread.hh>
#include <rte_flow.h>
#include "devicebalancer.hh"
#include "xdploader.hh"
#include <queue>
#include <vector>
#include <limits.h>
#include <float.h>
//#include <glpk.h>             /* GNU GLPK linear/mixed integer solver */
#include <nlopt.hpp>

#include "../analysis/rssverifier.hh"
#include "../flow/flowipmanager.hh"

#if RTE_VERSION >= RTE_VERSION_NUM(17,5,0,0)
    #include <click/flowdirector.hh>
#endif

CLICK_DECLS

String rewrite_id(String rule, int core) {
    int pos = rule.find_left("queue index");
    int end = pos + 12 + rule.substring(pos+12).find_left(' ');
    return rule.substring(0,pos)+ "queue index "+String(core)+ rule.substring(end);

}

int LoadTracker::load_tracker_initialize(ErrorHandler* errh) {
    _past_load.resize(click_max_cpu_ids());
    _moved.resize(click_max_cpu_ids(),false);

    _last_movement.resize(click_max_cpu_ids());
    click_jiffies_t now = click_jiffies();
    for (int i =0; i < _last_movement.size(); i++) {
        _last_movement[i] = now;
    }

    return 0;
}

int BalanceMethod::configure(Vector<String> &, ErrorHandler *)  {
	return 0;
}


int BalanceMethod::initialize(ErrorHandler *errh, int startwith) {
    return 0;
};

BalanceMethodDevice::BalanceMethodDevice(DeviceBalancer* b, Element* fd) : BalanceMethod(b) {
    _fd = (EthernetDevice*)fd->cast("EthernetDevice");
    if (!_fd) {
        click_chatter("Not an Ethernet Device");
    }
    _is_dpdk = fd->cast("DPDKDevice");
    if (_is_dpdk) {
        click_chatter("DPDK mode");
    } else
        click_chatter("Kernel mode");
    assert(_fd);
    assert(_fd->set_rss_reta);
    assert(_fd->get_rss_reta_size);
}


/**
 * Metron
 */

int MethodMetron::initialize(ErrorHandler *errh, int startwith) {
	if (!_is_dpdk)
		return errh->error("Metron only works with DPDK");
    FlowDirector *flow_dir = FlowDirector::get_flow_director(((DPDKDevice*)_fd)->port_id);
    assert(flow_dir);

    // Invoke Flow Director only if active
    if (flow_dir->active()) {
        // There is a file with (user-defined) rules
        if (!_rules_file.empty()) {
            HashMap<long, String> rules_map;
            const String rules_str = (const String) flow_dir->load_rules_from_file_to_string(_rules_file);

            if (rules_str.empty()) {
                return errh->error("Failed to add rules due to empty input from file");
            }

            // Tokenize them to facilitate the insertion in the flow cache
            Vector<String> rules_vec = rules_str.trim_space().split('\n');

            for (uint32_t i = 0; i < (unsigned)rules_vec.size(); i++) {
                String rule = rules_vec[i] + "\n";

                // Add rule to the map
                rules_map.insert((long) i, rewrite_id(rule,i % startwith));
            }

            int ret = flow_dir->update_rules(rules_map, true);

            if (ret < 0)
                return errh->error("Could not install rules.");
            else
                click_chatter("Installed %d rules", ret);
        } else
            return errh->error("No rule file !");
    } else {
        return errh->error("Flow director is not active.");
    }


    ((DPDKDevice*)_fd)->dpdk_set_rss_max(startwith);

    load_tracker_initialize(errh);

    auto cache = flow_dir->get_flow_cache();
    if (!cache) {
        return errh->error("Flow dir is not initialized");
    }

    if (!cache->has_rules())
        return errh->error("Cache error !");

    return 0;
}


int MethodMetron::configure(Vector<String> &conf, ErrorHandler *errh)  {
	if (Args(balancer, errh).bind(conf)
			.read_or_set("MIN_MOVEMENT", _min_movement, 3)
			.read_or_set("DEFLATION_FACTOR", _deflation_factor, 2)
			.consume() < 0)
		return -1;

	return 0;
}

void MethodMetron::rebalance(Vector<Pair<int,float>> load) {
    FlowDirector *flow_dir = FlowDirector::get_flow_director(((DPDKDevice*)_fd)->port_id);
    click_jiffies_t now = click_jiffies();
    std::vector<Pair<int,float>> underloaded;
    std::vector<Pair<int,float>> overloaded;
    for (int j = 0; j < load.size(); j++) {
        int cpuid = load[j].first;
        float load_current = load[j].second;

        float load_past = _past_load[cpuid];
        float load_diff = load_current-load_past;
        float load_future = load_current + load_diff;
        _past_load[cpuid] = load_current;

        click_chatter("Load of core %d : %f. Past : %f, Future : %f, Last movement : %f (seconds)", cpuid, load_current, load_past, load_future, (float)(now - _last_movement[cpuid]) / CLICK_HZ);
        if (now - _last_movement[cpuid] < CLICK_HZ * _min_movement)
            continue;
        if (load_future > balancer->_overloaded_thresh || load_current > balancer->_overloaded_thresh) {
            overloaded.push_back(Pair<int,float>{cpuid,load_current});
        } else if (load_future < balancer->_underloaded_thresh && load_current < balancer->_underloaded_thresh && load_past < balancer->_underloaded_thresh) {
            underloaded.push_back(Pair<int,float>{cpuid,load_current});
        }
    }


    click_chatter("Have %d overloaded and %d underloaded",overloaded.size(), underloaded.size());

    std::sort(overloaded.begin(), overloaded.end(), [](Pair<int,float> a, Pair<int,float> b) {return a.second < b.second; });

    std::sort(underloaded.begin(), underloaded.end(), [load](Pair<int,float> a, Pair<int,float> b) {return a.second > b.second; });

    while (overloaded.size() > 0) {
        int o = overloaded[overloaded.size() - 1].first;
        float load_o = overloaded[overloaded.size() - 1].second;
        overloaded.pop_back();
        int a_phys_id;
        float load_a = 0;
        if (balancer->_autoscale && balancer->_available_cpus.size() > 0) {
            a_phys_id = balancer->addCore();
        } else if (underloaded.size() > 0) {
		    a_phys_id = underloaded.back().first;
            load_a = underloaded.back().second;

            underloaded.pop_back();
        } else {
            click_chatter("Not enough CPUs available");
            break;
        }

        _last_movement[o] = now;
        _last_movement[a_phys_id] = now;

        click_chatter("Deflating cpu phys core %d to %d", o, a_phys_id);

        auto cache = flow_dir->get_flow_cache();
        if (!cache) {
            click_chatter("Flow dir is not initialized");
            assert(false);
        }


        HashMap<long, String> * rmap;
        rmap = cache->rules_map_by_core_id(o);
        int mig;
        if (_deflation_factor == 0) {
            mig = rmap->size() * ((load_o - load_a) / 2.0);
        } else {
            mig = rmap->size() / _deflation_factor;
        }

        if (mig == 0)
            mig = 1;
        //click_chatter("Rmap %p", rmap);
        click_chatter("Migrating %d flows / %d flows", mig , rmap->size());

        Vector<uint32_t> dlist;
        HashMap<long, String> nmap;
//      flow_dir->flow_cache()->delete_rule_by_global_id(rmap.keys());
//      flow_dir->flow_cache()->
        unsigned n = 0;
        auto it = rmap->begin();
        while (it != rmap->end()) {
            if (n++ < rmap->size() - mig) {
                it++;
                continue;
            }

            long rule_id = it.key();
            String rule = String(it.value());

//            int intid = cache->internal_from_global_rule_id(rule_id);
            //click_chatter("Deleting %d %d", intid,rule_id);
//            dlist.push_back(intid);
            String nrule = rewrite_id(rule, a_phys_id);
            nmap.insert(rule_id,nrule);
        //    click_chatter("New rule %s",nrule.c_str());
            it++;
        }

        // Delete the flow rules
  //      int status = flow_dir->flow_rules_delete(dlist.data(), dlist.size());

//        click_chatter("Deleted %d/%d", status, dlist.size(), true);
        int status = flow_dir->update_rules(nmap, true, a_phys_id);
        click_chatter("Update %d", status);

    }

    while (underloaded.size() >= 2 && balancer->_autoscale) {
        int remove_i = underloaded[underloaded.size() - 2].first;
        int with_i = underloaded[underloaded.size() - 1].first;
        click_chatter("Removing underloaded CPU %d!", remove_i);
        float load_r =  underloaded[underloaded.size() - 2].second;
        float load_i =  underloaded[underloaded.size() - 1].second;

        _last_movement[remove_i] = now;


        _last_movement[with_i] = now;
        underloaded.pop_back();
        underloaded.pop_back();
//        click_chatter("Core %d and %d", remove_i, with_i);
        if (load_r + load_i >= balancer->_overloaded_thresh) {
		click_chatter("ERROR : too much load, select another core");
                    continue;
        }
        balancer->removeCore(remove_i);
//        click_chatter("Core %d is now available", remove_i);

        HashMap<long, String> * rmap;

        HashMap<long, String> nmap;
        auto cache = flow_dir->get_flow_cache();

        rmap = cache->rules_map_by_core_id(remove_i);
//        click_chatter("Rmap %p", rmap);
        click_chatter("Migrating %d flows to core %d", rmap->size(), with_i);
//      flow_dir->flow_cache()->delete_rule_by_global_id(rmap.keys());
//      flow_dir->flow_cache()->

        auto it = rmap->begin();
        while (it != rmap->end()) {
            long rule_id = it.key();
            String rule = it.value();

            //int intid = cache->internal_from_global_rule_id(rule_id);
            nmap.insert(rule_id, rewrite_id(rule, with_i));
            it++;
        }
        int status = flow_dir->update_rules(nmap, true, with_i);
        if (status == 0)
            click_chatter("Failed to update rules");
//        click_chatter("Update %d", status);
    }
}


/**
 * RSS base
 */
BalanceMethodRSS::BalanceMethodRSS(DeviceBalancer* b, Element* fd) : BalanceMethodDevice(b,fd), _verifier(0) {
}

int BalanceMethodRSS::initialize(ErrorHandler *errh, int startwith) {
	int reta = _fd->get_rss_reta_size(_fd);
	click_chatter("Reta size %d", reta);
	if (reta <= 0)
		return errh->error("Device not initialized or RSS is misconfigured");
	if (_fd->get_rss_reta)
		_table = _fd->get_rss_reta(_fd);

	_table.resize(_reta_size);
	/*
    if (_table.size() < 128) {
	return errh->error("RSS reta table is %d long. It should be at least 128.", _table.size());
    }*/
	//We update the default reta to 0 to be sure it works
    for (int i = 0; i < _table.size(); i++) {
		_table[i] = 0;
	}
    _fd->set_rss_reta(_fd, _table);

    if (_is_dpdk) {

		int port_id = ((DPDKDevice*)_fd)->port_id;

		_rss_conf.rss_key = (uint8_t*)CLICK_LALLOC(128);
		_rss_conf.rss_key_len = 128; //This is only a max
		rte_eth_dev_rss_hash_conf_get(port_id, &_rss_conf);


		struct rte_flow_error error;
		rte_eth_dev_stop(port_id);
		//rte_eth_promiscuous_disable(port_id);
		int res = rte_flow_isolate(port_id, 1, &error);
		if (res != 0)
			errh->warning("Warning %d : Could not set isolated mode because %s !",res,error.message);

		rte_eth_dev_start(port_id);
    }

	for (int i = 0; i < _table.size(); i++) {

		_table[i] = i % startwith;
	}

    click_chatter("RSS initialized with %d CPUs", startwith);
    int err = BalanceMethodDevice::initialize(errh, startwith);
    if (err != 0)
        return err;

    _update_reta_flow = true;
    if (_is_dpdk) {
		if (!update_reta_flow(true)) {
			_update_reta_flow = false;
			if (_fd->set_rss_reta(_fd, _table) != 0)
                return errh->error("Neither flow RSS or global RSS works to program the RSS table.");
		} else
	        click_chatter("RETA update method is flow");
    } else {
	_update_reta_flow = false;

        if (_fd->set_rss_reta(_fd, _table) != 0)
            return errh->error("Cannot program the RSS table.");
    }
    if (!_update_reta_flow)  {
	click_chatter("RETA update method is global");
    }

    if (balancer->_manager) {
	balancer->_manager->init_assignment(_table);
    }
    return err;
}

void BalanceMethodRSS::rebalance(Vector<Pair<int,float>> load) {
    update_reta();
}


int BalanceMethodRSS::configure(Vector<String> &conf, ErrorHandler *errh)  {
	Element* e = 0;
	if (Args(balancer, errh).bind(conf)
			.read("VERIFIER", e)
			.read_or_set("RETA_SIZE", _reta_size, 128)
			.consume() < 0)
		return -1;
	if (e) {
		_verifier = (RSSVerifier*)e->cast("RSSVerifier");
		if (!_verifier)
			return errh->error("Verifier must be of the type RSSVerifier");
	}
	return 0;
}



/**
 * RSS - RoundRobin
 */
void MethodRSSRR::rebalance(Vector<Pair<int,float>> load) {
    for (int r = 0; r < _table.size(); r++) {
		_table[r] = (_table[r] + 1) % load.size();
    }
    update_reta();
};

/**
 * Piano RSS
 */
MethodPianoRSS::MethodPianoRSS(DeviceBalancer* b, Element* fd, String config) : BalanceMethodRSS(b,fd) {
}


int MethodPianoRSS::initialize(ErrorHandler* errh, int startwith) {
    if (BalanceMethodRSS::initialize(errh, startwith) != 0)
        return -1;
    load_tracker_initialize(errh);

    if (_counter_is_xdp) {
	click_chatter("Resizing counte to %d",_reta_size);
		_count.resize(_reta_size);
		_xdp_table_fd = ((XDPLoader*)_counter)->get_map_fd("count_map");
		if (!_xdp_table_fd)
			return errh->error("Could not find map !");
    }
    return 0;
}
#define EPSILON 0.0001f



class Problem
{ public:
	Vector<int> oid;
	Vector<int> uid;
	Vector<int> transfer;
	float min_cost;
	Vector<float> imbalance;
	float target;
	int N;

	Problem() : oid(), uid(), min_cost(FLT_MAX) {

	}

	bool computeSol() {
		float newload[N] = {0.0f};
		for (int c = 0; c < N; c++) { //Sum of imbalance for all cores
			newload[transfer[c]] += imbalance[c];
		}
		float imb = 0;
		for (int c = 0; c < N; c++) {
			imb += newload[c] * newload[c];
		}
		if (imb < min_cost) {
		    min_cost = imb;
			return true;
		}
		return false;
	}

	bool solve() {
		tryK(0);
	}

    /**
     * Problem of this method:
     * If all sending cores have quite big buckets, then correction may lead to all
     * cores having too low corrected transfer so no bucket will move. Moving some
     * of them would probably be best.
     *
     * One solution is to consider all buckets that fit the correction, then select a subset
     * that min the invariance of the moving set.
     *
     * But that's again a bad minimization
     *
     * So if no buckets were moved, we do a second pass without the correction
     */
	Vector<float> fixedTransfert() {
        //See above
		float newload[N] = {0.0f};
		for (int c = 0; c < N; c++) { //Sum of imbalance for all cores
			newload[transfer[c]] += imbalance[c];
		}
		Vector<float> fT;
		fT.resize(N,1.0);

		//Eg c0 had  0.15 of imbalance //underloaded, newload -0.03
		//   c1 had -0.10 of imbalance //overloaded, newload 0
		//   c2 has -0.08 of imbalance  //overloaded, newload 0

		// c1 and c2 will give all to c0, but c0 will become overloaded by -0.03.
		//Let's share the final imbalance to -0.03/3, so every of those are at -0.01
						//
		for (int c = 0; c < N; c++) { //C is the receiving core, ie, for each receiving core C
			if (newload[c] < 0) {

				int nsources = 0;
				float imb = 0; //Imbalance of all senders

				for (int j = 0; j < N; j++) {
					if (j == c)
						continue;
					if (transfer[j] == c) {
						nsources++;
						imb += imbalance[j]; // += -0.10 += -0.08  -> -0.18
					}
				}
				if (nsources == 0)
					continue;
				float avgimb = (imbalance[c] + imb) / (nsources + 1); //Average imbalance Eg 0.15 - 0.18 / 3 = -0.03 / 3 = -0.01



				//If we just give the factor to all sources, then some would be left with more imbalance than others
				// We want all the sources to have exactly the same amount of final imbalance ->
				//   avgimb
				if (avgimb < -0.001) {
					for (int j = 0; j < N; j++) {
						if (j == c)
							continue;
						if (transfer[j] == c) {

							// c1 :: -0.01 / -0.10 -> 0.1 -> 0.9
							// c2 :: -0.01 / -0.08 -> 0.125 -> 0.875
							fT[j] = 1 - avgimb / imbalance[j];
						}
					}
				}

				//Final example : j gives 0.875*-0.2 == -0.175 --> n = -0.025
				//c reveives -0.175 + 0.15 = -0.025

			}
		}
		return fT;

	}

private:
	bool tryK(int i)
	{
		if (i == oid.size()) {
			return computeSol();
		}
		int best = -1;
		for (int j = 0; j < uid.size(); j++) {
				transfer[oid[i]] = uid[j];
				if (tryK(i + 1)) {
					best = j;
				}
		}
		if (best > -1) {
			transfer[oid[i]] = uid[best];

			return true;
		}
		return false;
	}
};

class BucketMapProblem
{ public:
	Vector<int> transfer; //existing core id for each buckets
	Vector<float> imbalance; //Imbalance for each existing cores (no holes)
	Vector<float> buckets_load; //Load for each bucket


	BucketMapProblem(int nbuckets, int ncpu) : min_cost(FLT_MAX) {
		transfer.resize(nbuckets);
		buckets_load.resize(nbuckets);
		imbalance.resize(ncpu);
	}

	/*This is WAY too long
	bool tryK(int i)
	{

		if (i == buckets_load.size()) { //All buckets have been mapped
			float newload[imbalance.size()] = {0.0f};
			for (int c = 0; c < buckets_load.size(); c++) { //Sum of imbalance for all cores
				newload[transfer[c]] += buckets_load[c];
			}
			float imb = 0;
			for (int c = 0; c < imbalance.size(); c++) {
				float coreload = imbalance[c] + newload[c];
				imb += (coreload) * (coreload);
			}
			if (imb < min_cost) {
			    min_cost = imb;
				return true;
			}
			return false;
		}
		int best = -1;
		for (int j = 0; j < imbalance.size(); j++) {
//				click_chatter("Level %d:%d",i,j);
				transfer[i] = j;
				if (tryK(i + 1)) {
					best = j;
				}
		}
		if (best > -1) {
			transfer[i] = best;

			return true;
		}
		return false;
	}

	void solve() {
		return tryK(0);
	}
	*/

	void solve() {
		typedef struct {
			int id;
			float load;
		} bref;
		auto cmp = [](bref left, bref right) { return left.load > right.load; };
		std::priority_queue<bref, std::vector<bref>, decltype(cmp)> q(cmp);
		for(int i = 0; i < buckets_load.size(); i++) {
			float f = buckets_load[i];

			q.push(bref{.id = i,.load =   f});
		}

		auto cmpc = [](bref left, bref right) { return left.load < right.load; };
		std::priority_queue<bref, std::vector<bref>, decltype(cmp)> cores(cmp);
		for(int i = 0; i < imbalance.size(); i++) {
			float f = imbalance[i];
			//click_chatter("Core %d should receive %f load",i,f);
			cores.push(bref{.id = i,.load = - f});//negative of imbalance, so we should reach a nice 0 everywhere by adding some load
		}


		while (!q.empty()) {
			bref t = q.top();
			q.pop();
			bref c = cores.top();
			cores.pop();
			transfer[t.id] = c.id;
			c.load += t.load;
			cores.push(c);
		}

	}

private:
	float min_cost;
};

double myvfunc(const std::vector<double> &x, std::vector<double> &grad, void *my_func_data)
{
    if (!grad.empty()) {
        grad[0] = 0.0;
        grad[1] = 0.5 / sqrt(x[1]);
    }
    return sqrt(x[1]);
}

struct cref {
	float load;
	int id;
};
class Compare
{
public:
    bool operator() (cref a, cref b)
    {
        return a.load < b.load;
    }
};

/**
 * Take some buckets of a list to make the load of a set of CPU reach target
 */
class BucketMapTargetProblem
{ public:
	Vector<int> transfer; //existing core id for each buckets
	Vector<float> target; //Target for destination
	Vector<float> max; //Max load to take from overloaded
	Vector<int> buckets_max_idx; //Load for each bucket
	Vector<float> buckets_load; //Load for each bucket


	BucketMapTargetProblem(int nbuckets, int nunderloaded, int noverloaded ) : min_cost(FLT_MAX) {
		transfer.resize(nbuckets,-1);
		buckets_load.resize(nbuckets);
		buckets_max_idx.resize(nbuckets);
		target.resize(nunderloaded);
		max.resize(noverloaded);
	}

	void solve(DeviceBalancer* balancer) {
		click_chatter("Assigning %d buckets to %d targets :", buckets_load.size(), target.size());
		if (unlikely(balancer->_verbose)) {
			for (int i = 0; i < max.size(); i++) {
				click_chatter("Overloaded %d , %f",i, max[i]);
			}
			for (int i = 0; i < target.size(); i++) {
				click_chatter("Underloaded %d , %f",i, target[i]);
			}
			for (int i = 0; i < buckets_load.size(); i++) {
				click_chatter("Bucket %d , %f, oid %d",i, buckets_load[i], buckets_max_idx[i]);
			}
		}

		float target_imbalance = 0.01;//(o_load - u_loadà ;
		//Build priority for each overloaded
		//

		auto cmplt = [](cref left, cref right) { return left.load > right.load; };
		auto cmpgt = [](cref left, cref right) { return left.load < right.load; };

		float overload_allowed =  -0.01; //Take just not enough load of overloaded
		float underload_allowed = -0.01; //Give just not enough to underloaded
		float imbalance_u = 0;
		float imbalance_o = 0;
		float oa_min;
		float ua_min;
		float bottom_oa = overload_allowed;
		float top_oa;
		float bottom_ua = underload_allowed;
		float top_ua;
		float square_imbalance = 0;
		float last_sq = FLT_MAX;
		float min_sq = FLT_MAX;
		float bottom_sq = FLT_MAX;
		float top_sq = FLT_MAX;
		int run = 1;
		float m = -1;
		int phase;
#define max_runs 10

		while(true) {
			imbalance_u = 0;
			imbalance_o = 0;
			square_imbalance = 0;
			typedef std::priority_queue<cref, std::vector<cref>, Compare> bstack;
			std::vector<bstack> bstacks;
			bstacks.resize(max.size(),bstack());
			for (int i = 0; i < buckets_load.size(); i++) {
				bstacks[buckets_max_idx[i]].push(cref{buckets_load[i],i});
				transfer[i] = -1;
			}

			std::priority_queue<cref, std::vector<cref>, decltype(cmpgt)> overloaded(cmpgt);
			for (int i = 0; i < max.size(); i++) {
				overloaded.push(cref{max[i],i});
			}

			std::priority_queue<cref, std::vector<cref>, decltype(cmpgt)> underloaded(cmpgt);
			for (int i = 0; i < target.size(); i++) {
				underloaded.push(cref{target[i],i});
			}



			while (!overloaded.empty() && !underloaded.empty()) {
				next_core:
				//Select most overloaded core
				cref o = overloaded.top();
				overloaded.pop(); //Will be added bacj

				//Select biggest bucket
				bstack& buckets = bstacks[o.id];
				cref bucket = buckets.top();
				buckets.pop();//Bucket is removed forever

				//Assign its most used buckets
				std::vector<cref> save;
				while (!underloaded.empty()) {
					cref u = underloaded.top();
					underloaded.pop();

					if (unlikely(balancer->_verbose))
						click_chatter("U%d load %f",u.id, u.load);
					if (bucket.load < u.load + underload_allowed) {
						u.load -= bucket.load;
						o.load -= bucket.load;
						transfer[bucket.id] = u.id;

						if (unlikely(balancer->_verbose))
							click_chatter("Bucket %d to ucore %d, bucket load %f", bucket.id, u.id, bucket.load);
						if (u.load > target_imbalance) {
							underloaded.push(u);
						} else {
							imbalance_u += abs(u.load);
							square_imbalance += u.load * u.load;
							if (unlikely(balancer->_verbose))
									click_chatter("Underloaded core %d is now okay with %f load", u.id, u.load);
						}
						goto bucket_assigned;
					} else {
						save.push_back(u);
					}
				}

				if (unlikely(balancer->_verbose))
					click_chatter("Bucket %d UNMOVED, load %f", bucket.id, o.id, bucket.load);

				bucket_assigned:
				while (!save.empty()) {
					underloaded.push(save.back());
					save.pop_back();
				}

				if (o.load > - overload_allowed
						&& !buckets.empty()) {
					overloaded.push(o);
				} else {
					imbalance_o += abs(o.load);
					square_imbalance += o.load * o.load;
					if (unlikely(balancer->_verbose))
						click_chatter("Overloaded core %d is now okay with %f load. Empty : %d", o.id, o.load,buckets.empty());
				}
			}
			while (!overloaded.empty()) {
				auto o = overloaded.top();
				imbalance_o += abs(o.load);
				square_imbalance += o.load * o.load;
				overloaded.pop();
			}
			while (!underloaded.empty()) {
				auto u = underloaded.top();
				imbalance_u += abs(u.load);
				square_imbalance += u.load * u.load;
				underloaded.pop();
			}

			click_chatter("Imbalance at run %d : %f-%f %f-%f, square %f, m %f",run,imbalance_o,imbalance_u,overload_allowed,underload_allowed,square_imbalance, m);

			if (run == max_runs || square_imbalance < target_imbalance) break;

			if (run == 1) {
				overload_allowed = imbalance_o;
				underload_allowed = imbalance_u;
				//bottom_square_imbalance = square_imbalance;
				//min_sq = square_imbalance; //min is the sq for the min allowed, not the min seen sq
				phase = 1; //searching top
			} else {

				if (square_imbalance < min_sq) {
					oa_min = overload_allowed;
					ua_min = underload_allowed;
				}


				if (phase == 1) {
					if (square_imbalance <= last_sq) { //Continue finding top
						overload_allowed = overload_allowed + overload_allowed / 2;
						underload_allowed = underload_allowed + underload_allowed / 2;
					} else {
						phase = 2;
						top_oa = overload_allowed;
						top_ua = underload_allowed;
						m = 0.5;
						overload_allowed = overload_allowed / 2;
						underload_allowed = underload_allowed / 2;
					}
				} else if (phase == 2) { //Searching left inflation
					if (square_imbalance <= last_sq) {
						m = (0 + m) / 2; //continue left;
						if (m < 0.01)
							break; //Border is the max
						overload_allowed = bottom_oa + (top_oa - bottom_oa) * m;
						underload_allowed = bottom_ua + (top_ua - bottom_ua) * m;
						//Either we still need to descend left
						//Or we hit the left inflation and will need to descend right afterwards
					} else if (square_imbalance > last_sq) { //we found a new bottom
						bottom_oa = overload_allowed;
						bottom_ua = underload_allowed;
						bottom_sq = square_imbalance;
						phase = 3;
						m = 0.5;
						overload_allowed = (bottom_oa + top_oa) / 2;
						underload_allowed = (bottom_ua + top_ua) / 2;
					}
				} else if (phase == 3) { //Searching right inflation
					if (square_imbalance <= last_sq) {
						m = (1 + m) / 2; // continue right
						if (m > 0.99)
								break; //Border is the min
						overload_allowed = bottom_oa + (top_oa - bottom_oa) * m;
						underload_allowed = bottom_ua + (top_ua - bottom_ua) * m;
					} else 	if (square_imbalance > last_sq) { //we found a new top
						phase = 2;
						top_oa = overload_allowed;
						top_ua = underload_allowed;
						top_sq = square_imbalance;
						m = 0.5;
						overload_allowed = bottom_oa + (top_oa - bottom_oa) * m;
						underload_allowed = bottom_ua + (top_ua - bottom_ua) * m;
					}
				}
				if (top_sq == bottom_sq && bottom_sq == square_imbalance && square_imbalance == min_sq) {
					break;
				}


				/*
				if (square_imbalance <= last_sq) { // Continue in this direction

				} else {
					//invert_direction, set max according to dir
					if (dir > 0) {
						top_oa = overload_allowed;
						top_ua = underload_allowed;
						dir = -1;
						overload_allowed = (bottom_oa + overload_allowed) / 2;
						underload_allowed = (bottom_ua + underload_allowed) / 2;
					} else {
						bottom_oa = overload_allowed;
						bottom_ua = underload_allowed;
						dir = 1;
						overload_allowed = (top_oa + overload_allowed) / 2;
						underload_allowed = (top_ua + underload_allowed) / 2;
						bottom_square_imbalance = square_imbalance;
					}
					if (bottom_square_imbalance)
				}*/
				/*if (square_imbalance == last_sq)
					break;*/


				/*if (square_imbalance > last_sq) {
					dir = -dir;
					n_change++;
					if (nchange > 2)
				} else if (square_imbalance < last_sq) {
					n_change = 0;
				}
				overload_allowed = overload_allowed + dir * (overload_allowed / 2);
				underload_allowed = underload_allowed + dir * (underload_allowed / 2);*/
			}
			click_chatter("Phase %d", phase);
			run++;
			last_sq = square_imbalance;
			if (run == max_runs) {
				if (square_imbalance != min_sq) {
					overload_allowed = oa_min;
					underload_allowed = ua_min;
				} else
					break;
			}

		}


/*
 * Not sure about the gradient
		nlopt::opt opt(nlopt::LD_MMA, 2);
		std::vector<double> lb(2);
		lb[0] = -HUGE_VAL; lb[1] = 0;
		opt.set_lower_bounds(lb);
		opt.set_min_objective(myfunc, NULL);
		my_constraint_data data[2] = { {2,0}, {-1,1} };
		opt.add_inequality_constraint(myconstraint, &data[0], 1e-8);
		opt.add_inequality_constraint(myconstraint, &data[1], 1e-8);
		opt.set_xtol_rel(1e-4);
		std::vector<double> x(2);
		x[0] = 1.234; x[1] = 5.678;
		double minf;

		try{
		    nlopt::result result = opt.optimize(x, minf);
		    std::cout << "found minimum at f(" << x[0] << "," << x[1] << ") = "
		        << std::setprecision(10) << minf << std::endl;
		}
		catch(std::exception &e) {
		    std::cout << "nlopt failed: " << e.what() << std::endl;
		}*/
		/*
        lsint weights[] = {10, 60, 30, 40, 30, 20, 20, 2};
        lsint values[] = {1, 10, 15, 40, 60, 90, 100, 15};
        lsint knapsackBound = 102;

        // Declares the optimization model.
        LocalSolver localsolver;
        LSModel model = localsolver.getModel();

        // 0-1 decisions
        LSExpression x[8];
        for (int i = 0; i < 8; i++)
            x[i] = model.boolVar();

        // knapsackWeight <- 10*x0 + 60*x1 + 30*x2 + 40*x3 + 30*x4 + 20*x5 + 20*x6 + 2*x7;
        LSExpression knapsackWeight = model.sum();
        for (int i = 0; i < 8; i++)
            knapsackWeight += weights[i]*x[i];

        // knapsackWeight <= knapsackBound;
        model.constraint(knapsackWeight <= knapsackBound);

        // knapsackValue <- 1*x0 + 10*x1 + 15*x2 + 40*x3 + 60*x4 + 90*x5 + 100*x6 + 15*x7;
        LSExpression knapsackValue = model.sum();
        for (int i = 0; i < 8; i++)
            knapsackValue += values[i]*x[i];

        // maximize knapsackValue;
        model.maximize(knapsackValue);

        // close model, then solve
        model.close();

        // Parameterizes the solver.
        localsolver.getParam().setTimeLimit(1);
        localsolver.solve();*/

	/*	glp_prob *lp;
		lp = glp_create_prob();
		glp_set_prob_name(mip, "sample");
		//int x[target.size()]
		int* x = transfer.data();
		glp_set_obj_dir(lp, GLP_MIN);

		//glp_add_rows(lp,buckets_load.size());

		glp_add_cols(lp,buckets_load.size());
		for (int i = 0; i < buckets_load.size(); i++) {
			//glp_set_col_name(mip, 1, "x1");
			// glp_set_col_bnds(mip, 1, GLP_DB, 0.0, 40.0);
			 // glp_set_obj_coef(mip, 1, 1.0);
			  //glp_set_col_name(lp, 1, "b1");
			  glp_set_col_bnds(lp, i, GLP_DB, 0.0, 1.0);
			  glp_set_obj_coef(lp, i, buckets_load.size());
			  glp_set_col_kind(lp, i, GLP_IV);
		}

		glp_iocp parm;
		  glp_init_iocp(&parm);
		  parm.presolve = GLP_ON;
		  int err = glp_intopt(mip, &parm);*/

	}

private:
	float min_cost;
};



class Load { public:

	Load() : Load(-1) {

	}

	Load(int phys_id) : cpu_phys_id(phys_id), load(0), high(false), npackets(0), nbuckets(0), nbuckets_nz(0)  {

	}
	int cpu_phys_id;
	float load;
	bool high;
	unsigned long long npackets;
    unsigned nbuckets;
    unsigned nbuckets_nz;
};

void MethodPianoRSS::rebalance(Vector<Pair<int,float>> rload) {
    click_jiffies_t now = click_jiffies();

    if (_counter_is_xdp) {
	click_chatter("Reading XDP table");
	int cpus = balancer->_max_cpus;
	unsigned int nr_cpus = bpf_num_possible_cpus();
	    uint64_t values[nr_cpus];
	    for (uint32_t key = 0; key < _count.size(); key++) {
	        if (bpf_map_lookup_elem(_xdp_table_fd, &key, values)) {
			click_chatter("XDP lookup failed");
			continue;
	        }
			uint64_t tot = 0;
//		    for (int i = 0; i < nr_cpus; i++) {
				tot += values[_table[key]];
//			if (values[i] && i != _table[key])
			//if (values[i] != 0)
//				click_chatter("BPF map ha value on core %d, but RSS is assigned to %d. Val %d",i,_table[key],values[i]);
  //              tot += values[i];
//			}
//            tot -= _count[key].count;

			_count[key].count  = tot;
			uint64_t var =  (_count[key].variance  / 3) + (2 * tot / 3);
			_count[key].variance  = var;
			if (unlikely(tot != 0 && balancer->_verbose))
                click_chatter("Key %d core %d val %d, var %f",key,_table[key], tot,(float)min(tot,var) / (float)max(tot,var));
	    }
    }
//    return;

#define get_node_count(i) ((_counter_is_xdp)?_count.unchecked_at(i).count:((AggregateCounterVector*)_counter)->find_node_nocheck(i).count)
#define get_node_variance(i) ((_counter_is_xdp)?_count.unchecked_at(i).variance:((AggregateCounterVector*)_counter)->find_node_nocheck(i).variance)
#define get_node_moved(i) (_count.unchecked_at(i).moved)

/*
    uint64_t total_packets = 0;
    for (int i = 0; i < _table.size(); i++) {
	total_packets += _counter->find_node(i)->count;
    }
*/
    Timestamp begin = Timestamp::now_steady();
    float _min_load = 0.01;
    float _threshold_force_overload = 0.90;
    float _load_alpha = 1;
    //float _high_load_threshold = 0.;
    const float _imbalance_threshold = 0.01;
    assert(_imbalance_threshold <= (_threshold / 2) + EPSILON); //Or we'll always miss
    //Vector<float> corrections;
    //corrections.resize(click_max_cpu_ids(), 0);

    float suppload = 0;
    Problem p;
    float totalload = 0;
    int min_core = -1;
    float min_core_load = 1;
    float max_core_load = 0;
    bool has_high_load = false;
    Vector<Load> load(rload.size(),Load());
    Vector<unsigned> map_phys_to_id(click_max_cpu_ids(), -1);
    for (int j = 0; j < rload.size(); j++) {
        int cpuid = rload[j].first;
        float load_current = rload[j].second;
        //assert(load_current <= 1);
        if (load_current < min_core_load)
		min_core = j;
        load[j].cpu_phys_id = cpuid;
        if (_past_load[cpuid] == 0)
		load[j].load = load_current;
        else
		load[j].load = load_current * _load_alpha + _past_load[cpuid] * (1.0-_load_alpha);
        _past_load[cpuid] = load[j].load;
        float diff = _target_load - load[j].load; // >0 if underloaded diff->quantity to add
        suppload += diff;
        if (abs(diff) <= _threshold)
		    diff = 0;
        if (load[j].load > _target_load) {
			has_high_load = true;
			load[j].high = true;
        }
        if (load[j].load > max_core_load) {
		max_core_load = load[j].load;
        }
        map_phys_to_id[load[j].cpu_phys_id] = j;

        //corrections[cpuid] = diff;
        //suppload += diff;
        totalload+=load[j].load;
    }

    p.N = load.size();
    p.target = totalload / (float)p.N;

    //suppload = (p.N *_target_load) - totalload; //Total power that we have in excess

    float var_protect = (p.N + 1) * (max(0.05f,(max_core_load-_target_load)));
    if (unlikely(balancer->_verbose))
	    click_chatter("Target %f. Total load %f. %d cores. Suppload %f, Var %f", p.target, totalload, p.N, suppload, var_protect);

    p.imbalance.resize(p.N);

    //Count the number of packets for this core
	for (int j = 0; j < _table.size(); j++) {
		unsigned long long c = get_node_count(j);
		load[map_phys_to_id[_table[j]]].npackets += c;
		load[map_phys_to_id[_table[j]]].nbuckets += 1;
        if (c > 0)
		    load[map_phys_to_id[_table[j]]].nbuckets_nz += 1;
	}

    /**
     * Scaling. If we need more cores, we just add it and the imbalance mapping will do the rest
     * If we remove a core, we minimize the selection of buckets of the removed core -> existing cores to
     * diminish the imbalance already
     */
    if (balancer->_autoscale) {
        if (suppload > (1 + (1 - _target_load) + var_protect)  && p.N > 1) { // we can remove a core
            if (unlikely(balancer->_verbose)) {
                click_chatter("Removing a core (target %f)",_target_load);
            }

            click_chatter("Removing core %d", min_core);
            int min_core_phys = load[min_core].cpu_phys_id;
            unsigned long long totcount = load[min_core].npackets;
            load[min_core] = load[load.size() - 1];
            load.pop_back();
            balancer->removeCore(min_core_phys);


            //Count the number of packets for this core
            Vector<int> buckets_indexes;
			for (int j = 0; j < _table.size(); j++) {
				if (_table[j] == min_core_phys) {
					buckets_indexes.push_back(j);
				}
			}
			if (unlikely(balancer->_verbose))
				click_chatter("Removing %d buckets of core %d", buckets_indexes.size(), min_core);


            BucketMapProblem cp(buckets_indexes.size(), load.size());

			//Compute load for each buckets
			for (int i = 0; i < buckets_indexes.size(); i++) {
				int j = buckets_indexes[i];
				double c = get_node_count(j);
				cp.buckets_load[i] = (((double) c / (double)totcount)*min_core_load);
				//click_chatter("Bucket %d (%d) load %f", i ,j, cp.buckets_load[i]);
			}

			//Add imbalance of all cores
			p.imbalance.resize(p.N-1);
            p.N = p.N-1;
            p.target = totalload / (float)p.N;

            //Fix imbalance without the removed core
            for (int i = 0; i < p.N; i++) {
				//Imbalance is positive if the core should loose some load
				p.imbalance[i] = 0 * (1.0-_imbalance_alpha) + ( (p.target - load[i].load) * _imbalance_alpha); //(_last_imbalance[load[i].first] / 2) + ((p.target - load[i].second) / 2.0f);
            }
            cp.imbalance = p.imbalance;

            click_chatter("Solving problem...");
            //Solve problem ; map all buckets of removed core to other cores
            cp.solve();

            //assert(false); //unfinished

			for (int i = 0; i < buckets_indexes.size(); i++) {
				int j = buckets_indexes[i];
				_table[j] = cp.transfer[i];
			//	click_chatter("Bucket %d (%d) -> core %d", i ,j, cp.transfer[i]);
			}

			Timestamp t = Timestamp::now_steady();
			click_chatter("Solution computed in %d usec", (t-begin).usecval());
			update_reta();
			if (!_counter_is_xdp)
				((AggregateCounterVector*)_counter)->advance_epoch();
			return;

        } else if (suppload < -0.1) { //We need a new core because the total load even with perfect spread incurs 10% overload
            if (unlikely(balancer->_verbose))
                click_chatter("Adding a core as load is %f");
            int a_phys_id = balancer->addCore();
            if (a_phys_id == -1) {
                if (unlikely(balancer->_verbose))
                    click_chatter("Not enough cores...");
            } else {
				p.imbalance.resize(p.imbalance.size() + 1);
				p.N = p.N+1;
				p.target = totalload / (float)p.N;
				int aid = load.size();
				load.push_back(Load(a_phys_id));
				map_phys_to_id[a_phys_id] = aid;
            }
        }
    } else {
	    if (p.target <  _min_load && !_threshold_force_overload) {
		    click_chatter("Underloaded, skipping balancing");
	}
    }

    /*
     * Dancers
     * We move the bucket that have more load than XXX 50% to other cores
     */
     if (has_high_load) {
	    if (unlikely(balancer->_verbose))
		    click_chatter("Has high load !");
		for (int j = 0; j < _table.size(); j++) {
			unsigned long long c = get_node_count(j);
			unsigned phys_id = _table[j];
			unsigned id = map_phys_to_id[phys_id];
			if (!load[id].high || load[id].npackets == 0 || load[id].nbuckets <= 1)
				continue;
			unsigned long long pc = (c * 1024) / load[id].npackets;
			float l = ((float)pc / 1024.0) * (load[id].load);
			if (l > 0.5) {
				click_chatter("Bucket %d (cpu id %d) is a dancer ! %f%% of the load", j, id, l);

				float min_load = 1;
				int least = -1;
				for (int i = 0; i < load.size(); i++) {
					if (load[i].load < min_load) {
						min_load = load[i].load;
						least = i;
					}
				}

				if (unlikely(balancer->_verbose))
					click_chatter("Moving to %d", least);

                //We fix the load here. So the next step of the algo will rebalance the other flows
                // as if nothing happened
				load[least].load += l;
				load[least].npackets += c;
				load[least].nbuckets+=1;
				load[least].nbuckets_nz+=1;
				load[id].load -= l;
				load[id].nbuckets-=1;
				load[id].nbuckets_nz-=1;
				load[id].npackets-= c;
				_table[j] = load[least].cpu_phys_id;
			}
		}
    }

	float total_imbalance = 0;

	//Compute the imbalance
	for (unsigned i = 0; i < p.N; i++) {
		int cpuid = load[i].cpu_phys_id;
		if (_moved[cpuid]) {
			p.imbalance[i] = 0;
			_moved[cpuid] = false;
			continue;
		}

		p.imbalance[i] = 0 * (1.0-_imbalance_alpha) + ( (p.target - load[i].load) * _imbalance_alpha); //(_last_imbalance[load[i].first] / 2) + ((p.target - load[i].load) / 2.0f);i
        total_imbalance += abs(p.imbalance[i]);

		if (p.imbalance[i] > _threshold) {
            if (unlikely(balancer->_verbose))
		click_chatter("Underloaded %d is cpu %d, imb %f, buckets %d",p.uid.size(),i, p.imbalance[i],load[i].nbuckets_nz);
			p.uid.push_back(i);
		}
		else if (p.imbalance[i] < - _threshold) {
            if (load[i].nbuckets_nz == 0) {
                click_chatter("WARNING : A core is overloaded but has no buckets !");
            } else if (load[i].nbuckets_nz > 1) { //Else there is nothing we can do
		if (unlikely(balancer->_verbose))
			click_chatter("Overloaded %d is cpu %d, imb %f, buckets %d",p.oid.size(),i, p.imbalance[i],load[i].nbuckets_nz);
			    p.oid.push_back(i);
            } else
		if (unlikely(balancer->_verbose))
		    click_chatter("Overloaded cpu %d, imb %f, buckets %d",i, p.imbalance[i],load[i].nbuckets_nz);
		}
    }


	bool moved = false;
    /**
     * Re-balancing
     * We minimize the overall imbalance by moving some buckets from each cores to other cores
     */
    if (p.oid.size() > 0) {

#ifndef BALANCE_TWO_PASS

        //Count the number of packets for this core
        /*
         * We rebalance all buckets of overloaded to the set of underloaded
         */
        unsigned long long all_overloaded_count = 0;
        for (int u = 0; u < p.oid.size(); u++) {
            all_overloaded_count += load[p.oid[u]].npackets;
        }
        struct Bucket{
            int oid_id;
            int bucket_id;
            int cpu_id;
        };
        Vector<Bucket> buckets_indexes;
		for (int j = 0; j < _table.size(); j++) {
			if (get_node_count(j) == 0) continue;
			for (int u = 0; u < p.oid.size(); u++) {
				if (map_phys_to_id[_table[j]] == p.oid[u]) {
					buckets_indexes.push_back(Bucket{.oid_id =u,.bucket_id =j, .cpu_id = p.oid[u]});
					break;
				}
			}
		}

		BucketMapTargetProblem pm(buckets_indexes.size(), p.uid.size(), p.oid.size());

		//Compute load for each buckets
		for (int i = 0; i < buckets_indexes.size(); i++) {
			Bucket& b = buckets_indexes[i];
			int j = b.bucket_id;
			double c = get_node_count(j);
			pm.buckets_load[i] = (((double) c / (double)load[buckets_indexes[i].cpu_id].npackets)*load[b.cpu_id].load);
			pm.buckets_max_idx[i] = buckets_indexes[i].oid_id;
			//click_chatter("Bucket %d (%d) load %f", i ,j, cp.buckets_load[i]);
		}

		for (int i = 0; i < p.uid.size(); i++) {
			pm.target[i] = p.imbalance[p.uid[i]];
		}

		for (int i = 0; i < p.oid.size(); i++) {
			pm.max[i] = -p.imbalance[p.oid[i]];
		}


		pm.solve(balancer);

		for (int i = 0; i < pm.buckets_load.size(); i++) {
			int to_uid = pm.transfer[i];
			if (to_uid == -1) continue;

            if (unlikely(balancer->_verbose))
            click_chatter("B idx %d to ucpu %d",i,to_uid);
			int to_cpu = p.uid[to_uid];
			int from_cpu = p.oid[pm.buckets_max_idx[i]];
			p.imbalance[from_cpu] += pm.buckets_load[i];
			p.imbalance[to_cpu] -= pm.buckets_load[i];

            if (unlikely(balancer->_verbose))
			click_chatter("Move bucket %d to core %d",buckets_indexes[i].bucket_id,load[to_cpu].cpu_phys_id);
			_table[buckets_indexes[i].bucket_id] = load[to_cpu].cpu_phys_id;
			_moved[load[to_cpu].cpu_phys_id] = true;
			moved = true;

		}

#else

        //Set problem parameters
        float min_cost = p.N * p.N;
		p.transfer.resize(p.N);
		for (unsigned i = 0; i < p.N; i++) {
			p.transfer[i] = i;//Default is to transfer load to itself
		}
		p.min_cost = p.N;

        //Solve assignment of imbalance of overloaded -> underloaded
        p.solve();

        //We "correct" the imbalance. Indeed even from the best solution, moving overloaded
        // to underloaded may lead to a unoptimal local scenario, we could move just
        // a little less of every overloaded core to the single underloaded
        Vector<float> fT = p.fixedTransfert();


		if (unlikely(balancer->_verbose))
			click_chatter("Transfer solution for %d uid, %d oid:",p.uid.size(), p.oid.size());

        //We start at a random core index, as the first cores to be balanced may be advantaged
	    int iRand = click_random();

		for (unsigned iOffset = 0; iOffset < p.N; iOffset++) {
            unsigned i = (iRand + iOffset) % p.N;
			if (unlikely(balancer->_verbose))
				click_chatter("Core %d (load %f, imbalance %f (corr %f)) -> %d (load %f)", i, load[i].load, p.imbalance[i], fT[i], p.transfer[i], load[p.transfer[i]].load);
			if ((unsigned)p.transfer[i] == i)
				continue;

			/*
		    uint64_t total_packets = 0;
		    for (int i = 0; i < _table.size(); i++) {
			total_packets += _counter->find_node(i)->count;
		    }*/
				/**
				 * Minimize the number of state transfer, or bucket transfer
				 *for now, just add up the number of packets
				 */
			uint64_t core_tot = 0;
			double core_tot_min = 0;
			int cn = 0;
			/*auto cmp = [](Pair<int,float> left, Pair<int,float> right) { return left.second > right.second;};
			std::priority_queue<Pair<int,float>, std::vector<Pair<int,float> >, decltype(cmp)> q(cmp);*/
			/*struct BucketInfo {
				int index;
				float var;
				int flows;
			}*/
			//Vector<BucketInfo> q;

            //We find all buckets of this core, starting at a random point, so we do not
            //always move the same ones
			Vector<Pair<int,float> > q; //List of (i core) bucket->load
			int j = click_random() % _table.size();
			for (int r = 0; r < _table.size(); r++) {
				if (++j == _table.size())
					j = 0;
				if (map_phys_to_id[_table[j]] == i) {

					uint64_t cu = get_node_count(j);
					core_tot += cu;
                    double c = cu;
                    double p = get_node_variance(j);
					double var;
					double m = max(c,p);
					if (m == 0) {
						var = 0.2;
					} else
						var = min(c,p) / max(c,p);
					if (var < 0.2)
						var = 0.2;

					core_tot_min += c * var;
					//if (var < 1) SEE BELOW
					q.push_back(Pair<int,float>(j, var));
					cn++;
				}
			}

reagain:
			int n = 0;
            int nSkipped = 0;
			float tot_bload = 0;
			float tot_bpc = 0;
			int qsz = q.size();

			/**
			 * We cannot leave the core with only a high variance bucket
			 * We want :
			 * - The splits to have an equal amount of variance -> we just take bucket in order, it's like random
			 * - Not overflow the imbalance
			 */
			bool miss = false;
			if (!q.empty())
                if (core_tot_min == 0) {
                    click_chatter("WARNING : Core has not seen any packet... But has %f load.", load[i].load);
                    continue;
                }
			while (!q.empty()) {
				int j = q.back().first;
				float var = q.back().second;
				q.pop_back();

				double c = (double)get_node_count(j);

				double bpc = ((c * var) / core_tot_min); //Percentage of all flows seen (minimum load observed)
                //click_chatter("ctm %f bpc %f, load %f",core_tot_min,bpc,load[i].load);
				assert( bpc != NAN);
				double bload =  bpc * load[i].load; //How much CPU load this bucket represents
				if (bload > -p.imbalance[i] * fT[i] || bload > 0.5 ) { //Never over balance! Also dancers are handled separately
					//if (bload > _threshold)
					if (!miss) {
                        //if (unlikely(balancer->_verbose))
						//click_chatter("Trying to fill %f of imbalance", p.imbalance[i]*fT[i]);
						miss =true;
					}
					if (unlikely(balancer->_verbose > 1)) {
						click_chatter("Skipping bucket %d var %f, load %f%% (min %f%% of bucket, %f%% real), imbalance %f, corrected imbalance %f, at least %d flows",j, var, bload*100, bpc*100, ((c * 100.0)/core_tot),p.imbalance[i], p.imbalance[i]*fT[i],-1);//get_node(j).flows);
                    }
                    if (bload < -p.imbalance[i] && bload < 0.5)
                        nSkipped ++;
                        /*
						left_max += bload * (1 + var);
						left += bload;
						left_min += bload / (1 + var);*/

					continue;
				} else
					miss = false;
				if (p.imbalance[i]*fT[i] > -_imbalance_threshold)
					break;


				tot_bload += bload;
				tot_bpc += bpc;
				//click_chatter("Bucket %d var %f, load %f, imbalance %f",j, var, bload,p.imbalance[i]);

				moved = true;
				_table[j] = load[p.transfer[i]].cpu_phys_id;

				p.imbalance[i] += bload;
				//p.imbalance[p.transfer[i]] -= bload;
				n++;
				/*if (n == cn) {
					click_chatter("All buckets moved... This should not happen");
					assert(false);
				}*/
			}
            if (n == 0 && nSkipped > 0) {
                fT[i] = 1;
                goto reagain;
                //We did not move any bucket, but we skipped some that were in without the correction
                //we re-do a pass and move those ones.
                //TODO : we should do this on the highly overloaded first, and re-compute the correction each time we force a move
            }
			if (unlikely(balancer->_verbose))
				click_chatter("Moving %d/%d/%d buckets (%f%% of load, %f%% of buckets). Core has seen %llu packets.", n, qsz, cn,  tot_bload *100, tot_bpc*100, core_tot);
		} //For each cores

#endif
        total_imbalance = 0;
        for (unsigned i = 0; i < p.N; i++) {
            total_imbalance += abs(p.imbalance[i]);
			if (abs(p.imbalance[i]) > _threshold * 2 ) {
				if (unlikely(balancer->_verbose))
					click_chatter("Imbalance of core %d left to %f, that's quite a MISS.", i, p.imbalance[i]);
			}
		}


		if (moved) {
			Timestamp t = Timestamp::now_steady();
			auto v = (t-begin).usecval();
			if (unlikely(balancer->_verbose || v > 100))
				click_chatter("Solution computed in %d usec", v);
			update_reta();
		}
    }

	if (!_counter_is_xdp)
		((AggregateCounterVector*)_counter)->advance_epoch();
	else {
	click_chatter("Reseting XDP table");
	    int cpus = balancer->_max_cpus;
	    unsigned int nr_cpus = bpf_num_possible_cpus();
	    uint64_t values[nr_cpus] = {0};
	    for (uint32_t key = 0; key < _count.size(); key++) {
	        if (bpf_map_update_elem(_xdp_table_fd, &key, values, BPF_ANY)) {
				click_chatter("XDP set failed");
				continue;
	        }
	    }
	}

    //Change the speed of ticks if necessary
    if (total_imbalance > 0.2) {
        if (total_imbalance > 0.4)
            balancer->_current_tick = balancer->_tick;
        else
            balancer->_current_tick /= 2;
    } else if (!moved)
        balancer->_current_tick *= 2;

    if (balancer->_current_tick > balancer->_tick_max)
        balancer->_current_tick = balancer->_tick_max;
    if (balancer->_current_tick < balancer->_tick)
        balancer->_current_tick = balancer->_tick;

}

bool BalanceMethodRSS::update_reta_flow(bool validate) {
		int port_id = ((DPDKDevice*)_fd)->port_id;
		struct rte_flow_attr attr;
		    Vector<rte_flow*> newflows;
             int tot;
             if (_flows.size() == 1)
                 tot = 2;
             else tot = 1;
             for (int i = 0; i < tot; i++) {
                memset(&attr, 0, sizeof(struct rte_flow_attr));
                attr.ingress = 1;

                struct rte_flow_action action[2];
                struct rte_flow_action_mark mark;
                struct rte_flow_action_rss rss;

                memset(action, 0, sizeof(action));
                memset(&rss, 0, sizeof(rss));
    /*
		        action[0].type = RTE_FLOW_ACTION_TYPE_MARK;
		        mark.id = _matches.size();
		        action[0].conf = &mark;
*/
int aid = 0;
		            action[aid].type = RTE_FLOW_ACTION_TYPE_RSS;
		            assert(_table.size() > 0);
		            uint16_t queue[_table.size()];
		            for (int i = 0; i < _table.size(); i++) {
		                queue[i] = _table[i];
		            }
		            rss.types = _rss_conf.rss_hf;
		            rss.key_len = _rss_conf.rss_key_len;
		            rss.queue_num = _table.size();
		            rss.key = _rss_conf.rss_key;
		            rss.queue = queue;
		            rss.level = 0;
		            rss.func = RTE_ETH_HASH_FUNCTION_DEFAULT;
		            action[aid].conf = &rss;
++aid;
		        action[aid].type = RTE_FLOW_ACTION_TYPE_END;
++aid;

		        Vector<rte_flow_item> pattern;
		        //Ethernet
		        /*
		        struct rte_flow_item_eth* eth = (struct rte_flow_item_eth*) malloc(sizeof(rte_flow_item_eth));
		        struct rte_flow_item_eth* mask = (struct rte_flow_item_eth*) malloc(sizeof(rte_flow_item_eth));
		        bzero(eth, sizeof(rte_flow_item_eth));
		        bzero(mask, sizeof(rte_flow_item_eth));*/
		        rte_flow_item pat;
		        pat.type = RTE_FLOW_ITEM_TYPE_ETH;
			   pat.spec = 0;
			   pat.mask = 0;
			   pat.last = 0;
			  pattern.push_back(pat);

		       pat.type = RTE_FLOW_ITEM_TYPE_IPV4;

		       if (tot == 2) {
				   struct rte_flow_item_ipv4* spec = (struct rte_flow_item_ipv4*) malloc(sizeof(rte_flow_item_ipv4));
				   struct rte_flow_item_ipv4* mask = (struct rte_flow_item_ipv4*) malloc(sizeof(rte_flow_item_ipv4));
				   bzero(spec, sizeof(rte_flow_item_ipv4));
				   bzero(mask, sizeof(rte_flow_item_ipv4));
				   spec->hdr.dst_addr = i;
				   mask->hdr.dst_addr = 1;
			   pat.spec = spec;
			  pat.mask = mask;
		       } else {
		                           pat.spec = 0;
		                           pat.mask = 0;
		       }

		                           pat.last = 0;
		                           pattern.push_back(pat);

		    rte_flow_item end;
		    memset(&end, 0, sizeof(struct rte_flow_item));
		    end.type =  RTE_FLOW_ITEM_TYPE_END;
		    pattern.push_back(end);

		    struct rte_flow_error error;
		    int res = 0;
		    if (validate)
			res = rte_flow_validate(port_id, &attr, pattern.data(), action, &error);
		    if (!res) {

		        struct rte_flow *flow = rte_flow_create(port_id, &attr, pattern.data(), action, &error);
		        if (flow) {
				if (unlikely(balancer->_verbose))
					click_chatter("Flow added succesfully with %d patterns!", pattern.size());
		        } else {
				if (unlikely(balancer->_verbose))
					click_chatter("Could not add pattern with %d patterns, error %d : %s", pattern.size(),  res, error.message);
                    return false;
		        }

		        newflows.push_back(flow);
		    } else {
			if (unlikely(balancer->_verbose))
				click_chatter("Could not validate pattern with %d patterns, error %d : %s", pattern.size(),  res, error.message);
                return false;
		    }
		 }
	        while (!_flows.empty()) {
			struct rte_flow_error error;
			rte_flow_destroy(port_id,_flows.back(), &error);
			_flows.pop_back();
	        }
		 _flows = newflows;
         return true;

}
bool BalanceMethodRSS::update_reta(bool validate) {
    Timestamp t = Timestamp::now_steady();
	if (_verifier) {
		_verifier->_table = _table;
	}

	if (_update_reta_flow) {
        if (!update_reta_flow(validate))
            return false;
    } else {
        if (!_fd->set_rss_reta(_fd, _table))
            return false;
    }

	Timestamp s = Timestamp::now_steady();
    if (validate || balancer->_verbose)
	click_chatter("Reta updated in %d usec",(s-t).usecval());
    return true;
}

int MethodPianoRSS::configure(Vector<String> &conf, ErrorHandler *errh)  {
	Element* e = 0;
	double t;
	double i;
    double threshold;
	if (Args(balancer, errh).bind(conf)
            .read_or_set("TARGET_LOAD", t, 0.8)
			.read("RSSCOUNTER", e)
			.read_or_set("IMBALANCE_ALPHA", i, 1)
            .read_or_set("IMBALANCE_THRESHOLD", threshold, 0.02) //Do not scale core underloaded or overloaded by this threshold
			.consume() < 0)
		return -1;
	_target_load = t;
	_imbalance_alpha = i;
    _threshold = threshold;

	int err = BalanceMethodRSS::configure(conf, errh);
	if (err != 0)
		return err;

	//Reta_size must be set before this
	if (e) {
		_counter = (AggregateCounterVector*)e->cast("AggregateCounterVector");
		if (!_counter) {
			_counter = (XDPLoader*)e->cast("XDPLoader");
			if (!_counter) {
				return errh->error("Counter must be of the type AggregateCounterVector or XDPLoader");
			}
			_counter_is_xdp = true;
		} else {
            _counter_is_xdp = false;
        }
	} else {
        return errh->error("You must set a RSSCOUNTER element");
    }
	return 0;
}



DeviceBalancer::DeviceBalancer() : _timer(this), _verbose(false) {
}

DeviceBalancer::~DeviceBalancer() {
}


int DeviceBalancer::addCore() {
	if (_available_cpus.size() < 1)
		return -1;
	int a_id = _available_cpus.back();
	_available_cpus.pop_back();
	_used_cpus.push_back(make_info(a_id));
	return a_id;
}

void DeviceBalancer::removeCore(int remove_phys_id) {
	_available_cpus.push_back(remove_phys_id);
	for (int uidx = 0; uidx < _used_cpus.size(); uidx++) {
		if (_used_cpus[uidx].id == remove_phys_id) {
				_used_cpus[uidx] = _used_cpus.back();
				_used_cpus.pop_back();
				return;
		}
	}
	assert(false);
}

int
DeviceBalancer::configure(Vector<String> &conf, ErrorHandler *errh) {
    Element* dev = 0;
    String config;
    String method;
    String target;
    String load;
    String source;
    String cycles;
    int startcpu;
    bool havemax;
    Element* manager = 0;
    if (Args(this, errh).bind(conf)
        .read_mp("METHOD", method)
        .read_mp("DEV", dev)
        .read_or_set("CONFIG", config, "")
        .read_or_set("CORE_OFFSET", _core_offset, 0)
        .read_or_set("TIMER", _tick, 100)
        .read("TIMER_MAX", _tick_max).read_status(havemax)
        .read_or_set("CPUS", _max_cpus, click_max_cpu_ids())
        .read_or_set("TARGET", target, "load")
        .read_or_set("STARTCPU", startcpu, -1)
        .read_or_set("UNDERLOAD", _underloaded_thresh, 0.25)
        .read_or_set("OVERLOAD", _overloaded_thresh, 0.75)
        .read_or_set("LOAD", load, "cpu")
		.read_or_set("CYCLES", cycles, "cycles")
		.read_or_set("AUTOSCALE", _autoscale, false)
		.read_or_set("ACTIVE", _active, true)
		.read_or_set("VERBOSE", _verbose, true)
		.read("MANAGER", manager)
        .consume() < 0)
        return -1;

    if (!havemax)
        _tick_max = _tick;

    if (cycles == "cycles") {
	_load = LOAD_CYCLES;
    } else  if (cycles == "cyclesqueue") {
	_load = LOAD_CYCLES_THEN_QUEUE;
    } else  if (cycles == "cpu") {
	_load = LOAD_CPU;
    } else  if (cycles == "queue") {
		_load = LOAD_QUEUE;
    } else  if (cycles == "realcpu") {
	_cpustats.resize(_max_cpus);
		_load = LOAD_REALCPU;
	} else {
		return errh->error("Unknown cycle method !");
	}

    if (startcpu == -1) {
        startcpu = _max_cpus;
    }

    _startwith = startcpu;

    if (method == "metron") {
        _method = new MethodMetron(this, dev, config);
    } else if (method == "pianorss") {
        _method = new MethodPianoRSS(this, dev, config);
    } else if (method == "rss") {
        _method = new BalanceMethodRSS(this, dev);
    } else if (method == "rssrr") {
        _method = new MethodRSSRR(this, dev);
    } else {
        return errh->error("Unknown method %s", method.c_str());
    }

    if (manager) {
	_manager = dynamic_cast<FlowIPManager*>(manager);
    }

    if (target=="load")
        _target = TARGET_LOAD;
    else if (target == "balance")
        _target = TARGET_BALANCE;
    else
        return errh->error("Unknown target %s", target.c_str());

    if (_method->configure(conf,errh) !=0)
	return -1;

    if (Args(this, errh).bind(conf)
            .complete() < 0)
        return -1;


    return 0;
}


int
DeviceBalancer::initialize(ErrorHandler *errh) {
    int startwith = _startwith;
    if (_method->initialize(errh, startwith) != 0)
        return -1;
    for (int i = 0; i < startwith; i++) {
       _used_cpus.push_back(CpuInfo{.id= i,.last_cycles=0});
    }
    for (int i = _max_cpus - 1; i >=startwith; i--) {
       _available_cpus.push_back(i);
    }
    _timer.initialize(this);
    _current_tick = _tick;
    if (_active)
        _timer.schedule_after_msec(_current_tick);
    return 0;
}



void
DeviceBalancer::run_timer(Timer* t) {
    Vector<Pair<int,float>> load;
    float totload = 0;
    if (_load == LOAD_CPU) {
        for (int u = 0; u < _used_cpus.size(); u++) {
            int i = _used_cpus[u].id;
            float l = master()->thread(i)->load();
            load.push_back(Pair<int,float>{i,l});
            totload += l;
        }
    } else if (_load == LOAD_CYCLES || _load == LOAD_CYCLES_THEN_QUEUE) {
		/**
		 * Use the amount of cycles since last tick, more precise than LOAD_CPU.
		 * We use the raw amount of cycles, divided by the total amount of cycles for all CPUs
		 * This will give a number between 0 and 1, 1 being the total for all CPUs
		 * We therefore multiply the load by the total (unprecise) load, to give a realistic
		 * scale in term of "amount of cores" load but giving a relative better precision.
		 */
		unsigned long long utotload = 0;
		Vector<unsigned long long> uload;
		for (int u = 0; u < _used_cpus.size(); u++) {
			int phys_id = _used_cpus[u].id;
			unsigned long long ul = master()->thread(phys_id)->useful_kcycles();
			unsigned long long pl = _used_cpus[u].last_cycles;
			unsigned long long dl = ul - pl;
			_used_cpus[u].last_cycles = ul;
			uload.push_back(dl);
			utotload += dl;
			//click_chatter("core %d kcycles %llu %llu load %f", i, ul, dl, master()->thread(i)->load());
			totload += master()->thread(phys_id)->load();
		}
		if (utotload > 0) {
				for (int u = 0; u < _used_cpus.size(); u++) {
					double pc = (double)uload[u] / (double)utotload;
					if (totload * pc > 1.0)
						totload = 1.0 / pc;
				}
		}

		int overloaded = 0;

		for (int u = 0; u < _used_cpus.size(); u++) {
			int i = _used_cpus[u].id;
			float l;
			if (utotload == 0)
				l = 0;
			else {
				double pc = (double)uload[u] / (double)utotload;
				l = pc * totload;

				//click_chatter("core %d load %f -> %f", i, pc, l);
			}
			assert(l <= 1 && l>=0);
			if (l > 0.98)
				overloaded ++;
			load.push_back(Pair<int,float>{i,l});
		}

		if (_load == LOAD_CYCLES_THEN_QUEUE && overloaded > 1) {
	        DPDKDevice* fd = (DPDKDevice*)((BalanceMethodDevice*)_method)->_fd;
	        int port_id = fd->port_id;
	        float rxdesc = fd->get_nb_rxdesc();
	        for (int u = 0; u < _used_cpus.size(); u++) {
	            int i = _used_cpus[u].id;
	            int v = rte_eth_rx_queue_count(port_id, i);
	            if (v < 0) {
			click_chatter("WARNING : unsupported rte_eth_rx_queue_count for queue %d, error %d", i, v);
			continue;
	            }
	            float l = (float)v / rxdesc;
	            //click_chatter("Core %d %f %f",u,load[u].second,l);

	            load[u].second = load[u].second * 0.90 + 0.10 * l;
	        }
		}
    } else if (_load == LOAD_REALCPU) {
        Vector<float> l(_max_cpus, 0);
        unsigned long long totalUser, totalUserLow, totalSys, totalIdle, totalIoWait, totalIrq, totalSoftIrq;
        int cpuId;
        FILE* file = fopen("/proc/stat", "r");
        char buffer[1024];
        char *res = fgets(buffer, 1024, file);
        assert(res);
        while (fscanf(file, "cpu%d %llu %llu %llu %llu %llu %llu %llu", &cpuId, &totalUser, &totalUserLow, &totalSys, &totalIdle, &totalIoWait, &totalIrq, &totalSoftIrq) > 0) {
		if (cpuId < l.size()) {
                unsigned long long newTotal = totalUser + totalUserLow + totalSys + totalIrq + totalSoftIrq;
                unsigned long long tdiff =  (newTotal - _cpustats[cpuId].lastTotal);
                unsigned long long idiff =  (totalIdle - _cpustats[cpuId].lastIdle);
                if (tdiff + idiff > 0)
                    l[cpuId] =  (float)tdiff / (float)(tdiff + idiff);
                _cpustats[cpuId].lastTotal = newTotal;
                _cpustats[cpuId].lastIdle = totalIdle;
                //click_chatter("C %d total %d %d %d",cpuId,newTotal,tdiff, idiff);
                res = fgets(buffer, 1024, file);
            }
        }
        fclose(file);
        for (int u = 0; u < _used_cpus.size(); u++) {
		//click_chatter("Used %d load %f",u,l[u]);
		int i = _used_cpus[u].id;
			float cl = l[i];
			load.push_back(Pair<int,float>{i,cl});
			totload += cl;
        }
    } else { //_load == LOAD_QUEUE
        DPDKDevice* fd = (DPDKDevice*)((BalanceMethodDevice*)_method)->_fd;
        int port_id = fd->port_id;
        float rxdesc = fd->get_nb_rxdesc();
        for (int u = 0; u < _used_cpus.size(); u++) {
            int i = _used_cpus[u].id;
            int v = rte_eth_rx_queue_count(port_id, i);
            float l = (float)v / rxdesc;
            load.push_back(Pair<int,float>{i,l});
            totload += l;
        }
    }

    if (unlikely(_verbose > 1)) {
		String s = "load ";
		for (int u = 0; u < load.size(); u++) {
			s += String(load[u].second) + " ";
		}
		s += "\n";
		click_chatter("%s",s.c_str());
    }

    if (unlikely(_target == TARGET_BALANCE)) {
        float target = totload / _max_cpus;
        for (int i = 0; i < load.size(); i ++) {
            if (target < 0.1)
                load[i].second = 0.5;
            else {
            // 0  0.2  -> 0
            // 0.1 0.2 -> 0.25
            // 0.2 0.2 -> 0.5
            // 0.3 0.2 -> 0.75
            // 0.4 0.2 -> 1
//                click_chatter("Load %f target %f -> ",load[i].second, target);
                load[i].second =  load[i].second / target / 2;
            }

            if (load[i].second > 1)
                load[i].second = 1;
            if (load[i].second < 0)
                load[i].second = 0;
            //click_chatter("%f",load[i].second);
        }

    }

 /*   assert(_method);
    assert(load.size() > 0);
    for (int i = 0; i < load.size(); i++) {
	click_chatter("Load of core %d is %f",load[i].first,load[i].second);
    }*/

    _method->rebalance(load);
    if (_active)
        _timer.schedule_after_msec(_current_tick);
}

DeviceBalancer::CpuInfo
DeviceBalancer::make_info(int _id) {
	CpuInfo i;
	i.id = _id;
	i.last_cycles = master()->thread(_id)->useful_kcycles();
	return i;
}


enum {h_active, h_autoscale};


String
DeviceBalancer::read_param(Element *e, void *thunk)
{
	DeviceBalancer *td = (DeviceBalancer *)e;
    switch((uintptr_t) thunk) {
    case h_active:
        return String(td->_active);
    case h_autoscale:
        return String(td->_autoscale);
    default:
        return String();
    }
}

int
DeviceBalancer::write_param(const String &in_s, Element *e, void *vparam,
                     ErrorHandler *errh)
{
	DeviceBalancer *td = (DeviceBalancer *)e;
    String s = cp_uncomment(in_s);
    switch ((intptr_t)vparam) {
    case h_active: {
        bool active;
        if (!BoolArg().parse(s, active))
            return errh->error("type mismatch");
        if (active && !td->_active)
		    td->_timer.schedule_after_msec(td->_current_tick);

        td->_active = active;
        break;
    }
    case h_autoscale: {
        bool active;
        if (!BoolArg().parse(s, active))
            return errh->error("type mismatch");
        td->_autoscale = active;
        break;
    }
    }
    return 0;
}

void
DeviceBalancer::add_handlers()
{
    add_read_handler("active", read_param, h_active, Handler::CHECKBOX);
    add_read_handler("autoscale", read_param, h_autoscale, Handler::CHECKBOX);
    add_write_handler("active", write_param, h_active);
    add_write_handler("autoscale", write_param, h_autoscale);
}


CLICK_ENDDECLS

EXPORT_ELEMENT(DeviceBalancer)
