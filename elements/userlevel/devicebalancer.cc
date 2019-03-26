// -*- c-basic-offset: 4; related-file-name: "fromdpdkdevice.hh" -*-
/*
 * fromdpdkdevice.{cc,hh} -- element reads packets live from network via
 * the DPDK. Configures DPDK-based NICs via DPDK's Flow API.
 *
 * Copyright (c) 2014-2015 Cyril Soldani, University of Liège
 * Copyright (c) 2016-2017 Tom Barbette, University of Liège
 * Copyright (c) 2017 Georgios Katsikas, RISE SICS
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
#include <rte_flow.h>
#include "devicebalancer.hh"
#include <queue>
#include <vector>

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
int MethodMetron::initialize(ErrorHandler *errh, int startwith) {
    FlowDirector *flow_dir = FlowDirector::get_flow_director(_fd->get_device()->port_id);
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

            for (uint32_t i = 0; i < rules_vec.size(); i++) {
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

    load_tracker_initialize(errh);

    auto cache = flow_dir->get_flow_cache();
    if (!cache) {
        return errh->error("Flow dir is not initialized");
    }

    if (!cache->has_rules())
        return errh->error("Cache error !");

    return 0;
}

void MethodMetron::rebalance(Vector<Pair<int,float>> load) {
    FlowDirector *flow_dir = FlowDirector::get_flow_director(_fd->get_device()->port_id);
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
        if (now - _last_movement[cpuid] < CLICK_HZ * 3)
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
        int a;
        float load_a;
        if (balancer->_available_cpus.size() > 0) {
            a = balancer->_available_cpus.back();
            load_a = 0;
            balancer->_available_cpus.pop_back();
            balancer->_used_cpus.push_back(balancer->make_info(a));
        } else if (underloaded.size() > 0) {
            a = underloaded.back().first;
            load_a = underloaded.back().second;

            underloaded.pop_back();
        } else {
            click_chatter("Not enough CPUs available");
            break;
        }

        _last_movement[o] = now;
        _last_movement[a] = now;

        click_chatter("Deflating cpu core %d to %d", o, a);

        auto cache = flow_dir->get_flow_cache();
        if (!cache) {
            click_chatter("Flow dir is not initialized");
            assert(false);
        }


        HashMap<long, String> * rmap;
        rmap = cache->rules_map_by_core_id(o);
        int mig = rmap->size() * ((load_o - load_a) / 2.0);
        //click_chatter("Rmap %p", rmap);
        click_chatter("Migrating %d flows / %d flows", mig , rmap->size());

        Vector<uint32_t> dlist;
        HashMap<long, String> nmap;
//      flow_dir->flow_cache()->delete_rule_by_global_id(rmap.keys());
//      flow_dir->flow_cache()->
        int n = 0;
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
            String nrule = rewrite_id(rule, a);
            nmap.insert(rule_id,nrule);
        //    click_chatter("New rule %s",nrule.c_str());
            it++;
        }

        // Delete the flow rules
  //      int status = flow_dir->flow_rules_delete(dlist.data(), dlist.size());

//        click_chatter("Deleted %d/%d", status, dlist.size(), true);
        int status = flow_dir->update_rules(nmap, true, a);
        click_chatter("Update %d", status);

    }

    while (underloaded.size() >= 2 && balancer->_target != TARGET_BALANCE) {
        int remove_i = underloaded[underloaded.size() - 2].first;
        int with_i = underloaded[underloaded.size() - 1].first;
        click_chatter("Removing underloaded CPU %d!", remove_i);
        float load_r =  underloaded[underloaded.size() - 2].second;
        float load_i =  underloaded[underloaded.size() - 1].second;

        _last_movement[remove_i] = now;

        for (int uidx = 0; uidx < balancer->_used_cpus.size(); uidx++) {
            if (balancer->_used_cpus[uidx].id == remove_i) {
                balancer->_used_cpus[uidx] = balancer->_used_cpus.back();
                balancer->_used_cpus.pop_back();
                break;
            }
        }
        _last_movement[with_i] = now;
        underloaded.pop_back();
        underloaded.pop_back();
//        click_chatter("Core %d and %d", remove_i, with_i);
        if (load_r + load_i >= balancer->_overloaded_thresh)
            continue;
        balancer->_available_cpus.push_back(remove_i);
//        click_chatter("Core %d is now available", remove_i);

        HashMap<long, String> * rmap;

        HashMap<long, String> nmap;
        auto cache = flow_dir->get_flow_cache();

        rmap = cache->rules_map_by_core_id(remove_i);
//        click_chatter("Rmap %p", rmap);
        click_chatter("Migrating %d flows", rmap->size());
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

//        click_chatter("Update %d", status);
    }
}


int MethodPianoRSS::initialize(ErrorHandler *errh, int startwith) {

    // Invoke Flow Director only if active
/*    if (flow_dir->active()) {

    } else {
        return errh->error("Flow director is not active.");
    }
*/
	int reta = _fd->get_device()->get_reta_size();

	if (reta <= 0)
		return errh->error("DPDK device not initialized or RSS is misconfigured");
	_table = _fd->get_device()->get_rss_reta();
	for (int i = 0; i < startwith; i++) {
		_table[i] = i % startwith;
	}
    load_tracker_initialize(errh);
    click_chatter("PIANO initialized");
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

	Problem() : oid(), uid() {

	}

	bool tryK(int i)
	{
		if (i == oid.size()) {

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


void MethodPianoRSS::rebalance(Vector<Pair<int,float>> rload) {
    click_jiffies_t now = click_jiffies();

/*
    uint64_t total_packets = 0;
    for (int i = 0; i < _table.size(); i++) {
	total_packets += _counter->find_node(i)->count;
    }
*/
    Timestamp begin = Timestamp::now_steady();
    float _threshold = 0.05; //Do not scale core underloaded or overloaded by this threshold
    float _imbalance_alpha = 2.0f/3.0f; //Percentage of imbalance to consider
    float _min_load = 0.15;
    float _threshold_force_overload = 0.90;
    float _load_alpha = 0.2;
    //Vector<float> corrections;
    //corrections.resize(click_max_cpu_ids(), 0);
    float target_load = _target_load;

    //float target_load = balancer->_used_cpus.size();
    float suppload = 0;
    Problem p;
    float totalload = 0;
    Vector<Pair<int,float>> load(rload.size(),Pair<int,float>(0,0));
    for (int j = 0; j < rload.size(); j++) {
        int cpuid = rload[j].first;
        float load_current = rload[j].second;
        load[j].first = cpuid;
        if (_past_load[cpuid] == 0)
		load[j].second = load_current;
        else
		load[j].second = load_current * _load_alpha + _past_load[cpuid] * (1-_load_alpha);
        _past_load[cpuid] = load[j].second;
        float diff = target_load - load_current; // >0 if underloaded diff->quantity to add
        if (abs(diff) <= _threshold)
		diff = 0;
        //corrections[cpuid] = diff;
        suppload += diff;
        totalload+=load_current;
    }

/*

    if (suppload > 1.1) { // we can remove a core
	if (unlikely(balancer->_verbose))
		click_chatter("Removing a core");
    } else if (suppload < EPSILON) { //We need a new core
	if (unlikely(balancer->_verbose))
		click_chatter("Adding a core");
    }*/

    const int N = load.size();
    p.N = N;
    p.target = totalload / (float)N;

    if (unlikely(balancer->_verbose))
	click_chatter("Target %f. Total load %f. %d cores", p.target, totalload, N);

    p.imbalance.resize(N);
	if (p.target <  _min_load && !_threshold_force_overload) {
		click_chatter("Underloaded, skipping balancing");
	}

    for (int i = 0; i < N; i++) {

	p.imbalance[i] = 0 * (1.0-_imbalance_alpha) + ( (p.target - load[i].second) * _imbalance_alpha); //(_last_imbalance[load[i].first] / 2) + ((p.target - load[i].second) / 2.0f);

	if (p.imbalance[i] > _threshold)
		p.uid.push_back(i);
	else if (p.imbalance[i] < - _threshold) {
		p.oid.push_back(i);
	}
    }

    if (p.oid.size() > 0) {
		float min_cost = N * N;
		p.transfer.resize(N);
		for (int i = 0; i < N; i++) {
			p.transfer[i] = i;//Default is to transfer load to itself
		}
		p.min_cost = N;
		p.tryK(0);
		if (unlikely(balancer->_verbose))
			click_chatter("Transfer solution for %d uid, %d oid:",p.uid.size(), p.oid.size());
		for (int i = 0; i < N; i++) {
			if (unlikely(balancer->_verbose))
				click_chatter("Core %d (load %f, imbalance %f) -> %d (load %f)", i, load[i].second, p.imbalance[i], p.transfer[i], load[p.transfer[i]].second);
			if (p.transfer[i] == i)
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
			Vector<Pair<int,float> > q;
			int j = click_random() % _table.size();
			for (int r = 0; r < _table.size(); r++) {
				if (++j == _table.size())
					j = 0;
				if (_table[j] == i) {

					double c = _counter->find_node(j).count;
					double p = _counter->find_node(j).variance;
					core_tot += c;
					//double mean = (c + p) / 2;
					double var = min(c,p) / max(c,p);
					core_tot_min += c * var;
					//if (var < 1) SEE BELOW
						q.push_back(Pair<int,float>(j, var));
					cn++;
				}
			}
			int n = 0;
			float tot_bload = 0;
			float tot_bpc = 0;
			int qsz = q.size();
			/**
			 * We cannot leave the core with only a high variance bucket
			 * We want :
			 * - The splits to have an equal amount of variance -> we just take bucket in order, it's like random
			 * - Not overflow the imbalance
			 */
			while (!q.empty()) {
				int j = q.back().first;
				float var = q.back().second;
				q.pop_back();

				double c = (double)_counter->find_node(j).count;

				double bpc = ((c * var) / core_tot_min);
				double bload =  bpc * load[i].second; //How much CPU load this bucket represents
				if (bload > -p.imbalance[i] || p.imbalance[i] > 0) { //Never over balance!
					//if (bload > _threshold)
					if (unlikely(balancer->_verbose))
						click_chatter("Skipping bucket %d var %f, load %f%% (min %f%% of bucket, %f%% real), imbalance %f, at least %d flows",j, var, bload*100, bpc*100, ((c * 100.0)/core_tot),p.imbalance[i],_counter->find_node(j).flows);
/*
						left_max += bload * (1 + var);
						left += bload;
						left_min += bload / (1 + var);*/

					continue;
				}


				tot_bload += bload;
				tot_bpc += bpc;
				//click_chatter("Bucket %d var %f, load %f, imbalance %f",j, var, bload,p.imbalance[i]);

				_table[j] = p.transfer[i];

				p.imbalance[i] += bload;
				p.imbalance[p.transfer[i]] -= bload;
				n++;
				/*if (n == cn) {
					click_chatter("All buckets moved... This should not happen");
					assert(false);
				}*/
			}
			if (unlikely(balancer->_verbose))
			click_chatter("Moving %d/%d/%d buckets (%f%% of load, %f%%). Core has seen %llu packets.", n, qsz, cn, tot_bload *100, tot_bpc*100, core_tot);
		}
		for (int i =0; i < N;i ++) {
			if (abs(p.imbalance[i]) > _threshold) {
				if (unlikely(balancer->_verbose))
					click_chatter("Imbalance of core %d left to %f, that's quite a MISS.", i, p.imbalance[i]);
			}
		}
		Timestamp t = Timestamp::now_steady();
		click_chatter("Solution computed in %d usec", (t-begin).usecval());
		//_fd->get_device()->set_rss_reta(_table);
		int port_id = _fd->get_device()->port_id;
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
		            uint16_t queue[_table.size()];
		            for (int i = 0; i < _table.size(); i++) {
		                queue[i] = _table[i];
		            }
		            uint8_t rss_key[40];
		            struct rte_eth_rss_conf rss_conf;
		            rss_conf.rss_key = rss_key;
		            rss_conf.rss_key_len = 40;
		            rte_eth_dev_rss_hash_conf_get(port_id, &rss_conf);
		            rss.types = rss_conf.rss_hf;
		            rss.key_len = rss_conf.rss_key_len;
		            rss.queue_num = _table.size();
		            rss.key = rss_key;
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
		    int res;
		    res = rte_flow_validate(port_id, &attr, pattern.data(), action, &error);
		    if (!res) {

		        struct rte_flow *flow = rte_flow_create(port_id, &attr, pattern.data(), action, &error);
		        if (flow) {
				if (unlikely(balancer->_verbose))
					click_chatter("Flow added succesfully with %d patterns!", pattern.size());
		        } else {
				if (unlikely(balancer->_verbose))
					click_chatter("Could not add pattern with %d patterns, error %d : %s", pattern.size(),  res, error.message);
		        }

		        newflows.push_back(flow);
		    } else {
			if (unlikely(balancer->_verbose))
				click_chatter("Could not validate pattern with %d patterns, error %d : %s", pattern.size(),  res, error.message);
		    }
		 }
	        while (!_flows.empty()) {
			struct rte_flow_error error;
			rte_flow_destroy(port_id,_flows.back(), &error);
			_flows.pop_back();
	        }
		 _flows = newflows;
		Timestamp s = Timestamp::now_steady();
		click_chatter("Reta updated in %d usec",(s-t).usecval());
    }
    assert(_counter);
    _counter->advance_epoch();
}

int MethodPianoRSS::configure(Vector<String> &conf, ErrorHandler *errh)  {
	Element* e;
	double t;
	if (Args(balancer, errh).bind(conf)
			.read("LOAD", t)
			.read("RSSCOUNTER", e)
			.consume() < 0)
		return -1;
	_target_load = t;
	_counter = (AggregateCounterVector*)e->cast("AggregateCounterVector");
	if (!_counter)
		return errh->error("Counter must be of the type AggregateCounterVector");
	return 0;
}



DeviceBalancer::DeviceBalancer() : _timer(this), _verbose(false) {
}

DeviceBalancer::~DeviceBalancer() {
}

int
DeviceBalancer::configure(Vector<String> &conf, ErrorHandler *errh) {
    Element* dev;
    String config;
    String method;
    String target;
    String load;
    String source;
    int startcpu;
    if (Args(this, errh).bind(conf)
        .read_mp("METHOD", method)
        .read_mp("DEV", dev)
        .read_or_set("CONFIG", config, "")
        .read_or_set("CORE_OFFSET", _core_offset, 0)
        .read_or_set("TIMER", _tick, 500)
        .read_or_set("CPUS", _max_cpus, click_max_cpu_ids())
        .read_or_set("TARGET", target, "load")
        .read_or_set("STARTCPU", startcpu, -1)
        .read_or_set("UNDERLOAD", _underloaded_thresh, 0.25)
        .read_or_set("OVERLOAD", _overloaded_thresh, 0.75)
        .read_or_set("LOAD", load, "cpu")
		.read_or_set("VERBOSE", _verbose, true)
        .consume() < 0)
        return -1;


    _load = LOAD_CYCLES;
    if (startcpu == -1) {
        startcpu = _max_cpus;
    }

    _startwith = startcpu;

    if (method == "metron") {
        _method = new MethodMetron(this, dev, config);
    } else if (method == "pianorss") {
        _method = new MethodPianoRSS(this, dev, config);
    } else {
        return errh->error("Unknown method %s", method.c_str());
    }

    if (target=="load")
        _target = TARGET_LOAD;
    else if (target == "balance")
        _target = TARGET_BALANCE;
    else
        return errh->error("Unknown target %s", target.c_str());

    if (_method->configure(conf,errh) !=0)
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
    _timer.schedule_after_msec(_tick);
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
    } else if (_load == LOAD_CYCLES) {
	/**
	 * Use the amount of cycles since last tick, more precise than LOAD_CPU.
	 * We use the raw amount of cycles, divided by the total amount of cycles for all CPUs
	 * This will give a number between 0 and 1, 1 being the total for all CPUs
	 * We therefore multiply the load by the total (unprecise) load, to give a realistic
	 * scale in term of "amount of cores" load but giving a relative better precision
	 */
	unsigned long long utotload = 0;
	Vector<unsigned long long> uload;
	for (int u = 0; u < _used_cpus.size(); u++) {
		int i = _used_cpus[u].id;
			unsigned long long ul = master()->thread(i)->useful_kcycles();
			unsigned long long pl = _used_cpus[u].last_cycles;
			unsigned long long dl = ul - pl;
			_used_cpus[i].last_cycles = ul;
			uload.push_back(dl);
			utotload += dl;
            click_chatter("core %d kcycles %llu %llu load %f", i, ul, dl, master()->thread(i)->load());
			totload += master()->thread(i)->load();
		}
	for (int u = 0; u < _used_cpus.size(); u++) {
		int i = _used_cpus[u].id;
		float l;
		if (utotload == 0)
			l = 0;
		else
			l = ((double)uload[u] / (double)utotload) * totload;
		load.push_back(Pair<int,float>{i,l});
	}
	//
    } else { //_load == LOAD_QUEUE
        FromDPDKDevice* fd = ((BalanceMethodDPDK*)_method)->_fd;
        int port_id = fd->get_device()->port_id;
        float rxdesc = fd->get_nb_desc();
        for (int u = 0; u < _used_cpus.size(); u++) {
            int i = _used_cpus[u].id;
            int v = rte_eth_rx_queue_count(port_id, i);
            float l = (float)v / rxdesc;
            load.push_back(Pair<int,float>{i,l});
            totload += l;
        }
    }
    if (_target == TARGET_BALANCE) {
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

    _method->rebalance(load);
    _timer.reschedule_after_msec(_tick);
}

DeviceBalancer::CpuInfo
DeviceBalancer::make_info(int _id) {
	CpuInfo i;
	i.id = _id;
	i.last_cycles = master()->thread(_id)->useful_kcycles();
	return i;
}


CLICK_ENDDECLS

EXPORT_ELEMENT(DeviceBalancer)
