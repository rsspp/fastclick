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
#include <click/multithread.hh>
#include <click/straccum.hh>

#include <click/numa.hh>
#include <rte_flow.h>
#include "devicebalancer.hh"
#ifdef HAVE_BPF
#include "xdploader.hh"
#endif
#include <click/solver.hh>
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


void BalanceMethod::cpu_changed() {
        balancer->_timer.schedule_now();
};

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
    click_chatter("Actual reta size %d, target %d", reta, _reta_size);
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

    click_chatter("RSS initialized with %d CPUs and %d buckets", startwith, _table.size());
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
    //update_reta();
}

void BalanceMethodRSS::cpu_changed() {
    int m =  balancer->_used_cpus.size();
    Vector<Vector<Pair<int,int>>> omoves(balancer->max_cpus(), Vector<Pair<int,int>>());
    /*Vector<int> epochs;
    epochs.resize(max_cpus());*/


    for (int i = 0; i < _table.size(); i++) {

    int newcpu = balancer->_used_cpus[i % m].id;
    if (balancer->_manager && newcpu!= _table[i]) {
            omoves[_table[i]].push_back(Pair<int,int>(i, newcpu));
        }
    //epochs(_table[i]) =
        _table[i] = newcpu;
    }

    if (balancer->_manager) {
        for (int i = 0; i < m; i++) {
            if (omoves[i].size() > 0) {
                balancer->_manager->pre_migrate((DPDKDevice*)_fd, i, omoves[i]);
            }
        }
    }
    click_chatter("Migration info written. Updating reta.");
    update_reta();
    click_chatter("Post migration");
    if (balancer->_manager) {
        for (int i = 0; i < m; i++) {
            if (omoves[i].size() > 0) {
                balancer->_manager->post_migrate((DPDKDevice*)_fd, i);
            }
        }
    }
    click_chatter("Post migration finished");
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

#ifdef HAVE_BPF
    if (_counter_is_xdp) {
        click_chatter("Resizing count to %d",_reta_size);
        _count.resize(_reta_size);
        _xdp_table_fd = ((XDPLoader*)_counter)->get_map_fd("count_map");
        if (!_xdp_table_fd)
            return errh->error("Could not find map !");
    }
#endif
    return 0;
}
#define EPSILON 0.0001f





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

struct MachineLoad {
    MachineLoad() : N(0), total_load(0), target(0) {

    }
    int N;
    float total_load;
    float target;
};

struct SocketLoad : MachineLoad {
    SocketLoad() : imbalance(),uid(),oid() {
    }

    Vector<float> imbalance;
    Vector<int> uid;
    Vector<int> oid;
};

#define print_cpu_assign() {for (int i = 0; i < rload.size(); i++) {click_chatter("CPU %d -> %d/%d",i, rload[i].first,load[i].cpu_phys_id);}}

void MethodPianoRSS::rebalance(Vector<Pair<int,float>> rload) {
    click_jiffies_t now = click_jiffies();
#ifdef HAVE_BPF
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
//            for (int i = 0; i < nr_cpus; i++) {
                tot += values[_table[key]];
//            if (values[i] && i != _table[key])
            //if (values[i] != 0)
//                click_chatter("BPF map ha value on core %d, but RSS is assigned to %d. Val %d",i,_table[key],values[i]);
  //              tot += values[i];
//            }
//            tot -= _count[key].count;

            _count[key].count  = tot;
            uint64_t var =  (_count[key].variance  / 3) + (2 * tot / 3);
            _count[key].variance  = var;
            if (unlikely(tot != 0 && balancer->_verbose))
                click_chatter("Key %d core %d val %d, var %f",key,_table[key], tot,(float)min(tot,var) / (float)max(tot,var));
        }
    }
#endif

#define get_node_count(i) ((_counter_is_xdp)?_count.unchecked_at(i).count:((AggregateCounterVector*)_counter)->find_node_nocheck(i).count)
#define get_node_variance(i) ((_counter_is_xdp)?_count.unchecked_at(i).variance:((AggregateCounterVector*)_counter)->find_node_nocheck(i).variance)
#define get_node_moved(i) (_count.unchecked_at(i).moved)


    Timestamp begin = Timestamp::now_steady();
    float _min_load = 0.01;
    float _threshold_force_overload = 0.90;
    float _load_alpha = 1;
    //float _high_load_threshold = 0.;
    const float _imbalance_threshold = _threshold / 2;
    assert(_imbalance_threshold <= (_threshold / 2) + EPSILON); //Or we'll always miss

    //Various flags
    float suppload = 0;
    int min_core = -1;
    float min_core_load = 1;
    float max_core_load = 0;
    bool has_high_load = false;

    //Per-core load statistic
    Vector<Load> load(rload.size(),Load());

    //Track physical CPU id to assigned ids
    Vector<unsigned> map_phys_to_id(click_max_cpu_ids(), -1);

    //Per NUMA-socket load. When do_numa is false, we consider a single NUMA node
    Vector<SocketLoad> sockets_load;
    bool do_numa = _numa;
    int numamax = do_numa?_numa_num:1;
    sockets_load.resize(numamax);

    //Keeps track of the whole machine load stats
    MachineLoad machine_load;

    //For each assigned cores, we compute the load and its smothed average
    //We check if some cores are completely overloaded, and report all that per-NUMA socket
    for (int j = 0; j < rload.size(); j++) {
        int physcpuid = rload[j].first;
        float load_current = rload[j].second;
        //assert(load_current <= 1);
        if (load_current < min_core_load)
            min_core = j;

        load[j].cpu_phys_id = physcpuid;

        if (_past_load[physcpuid] == 0)
            load[j].load = load_current;
        else
            load[j].load = load_current * _load_alpha + _past_load[physcpuid] * (1.0-_load_alpha);

        _past_load[physcpuid] = load[j].load;

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
        map_phys_to_id[physcpuid] = j;

        machine_load.total_load += load[j].load;

        if (do_numa) {
            int numaid =  Numa::get_numa_node_of_cpu(load[j].cpu_phys_id);
            sockets_load[numaid].N++;
            sockets_load[numaid].total_load += load[j].load;
        }
    }


    machine_load.N = load.size();
    machine_load.target = machine_load.total_load / (float)machine_load.N;

    for (int i = 0; i < sockets_load.size(); i++) {
        sockets_load[i].target = sockets_load[i].total_load / (float)sockets_load[i].N;
        if (do_numa && abs(sockets_load[i].target -  machine_load.target) > _threshold * 2) {
            click_chatter("Non-numa forced");
            do_numa=false;
        }
    }

    //suppload = (machine_load.N *_target_load) - totalload; //Total power that we have in excess

    float var_protect = (machine_load.N + 1) * (max(0.05f,(max_core_load-_target_load)));
    if (unlikely(balancer->_verbose))
        click_chatter("Target %f. Total load %f. %d cores. Suppload %f, Var %f", machine_load.target, machine_load.total_load, machine_load.N, suppload, var_protect);


    //Count the number of packets for each core
    for (int j = 0; j < _table.size(); j++) {
        unsigned long long c = get_node_count(j);
        unsigned cpu_phys_id =_table[j];
        if (unlikely(cpu_phys_id >= map_phys_to_id.size())) {
            click_chatter("ERROR : invalid phys id %d", cpu_phys_id);
            print_cpu_assign();
            abort();
        }
        unsigned cpuid = map_phys_to_id.unchecked_at(cpu_phys_id);
        if (unlikely(cpuid >= load.size())) {
            click_chatter("ERROR : invalid cpu id %d for physical id %d, table index %d", cpuid, cpu_phys_id, j);
            print_cpu_assign();
            abort();
        }
        auto &l = load.unchecked_at(cpuid);
        l.npackets += c;
        l.nbuckets += 1;
        if (c > 0)
            load[map_phys_to_id[_table[j]]].nbuckets_nz += 1;
    }


    /**
     * Scaling. If we need more cores, we just add it and the imbalance mapping will do the rest
     * If we remove a core, we minimize the selection of buckets of the removed core -> existing cores to
     * diminish the imbalance already
     */
    if (balancer->_autoscale) {
        if (suppload > (1 + (1 - _target_load) + var_protect)  && machine_load.N > 1) { // we can remove a core
            if (unlikely(balancer->_verbose)) {
                click_chatter("Removing a core (target %f)",_target_load);
            }

            int min_core_phys = load[min_core].cpu_phys_id;

            click_chatter("Removing core %d (phys %d)", min_core, min_core_phys);
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


            BucketMapProblem cp(buckets_indexes.size(), load.size()); //load.size() is already fixed

            //Compute load for each buckets
            for (int i = 0; i < buckets_indexes.size(); i++) {
                int j = buckets_indexes[i];
                double c = get_node_count(j);
                cp.buckets_load[i] = (((double) c / (double)totcount)*min_core_load);
            }

            //Add imbalance of all cores
            Vector<float> imbalance(machine_load.N-1,0);
            machine_load.N = machine_load.N-1;
            machine_load.target = machine_load.total_load / (float)machine_load.N;

            //Fix imbalance without the removed core
            for (int i = 0; i < machine_load.N; i++) {
                //Imbalance is positive if the core should loose some load
                imbalance[i] = 0 * (1.0-_imbalance_alpha) + ( (machine_load.target - load[i].load) * _imbalance_alpha); //(_last_imbalance[load[i].first] / 2) + ((p.target - load[i].second) / 2.0f);
            }

            click_chatter("Solving problem...");

            //Solve problem ; map all buckets of removed core to other cores
            cp.solve();

            for (int i = 0; i < buckets_indexes.size(); i++) {
                unsigned j = buckets_indexes[i];
                unsigned raw_id = cp.transfer[i];
                //assert(raw_id != -1);
                unsigned cpu_phys_id = load[raw_id].cpu_phys_id;
                _table[j] = cpu_phys_id;
                //assert(cpu_phys_id != -1);
                //assert(cpu_phys_id != min_core_phys);
                //click_chatter("Bucket %d (%d) -> new id %d phys core %d", i ,j, raw_id, cpu_phys_id);
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
                //imbalance.resize(imbalance.size() + 1);
                machine_load.N = machine_load.N+1;
                machine_load.target = machine_load.total_load / (float)machine_load.N;
                int aid = load.size();
                load.push_back(Load(a_phys_id));
                map_phys_to_id[a_phys_id] = aid;
                if (do_numa) {
                    //We disable numa in any case, to allow inter-numa shuffling when we add a core
                    //int numaid =  Numa::get_numa_node_of_cpu(a_phys_id);
                    //sockets_load[numaid].N++;
                    //sockets_load[numaid].target = sockets_load[numaid].total_load / (float)sockets_load[numaid].N;;
                    ////sockets_load[numaid].imbalance.resize(sockets_load[numaid].imbalance.size() + 1);
                    do_numa = false;
                }
            }
        }
    } else {
        if (machine_load.target <  _min_load && !_threshold_force_overload) {
            click_chatter("Underloaded, skipping balancing");
        }
    }

    //We need to fix the first numa socket if numa awareness has been disabled
    if (!do_numa) {
        sockets_load.resize(1);
        sockets_load[0].N = machine_load.N;
        sockets_load[0].target = machine_load.target;
        sockets_load[0].total_load = machine_load.total_load;
        numamax = 1;
    }


    /*
     * Dancers
     * We move the bucket that have more load than XXX 50% to other cores
     * Not shown in paper. Kept for later. Study of doing this VS splitting TC
     */
     if (has_high_load && _dancer) {
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
    //Vector<int> machine_oid;
    //Vector<int> machine_uid;

    //imbalance.resize(machine_load.N);
    //Prepare per-numa socket data structures
    for (int i = 0; i < numamax; i++) {
        SocketLoad &socket = sockets_load[i];
        Vector<float> &imbalance = socket.imbalance;
        imbalance.resize(machine_load.N); //imbalance have "load" ids, not per-socket ids
        //click_chatter("Socket %d/%d has %d cores", i,numamax,socket.N);
    }

    int noverloaded = 0;
    int nunderloaded = 0;
    //Compute the imbalance
    for (unsigned i = 0; i < machine_load.N; i++) {
        int cpuid = load[i].cpu_phys_id;
        int numa = 0;
        if (do_numa) {
            numa = Numa::get_numa_node_of_cpu(cpuid);
        }
        SocketLoad &socket = sockets_load[numa];
        Vector<float> &imbalance = socket.imbalance;

        if (_moved[cpuid]) {
            imbalance[i] = 0;
            _moved[cpuid] = false;
            //click_chatter("Cpu %d just moved",i);
            continue;
        }

        imbalance[i] = 0 * (1.0-_imbalance_alpha) + ( (socket.target - load[i].load) * _imbalance_alpha); //(_last_imbalance[load[i].first] / 2) + ((p.target - load[i].load) / 2.0f);i
        total_imbalance += abs(imbalance[i]);

        if (imbalance[i] > _threshold) {
            if (unlikely(balancer->_verbose))
                click_chatter("Underloaded %d is cpu %d, imb %f, buckets %d",socket.uid.size(),i, imbalance[i],load[i].nbuckets_nz);
            socket.uid.push_back(i);
            nunderloaded++;
        } else if (imbalance[i] < - _threshold) {
            if (load[i].nbuckets_nz == 0) {
                click_chatter("WARNING : A core is overloaded but has no buckets !");
            } else if (load[i].nbuckets_nz > 1) { //Else there is nothing we can do
               if (unlikely(balancer->_verbose))
                    click_chatter("Overloaded %d is cpu %d, imb %f, buckets %d",socket.oid.size(),i, imbalance[i],load[i].nbuckets_nz);
                socket.oid.push_back(i);
                noverloaded++;
            } else {
                if (unlikely(balancer->_verbose))
                    click_chatter("Overloaded cpu %d, imb %f, buckets %d",i, imbalance[i],load[i].nbuckets_nz);
            }
        } else {
            if (unlikely(balancer->_verbose))
                click_chatter("Ignored cpu %d, imb %f, buckets %d",i, imbalance[i],load[i].nbuckets_nz);

        }
    }


    bool moved = false;
    /**
     * Re-balancing
     * We minimize the overall imbalance by moving some buckets from overloaded cores to underloaded cores
     */


#ifndef BALANCE_TWO_PASS

    Vector<Vector<Pair<int,int>>> omoves(noverloaded, Vector<Pair<int,int>>());

    for (int nid = 0; nid < numamax; nid++) {
        SocketLoad &socket = sockets_load[nid];
        Vector<float> &imbalance = socket.imbalance;

        if (!(socket.oid.size() > 0 && socket.uid.size() > 0))
            continue;
            //Count the number of packets for this core
            /*
             * We rebalance all buckets of overloaded to the set of underloaded
             */
            unsigned long long all_overloaded_count = 0;
            for (int u = 0; u < socket.oid.size(); u++) {
                all_overloaded_count += load[socket.oid[u]].npackets;
            }

            struct Bucket{
                int oid_id;
                int bucket_id;
                int cpu_id;
            };

        Vector<Bucket> buckets_indexes;
        Vector<int> oid;
        for (int i = 0; i < socket.oid.size(); i++) {
            if (!do_numa || Numa::get_numa_node_of_cpu(socket.oid[i]) == nid) {
                oid.push_back(socket.oid[i]);
            }
        }
        Vector<int> uid;
        for (int i = 0; i < socket.uid.size(); i++) {
            if (!do_numa || Numa::get_numa_node_of_cpu(socket.uid[i]) == nid) {
                uid.push_back(socket.uid[i]);
            }
        }
        for (int j = 0; j < _table.size(); j++) {
            if (get_node_count(j) == 0) continue;
            for (int u = 0; u < oid.size(); u++) {
                if (map_phys_to_id[_table[j]] == oid[u]) {
                    int numa = 0;
                    if (do_numa) {
                        numa = Numa::get_numa_node_of_cpu(_table[j]);
                    }
                    buckets_indexes.push_back(Bucket{.oid_id =u,.bucket_id =j, .cpu_id = oid[u]});
                    break;
                }
            }
        }

        if (buckets_indexes.size() == 0)
            continue;
        BucketMapTargetProblem pm(buckets_indexes.size(), uid.size(), oid.size());

        //Compute load for each buckets
        for (int i = 0; i < buckets_indexes.size(); i++) {
            Bucket& b = buckets_indexes[i];
            int j = b.bucket_id;
            double c = get_node_count(j);
            pm.buckets_load[i] = (((double) c / (double)load[b.cpu_id].npackets)*load[b.cpu_id].load);
            pm.buckets_max_idx[i] = b.oid_id;
            //click_chatter("Bucket %d (%d) load %f", i ,j, b);
        }

        for (int i = 0; i < uid.size(); i++) {
            pm.target[i] = imbalance[uid[i]];
        }

        for (int i = 0; i < oid.size(); i++) {
            pm.max[i] = -imbalance[oid[i]];
        }

        pm.solve(balancer);

        for (int i = 0; i < pm.buckets_load.size(); i++) {
            Bucket& b = buckets_indexes[i];
            int to_uid = pm.transfer[i];
            if (to_uid == -1) continue;

            if (unlikely(balancer->_verbose > 2))
                click_chatter("B idx %d to ucpu %d",i,to_uid);
            int to_cpu = uid[to_uid];
            int from_cpu = oid[pm.buckets_max_idx[i]];
            imbalance[from_cpu] += pm.buckets_load[i];
            imbalance[to_cpu] -= pm.buckets_load[i];

            if (unlikely(balancer->_verbose))
                click_chatter("Move bucket %d from core %d to core %d",b.bucket_id,_table[b.bucket_id], load[to_cpu].cpu_phys_id);

            if (balancer->_manager) {
                omoves[pm.buckets_max_idx[i]].push_back(Pair<int,int>(b.bucket_id, load[to_cpu].cpu_phys_id));
            }
            _table[b.bucket_id] = load[to_cpu].cpu_phys_id;

            _moved[load[to_cpu].cpu_phys_id] = true;
            moved = true;
        }
        if (balancer->_manager) {
            for (int i = 0; i < oid.size(); i++) {
                if (omoves[i].size() > 0) {
                    int from_phys_id = load[oid[i]].cpu_phys_id;
                    balancer->_manager->pre_migrate((DPDKDevice*)_fd,from_phys_id,omoves[i]);
                }
            }
        }
    }

#else
        if (!(socket.oid.size() > 0 && socket.uid.size() > 0))
        //Set problem parameters
        float min_cost = socket.N * socket.N;
        p.transfer.resize(socket.N);
        for (unsigned i = 0; i < socket.N; i++) {
            p.transfer[i] = i;//Default is to transfer load to itself
        }
        p.min_cost = socket.N;

        //Solve assignment of imbalance of overloaded -> underloaded
        p.solve();

        //We "correct" the imbalance. Indeed even from the best solution, moving overloaded
        // to underloaded may lead to a unoptimal local scenario, we could move just
        // a little less of every overloaded core to the single underloaded
        Vector<float> fT = p.fixedTransfert();


        if (unlikely(balancer->_verbose))
            click_chatter("Transfer solution for %d uid, %d oid:",socket.uid.size(), socket.oid.size());

        //We start at a random core index, as the first cores to be balanced may be advantaged
        int iRand = click_random();

        for (unsigned iOffset = 0; iOffset < socket.N; iOffset++) {
            unsigned i = (iRand + iOffset) % socket.N;
            if (unlikely(balancer->_verbose))
                click_chatter("Core %d (load %f, imbalance %f (corr %f)) -> %d (load %f)", i, load[i].load, imbalance[i], fT[i], p.transfer[i], load[p.transfer[i]].load);
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
                if (bload > -imbalance[i] * fT[i] || bload > 0.5 ) { //Never over balance! Also dancers are handled separately
                    //if (bload > _threshold)
                    if (!miss) {
                        //if (unlikely(balancer->_verbose))
                        //click_chatter("Trying to fill %f of imbalance", imbalance[i]*fT[i]);
                        miss =true;
                    }
                    if (unlikely(balancer->_verbose > 1)) {
                        click_chatter("Skipping bucket %d var %f, load %f%% (min %f%% of bucket, %f%% real), imbalance %f, corrected imbalance %f, at least %d flows",j, var, bload*100, bpc*100, ((c * 100.0)/core_tot),imbalance[i], imbalance[i]*fT[i],-1);//get_node(j).flows);
                    }
                    if (bload < -imbalance[i] && bload < 0.5)
                        nSkipped ++;
                        /*
                        left_max += bload * (1 + var);
                        left += bload;
                        left_min += bload / (1 + var);*/

                    continue;
                } else
                    miss = false;
                if (imbalance[i]*fT[i] > -_imbalance_threshold)
                    break;


                tot_bload += bload;
                tot_bpc += bpc;
                //click_chatter("Bucket %d var %f, load %f, imbalance %f",j, var, bload,imbalance[i]);

                moved = true;
                _table[j] = load[p.transfer[i]].cpu_phys_id;

                imbalance[i] += bload;
                //imbalance[p.transfer[i]] -= bload;
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

    if (moved) {
        for (int nid = 0; nid < numamax; nid++) {
            SocketLoad &socket = sockets_load[nid];
            for (unsigned i = 0; i < socket.N; i++) {
                total_imbalance += abs(socket.imbalance[i]);
                if (abs(socket.imbalance[i]) > _threshold * 2 ) {
                    if (unlikely(balancer->_verbose))
                        click_chatter("Imbalance of core %d left to %f, that's quite a MISS.", i, socket.imbalance[i]);
                }
            }
        }

        Timestamp t = Timestamp::now_steady();
        auto v = (t-begin).usecval();
        if (unlikely(balancer->_verbose || v > 100))
            click_chatter("Solution computed in %d usec", v);

        update_reta();
        if (balancer->_manager) {
            for (int nid = 0; nid < numamax; nid++) {
                SocketLoad &socket = sockets_load[nid];
                for (int i = 0; i < socket.oid.size(); i++) {
                    if (omoves[i].size() > 0) {
                        int from_phys_id =  load[socket.oid[i]].cpu_phys_id;
                        balancer->_manager->post_migrate((DPDKDevice*)_fd,from_phys_id);
                    }
                }
            }
        }
    }


    if (!_counter_is_xdp)
        ((AggregateCounterVector*)_counter)->advance_epoch();
    else {
#ifdef HAVE_BPF
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
#endif
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
             else
                 tot = 1;
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
                        assert(_table[i] >= 0);
                        //click_chatter("%d->%d",i,_table[i]);
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
            .read_or_set("DANCER", _dancer, false)
            .read_or_set("NUMA", _numa, false)
            .consume() < 0)
        return -1;
    _target_load = t;
    _imbalance_alpha = i;
    _threshold = threshold;
    if (_numa)
        _numa_num = Numa::get_max_numas();
    else
        _numa_num = 1;

    int err = BalanceMethodRSS::configure(conf, errh);
    if (err != 0)
        return err;

    //Reta_size must be set before this
    if (e) {
        _counter = (AggregateCounterVector*)e->cast("AggregateCounterVector");
        if (!_counter) {
#ifdef HAVE_BPF
            _counter = (XDPLoader*)e->cast("XDPLoader");
#endif
            if (!_counter) {
                return errh->error("COUNTER must be of the type AggregateCounterVector or XDPLoader");
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
        .read("MANAGER",ElementArg(), manager)
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
    for (int i = 0; i < _stats.weight(); i++) {
        _stats.get_value_for_thread(i).resize(10);
    }
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
} //end of balancing timer

DeviceBalancer::CpuInfo
DeviceBalancer::make_info(int _id) {
    CpuInfo i;
    i.id = _id;
    i.last_cycles = master()->thread(_id)->useful_kcycles();
    return i;
}


enum {h_active, h_autoscale, h_force_cpu, h_run_stats, h_run_stats_imbalance_first, h_run_stats_time_first};


String
DeviceBalancer::read_param(Element *e, void *thunk)
{
    DeviceBalancer *td = (DeviceBalancer *)e;
    switch((uintptr_t) thunk) {
    case h_active:
        return String(td->_active);
    case h_autoscale:
        return String(td->_autoscale);
    case h_run_stats: {
        StringAccum acc;
        for (int j = 0; j < 10; j++) {
            int count = 0;
            double imbalance = 0;
            uint64_t time = 0;
            for (int i = 0; i < td->_stats.weight(); i++) {
                Vector<DeviceBalancer::RunStat> &v = td->_stats.get_value_for_thread(i);
                count += v[j].count;
                imbalance += v[j].imbalance;
                time += v[j].time;
            }
            if (count > 0)
                acc << j << " " << count << " " << imbalance/count << " " << time/count << "\n";
            else
                acc << j << " 0 nan nan\n";
        }
        return acc.take_string();
       }
   case h_run_stats_imbalance_first: {
            int count = 0;
            double imbalance = 0;
            for (int i = 0; i < td->_stats.weight(); i++) {
                Vector<DeviceBalancer::RunStat> &v = td->_stats.get_value_for_thread(i);
                count += v[0].count;
                imbalance += v[0].imbalance;
            }
            return String(imbalance/count);
      }
   case h_run_stats_time_first: {
            int count = 0;
            uint64_t time = 0;
            for (int i = 0; i < td->_stats.weight(); i++) {
                Vector<DeviceBalancer::RunStat> &v = td->_stats.get_value_for_thread(i);
                count += v[0].count;
                time += v[0].time;
            }
            if (count > 0) {
                return String(time/count);
            } else
                return "nan";
      }
    default:
        return String();
    }
}

int
DeviceBalancer::write_param(const String &in_s, Element *e, void *vparam,
                     ErrorHandler *errh)
{
    DeviceBalancer *db = (DeviceBalancer *)e;
    String s = cp_uncomment(in_s);
    switch ((intptr_t)vparam) {
    case h_active: {
        bool active;
        if (!BoolArg().parse(s, active))
            return errh->error("type mismatch");
        if (active && !db->_active)
        db->_timer.schedule_after_msec(db->_current_tick);

        db->_active = active;
        break;
    }
    case h_autoscale: {
        bool active;
        if (!BoolArg().parse(s, active))
            return errh->error("type mismatch");
        db->_autoscale = active;
        break;
    }
    case h_force_cpu: {
    int cpus;
        if (!IntArg().parse(s, cpus))
            return errh->error("type mismatch");
        bool moved = false;
        if (cpus > db->_used_cpus.size())
        for (int i = 0; i < cpus-db->_used_cpus.size(); i++) {
            if (db->_available_cpus.size() > 0) {
                db->addCore();
                moved = true;
            }
        }
        else if (cpus < db->_used_cpus.size()) {
            //in_s = "UNSUPPORTED";
            return -1;
        //for (int i = 0; i < db->_used_cpus.size() - cpus; i++) {

            /*if (db->_used_cpus.size() > 0) {
                db->removeCore();
                moved = false;
            }*/
        //}
        }
        if (moved) {
            db->_method->cpu_changed();
        }
        //in_s = cpus-db->_used_cpus.size();
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
    add_write_handler("cpus", write_param, h_force_cpu);
    add_read_handler("run_stats", read_param, h_run_stats);
    add_read_handler("run_stats_imbalance_first", read_param, h_run_stats_imbalance_first);
    add_read_handler("run_stats_time_first", read_param, h_run_stats_time_first);
}


CLICK_ENDDECLS

EXPORT_ELEMENT(DeviceBalancer)
