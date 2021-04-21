/**
 * Metron
 */
#include <click/flowdispatcher.hh>

String rewrite_id(String rule, int core) {
    int pos = rule.find_left("queue index");
    int end = pos + 12 + rule.substring(pos+12).find_left(' ');
    return rule.substring(0,pos)+ "queue index "+String(core)+ rule.substring(end);

}

int MethodMetron::initialize(ErrorHandler *errh, int startwith) {
    if (!_is_dpdk)
        return errh->error("Metron only works with DPDK");
    FlowDispatcher *flow_dir = FlowDispatcher::get_flow_dispatcher(((DPDKDevice*)_fd)->port_id);
    assert(flow_dir);

    // Invoke Flow Dispatcher only if active
    if (flow_dir->active()) {
        // There is a file with (user-defined) rules
        if (!_rules_file.empty()) {
            HashMap<uint32_t, String> rules_map;
            const String rules_str = (const String) flow_dir->load_rules_from_file_to_string(_rules_file);

            if (rules_str.empty()) {
                return errh->error("Failed to add rules due to empty input from file");
            }

            // Tokenize them to facilitate the insertion in the flow cache
            Vector<String> rules_vec = rules_str.trim_space().split('\n');

            for (uint32_t i = 0; i < (unsigned)rules_vec.size(); i++) {
                String rule = rules_vec[i] + "\n";

                // Add rule to the map
                rules_map.insert(i, rewrite_id(rule,i % startwith));
            }

            int ret = flow_dir->update_rules(rules_map, true);

            if (ret < 0)
                return errh->error("Could not install rules.");
            else
                click_chatter("Installed %d rules", ret);
        } else
            return errh->error("No rule file !");
    } else {
        return errh->error("Flow dispatcher is not active.");
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

void MethodMetron::rebalance(std::vector<std::pair<int,float>> load) {
    FlowDispatcher *flow_dir = FlowDispatcher::get_flow_dispatcher(((DPDKDevice*)_fd)->port_id);
    click_jiffies_t now = click_jiffies();
    std::vector<std::pair<int,float>> underloaded;
    std::vector<std::pair<int,float>> overloaded;
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
        if (load_future > balancer->overloaded_thresh() || load_current > balancer->overloaded_thresh()) {
            overloaded.push_back(std::pair<int,float>{cpuid,load_current});
        } else if (load_future < balancer->underloaded_thresh() && load_current < balancer->underloaded_thresh() && load_past < balancer->underloaded_thresh()) {
            underloaded.push_back(std::pair<int,float>{cpuid,load_current});
        }
    }


    click_chatter("Have %d overloaded and %d underloaded",overloaded.size(), underloaded.size());

    std::sort(overloaded.begin(), overloaded.end(), [](std::pair<int,float> a, std::pair<int,float> b) {return a.second < b.second; });

    std::sort(underloaded.begin(), underloaded.end(), [load](std::pair<int,float> a, std::pair<int,float> b) {return a.second > b.second; });

    while (overloaded.size() > 0) {
        int o = overloaded[overloaded.size() - 1].first;
        float load_o = overloaded[overloaded.size() - 1].second;
        overloaded.pop_back();
        int a_phys_id;
        float load_a = 0;
        if (balancer->autoscale() && balancer->num_spare_cpus() > 0) {
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


        HashMap<uint32_t, String> * rmap;
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

        std::vector<uint32_t> dlist;
        HashMap<uint32_t, String> nmap;
//      flow_dir->flow_cache()->delete_rule_by_global_id(rmap.keys());
//      flow_dir->flow_cache()->
        unsigned n = 0;
        auto it = rmap->begin();
        while (it != rmap->end()) {
            if (n++ < rmap->size() - mig) {
                it++;
                continue;
            }

            uint32_t rule_id = it.key();
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

    while (underloaded.size() >= 2 && balancer->autoscale()) {
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
        if (load_r + load_i >= balancer->overloaded_thresh()) {
        click_chatter("ERROR : too much load, select another core");
                    continue;
        }
        balancer->removeCore(remove_i);
//        click_chatter("Core %d is now available", remove_i);

        HashMap<uint32_t, String> * rmap;

        HashMap<uint32_t, String> nmap;
        auto cache = flow_dir->get_flow_cache();

        rmap = cache->rules_map_by_core_id(remove_i);
//        click_chatter("Rmap %p", rmap);
        click_chatter("Migrating %d flows to core %d", rmap->size(), with_i);
//      flow_dir->flow_cache()->delete_rule_by_global_id(rmap.keys());
//      flow_dir->flow_cache()->

        auto it = rmap->begin();
        while (it != rmap->end()) {
            uint32_t rule_id = it.key();
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
