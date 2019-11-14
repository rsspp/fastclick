/**
 * RSS base
 */
MethodRSS::MethodRSS(NICScheduler* b, EthernetDevice* fd) :
    BalanceMethodDevice(b,fd),
    //_verifier(0),
    _isolate(0) {
}

int MethodRSS::initialize(ErrorHandler *errh, int startwith) {
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
    _fd->set_rss_reta(_fd, _table.data(), _table.size());

    if (_is_dpdk) {
        int port_id = ((DPDKEthernetDevice*)_fd)->get_port_id();

        _rss_conf.rss_key = (uint8_t*)CLICK_LALLOC(128);
        _rss_conf.rss_key_len = 128; //This is only a max
        rte_eth_dev_rss_hash_conf_get(port_id, &_rss_conf);


        struct rte_flow_error error;
        rte_eth_dev_stop(port_id);
        //rte_eth_promiscuous_disable(port_id);
        int res = rte_flow_isolate(port_id, _isolate, &error);
        if (res != 0)
            errh->warning("Warning %d : Could not set isolated mode because %s !",res,error.message);

        rte_eth_dev_start(port_id);
    }

    for (int i = 0; i < _table.size(); i++) {
        _table[i] = i % startwith;
    }

    _fd->set_rss_reta(_fd, _table.data(), _table.size());
    click_chatter("RSS initialized with %d CPUs and %d buckets", startwith, _table.size());
    int err = BalanceMethodDevice::initialize(errh, startwith);
    if (err != 0)
        return err;

    _update_reta_flow = true;
    if (_is_dpdk) {
       if (!update_reta_flow(true)) {
            _update_reta_flow = false;
            if (_fd->set_rss_reta(_fd, _table.data(), _table.size()) != 0)
                return errh->error("Neither flow RSS or global RSS works to program the RSS table.");
       } else
           click_chatter("RETA update method is flow");
    } else {
        _update_reta_flow = false;
        if (_fd->set_rss_reta(_fd, _table.data(), _table.size()) != 0)
            return errh->error("Cannot program the RSS table.");
    }
    if (!_update_reta_flow)  {
        click_chatter("RETA update method is global");
    }

    return err;
}

void MethodRSS::rebalance(std::vector<std::pair<int,float>> load) {
    //update_reta();
}

void MethodRSS::cpu_changed() {
    int m =  balancer->num_used_cpus();
    std::vector<std::vector<std::pair<int,int>>> omoves(balancer->num_used_cpus(), std::vector<std::pair<int,int>>());
    /*std::vector<int> epochs;
    epochs.resize(max_cpus());*/


    for (int i = 0; i < _table.size(); i++) {
        int newcpu = balancer->get_cpu_info(i % m).id;
        if (balancer->_manager && newcpu!= _table[i]) {
            omoves[_table[i]].push_back(std::pair<int,int>(i, newcpu));
        }
        //epochs(_table[i]) =
        _table[i] = newcpu;
    }

    if (balancer->_manager) {
        for (int i = 0; i < m; i++) {
            if (omoves[i].size() > 0) {
                balancer->_manager->pre_migrate((DPDKEthernetDevice*)_fd, i, omoves[i]);
            }
        }
    }
    click_chatter("Migration info written. Updating reta.");
    update_reta();
    click_chatter("Post migration");
    if (balancer->_manager) {
        for (int i = 0; i < m; i++) {
            if (omoves[i].size() > 0) {
                balancer->_manager->post_migrate((DPDKEthernetDevice*)_fd, i);
            }
        }
    }
    click_chatter("Post migration finished");
}
