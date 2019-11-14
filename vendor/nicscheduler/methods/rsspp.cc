/**
 * RSS++
 */

#include "solver.hh"

MethodRSSPP::MethodRSSPP(NICScheduler* b, EthernetDevice* fd) : MethodRSS(b,fd) {
}


int MethodRSSPP::initialize(ErrorHandler* errh, int startwith) {
    if (MethodRSS::initialize(errh, startwith) != 0)
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

    std::vector<float> imbalance;
    std::vector<int> uid;
    std::vector<int> oid;
};

#define print_cpu_assign() {for (int i = 0; i < rload.size(); i++) {click_chatter("CPU %d -> %d/%d",i, rload[i].first,load[i].cpu_phys_id);}}

void MethodRSSPP::rebalance(std::vector<std::pair<int,float>> rload) {
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
            if (unlikely(tot != 0 && balancer->verbose()))
                click_chatter("Key %d core %d val %d, var %f",key,_table[key], tot,(float)min(tot,var) / (float)max(tot,var));
        }
    }
#endif

#define get_node_count(i) ((_counter_is_xdp)?_count[i].count:((AggregateCounterVector*)_counter)->find_node_nocheck(i).count)
#define get_node_variance(i) ((_counter_is_xdp)?_count[i].variance:((AggregateCounterVector*)_counter)->find_node_nocheck(i).variance)
#define get_node_moved(i) (_count[i).moved]

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
    std::vector<Load> load(rload.size(),Load());

    //Track physical CPU id to assigned ids
    std::vector<unsigned> map_phys_to_id(click_max_cpu_ids(), -1);

    //Per NUMA-socket load. When do_numa is false, we consider a single NUMA node
    std::vector<SocketLoad> sockets_load;
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
#if HAVE_NUMA
            int numaid =  Numa::get_numa_node_of_cpu(load[j].cpu_phys_id);
#else
            int numaid = 0;
#endif
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
    if (unlikely(balancer->verbose())) {
        click_chatter("Target %f. Total load %f. %d cores. Suppload %f, Var %f", machine_load.target, machine_load.total_load, machine_load.N, suppload, var_protect);
    }


    //Count the number of packets for each core
    for (int j = 0; j < _table.size(); j++) {
        unsigned long long c = get_node_count(j);
        unsigned cpu_phys_id =_table[j];
        if (unlikely(cpu_phys_id >= map_phys_to_id.size())) {
            click_chatter("ERROR : invalid phys id %d", cpu_phys_id);
            print_cpu_assign();
            abort();
        }
        unsigned cpuid = map_phys_to_id[cpu_phys_id];
        if (unlikely(cpuid >= load.size())) {
            click_chatter("ERROR : invalid cpu id %d for physical id %d, table index %d", cpuid, cpu_phys_id, j);
            print_cpu_assign();
            abort();
        }
        auto &l = load[cpuid];
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
    if (balancer->autoscale()) {
        if (suppload > (1 + (1 - _target_load) + var_protect)  && machine_load.N > 1) { // we can remove a core
            if (unlikely(balancer->verbose())) {
                click_chatter("Removing a core (target %f)",_target_load);
            }

            int min_core_phys = load[min_core].cpu_phys_id;

            click_chatter("Removing core %d (phys %d)", min_core, min_core_phys);
            unsigned long long totcount = load[min_core].npackets;
            load[min_core] = load[load.size() - 1];
            load.pop_back();
            balancer->removeCore(min_core_phys);

            //Count the number of packets for this core
            std::vector<int> buckets_indexes;
            for (int j = 0; j < _table.size(); j++) {
                if (_table[j] == min_core_phys) {
                    buckets_indexes.push_back(j);
                }
            }
            if (unlikely(balancer->verbose()))
                click_chatter("Removing %d buckets of core %d", buckets_indexes.size(), min_core);


            BucketMapProblem cp(buckets_indexes.size(), load.size()); //load.size() is already fixed

            //Compute load for each buckets
            for (int i = 0; i < buckets_indexes.size(); i++) {
                int j = buckets_indexes[i];
                double c = get_node_count(j);
                cp.buckets_load[i] = (((double) c / (double)totcount)*min_core_load);
            }

            //Add imbalance of all cores
            std::vector<float> imbalance(machine_load.N-1,0);
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
            if (unlikely(balancer->verbose()))
                click_chatter("Adding a core as load is %f");
            int a_phys_id = balancer->addCore();
            if (a_phys_id == -1) {
                if (unlikely(balancer->verbose()))
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
        if (unlikely(balancer->verbose()))
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

                if (unlikely(balancer->verbose()))
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
    //std::vector<int> machine_oid;
    //std::vector<int> machine_uid;

    //imbalance.resize(machine_load.N);
    //Prepare per-numa socket data structures
    for (int i = 0; i < numamax; i++) {
        SocketLoad &socket = sockets_load[i];
        std::vector<float> &imbalance = socket.imbalance;
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
#if HAVE_NUMA
            numa = Numa::get_numa_node_of_cpu(cpuid);
#endif
        }
        SocketLoad &socket = sockets_load[numa];
        std::vector<float> &imbalance = socket.imbalance;

        if (_moved[cpuid]) {
            imbalance[i] = 0;
            _moved[cpuid] = false;
            //click_chatter("Cpu %d just moved",i);
            continue;
        }

        imbalance[i] = 0 * (1.0-_imbalance_alpha) + ( (socket.target - load[i].load) * _imbalance_alpha); //(_last_imbalance[load[i].first] / 2) + ((p.target - load[i].load) / 2.0f);i
        total_imbalance += abs(imbalance[i]);

        if (imbalance[i] > _threshold) {
            if (unlikely(balancer->verbose()))
                click_chatter("Underloaded %d is cpu %d, imb %f, buckets %d",socket.uid.size(),i, imbalance[i],load[i].nbuckets_nz);
            socket.uid.push_back(i);
            nunderloaded++;
        } else if (imbalance[i] < - _threshold) {
            if (load[i].nbuckets_nz == 0) {
                click_chatter("WARNING : A core is overloaded but has no buckets !");
            } else if (load[i].nbuckets_nz > 1) { //Else there is nothing we can do
               if (unlikely(balancer->verbose()))
                    click_chatter("Overloaded %d is cpu %d, imb %f, buckets %d",socket.oid.size(),i, imbalance[i],load[i].nbuckets_nz);
                socket.oid.push_back(i);
                noverloaded++;
            } else {
                if (unlikely(balancer->verbose()))
                    click_chatter("Overloaded cpu %d, imb %f, buckets %d",i, imbalance[i],load[i].nbuckets_nz);
            }
        } else {
            if (unlikely(balancer->verbose()))
                click_chatter("Ignored cpu %d, imb %f, buckets %d",i, imbalance[i],load[i].nbuckets_nz);

        }
    }


    bool moved = false;
    /**
     * Re-balancing
     * We minimize the overall imbalance by moving some buckets from overloaded cores to underloaded cores
     */


#ifndef BALANCE_TWO_PASS

    std::vector<std::vector<std::pair<int,int>>> omoves(noverloaded, std::vector<std::pair<int,int>>());

    for (int nid = 0; nid < numamax; nid++) {
        SocketLoad &socket = sockets_load[nid];
        std::vector<float> &imbalance = socket.imbalance;

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

        std::vector<Bucket> buckets_indexes;
        std::vector<int> oid;
        for (int i = 0; i < socket.oid.size(); i++) {
            if (!do_numa
#if HAVE_NUMA
                    || Numa::get_numa_node_of_cpu(socket.oid[i]) == nid
#endif
                    ) {
                oid.push_back(socket.oid[i]);
            }
        }
        std::vector<int> uid;
        for (int i = 0; i < socket.uid.size(); i++) {
            if (!do_numa
#if HAVE_NUMA
                    || Numa::get_numa_node_of_cpu(socket.uid[i]) == nid
#endif
                    ) {
                uid.push_back(socket.uid[i]);
            }
        }
        for (int j = 0; j < _table.size(); j++) {
            if (get_node_count(j) == 0) continue;
            for (int u = 0; u < oid.size(); u++) {
                if (map_phys_to_id[_table[j]] == oid[u]) {
                    int numa = 0;
                    if (do_numa) {
#if HAVE_NUMA
                        numa = Numa::get_numa_node_of_cpu(_table[j]);
#endif
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

            if (unlikely(balancer->verbose() > 2))
                click_chatter("B idx %d to ucpu %d",i,to_uid);
            int to_cpu = uid[to_uid];
            int from_cpu = oid[pm.buckets_max_idx[i]];
            imbalance[from_cpu] += pm.buckets_load[i];
            imbalance[to_cpu] -= pm.buckets_load[i];

            if (unlikely(balancer->verbose()))
                click_chatter("Move bucket %d from core %d to core %d",b.bucket_id,_table[b.bucket_id], load[to_cpu].cpu_phys_id);

            if (balancer->_manager) {
                omoves[pm.buckets_max_idx[i]].push_back(std::pair<int,int>(b.bucket_id, load[to_cpu].cpu_phys_id));
            }
            _table[b.bucket_id] = load[to_cpu].cpu_phys_id;

            _moved[load[to_cpu].cpu_phys_id] = true;
            moved = true;
        }
        if (balancer->_manager) {
            for (int i = 0; i < oid.size(); i++) {
                if (omoves[i].size() > 0) {
                    int from_phys_id = load[oid[i]].cpu_phys_id;
                    balancer->_manager->pre_migrate((DPDKEthernetDevice*)_fd,from_phys_id,omoves[i]);
                }
            }
        }
    }

#else

        //OLD VERSION, STOP LOOKING
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
        std::vector<float> fT = p.fixedTransfert();


        if (unlikely(balancer->verbose()))
            click_chatter("Transfer solution for %d uid, %d oid:",socket.uid.size(), socket.oid.size());

        //We start at a random core index, as the first cores to be balanced may be advantaged
        int iRand = click_random();

        for (unsigned iOffset = 0; iOffset < socket.N; iOffset++) {
            unsigned i = (iRand + iOffset) % socket.N;
            if (unlikely(balancer->verbose()))
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
            /*auto cmp = [](std::pair<int,float> left, std::pair<int,float> right) { return left.second > right.second;};
            std::priority_queue<std::pair<int,float>, std::vector<std::pair<int,float> >, decltype(cmp)> q(cmp);*/
            /*struct BucketInfo {
                int index;
                float var;
                int flows;
            }*/
            //std::vector<BucketInfo> q;

            //We find all buckets of this core, starting at a random point, so we do not
            //always move the same ones
            std::vector<std::pair<int,float> > q; //List of (i core) bucket->load
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
                    q.push_back(std::pair<int,float>(j, var));
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
                        //if (unlikely(balancer->verbose()))
                        //click_chatter("Trying to fill %f of imbalance", imbalance[i]*fT[i]);
                        miss =true;
                    }
                    if (unlikely(balancer->verbose() > 1)) {
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
            if (unlikely(balancer->verbose()))
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
                    if (unlikely(balancer->verbose()))
                        click_chatter("Imbalance of core %d left to %f, that's quite a MISS.", i, socket.imbalance[i]);
                }
            }
        }

        Timestamp t = Timestamp::now_steady();
        auto v = (t-begin).usecval();
        if (unlikely(balancer->verbose() || v > 100))
            click_chatter("Solution computed in %d usec", v);

        update_reta();
        if (balancer->_manager) {
            for (int nid = 0; nid < numamax; nid++) {
                SocketLoad &socket = sockets_load[nid];
                for (int i = 0; i < socket.oid.size(); i++) {
                    if (omoves[i].size() > 0) {
                        int from_phys_id =  load[socket.oid[i]].cpu_phys_id;
                        balancer->_manager->post_migrate((DPDKEthernetDevice*)_fd,from_phys_id);
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
            balancer->set_tick(balancer->tick_min());
        else
            balancer->set_tick(balancer->current_tick() / 2);
    } else if (!moved) {
        balancer->set_tick(balancer->current_tick() * 2);
    }

    if (balancer->current_tick() > balancer->tick_max())
        balancer->set_tick(balancer->tick_max());
    if (balancer->current_tick() < balancer->tick_min())
        balancer->set_tick(balancer->tick_min());

}

bool MethodRSS::update_reta_flow(bool validate) {
        int port_id = ((DPDKEthernetDevice*)_fd)->port_id;
        struct rte_flow_attr attr;
            std::vector<rte_flow*> newflows;
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

                std::vector<rte_flow_item> pattern;
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
                if (unlikely(balancer->verbose()))
                    click_chatter("Flow added succesfully with %d patterns!", pattern.size());
                } else {
                if (unlikely(balancer->verbose()))
                    click_chatter("Could not add pattern with %d patterns, error %d : %s", pattern.size(),  res, error.message);
                    return false;
                }

                newflows.push_back(flow);
            } else {
            if (unlikely(balancer->verbose()))
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

bool MethodRSS::update_reta(bool validate) {
    Timestamp t = Timestamp::now_steady();
    /*if (_verifier) {
        _verifier->_table = _table;
    }*/

    if (_update_reta_flow) {
        if (!update_reta_flow(validate))
            return false;
    } else {
        if (!_fd->set_rss_reta(_fd, _table.data(), _table.size()))
            return false;
    }

    Timestamp s = Timestamp::now_steady();
    if (validate || balancer->verbose())
    click_chatter("Reta updated in %d usec",(s-t).usecval());
    return true;
}
