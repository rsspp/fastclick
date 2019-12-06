/*
 * FlowIPManager.{cc,hh}
 */

#include <click/config.h>
#include <click/glue.hh>
#include <click/args.hh>
#include <click/ipflowid.hh>
#include <click/routervisitor.hh>
#include <click/straccum.hh>
#include "flowipmanager.hh"
#include <rte_hash.h>
#include <rte_hash_crc.h>
#include "../userlevel/fromdpdkdevice.hh"
#include <rte_ethdev.h>

CLICK_DECLS

FlowIPManager::FlowIPManager() : _verbose(1), _tables(0), _groups(0), _def_thread(0) {

}

FlowIPManager::~FlowIPManager() {

}

int
FlowIPManager::configure(Vector<String> &conf, ErrorHandler *errh)
{
    _def_thread = click_max_cpu_ids();

    if (Args(conf, this, errh)
    		.read_or_set_p("GROUPS", _groups, 512)
    		.read_or_set_p("CAPACITY", _table_size, 65536)
            .read_or_set("RESERVE",_reserve, 0)
            .read("DEF_THREAD", _def_thread)
            .read("VERBOSE", _verbose)
            .complete() < 0)
        return -1;

    find_children(_verbose);

    router()->get_root_init_future()->postOnce(&_fcb_builded_init_future);
    _fcb_builded_init_future.post(this);

    return 0;
}

static inline uint32_t
ipv4_hash_crc(const void *data, __rte_unused uint32_t data_len,
                uint32_t init_val)
{
        const IPFlow5ID *k;
        uint32_t t;
        const uint32_t *p;
        k = (const IPFlow5ID *)data;
        t = k->proto();
        p = ((const uint32_t *)k) + 2;
        init_val = rte_hash_crc_4byte(t, init_val);
        init_val = rte_hash_crc_4byte(k->saddr(), init_val);
        init_val = rte_hash_crc_4byte(k->daddr(), init_val);
        init_val = rte_hash_crc_4byte(*p, init_val);
        return init_val;
}

int FlowIPManager::solve_initialize(ErrorHandler *errh) {
    struct rte_hash_parameters hash_params = {0};
    char buf[32];
    hash_params.name = buf;
    hash_params.entries = _table_size;
    hash_params.key_len = sizeof(IPFlow5ID);
    hash_params.hash_func = ipv4_hash_crc;
    hash_params.hash_func_init_val = 0;
    hash_params.extra_flag = 0; //RTE_HASH_EXTRA_FLAGS_MULTI_WRITER_ADD | RTE_HASH_EXTRA_FLAGS_RW_CONCURRENCY; //| RTE_HASH_EXTRA_FLAGS_RW_CONCURRENCY_LF

    _flow_state_size_full = sizeof(FlowControlBlock) + _reserve;

    _tables = CLICK_ALIGNED_NEW(gtable,_groups);
    CLICK_ASSERT_ALIGNED(_tables);

    for (int i = 0; i < _groups; i++) {
        sprintf(buf, "flowipmanager%d", i);
        _tables[i].hash = rte_hash_create(&hash_params);
        if (!_tables[i].hash)
            return errh->error("Could not init flow table %d!", i);

        _tables[i].fcbs =  (FlowControlBlock*)CLICK_ALIGNED_ALLOC(_flow_state_size_full * _table_size);
        bzero(_tables[i].fcbs,_flow_state_size_full * _table_size);
        CLICK_ASSERT_ALIGNED(_tables[i].fcbs);
        if (!_tables[i].fcbs)
            return errh->error("Could not init data table %d!", i);
        if (_def_thread > 0)
            _tables[i].owner = i % _def_thread;
    }
    click_chatter("%p{element} initialized with %d groups and %d threads", this, _groups, _def_thread);

    return Router::InitFuture::solve_initialize(errh);
}

void FlowIPManager::cleanup(CleanupStage stage) {

}


void FlowIPManager::pre_migrate(DPDKEthernetDevice* dev, int from, std::vector<std::pair<int,int>> gids) {
	CoreInfo &coref = _cores.get_value_for_thread(from);
    coref.lock.acquire();
    for (int i = 0;i < gids.size(); i++) {
        coref.moves.push_back(gids[i]);
    }
    coref.pending = 1;
    coref.lock.release();
}


void FlowIPManager::post_migrate(DPDKEthernetDevice* dev, int from) {
    int port_id = dev->get_port_id();
    int v = rte_eth_rx_queue_count(port_id, from);

    CoreInfo &coref = _cores.get_value_for_thread(from);
    uint64_t w = coref.count + v;

    if (coref.watch < w) {
        coref.lock.acquire();
        if (coref.pending == 1) {
            coref.watch = w;
        } else {
            click_chatter("Useless post migration for cpu %d. Already done...", from);
        }
        coref.lock.release();
    }
    //TODO fire migration task so if v == 0 or very low we do not wait for the next packet to migrate
}

/**
 * Classify a packet in the given group
 * FCB may change!
 */
void FlowIPManager::process(int groupid, Packet* p, BatchBuilder& b) {
    IPFlow5ID fid = IPFlow5ID(p);
    bool first = false;
    gtable& t = _tables[groupid];
    rte_hash*& table = t.hash;

    int ret = rte_hash_lookup(table, &fid);

    if (ret < 0) { //new flow
        ret = rte_hash_add_key(table, &fid);
        if (ret < 0) {
            click_chatter("Cannot add key (have %d items)!", rte_hash_count(table));
            return;
        }
        first = true;
    }

    if (unlikely(_verbose > 2))
        click_chatter("Packet of flow group %d, id %d", groupid, ret);

    if (b.last == ret) {
        b.append(p);
    } else {
        PacketBatch* batch;
        batch = b.finish();
        if (batch) {
            fcb_stack->acquire(batch->count());
            output_push_batch(0, batch);
        }
        fcb_stack = (FlowControlBlock*)((unsigned char*)t.fcbs + _flow_state_size_full * ret);
        b.init();
        b.append(p);
    }
}

inline void FlowIPManager::flush_queue(int groupid, BatchBuilder &b) {
    if (unlikely(_tables[groupid].queue)) {
        if (_verbose > 0)
            click_chatter("Flush table %d from core %d", groupid, click_current_cpu_id());
        Packet* next = _tables[groupid].queue;
        while (next != 0) {
            Packet* r = next;
            next = r->next();
            process(groupid, r, b);
        }
        _tables[groupid].queue = 0;
    }
}

void FlowIPManager::init_assignment(unsigned* table, int sz) {
    click_chatter("Initializing flow table assignment with %d buckets", sz);
    assert(_tables && _groups);
    if (sz != _groups) {
        click_chatter("ERROR: Initializing %p{element} with %d buckets, but configured with %d buckets", this, sz, _groups);
        abort();
    }
    for (int i = 0; i < sz; i++) {
        assert(table[i] < 64);
        _tables[i].owner = table[i];
    }
}

void FlowIPManager::do_migrate(CoreInfo &core) {
    click_chatter("Core %d is releasing its migrated buckets (%d packets reached)", click_current_cpu_id(), core.watch);
    core.lock.acquire(); //No race condition here, only this core can release the pending flag
    for (int i = 0;i < core.moves.size(); i++) {
        if (unlikely(_verbose > 1)) {
            click_chatter("Table %d now owned by %d", core.moves[i].first, core.moves[i].second);
        }
        assert(core.moves[i].second < 64 && core.moves[i].second >= 0 );
        _tables[core.moves[i].first].owner = core.moves[i].second;
    }
    core.moves.clear();
    core.pending = 0;
    core.watch = 0;
    core.lock.release();
}

void FlowIPManager::push_batch(int, PacketBatch* batch) {
    CoreInfo& core = *_cores;
    BatchBuilder b;
    int count = batch->count();

    core.count += count;

    if (core.pending) {
        if (core.watch > 0 && core.watch < core.count) { //If core watch is 0, post_migration has not run yet and the watch is invalid

            do_migrate(core);
        }
    }

    int last = -1;
    FOR_EACH_PACKET_SAFE(batch, p) {
        int groupid = AGGREGATE_ANNO(p) % _groups;
            if (_tables[groupid].owner != click_current_cpu_id()) {
                //First : the owner is for the old CPU core
                //We enter here as the migration has not been done
                //Then, core "from" release the table and change owner
                //--> now we do not enter here anymore
                //
                //The table is not for us (it was moved in as only a core that actually saw a change can change the owner)

                if (core.pending) //a new assignment has been seen ! Move our CPU if it wasn't done
                    do_migrate(core);

                if (unlikely(_verbose > 2))
                    click_chatter("Packet of group %d pushed on core %d while %d still holds the lock", groupid, click_current_cpu_id(), _tables[groupid].owner);

                if (_tables[groupid].queue) {
                    _tables[groupid].queue->prev()->set_next(p);
                } else {
                    _tables[groupid].queue = p;
                }
                _tables[groupid].queue->set_prev(p);
                p->set_next(0);

            } else {
                if (!core.pending) //While waiting for the new balancing to be written, we may receive packets of a table that we have to give to someone else, that is CURRENTLY enqueing
                    flush_queue(groupid, b); //This is a table for us. If there is a trailing queue (set by this same core before, we flush it)

                process(groupid, p, b);
            }
            last = groupid;
    }

    batch = b.finish();
    if (batch) {
        fcb_stack->acquire(batch->count());
        output_push_batch(0, batch);
    }

    //fcb_table = &_table;
    //fcb_stack =
}

enum {h_count};
String FlowIPManager::read_handler(Element* e, void* thunk) {
   FlowIPManager* fc = static_cast<FlowIPManager*>(e);

   switch ((intptr_t)thunk) {
       case h_count: {
        StringAccum acc;
        for (int i = 0; i < fc->_groups; i++) {
            acc << rte_hash_count(fc->_tables[i].hash) << "\n";
        }
        return acc.take_string(); }
       default:
        return "<error>";

   };
}

void FlowIPManager::add_handlers() {
    add_read_handler("count", read_handler, h_count);
}





CLICK_ENDDECLS

EXPORT_ELEMENT(FlowIPManager)
ELEMENT_MT_SAFE(FlowIPManager)