/*
 * FlowIPManager.{cc,hh}
 */

#include <click/config.h>
#include <click/glue.hh>
#include <click/args.hh>
#include <click/ipflowid.hh>
#include <click/routervisitor.hh>
#include "flowipmanager.hh"
#include <rte_hash.h>
#include <rte_hash_crc.h>
#include <click/dpdkdevice.hh>
#include "../userlevel/fromdpdkdevice.hh"
#include <rte_ethdev.h>

CLICK_DECLS

FlowIPManager::FlowIPManager() : _verbose(1) {

}

FlowIPManager::~FlowIPManager() {

}

int
FlowIPManager::configure(Vector<String> &conf, ErrorHandler *errh)
{

    if (Args(conf, this, errh)
    		.read_or_set_p("GROUPS", _groups, 512)
    		.read_or_set_p("CAPACITY", _table_size, 65536)
            .read_or_set("RESERVE",_reserve, 0)
            .complete() < 0)
        return -1;


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

int FlowIPManager::initialize(ErrorHandler *errh) {
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
	    CLICK_ASSERT_ALIGNED(_tables[i].fcbs);
	    if (!_tables[i].fcbs)
	    	return errh->error("Could not init data table %d!", i);
	}



    return 0;
}

void FlowIPManager::cleanup(CleanupStage stage) {

}


void FlowIPManager::migrate(DPDKDevice* dev, int from, Vector<Pair<int,int>> gids) {
	int port_id = dev->port_id;
	int v = rte_eth_rx_queue_count(port_id, from);
//	struct rte_eth_stats stats;
 //   rte_eth_stats_get(port_id, &stats);
	CoreInfo &coref = _cores.get_value_for_thread(from);
    uint64_t w = coref.count + v;
    if (coref.watch < w) {
    	coref.lock.acquire();
    	for (int i = 0;i < gids.size(); i++) {
    		coref.moves.push_back(gids[i]);
    	}
    	coref.lock.release();
    	coref.watch = w;
    }
}

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

	if (unlikely(_verbose > 1))
		click_chatter("Packet of flow group %d, id %d", groupid, ret);

	if (b.last == ret) {
		b.append(p);
	} else {
		PacketBatch* batch;
		batch = b.finish();
		if (batch)
			output_push_batch(0, batch);
		fcb_stack = (FlowControlBlock*)((unsigned char*)t.fcbs + _flow_state_size_full * ret);
		b.init();
        b.append(p);
	}
}

inline void FlowIPManager::flush_queue(int groupid, BatchBuilder &b) {
	if (_tables[groupid].queue) {
		Packet* next = _tables[groupid].queue;
		while (next != 0) {
			Packet* r = next;
			next = r->next();
			process(groupid, r, b);
		}
		_tables[groupid].queue = 0;
	}
}

void FlowIPManager::init_assignment(Vector<unsigned> table) {
	click_chatter("Initializing flow table assignment");
	for (int i = 0; i < table.size(); i++) {
		_tables[i].owner = table[i];
	}

}

void FlowIPManager::push_batch(int, PacketBatch* batch) {
	CoreInfo& core = *_cores;
	BatchBuilder b;
	int count = batch->count();
	int last = -1;
	FOR_EACH_PACKET(batch, p) {
		int groupid = AGGREGATE_ANNO(p) % _groups;
		if (groupid != last) {
			if (_tables[groupid].owner != click_current_cpu_id()) {
				if (unlikely(_verbose > 0))
					click_chatter("Packet of group %d pushed on core %d while %d still holds the lock", groupid, click_current_cpu_id(), _tables[groupid].owner);
				if (_tables[groupid].queue)
					_tables[groupid].queue->prev()->set_next(p);
				else
					_tables[groupid].queue = p;
				_tables[groupid].queue->set_prev(p);
			} else {
				flush_queue(groupid, b);
			}
		}
		last = groupid;
		process(groupid, p, b);
	}

	batch = b.finish();
	if (batch)
		output_push_batch(0, batch);

	core.count += count;

	if (core.watch > 0) {
		if (core.watch < core.count) {
	    	core.lock.acquire();
	    	for (int i = 0;i < core.moves.size(); i++) {
	    		_tables[core.moves[i].first].owner = core.moves[i].second;
	    	}
	    	core.moves.clear();
	    	core.lock.release();
			core.watch = 0;
		}
	}

    //fcb_table = &_table;
	//fcb_stack =

}

/*
enum {h_leaves_count, h_active_count, h_print, h_timeout_count};
String FlowIPManager::read_handler(Element* e, void* thunk) {
    FlowIPManager* fc = static_cast<FlowIPManager*>(e);

    fcb_table = &fc->_table;
    switch ((intptr_t)thunk) {

};

void FlowIPManager::add_handlers() {


}*/





CLICK_ENDDECLS

EXPORT_ELEMENT(FlowIPManager)
ELEMENT_MT_SAFE(FlowIPManager)
