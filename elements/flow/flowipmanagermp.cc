/*
 * FlowIPManagerMP.{cc,hh}
 */

#include <click/config.h>
#include <click/glue.hh>
#include <click/args.hh>
#include <click/ipflowid.hh>
#include <click/routervisitor.hh>
#include <click/error.hh>
#include "flowipmanagermp.hh"
#include <rte_hash.h>
#include <rte_hash_crc.h>
#include <rte_ethdev.h>

CLICK_DECLS

FlowIPManagerMP::FlowIPManagerMP() : _verbose(1) {

}

FlowIPManagerMP::~FlowIPManagerMP() {

}

int
FlowIPManagerMP::configure(Vector<String> &conf, ErrorHandler *errh)
{

    if (Args(conf, this, errh)
    		.read_or_set_p("CAPACITY", _table_size, 65536)
            .read_or_set("RESERVE",_reserve, 0)
            .complete() < 0)
        return -1;

    find_children(_verbose);

    router()->get_root_init_future()->postOnce(&_fcb_builded_init_future);
    _fcb_builded_init_future.post(this);

    return 0;
}

inline uint32_t
ipv4_hash_crc(const void *data,  uint32_t data_len,
                uint32_t init_val)
{
	(void)data_len;
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

int FlowIPManagerMP::solve_initialize(ErrorHandler *errh) {
	struct rte_hash_parameters hash_params = {0};
	char buf[32];
	hash_params.name = buf;
	hash_params.entries = _table_size;
	hash_params.key_len = sizeof(IPFlow5ID);
	hash_params.hash_func = ipv4_hash_crc;
	hash_params.hash_func_init_val = 0;
	hash_params.extra_flag = RTE_HASH_EXTRA_FLAGS_MULTI_WRITER_ADD | RTE_HASH_EXTRA_FLAGS_RW_CONCURRENCY; //| RTE_HASH_EXTRA_FLAGS_RW_CONCURRENCY_LF

	_flow_state_size_full = sizeof(FlowControlBlock) + _reserve;

	sprintf(buf, "FlowIPManagerMP");
	hash = rte_hash_create(&hash_params);
	if (!hash)
		return errh->error("Could not init flow table !");

	fcbs =  (FlowControlBlock*)CLICK_ALIGNED_ALLOC(_flow_state_size_full * _table_size);
    bzero(fcbs,_flow_state_size_full * _table_size);
	CLICK_ASSERT_ALIGNED(fcbs);
	if (!fcbs)
		return errh->error("Could not init data table !");



    return 0;
}

void FlowIPManagerMP::cleanup(CleanupStage stage) {

}


void FlowIPManagerMP::pre_migrate(DPDKDevice* dev, int from, Vector<Pair<int,int>> gids) {

}


void FlowIPManagerMP::post_migrate(DPDKDevice* dev, int from) {

}

void FlowIPManagerMP::process(Packet* p, BatchBuilder& b) {
	IPFlow5ID fid = IPFlow5ID(p);
	bool first = false;
	rte_hash*& table = hash;

	int ret = rte_hash_lookup(table, &fid);

	if (ret < 0) { //new flow
		ret = rte_hash_add_key(table, &fid);
		if (ret < 0) {
			click_chatter("Cannot add key (have %d items)!", rte_hash_count(table));
			return;
		}
		first = true;
	}

	if (b.last == ret) {
		b.append(p);
	} else {
		PacketBatch* batch;
		batch = b.finish();
		if (batch)
			output_push_batch(0, batch);
		fcb_stack = (FlowControlBlock*)((unsigned char*)fcbs + _flow_state_size_full * ret);
		b.init();
        b.append(p);
	}
}


void FlowIPManagerMP::init_assignment(Vector<unsigned> table) {

}



void FlowIPManagerMP::push_batch(int, PacketBatch* batch) {
	BatchBuilder b;

	FOR_EACH_PACKET_SAFE(batch, p) {
		process(p, b);
	}

	batch = b.finish();
	if (batch)
		output_push_batch(0, batch);


    //fcb_table = &_table;
	//fcb_stack =

}

/*
enum {h_leaves_count, h_active_count, h_print, h_timeout_count};
String FlowIPManagerMP::read_handler(Element* e, void* thunk) {
    FlowIPManagerMP* fc = static_cast<FlowIPManagerMP*>(e);

    fcb_table = &fc->_table;
    switch ((intptr_t)thunk) {

};

void FlowIPManagerMP::add_handlers() {


}*/





CLICK_ENDDECLS

EXPORT_ELEMENT(FlowIPManagerMP)
ELEMENT_MT_SAFE(FlowIPManagerMP)
