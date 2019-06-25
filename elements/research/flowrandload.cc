/*
 * FlowRandLoad.{cc,hh}
 */

#include "flowrandload.hh"

#include <click/config.h>
#include <click/glue.hh>
#include <click/args.hh>
#include <clicknet/ip.h>
#include <clicknet/tcp.h>
#include <click/flow/flow.hh>

CLICK_DECLS

FlowRandLoad::FlowRandLoad() {

};

FlowRandLoad::~FlowRandLoad() {

}

int
FlowRandLoad::configure(Vector<String> &conf, ErrorHandler *errh)
{
    if (Args(conf, this, errh)
               .read_or_set("MIN", _min, 1)
               .read_or_set("MAX", _max, 100)
			   .read_mp("FCB_OFFSET",_flow_data_offset)
               .complete() < 0)
        return -1;

    return 0;
}


int FlowRandLoad::initialize(ErrorHandler *errh) {
    return 0;
}

void FlowRandLoad::push_batch(int port, RandLoadState* flowdata, PacketBatch* batch) {
	if (flowdata->w == 0) {
		flowdata->w =  _min + ((*_gens)() / (UINT_MAX / (_max - _min) ));  //click_random(_min, _max);
	}
	int r;
    auto fnt = [this,flowdata,&r](Packet* p) {
        for (int i = 0; i < flowdata->w; i ++) {
            r = (*_gens)();
        }
        return p;
    };
    EXECUTE_FOR_EACH_PACKET(fnt, batch);

    output_push_batch(0, batch);
}


CLICK_ENDDECLS

EXPORT_ELEMENT(FlowRandLoad)
ELEMENT_MT_SAFE(FlowRandLoad)
