/*
 * RandLoad.{cc,hh}
 */

#include "randload.hh"

#include <click/config.h>
#include <click/glue.hh>
#include <click/args.hh>

CLICK_DECLS

RandLoad::RandLoad() {

};

RandLoad::~RandLoad() {

}

int
RandLoad::configure(Vector<String> &conf, ErrorHandler *errh)
{
    if (Args(conf, this, errh)
               .read_or_set("MIN", _min, 1)
               .read_or_set("MAX", _max, 100)
               .complete() < 0)
        return -1;

    return 0;
}


int RandLoad::initialize(ErrorHandler *errh) {
    return 0;
}

inline void
RandLoad::load(Packet* p)  {
    int r;
    int w =  _min + ((*_gens)() / (UINT_MAX / (_max - _min) ));  //click_random(_min, _max);
    for (int i = 0; i < w - 1; i ++) {
        r = (*_gens)();
    }
};

void RandLoad::push(int port, Packet* p) {
    load(p);
    output_push(0, p);
}

#if HAVE_BATCH
void RandLoad::push_batch(int port, PacketBatch* batch) {
    FOR_EACH_PACKET(batch, p)
        load(p);

    output_push_batch(0, batch);
}
#endif


CLICK_ENDDECLS

EXPORT_ELEMENT(RandLoad)
ELEMENT_MT_SAFE(RandLoad)
