#ifndef CLICK_FlowRandLoad_HH
#define CLICK_FlowRandLoad_HH
#include <click/config.h>
#include <click/multithread.hh>
#include <click/hashtablemp.hh>
#include <click/glue.hh>
#include <click/vector.hh>
#include <click/flow/flowelement.hh>
#include <random>

CLICK_DECLS

struct RandLoadState {
	RandLoadState() : w(0) {

	}
	int w;
};

class FlowRandLoad : public FlowSpaceElement<RandLoadState> {

public:

    FlowRandLoad() CLICK_COLD;
    ~FlowRandLoad() CLICK_COLD;

    const char *class_name() const		{ return "FlowRandLoad"; }
    const char *port_count() const		{ return "1/1"; }
    const char *processing() const		{ return PUSH; }

    int configure(Vector<String> &, ErrorHandler *) override CLICK_COLD;
    int initialize(ErrorHandler *errh) override CLICK_COLD;

    void push_batch(int, RandLoadState*, PacketBatch *);

private:
    int _min;
    int _max;
    per_thread<std::mt19937> _gens;
};





CLICK_ENDDECLS
#endif
