#ifndef CLICK_RandLoad_HH
#define CLICK_RandLoad_HH
#include <click/config.h>
#include <click/glue.hh>
#include <click/vector.hh>
#include <click/batchelement.hh>
#include <random>

CLICK_DECLS

class RandLoad : public BatchElement {

public:

    RandLoad() CLICK_COLD;
    ~RandLoad() CLICK_COLD;

    const char *class_name() const		{ return "RandLoad"; }
    const char *port_count() const		{ return "1/1"; }
    const char *processing() const		{ return PUSH; }

    int configure(Vector<String> &, ErrorHandler *) CLICK_COLD;
    int initialize(ErrorHandler *errh);

    inline void load(Packet* p);
    void push(int, Packet *);
#if HAVE_BATCH
    void push_batch(int, PacketBatch *);
#endif
private:
    int _min;
    int _max;
    per_thread<std::mt19937> _gens;
};





CLICK_ENDDECLS
#endif
