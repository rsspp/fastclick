#ifndef CLICK_FlowIPManagerMPMP_HH
#define CLICK_FlowIPManagerMPMP_HH
#include <click/config.h>
#include <click/string.hh>
#include <click/timer.hh>
#include <click/vector.hh>
#include <click/multithread.hh>
#include <click/batchelement.hh>
#include <click/pair.hh>
#include <click/flow/common.hh>


#include "flowipmanager.hh"




CLICK_DECLS
class DPDKDevice;
struct rte_hash;


class FlowIPManagerMP: public BatchElement {
public:


    FlowIPManagerMP() CLICK_COLD;

	~FlowIPManagerMP() CLICK_COLD;

    const char *class_name() const		{ return "FlowIPManagerMP"; }
    const char *port_count() const		{ return "1/1"; }

    const char *processing() const		{ return PUSH; }
    int configure_phase() const     { return CONFIGURE_PHASE_PRIVILEGED + 1; }

    int configure(Vector<String> &, ErrorHandler *) CLICK_COLD;
    int initialize(ErrorHandler *errh) CLICK_COLD;
    void cleanup(CleanupStage stage) CLICK_COLD;

    //First : group id, second : destination cpu
    void pre_migrate(DPDKDevice* dev, int from, Vector<Pair<int,int>> gids);
    void post_migrate(DPDKDevice* dev, int from);


    void push_batch(int, PacketBatch* batch);

    void init_assignment(Vector<unsigned> table);

private:


    	volatile int owner;
    	Packet* queue;
    	rte_hash* hash;
    	FlowControlBlock *fcbs;


    int _reserve;
    int _table_size;
    int _flow_state_size_full;
    int _verbose;



    inline void process(Packet* p, BatchBuilder& b);

};



CLICK_ENDDECLS
#endif
