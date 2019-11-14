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
#include <click/flow/flowelement.hh>


#include "flowipmanager.hh"




CLICK_DECLS
class DPDKDevice;
struct rte_hash;

/**
 * FCB packet classifier - cuckoo shared-by-all-threads
 *
 * Initialize the FCB stack for every packets passing by.
 * The classification is done using a unique, thread safe cuckoo hash table.
 *
 * This element does not find automatically the FCB layout for FlowElement,
 * neither set the offsets for placement in the FCB automatically. Look at
 * the middleclick branch for alternatives.
 */
class FlowIPManagerMP: public VirtualFlowManager, Router::InitFuture, MigrationListener {
public:


    FlowIPManagerMP() CLICK_COLD;

	~FlowIPManagerMP() CLICK_COLD;

    const char *class_name() const		{ return "FlowIPManagerMP"; }
    const char *port_count() const		{ return "1/1"; }

    const char *processing() const		{ return PUSH; }
    int configure_phase() const     { return CONFIGURE_PHASE_PRIVILEGED + 1; }

    int configure(Vector<String> &, ErrorHandler *) override CLICK_COLD;
    int solve_initialize(ErrorHandler *errh) override CLICK_COLD;
    void cleanup(CleanupStage stage) override CLICK_COLD;

    //First : group id, second : destination cpu
    void pre_migrate(DPDKEthernetDevice* dev, int from, std::vector<std::pair<int,int>> gids) override;
    void post_migrate(DPDKEthernetDevice* dev, int from) override;

    void push_batch(int, PacketBatch* batch);

    void init_assignment(Vector<unsigned> table);

protected:


	volatile int owner;
	Packet* queue;
	rte_hash* hash;
	FlowControlBlock *fcbs;

    int _table_size;
    int _flow_state_size_full;
    int _verbose;

    inline void process(Packet* p, BatchBuilder& b);

};



CLICK_ENDDECLS
#endif
