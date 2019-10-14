#ifndef CLICK_FlowIPManager_HH
#define CLICK_FlowIPManager_HH
#include <click/config.h>
#include <click/string.hh>
#include <click/timer.hh>
#include <click/vector.hh>
#include <click/multithread.hh>
#include <click/batchelement.hh>
#include <click/pair.hh>
#include <click/flow/common.hh>
#include <click/flow/flowelement.hh>



CLICK_DECLS
class DPDKDevice;
struct rte_hash;

struct BatchBuilder {
	BatchBuilder() : first(0), last(-1), count(0) {

	};

	Packet* first;
	Packet* tail;
	int count;
	int last;

	inline void init() {
		count = 0;
		first = 0;
	}
	inline PacketBatch* finish() {
		if (!first)
			return 0;
		return PacketBatch::make_from_simple_list(first,tail,count);
	}

	inline void append(Packet* p) {
		count++;
		if (first) {
			tail->set_next(p);
			tail = p;
		} else {
			first = p;
			tail = p;
		}
	}

};

/**
 * FCB packet classifier - cuckoo per-thread
 *
 * Initialize the FCB stack for every packets passing by.
 * The classification is done using one cuckoo hashtable per threads. Hence,
 * when some flows are migratedthe do_migrate function must be called.
 *
 * This element does not find automatically the FCB layout for FlowElement,
 * neither set the offsets for placement in the FCB automatically. Look at
 * the middleclick branch for alternatives.
 */
class FlowIPManager: public VirtualFlowManager, public Router::InitFuture {
public:

    FlowIPManager() CLICK_COLD;
	~FlowIPManager() CLICK_COLD;

    const char *class_name() const		{ return "FlowIPManager"; }
    const char *port_count() const		{ return "1/1"; }

    const char *processing() const		{ return PUSH; }
    int configure_phase() const     { return CONFIGURE_PHASE_PRIVILEGED + 1; }

    int configure(Vector<String> &, ErrorHandler *) override CLICK_COLD;
    int solve_initialize(ErrorHandler *errh) override CLICK_COLD;
    void cleanup(CleanupStage stage) CLICK_COLD;

    //First : group id, second : destination cpu
    void pre_migrate(DPDKDevice* dev, int from, Vector<Pair<int,int>> gids);
    void post_migrate(DPDKDevice* dev, int from);


    void push_batch(int, PacketBatch* batch);

    void init_assignment(Vector<unsigned> table);

private:

    struct CoreInfo {
	CoreInfo() : watch(0), count(0), lock(), pending(false) {
    	};
    	uint64_t count;
    	uint64_t watch;
    	Spinlock lock;
	bool pending;
    	Vector<Pair<int,int> > moves; //First : group id, second : destination cpu
    } CLICK_ALIGNED(CLICK_CACHE_LINE_SIZE);

    struct gtable {
    	gtable() : queue(0) {

    	}
    	volatile int owner;
    	Packet* queue;
    	rte_hash* hash;
        FlowControlBlock *fcbs;
    } CLICK_ALIGNED(CLICK_CACHE_LINE_SIZE);


    int _groups;
    int _table_size;
    int _flow_state_size_full;
    int _verbose;
    int _def_thread;

    gtable* _tables;
    per_thread<CoreInfo> _cores;

    void do_migrate(CoreInfo &core);
    inline void flush_queue(int groupid, BatchBuilder &b);
    inline void process(int groupid, Packet* p, BatchBuilder& b);

};



CLICK_ENDDECLS
#endif
