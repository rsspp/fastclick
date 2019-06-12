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




CLICK_DECLS
class DPDKDevice;
struct rte_hash;

struct BatchBuilder {
	BatchBuilder() : count(0), first(0), last(-1) {

	}
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


class FlowIPManager: public BatchElement {
public:


    FlowIPManager() CLICK_COLD;

	~FlowIPManager() CLICK_COLD;

    const char *class_name() const		{ return "FlowIPManager"; }
    const char *port_count() const		{ return "1/1"; }

    const char *processing() const		{ return PUSH; }
    int configure_phase() const     { return CONFIGURE_PHASE_PRIVILEGED + 1; }

    int configure(Vector<String> &, ErrorHandler *) CLICK_COLD;
    int initialize(ErrorHandler *errh) CLICK_COLD;
    void cleanup(CleanupStage stage) CLICK_COLD;

    //First : group id, second : destination cpu
    void migrate(DPDKDevice* dev, int from, Vector<Pair<int,int>> gids);


    void push_batch(int, PacketBatch* batch);

    void init_assignment(Vector<unsigned> table);

private:
    inline void flush_queue(int groupid, BatchBuilder &b);
    inline void process(int groupid, Packet* p, BatchBuilder& b);
    struct CoreInfo {
    	CoreInfo() : watch(0), count(0), lock() {
    	};
    	uint64_t count;
    	uint64_t watch;
    	Spinlock lock;
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
    int _reserve;
    int _table_size;
    int _flow_state_size_full;
    int _verbose;


    gtable* _tables;
    per_thread<CoreInfo> _cores;

};



CLICK_ENDDECLS
#endif
