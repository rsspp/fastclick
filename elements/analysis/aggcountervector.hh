#ifndef CLICK_AGGCOUNTERVECTOR_HH
#define CLICK_AGGCOUNTERVECTOR_HH
#include <click/batchelement.hh>
#include <click/multithread.hh>
CLICK_DECLS
class HandlerCall;

/*
=c

AggregateCounter([I<KEYWORDS>])

=s aggregates

counts packets per aggregate annotation

=d

 */


class AggregateCounterVector : public BatchElement { public:


	struct Node {
	    uint64_t count;
	    uint64_t variance;
	    uint32_t epoch;
	    uint16_t flows; //Max 42*8
	    uint8_t map[42];

	    void add_flow(uint32_t agg) {
		uint16_t n = ((agg >> 24) ^ (agg >> 15)) % (42*8); //Number between 0 and 8*42
		//Set the bit n in map
		if (map[n / 8] & (1 << n % 8)) {

		} else {
			map[n / 8] |= (1 << n % 8);
			flows++;
		}

	    }
	} CLICK_CACHE_ALIGN;


    AggregateCounterVector() CLICK_COLD;
    ~AggregateCounterVector() CLICK_COLD;

    const char *class_name() const  { return "AggregateCounterVector"; }
    const char *port_count() const	{ return "1/1"; }

    int configure(Vector<String> &, ErrorHandler *) CLICK_COLD;
    int initialize(ErrorHandler *) CLICK_COLD;
    void cleanup(CleanupStage) CLICK_COLD;

    static String read_handler(Element *e, void *thunk);
    static int write_handler(const String &data, Element *e, void *thunk, ErrorHandler *errh);
    void add_handlers() CLICK_COLD;

    inline Node& find_node(uint32_t agg);
    inline bool update(Packet *);
    inline bool update_batch(PacketBatch *);
    void push(int, Packet *);
    Packet *pull(int);
#if HAVE_BATCH
    void push_batch(int, PacketBatch *) override;
    PacketBatch *pull_batch(int, unsigned) override;
#endif


    inline void advance_epoch() {
	_epoch++;
    }
  private:

    bool _bytes : 1;
    bool _ip_bytes : 1;
    bool _use_packet_count : 1;
    bool _use_extra_length : 1;
    bool _active;

    uint32_t _mask;

    uint32_t _epoch;

    Vector<Node> _nodes;



};

/**
 * @pre a is masked
 */
inline AggregateCounterVector::Node&
AggregateCounterVector::find_node(uint32_t a)
{
    Node& n = _nodes.unchecked_at(a);
    if (n.epoch != _epoch) {
	n.variance = n.variance / 3 + (n.count * 2 / 3);
	n.count = 0;
	n.epoch = _epoch;
	n.flows = 0;
	bzero(n.map,sizeof(n.map));
    }
    return n;
}

CLICK_ENDDECLS
#endif
