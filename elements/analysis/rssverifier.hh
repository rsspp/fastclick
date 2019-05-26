#ifndef CLICK_RSSVERIFIER_HH
#define CLICK_RSSVERIFIER_HH
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


class RSSVerifier : public BatchElement { public:

    RSSVerifier() CLICK_COLD;
    ~RSSVerifier() CLICK_COLD;

    const char *class_name() const  { return "RSSVerifier"; }
    const char *port_count() const	{ return "1/1"; }

    int configure(Vector<String> &, ErrorHandler *) CLICK_COLD;
    int initialize(ErrorHandler *) CLICK_COLD;
    void cleanup(CleanupStage) CLICK_COLD;

    static String read_handler(Element *e, void *thunk);
    static int write_handler(const String &data, Element *e, void *thunk, ErrorHandler *errh);
    void add_handlers() CLICK_COLD;

    inline bool update(Packet *);
    void push(int, Packet *);
    Packet *pull(int);
#if HAVE_BATCH
    void push_batch(int, PacketBatch *) override;
    PacketBatch *pull_batch(int, unsigned) override;
#endif


    Vector<unsigned> _table;

  private:

    bool _active;

    uint32_t _mask;

    atomic_uint64_t _bad;
    atomic_uint64_t _count;



};


CLICK_ENDDECLS
#endif
