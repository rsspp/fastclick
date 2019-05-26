/*
 * rssverifier.{cc,hh} -- Verify RSS mapping
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the "Software"),
 * to deal in the Software without restriction, subject to the conditions
 * listed in the Click LICENSE file. These conditions include: you must
 * preserve this copyright notice, and you cannot mention the copyright
 * holders in advertising related to the Software without their permission.
 * The Software is provided WITHOUT ANY WARRANTY, EXPRESS OR IMPLIED. This
 * notice is a summary of the Click LICENSE file; the license in that file is
 * legally binding.
 */

#include <click/config.h>
#include "rssverifier.hh"
#include <click/handlercall.hh>
#include <click/args.hh>
#include <click/error.hh>
#include <click/packet_anno.hh>
#include <click/integers.hh>	// for first_bit_set
#include <click/router.hh>
CLICK_DECLS



RSSVerifier::RSSVerifier()
{
	_bad = 0;
	_count = 0;
}


RSSVerifier::~RSSVerifier()
{
}


int
RSSVerifier::configure(Vector<String> &conf, ErrorHandler *errh)
{
    uint32_t mask = 511;
    int n = 1;

    if (Args(conf, this, errh)
    .read_mp("MASK", mask)
	.read_mp("CPU", n)
	.complete() < 0)
	return -1;

    _mask = mask;

    _table.resize(mask + 1);
    for (int i = 0; i < mask + 1; i++) {
    	_table[i] = i % n;
    }


    return 0;
}

int
RSSVerifier::initialize(ErrorHandler *errh)
{
    _active = true;
    return 0;
}

void
RSSVerifier::cleanup(CleanupStage)
{

}

inline bool
RSSVerifier::update(Packet *p)
{
    if (!_active)
	return false;

    uint32_t agg = AGGREGATE_ANNO(p) & _mask;
_count++;
    if (_table[agg] != click_current_cpu_id()) {
    	_bad++;
    	click_chatter("ERROR : AGG %u, idx %u does not map to %d but %d", AGGREGATE_ANNO(p),agg,_table[agg],click_current_cpu_id());

    }

}

void
RSSVerifier::push(int port, Packet *p)
{
    port = !update(p);
    output(noutputs() == 1 ? 0 : port).push(p);
}

Packet *
RSSVerifier::pull(int)
{
    Packet *p = input(0).pull();
    if (p && _active)
	update(p);
    return p;
}

#if HAVE_BATCH
void
RSSVerifier::push_batch(int port, PacketBatch *batch)
{
    auto fnt = [this,port](Packet*p){return !update(p);};
    CLASSIFY_EACH_PACKET(2,fnt,batch,[this](int port, PacketBatch* batch){ output(0).push_batch(batch);});
}

PacketBatch *
RSSVerifier::pull_batch(int port,unsigned max)
{
    PacketBatch *batch = input(0).pull_batch(max);
    if (batch && _active) {
        FOR_EACH_PACKET(batch,p) {
		update(p);
        }
    }
    return batch;
}
#endif


enum {
    AC_ACTIVE, AC_STOP,
};

String
RSSVerifier::read_handler(Element *e, void *thunk)
{
    RSSVerifier *ac = static_cast<RSSVerifier *>(e);
    switch ((intptr_t)thunk) {
      default:
	return "<error>";
    }
}

int
RSSVerifier::write_handler(const String &data, Element *e, void *thunk, ErrorHandler *errh)
{
    RSSVerifier *ac = static_cast<RSSVerifier *>(e);
    String s = cp_uncomment(data);
    switch ((intptr_t)thunk) {
      case AC_ACTIVE: {
	  bool val;
	  if (!BoolArg().parse(s, val))
	      return errh->error("type mismatch");
	  ac->_active = val;
	  return 0;
      }
      case AC_STOP:
	ac->_active = false;
	ac->router()->please_stop_driver();
	return 0;
      default:
	return errh->error("internal error");
    }
}


void
RSSVerifier::add_handlers()
{
    add_data_handlers("active", Handler::f_read | Handler::f_checkbox, &_active);
    add_write_handler("active", write_handler, AC_ACTIVE);
    add_write_handler("stop", write_handler, AC_STOP, Handler::f_button);
}


ELEMENT_REQUIRES(userlevel int64)
EXPORT_ELEMENT(RSSVerifier)
ELEMENT_MT_SAFE(RSSVerifier)
CLICK_ENDDECLS
