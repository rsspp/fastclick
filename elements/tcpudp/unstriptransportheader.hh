#ifndef CLICK_UnstripTransportHeader_HH
#define CLICK_UnstripTransportHeader_HH
#include <click/batchelement.hh>
CLICK_DECLS

/*
 * =c
 * UnstripTransportHeader()
 * =s ip
 * restores outermost transport header
 * =d
 *
 * Put outermost IP header back onto a stripped packet, based on the IP Header type,
 * and the IP header annotation from MarkIPHeader or CheckIPHeader. If IP header already on,
 * forwards packet unmodified.
 *
 * =a CheckIPHeader, MarkIPHeader, StripIPHeader */

class UnstripTransportHeader : public BatchElement {

public:
  UnstripTransportHeader() CLICK_COLD;
  ~UnstripTransportHeader() CLICK_COLD;

  const char *class_name() const		{ return "UnstripTransportHeader"; }
  const char *port_count() const		{ return PORTS_1_1; }

  Packet      *simple_action      (Packet      *p);
#if HAVE_BATCH
  PacketBatch *simple_action_batch(PacketBatch *batch);
#endif
};

CLICK_ENDDECLS
#endif
