/*
 * ipratemon.{cc,hh} -- counts packets clustered by src/dst addr.
 * Thomer M. Gil
 * Benjie Chen, Eddie Kohler (minor changes)
 *
 * Copyright (c) 1999-2000 Massachusetts Institute of Technology.
 *
 * This software is being provided by the copyright holders under the GNU
 * General Public License, either version 2 or, at your discretion, any later
 * version. For more information, see the `COPYRIGHT' file in the source
 * distribution.
 */

#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
#include "ipratemon.hh"
#include "confparse.hh"
#include "straccum.hh"
#include "click_ip.h"
#include "error.hh"
#include "glue.hh"

IPRateMonitor::IPRateMonitor()
  : _pb(COUNT_PACKETS), _offset(0), _thresh(1), _base(NULL)
{
}

IPRateMonitor::~IPRateMonitor()
{
}

IPRateMonitor *
IPRateMonitor::clone() const
{
  return new IPRateMonitor;
}

void
IPRateMonitor::notify_ninputs(int n)
{
  set_ninputs(n == 1 ? 1 : 2);
  set_noutputs(n == 1 ? 1 : 2);
}

int
IPRateMonitor::configure(const String &conf, ErrorHandler *errh)
{
#if IPVERSION == 4
  Vector<String> args;
  cp_argvec(conf, args);

  // Enough args?
  if(args.size() != 3)
    return errh->error("too few arguments.");

  // PACKETS/BYTES
  if(args[0] == "PACKETS")
    _pb = COUNT_PACKETS;
  else if(args[0] == "BYTES")
    _pb = COUNT_BYTES;
  else
    return errh->error("second argument should be \"PACKETS\" or \"BYTES\".");

  // OFFSET
  if(!cp_integer(args[1], _offset) || _offset < 0)
    return errh->error
      ("offset should be a non-negative integer.");

  // THRESH
  if(!cp_integer(args[2], _thresh) || _thresh < 0)
    return errh->error
      ("thresh should be non-negative integer.");
 
  set_resettime();

  // Make _base
  _base = new Stats();
  if(!_base)
    return errh->error("cannot allocate data structure.");

  return 0;
#else
  click_chatter("IPRateMonitor doesn't know how to handle non-IPv4!");
  return -1;
#endif
}

void
IPRateMonitor::push(int port, Packet *p)
{
  update_rates(p, port == 0);
  output(port).push(p);
}

Packet *
IPRateMonitor::pull(int port)
{
  Packet *p = input(port).pull();
  if (p)
    update_rates(p, port == 0);
  return p;
}


//
// Recursively destroys tables.
//

IPRateMonitor::Stats::Stats()
{
  for (int i = 0; i < MAX_COUNTERS; i++) {
    counter[i].fwd_rate.initialize();
    counter[i].rev_rate.initialize();
    counter[i].next_level = 0;
  }
}

IPRateMonitor::Stats::Stats(const MyEWMA &fwd, const MyEWMA &rev)
{
  for (int i = 0; i < MAX_COUNTERS; i++) {
    counter[i].fwd_rate = fwd;
    counter[i].rev_rate = rev;
    counter[i].next_level = 0;
  }
}

IPRateMonitor::Stats::~Stats()
{
  for (int i = 0; i < MAX_COUNTERS; i++)
    delete counter[i].next_level;
}

void
IPRateMonitor::Stats::clear()
{
  for (int i = 0; i < MAX_COUNTERS; i++) {
    delete counter[i].next_level;
    counter[i].next_level = 0;
    counter[i].rev_rate.initialize();
    counter[i].fwd_rate.initialize();
  }
}

//
// Prints out nice data.
//
String
IPRateMonitor::print(Stats *s, String ip = "")
{
  int jiffs = MyEWMA::now();
  String ret = "";
  for(int i = 0; i < MAX_COUNTERS; i++) {
    Counter &c = s->counter[i];
    if (c.rev_rate.average() > 0 || c.fwd_rate.average() > 0) {
      String this_ip;
      if (ip)
        this_ip = ip + "." + String(i);
      else
        this_ip = String(i);
      ret += this_ip;

      c.fwd_rate.update(jiffs, 0);
      c.rev_rate.update(jiffs, 0);
      ret += "\t"; 
      ret += cp_unparse_real(c.fwd_rate.average()*CLICK_HZ, c.fwd_rate.scale);
      ret += "\t"; 
      ret += cp_unparse_real(c.rev_rate.average()*CLICK_HZ, c.rev_rate.scale);
      
      ret += "\n";
      if (c.next_level) 
        ret += print(c.next_level, "\t" + this_ip);
    }
  }
  return ret;
}


inline void
IPRateMonitor::set_resettime()
{
  _resettime = MyEWMA::now();
}

String
IPRateMonitor::look_read_handler(Element *e, void *)
{
  IPRateMonitor *me = (IPRateMonitor*) e;

  String ret = String(MyEWMA::now() - me->_resettime) + "\n";
  return ret + me->print(me->_base);
}

String
IPRateMonitor::thresh_read_handler(Element *e, void *)
{
  IPRateMonitor *me = (IPRateMonitor *) e;
  return String(me->_thresh);
}

int
IPRateMonitor::reset_write_handler
(const String &, Element *e, void *, ErrorHandler *)
{
  IPRateMonitor* me = (IPRateMonitor *) e;
  me->_base->clear();
  me->set_resettime();
  return 0;
}

void
IPRateMonitor::add_handlers()
{
  add_read_handler("thresh", thresh_read_handler, 0);
  add_read_handler("look", look_read_handler, 0);

  add_write_handler("reset", reset_write_handler, 0);
}

EXPORT_ELEMENT(IPRateMonitor)

// template instances
#include "ewma.cc"
