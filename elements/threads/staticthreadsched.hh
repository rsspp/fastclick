// -*- c-basic-offset: 4 -*-
#ifndef STATICTHREADSCHED_HH
#define STATICTHREADSCHED_HH
#include <click/element.hh>
#include <click/standard/threadsched.hh>
CLICK_DECLS

/*
 * =c
 * StaticThreadSched(ELEMENT THREAD, ...)
 * =s threads
 * specifies element and thread scheduling parameters
 * =d
 * Statically binds elements to threads. If more than one StaticThreadSched
 * is specified, they will all run. The one that runs later may override an
 * earlier run.
 * Thread values may be negative, in which case the value is substracted from
 * the total number of threads used by Click. E.g. if click is launched with
 * 8 threads (-j 8), and -1 is passed, using -1 or 7 refer to the same thread.
 * =a
 * ThreadMonitor, BalancedThreadSched
 */

class StaticThreadSched : public Element, public ThreadSched { public:

    StaticThreadSched() CLICK_COLD;
    ~StaticThreadSched() CLICK_COLD;

    const char *class_name() const	{ return "StaticThreadSched"; }

    int configure(Vector<String> &, ErrorHandler *) CLICK_COLD;

    int initial_home_thread_id(const Element *e);

    Bitvector assigned_thread();

  private:
    Vector<int> _thread_preferences;
    ThreadSched *_next_thread_sched;

    bool set_preference(int eindex, int preference);
};

CLICK_ENDDECLS
#endif
