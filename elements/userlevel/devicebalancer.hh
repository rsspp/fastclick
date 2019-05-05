#ifndef CLICK_DEVICEBALANCER_HH
#define CLICK_DEVICEBALANCER_HH

#include <click/batchelement.hh>
#include "fromdpdkdevice.hh"
#include "../analysis/aggcountervector.hh"

CLICK_DECLS

enum target_method {
    TARGET_LOAD,
    TARGET_BALANCE
};

enum load_method {
    LOAD_CPU,
	LOAD_CYCLES,
    LOAD_QUEUE
};


class LoadTracker {
    public:

    protected:
    int load_tracker_initialize(ErrorHandler* errh);

    Vector<float> _past_load;
    Vector<click_jiffies_t> _last_movement;

};

class DeviceBalancer;

class BalanceMethod { public:
    BalanceMethod(DeviceBalancer* b) : balancer(b) {
    }

    virtual int configure(Vector<String> &, ErrorHandler *) CLICK_COLD;
    virtual void rebalance(Vector<Pair<int,float>> load) = 0;
    virtual int initialize(ErrorHandler *errh, int startwith);
protected:

    DeviceBalancer* balancer;
};

class BalanceMethodDPDK : public BalanceMethod { public:
    BalanceMethodDPDK(DeviceBalancer* b, Element* fd) : BalanceMethod(b) {
        _fd = dynamic_cast<FromDPDKDevice*>(fd);
        if (!_fd) {
            click_chatter("Not a DPDK device");
        }
    }
protected:
    FromDPDKDevice* _fd;
    friend class DeviceBalancer;
};

class MethodMetron : public BalanceMethodDPDK, public LoadTracker { public:

    MethodMetron(DeviceBalancer* b, Element* fd, String config) : BalanceMethodDPDK(b,fd) {
        _rules_file = config;
    }


    int initialize(ErrorHandler *errh, int startwith) override;

    void rebalance(Vector<Pair<int,float>> load) override;
private:
    String _rules_file;
};

class BalanceMethodRSS : public BalanceMethodDPDK { public:

    BalanceMethodRSS(DeviceBalancer* b, Element* fd) : BalanceMethodDPDK(b,fd) {
    }
    int initialize(ErrorHandler *errh, int startwith) override;

    void rebalance(Vector<Pair<int,float>> load) override;
    Vector<unsigned> _table;

    struct rte_eth_rss_conf _rss_conf;
    Vector<rte_flow*> _flows;

    bool update_reta_flow(bool validate = false);
    void update_reta(bool validate = false);
protected:
    bool _update_reta_flow;
};


class MethodRSSRR : public BalanceMethodRSS { public:

    MethodRSSRR(DeviceBalancer* b, Element* fd) : BalanceMethodRSS(b,fd) {
    }

    void rebalance(Vector<Pair<int,float>> load) override;
};


class MethodPianoRSS : public BalanceMethodRSS, public LoadTracker { public:

    MethodPianoRSS(DeviceBalancer* b, Element* fd, String config) : BalanceMethodRSS(b,fd) {
    }

    int configure(Vector<String> &, ErrorHandler *) override CLICK_COLD;
    int initialize(ErrorHandler *errh, int startwith) override;

    void rebalance(Vector<Pair<int,float>> load) override;
private:

    AggregateCounterVector* _counter;
    float _target_load;
    float _imbalance_alpha;
};

/*
=title DeviceBalancer

=c

DeviceBalancer()
*/

class DeviceBalancer : public Element {
public:

    DeviceBalancer() CLICK_COLD;
    ~DeviceBalancer() CLICK_COLD;

    const char *class_name() const { return "DeviceBalancer"; }
    const char *port_count() const { return "0/0"; }
    const char *processing() const { return AGNOSTIC; }

    bool can_live_reconfigure() const { return false; }

    int configure(Vector<String> &, ErrorHandler *) CLICK_COLD;
    int initialize(ErrorHandler *) CLICK_COLD;
/*    void add_handlers() CLICK_COLD;
    void cleanup(CleanupStage) CLICK_COLD;*/

    void run_timer(Timer* t) override;



    BalanceMethod* _method;
    int _core_offset;
    int _tick;
    Timer _timer;
    int _max_cpus;
    target_method _target;
    load_method _load;
    int _startwith;
    bool _autoscale;

    bool _verbose;

    struct CpuInfo {
	int id;
	unsigned long long last_cycles;
    };
    CpuInfo make_info(int _id);
	int addCore();
	void removeCore(int phys_id);

    Vector<CpuInfo> _used_cpus;
    Vector<int> _available_cpus;
    double _overloaded_thresh = 0.75;
    double _underloaded_thresh = 0.25;
};

CLICK_ENDDECLS

#endif
