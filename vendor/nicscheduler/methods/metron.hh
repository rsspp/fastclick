#ifndef LIBNICSCHEDULER_METHODMETRON_HH
#define LIBNICSCHEDULER_METHODMETRON_HH 1
/**
 * Metron's reimplemented balancing method
 */
class MethodMetron : public BalanceMethodDevice, public LoadTracker { public:

    MethodMetron(NICScheduler* b, EthernetDevice* fd) : BalanceMethodDevice(b,fd) {

    }

    virtual std::string name() override CLICK_COLD { return "metron"; }

    int initialize(ErrorHandler *errh, int startwith) override CLICK_COLD;

    void rebalance(std::vector<std::pair<int,float>> load) override;

    String _rules_file;
    int _min_movement;
    int _deflation_factor;
};
#endif
