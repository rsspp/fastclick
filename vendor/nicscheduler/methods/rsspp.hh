#ifndef LIBNICSCHEDULER_METHODRSSPP_HH
#define LIBNICSCHEDULER_METHODRSSPP_HH 1

class MethodRSSPP : public MethodRSS, public LoadTracker { public:

    MethodRSSPP(NICScheduler* b, EthernetDevice* fd);

    int initialize(ErrorHandler *errh, int startwith) override;

    virtual std::string name() override CLICK_COLD { return "rsspp"; }

    void rebalance(std::vector<std::pair<int,float>> load) override;

    //We keep the space here to avoid reallocation
    struct Node {
        uint64_t count;
        uint64_t variance;
        bool moved;
    };
    std::vector<Node> _count;

    bool _counter_is_xdp;
    int _xdp_table_fd;
    void* _counter;

    float _target_load;
    float _imbalance_alpha;
    float _threshold;

    bool _dancer;
    bool _numa;
    int _numa_num;
};

#endif
