# RSS++

This repository is a modified version of FastClick that implements RSS++, a load and state-aware intra-server load-balancer. RSS++ works by tweaking NICs' RSS indirection tables. To do so, RSS++ monitors the load of each RSS bucket and solves an optimization problem to re-assign RSS buckets to different CPU cores. Moreover, RSS++ proposes a state migration algorithm to avoid synchronization problems when rebalancing.
RSS++ can load-balance either FastClick applications or any socket application by attaching to XDP using BPF code and ethtool to change the indirection table. This is *NOT* Click in Kernel mode, it only uses standard APIs to communicate with the Kernel.

## Building

Note that if you only want to reproduce experiments, NPF can build RSS++ for you, so directly check the [experiments repo](https://github.com/rsspp/experiments) ;)

One compiles RSS++ as FastClick, only with the 3 supplementary flags (at the end) :
```
./configure --enable-multithread --disable-linuxmodule --enable-intel-cpu --enable-user-multithread --verbose CFLAGS="-g -O3" CXXFLAGS="-g -std=gnu++11 -O3" --disable-dynamic-linking --enable-poll --enable-bound-port-transfer --enable-dpdk --enable-batch --with-netmap=no --enable-zerocopy --disable-dpdk-pool --disable-dpdk-packet --enable-flow --disable-task-stats --enable-cpu-load
```
If you're only interested in Kernel mode (balancing your existing socket application), just remove all the "dpdk" stuffs in the line above (there are 3 options).

## Usage
This section explains how to use RSS++.

### Daemon mode
In daemon mode (sometimes refered to as  Kernel mode), no packets flow through Click, which is used as a daemon, much like irqbalance. The DeviceBalancer element, the core logic of RSS++, use a BPF program to count packets and reads them through BPF maps. According to the result of the optimization, the element will re-program the indirection table using the ethtool API.

```
click -e 'define ($IF eth0,
                  $CPU 18,
                  $RETA 128,
                  $VERBOSE 0,
                  $TIMER 10,
                  $BPF_XDP_COUNT "conf/bpf/xdp_count.o");
        
fd :: FromDevice($IF, ACTIVE false) -> Discard;

db :: DeviceBalancer(
    METHOD pianorss,
    DEV fd,
    CPUS $CPU,
    RSSCOUNTER xdp,
    CYCLES realcpu,
    RETA_SIZE $RETA,
    VERBOSE $VERBOSE,
    TARGET_LOAD 0.75,
    TIMER $TIMER,
    TIMER_MAX 1000,
    IMBALANCE_THRESHOLD 0.02,
    AUTOSCALE false);

StaticThreadSched(db 17) //Device balancer will run on the last core

xdp :: XDPLoader(PATH $BPF_XDP_COUNT, DEV $IF, CLEAN false)';

```
Note that "ACTIVE false" is passed to FromDevice, as Click is not reading packets to process them. The element is only used as a proxy for configuration.

### DPDK standalone mode
RSS++ also works on top of Click and DPDK (it could also work with Netmap with a few modifications). In this case, RSS++ counts packets using an AggregateCounteVector element and changes the NIC indirection table directly through the DPDK API.

An example of involved elements would be:
```
FromDPDKDevice(..., RSS_AGGREGATE 1)
          -> agg :: AggregateCounterVector(MASK 511)
          -> ...
          -> ToDPDKDevice(...)
    
balancer :: DeviceBalancer(DEV fd0, METHOD pianorss, VERBOSE 0, TIMER 100, CPUS $CPU, TARGET 0.75, STARTCPU -1, LOAD 0.90, RSSCOUNTER agg, AUTOSCALE 1, CYCLES realcpu, RETA_SIZE $RETA_SIZE, IMBALANCE_THRESHOLD 0.02);
```
Basically, we use AggregateCounterVector to count packets going through a given element instead of XDP.

Given that DPDK supports applications running as secondary processes (i.e., those attached to a primary process), one could use RSS++ to tune a NIC's indirection table, while another application (e.g., mTCP) is spawned. In this case, one has to implement his/her own mechanisms for reporting packet counters.

## Implementation
This section goes deeper into the implementation of those elements.

### Load balancing
The heart of RSS++ is the DeviceBalancer element ([elements/userlevel/devicebalancer.cc](https://github.com/rsspp/fastclick/blob/master/elements/userlevel/devicebalancer.hh)). This element periodically reads the load of a set of CPU cores using a method described by the "CYCLES" parameter, and then pass the load information to a balancing module. The balancing module is chosen using the "METHOD" parameter. The balancing module then programs the NIC through a common interface defined by the "DEV" argument.

#### CPU load (CYCLES)
The "realcpu" CYCLES parameter will parse /proc/stat to get the current CPU load. "cycles" will use the proportion of useful cycles over useless cycles as described in the paper. "cyclesqueue" will do the same but when a CPU is overloaded, will use the size of the dedicated NIC queue on top of the number of cycles. The kernelmode should use "realcpu", while DPDK mode should use "cycles" or "cyclesqueue". The code is available in [elements/userlevel/devicebalancer.cc#L1578](https://github.com/rsspp/fastclick/blob/master/elements/userlevel/devicebalancer.cc#L1578).

#### Balancing module (METHOD)
Each module is backed by a C++ class (metron is MethodMetron, rss++/rsspp is MethodPianoRSS, "rss" is BalanceMethodRSS, ...).

##### BalanceMethodRSS
This class set up RSS, the balancing periodic call is ignored.

##### MethodPianoRSS
This class implements RSS++'s logic. It inherits BalanceMethodRSS but do not ignore the balancing call.
The code is documented and available at [elements/userlevel/devicebalancer.cc#L530](https://github.com/rsspp/fastclick/blob/master/elements/userlevel/devicebalancer.cc#L530).
In a nutshell there are steps:
 - Copy the packet counters if in BPF mode (value is directly accessed in DPDK mode).
 - Do a first pass over the load to counts the overloaded cores and compute the average load.
 - Count the number of packets per core
 - If autoscale is enabled, remove or add a core according to the logic described in the paper. If a core is to be removed, an instance of BucketMapProblem is built and solved to move all buckets of the removed core to other cores. In that case the RSS table is rewritten and the function returns.
 - Then the CPUs are separated in two sets : overloaded and underloaded.
 - An instance of BucketMapTargetProblem  is built and solved, that realize the optimization described in the paper.
 - If a move was performed, the table is written to the NIC.

#### NIC programing (DEV)
The class is provided with an interface to the NIC through the "DEV" argument. It must be the name of another element that inherit the EthernetDevice class, currently only FromDevice or FromDPDKDevice. Both of these elements will implement two functions to program the NIC indirection table (set_rss_reta) and read its size for development purpose (get_rss_reta_size). FromDevice will use the standard ioctls to set it using the ethtool API ([elements/userlevel/fromdevice.cc#L768](https://github.com/rsspp/fastclick/blob/master/elements/userlevel/fromdevice.cc#L768)). The FromDPDKDevice elements will do the same using directly the DPDK interface ([lib/dpdkdevice.cc#L117](https://github.com/rsspp/fastclick/blob/master/lib/dpdkdevice.cc#L117))
Note that the RSS balancing elements may directly program the RSS table using the DPDK rte_flow API instead to overcome the limitation of Mellanox NICs mentioned in the paper.

### State migration
When DeviceBalancer decides to migrate a group of flows (i.e., an RSS indirection bucket), it sends a message to an element of choice. For example, the FlowIPManager can receive such a message and start the migration process of the corresponding group table as explained in the paper. The previous example is augmented below to include such flow processing element:

```
FromDPDKDevice(..., RSS_AGGREGATE 1)
          -> agg :: AggregateCounterVector(MASK 511)
          -> flow :: FlowIPManager(GROUPS 512)
          -> nat :: FlowIPNAT(SIP 10.0.0.1)
          -> ...
          -> ToDPDKDevice(...)
    
balancer :: DeviceBalancer(DEV fd0, METHOD pianorss, VERBOSE 0, TIMER 100, CPUS $CPU, TARGET 0.75, STARTCPU -1, LOAD 0.90, RSSCOUNTER agg, AUTOSCALE 1, CYCLES realcpu, RETA_SIZE $RETA_SIZE, IMBALANCE_THRESHOLD 0.02, MANAGER flow);
```


## Parent's code
Technically this repository is based on Metron's FastClick branch, supporting the paper of the same name (see README.metron.md) but this is only to compare against Metron, by using its flow parsing methods.

## Future plan
All RSS++ specific code will be merged into mainline FastClick, at which point this repository will be closed.
