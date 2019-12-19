# RSS++

This repository is a modified version of FastClick that includes support for [RSS++][rsspp-paper], a load and state-aware intra-server load-balancer.
RSS++ works by tweaking NICs' RSS indirection tables. To do so, RSS++ monitors the load of each RSS bucket and solves an optimization problem to re-assign RSS buckets to different CPU cores. Moreover, RSS++ proposes a state migration algorithm to avoid synchronization problems while rebalancing.
RSS++ can load-balance either FastClick applications or any socket application by attaching to XDP using BPF code and ethtool to change the indirection table. This is *NOT* Click in Kernel mode, it only uses standard APIs to communicate with the Kernel.

## NICScheduler library
The heart of RSS++ is the NICScheduler library (vendor/nicscheduler/). This library implements the logic to "balance" packets from an input device according to the load of CPUs. The NICScheduler library provides multiple scheduling strategies, RSS++ is one of them.

The library is included in FastClick's DeviceBalancer element, that uses another element to read packet counters and implements the logic to read the CPU load and set NICScheduler's parameters, such as the balancing method. FastClick is also used as a proxy to reconfigure either a Linux interface or a DPDK interface. Most dependencies to Click have been removed from the NICScheduler code, only a few utilities are left. Particularly, the Metron method is highly dependent to Click, and will only work under FastClick's directory.

## Building

Note that if you only want to reproduce experiments, NPF can build RSS++ for you, so directly check the [experiments repo](https://github.com/rsspp/experiments) ;)

One compiles RSS++ as FastClick, only with the 3 supplementary flags (at the end):
```
./configure --enable-multithread --disable-linuxmodule --enable-intel-cpu --enable-user-multithread --verbose CFLAGS="-g -O3" CXXFLAGS="-g -std=gnu++11 -O3" --disable-dynamic-linking --enable-poll --enable-bound-port-transfer --enable-dpdk --enable-batch --with-netmap=no --enable-zerocopy --disable-dpdk-pool --disable-dpdk-packet --enable-flow --disable-task-stats --enable-cpu-load
```
If you're only interested in Kernel mode (balancing your existing socket application), just remove all the "dpdk" stuff in the line above (there are 3 options).

## Usage
This section explains how to use RSS++.

### Daemon mode
In daemon mode (sometimes refered to as  Kernel mode), no packets flow through Click, which is used as a daemon, much like irqbalance. The DeviceBalancer element, the core logic of RSS++, uses a BPF program to count packets and reads them through BPF maps. According to the result of the optimization, the element will re-program the NIC indirection table using the ethtool API.

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

## NICScheduler library
This section goes deeper into the implementation of NICScheduler.

The heart of the NIC-driven scheduling library is the NICScheduler class ([vendor/nicscheduler/nicscheduler.hh](https://github.com/rsspp/fastclick/blob/mvendor/nicscheduler/nicscheduler.hh)).
The NICScheduler class retains some parameters, such as the underload and overloaded thresholds, the chosen balancing module (one of the class implementing BalanceMethod), the set of available and currently used CPU cores, etc.

One must call set_method(method, device), where method is one of the methods below. Device is an implementation of the EthernetDevice that allows to access some device parameters and reprogram the redirection table (RETA) table for the RSS based methods.

### Balancing module
Each module is backed by a C++ class (metron is MethodMetron, RSS++/rsspp is MethodRSSPP, rss is MethodRSS, etc.).

#### MethodRSS
This class sets up RSS indirection table as a normal driver would do and the balancing periodic call is ignored. Therefore this is the classical RSS sharded approach.

#### MethodRSSPP
This class implements RSS++'s logic. It inherits MethodRSS but does not ignore the balancing call.
The code is documented and available at [vendor/nicscheduler/methods/rsspp.hh](https://github.com/rsspp/fastclick/blob/master/vendor/nicscheduler/methods/rsspp.hh).
The various phases of the function work as follow:
 - Copy the packet counters if in BPF mode (value is directly accessed in DPDK mode).
 - Do a first pass over the load to count the overloaded cores and compute the average load.
 - Count the number of packets per core
 - If autoscale is enabled, remove or add a core according to the logic described in the paper. If a core is to be removed, an instance of BucketMapProblem is built and solved to move all buckets of the removed core to other cores. In that case the RSS table is rewritten and the function returns.
 - Then, the CPUs are separated in two sets : overloaded and underloaded.
 - An instance of BucketMapTargetProblem is built and solved in order to realize the optimization described in the paper.
 - If a move was performed, the table is written to the NIC.

### NIC programing
The class is provided with an interface to the NIC through an implementation of the EthernetDevice class.
Note that the RSS balancing elements may directly program the RSS table using the DPDK rte_flow API to overcome the limitation of Mellanox NICs mentioned in the paper.

### State migration
When the RSS++ method decides to migrate a group of flows (i.e., an RSS indirection bucket), it calls a function from a listener before and after migration, an implementation of the MigrationListener implementing some virtual functions. This is set up by calling the set_migration_listener() function.

## NICScheduler integration in FastClick
The DeviceBalancer element ([elements/userlevel/devicebalancer.cc].
This element periodically reads the load of a set of CPU cores using a method described by the "CYCLES" parameter and then passes the load information to a balancing module. The balancing module is chosen using the "METHOD" parameter. The balancing module then programs the NIC through a common interface defined by the "DEV" argument.

### Balancing methode (METHOD)
The value is simply passed to the NICScheduler library.

### CPU load (CYCLES)
The "realcpu" CYCLES parameter will parse /proc/stat to get the current CPU load. "cycles" will use the proportion of useful cycles over useless cycles as described in the paper. "cyclesqueue" will do the same, but when a CPU is overloaded it will use the size of the dedicated NIC queue on top of the number of cycles. The kernel mode should use "realcpu", while DPDK mode should use "cycles" or "cyclesqueue". The code is available at [elements/userlevel/devicebalancer.cc#L1578](https://github.com/rsspp/fastclick/blob/master/elements/userlevel/devicebalancer.cc#L1578).

### Device proxy (DEV)
It must be the name of another element that inherits the EthernetDevice class, currently only FromDevice or FromDPDKDevice. Both of these elements will implement two functions to program the NIC indirection table (set_rss_reta) and read its size for development purpose (get_rss_reta_size). FromDevice will use the standard ioctls to set it using the ethtool API ([elements/userlevel/fromdevice.cc#L768](https://github.com/rsspp/fastclick/blob/master/elements/userlevel/fromdevice.cc#L768)). The FromDPDKDevice elements will do the same directly using the DPDK interface ([lib/dpdkdevice.cc#L117](https://github.com/rsspp/fastclick/blob/master/lib/dpdkdevice.cc#L117)).

### Flow manager (MANAGER)
The FlowIPManager implements the MigrationManager. It uses the pre-migration and post-migration calls to receive such a message and start the migration process of the corresponding group table as explained in the paper. The previous example is augmented below to include such flow processing element:

```
FromDPDKDevice(..., RSS_AGGREGATE 1)
          -> agg :: AggregateCounterVector(MASK 511)
          -> flow :: FlowIPManager(GROUPS 512)
          -> nat :: FlowIPNAT(SIP 10.0.0.1)
          -> ...
          -> ToDPDKDevice(...)
    
balancer :: DeviceBalancer(DEV fd0, METHOD pianorss, VERBOSE 0, TIMER 100, CPUS $CPU, TARGET 0.75, STARTCPU -1, LOAD 0.90, RSSCOUNTER agg, AUTOSCALE 1, CYCLES realcpu, RETA_SIZE $RETA_SIZE, IMBALANCE_THRESHOLD 0.02, MANAGER flow);
```
If a "MANAGER" (in this case FlowIPManager) is set up in DeviceBalancer, the RSS++ balancing module will call the `pre_migrate()` function of that element before writing the indirection table. And the `post_migrate()` function will be called afterwards. The pre_migrate function saves the moves to be applied (migration of bucket tables from one core to another) to be operated when a new assignment is observed on any given core.
The post_migrate function will read the number of packets in the queue and set up the per-bucket trigger at which a migration is to be considered as done by the NIC. The code of the flow table itself strictly follows the migration process explained in the paper.

## Parent's code
Technically this repository is based on Metron's FastClick branch, supporting the paper of the same name (see README.metron.md) but this is only to compare against Metron, by using its flow parsing methods.

## Future plan
All RSS++ specific code will be merged into mainline FastClick, at which point this repository will be closed.

## Citing RSS++
If you use RSS++ in your work, please cite our [paper][rsspp-paper]:
```
@inproceedings{barbette-rsspp.conext19,
	author    = {Barbette, Tom and Katsikas, Georgios P. and Maguire,Jr., Gerald Q. and Kosti\'{c}, Dejan},
	title     = {{RSS++: load and state-aware receive side scaling}},
	booktitle = {Proceedings of the 15th International Conference on Emerging Networking Experiments And Technologies},
	series    = {CoNEXT '19},
	year      = {2019},
	isbn      = {978-1-4503-6998-5},
	location  = {Orlando, Florida},
	pages     = {318--333},
	numpages  = {16},
	url       = {http://doi.acm.org/10.1145/3359989.3365412},
	doi       = {10.1145/3359989.3365412},
	acmid     = {3365412},
	publisher = {ACM},
	address   = {New York, NY, USA},
	keywords  = {NIC indirection table, intra server, load-balancing, state-aware},
}
```

[rsspp-paper]:https://dl.acm.org/ft_gateway.cfm?id=3365412&ftid=2099399&dwn=1&CFID=117388292&CFTOKEN=f70f1ec24d3089b2-CBE9F99A-0D46-91B4-51B4F56D84A4B4D2
