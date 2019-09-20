RSS++
=====

This repository is a modified version of FastClick that implements RSS++, an intra-server load-balancer. RSS++ works by tweaking the indirection table. RSS++ will watch the load of each RSS buckets and solve an optimization problem to re-assign buckets to cores. RSS++ proposes a state migration algorithm to avoid synchronization problems when doing so.
RSS++ works with a FastClick application, or can load-balance any socket application by attaching to XDP using BPF code and ethtool to change the indirection table. This is *NOT* Click in Kernel mode.

Load-balancing
--------------
The heart of RSS++ is the DeviceBalancer element. It can load-balance the a "real" NIC flow table inside the kernel, using ethtool and counting entries with a BPF program. In that case no packets flow through Click, and it is used mainly as a library. We will be looking on packaging to remove the need of Click for this mode. Stay tuned.

TODO : add example of kernel config

We also have a real in-Click mode that uses DPDK (could use also Netmap with a few modifications) to flow packets through click, counting them using AggregateCounterMP and changing the indirection table through the DPDK apis directly. This one is of course much more high-speed.
Given that DPDK supports secondary application, one could use RSS++ to tune the indirection table while another application (like mTCP based). The problem is then reporting the counting, you have to implement your own mechanism.

TODO : add important bits of user code.

State migration
---------------
When DeviceBalancer decides it's time to migrate a group of flows (an RSS indirection bucket), it will send a message to an optional element. The FlowIPManager can receive such a message and start the migration process of the corresponding group table as explained in the paper.


Parent's code
-------------

Technically this repository is based on Metron's FastClick branch, supporting the paper of the same name (see README.metron.md) but this is only to compare against Metron, by using its flow parsing methods.


Future plan
-----------
All RSS++ specific code should be merged in mainline FastClick, at which point this repository will be closed.
