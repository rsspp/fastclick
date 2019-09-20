RSS++
=====

This repository is a modified version of FastClick that implements RSS++, a stateful intra-server load-balancer. RSS++ works by tweaking NICs' RSS indirection tables. To do so, RSS++ monitors the load of each RSS bucket and solves an optimization problem to re-assign RSS buckets to CPU cores. Moreover, RSS++ proposes a state migration algorithm to avoid synchronization problems when doing so.
RSS++ can load-balance either FastClick applications or any socket application by attaching to XDP using BPF code and ethtool to change the indirection table. This is *NOT* Click in Kernel mode.

Load-balancing
--------------
The heart of RSS++ is the DeviceBalancer element. This element can load-balance NIC flow tables inside the kernel using ethtool, while counting entries with a BPF program. In that case, no packets flow through Click, which is used mainly as a library. We are looking on packaging to remove the need for Click in this case, stay tuned.

TODO : add example of kernel config

RSS++ also works on top of Click and DPDK (could use also Netmap with a few modifications). In this case, RSS++ counts packets using an AggregateCounterMP element and changes the NIC indirection table directly through the DPDK API.
Given that DPDK supports applications running as secondary processes (i.e., those attached to a primary process), one could use RSS++ to tune a NIC's indirection table, while another application (e.g., mTCP) is spawned. In this case, one has to implement his/her own mechanisns for reporting packet counters.

TODO : add important bits of user code.

State migration
---------------
When DeviceBalancer decides to migrate a group of flows (i.e., an RSS indirection bucket), it sends a message to an element of choice. For example, the FlowIPManager can receive such a message and start the migration process of the corresponding group table as explained in the paper.

Parent's code
-------------
Technically this repository is based on Metron's FastClick branch, supporting the paper of the same name (see README.metron.md) but this is only to compare against Metron, by using its flow parsing methods.

Future plan
-----------
All RSS++ specific code will be merged into mainline FastClick, at which point this repository will be closed.
