/**
 *                      Module: Metron agent
 *                    Platform: FastClick
 *                 Network I/O: DPDK
 *         Traffic dispatching: MAC-based VMDq (mapping between MAC address values to CPU cores)
 *  Header field acting as tag: Destination MAC address
 *        Who tags the packets: A node that preceeds this one (e.g., an OpenFlow switch)
 *               Service chain: Sent by the controller on demand
 *
 * In this configuration, DPDK's VMDq is used to dispatch input traffic
 * to multiple NIC hardware queues. Metron associates each hardware
 * queue with a CPU core, on which the Metron controller deploys software
 * pipelines. Each software pipeline interacts with the underlying NICs
 * using DPDK.
 *
 * Each pipeline created by the controller runs as a secondary (slave)
 * process and can be represented as an individual Click configuration
 * as follows:
 *    FromDPDKDevice(QUEUE X) -> .... -> ToDPDKDevice(QUEUE X)
 * 
 * Each pipeline uses an Rx/Tx queue pair associated with a CPU core,
 * which runs to completion.
 * 
 * For VMDq to work properly, packets must arrive with their destination MAC
 * address already set to a known value. This value is matched by VMDq and an
 * associated CPU core index is found. This is where the packet will be dispatched.
 * Moreover, the NIC must support VMDq (Driver configuration might be required).
 */

/**
 * Assuming a server with 8 cores, deploy as follows:
 * sudo ../../bin/click --dpdk -l 0-7 -w 0000:01:00.0 -w 0000:01:00.1 -v -- metron-dispatcher-vmdq.click queues=8
 */

/* Configuration arguments */
define(
	$ifacePCI0      0000:01:00.0,
	$ifacePCI1      0000:01:00.1,

	$queues         8,

	$dpdkRxMode     vmdq,          // DPDK's VMDq
	$metronRxMode   mac,           // Using MAC address as a tag
	$vfPools        8,             // A pool of 8 MAC addresses, each mapped to a CPU core

	$agentIp        127.0.0.1,
	$agentPort      80,
	$agentVerbosity true,

	$discoverIp     127.0.0.1,
	$discoverPort   8181,
	$discoverUser   onos,
	$discoverPass   rocks,

	$monitoringMode false
);

/* Instantiate the http server to access Metron's handlers */
http :: HTTPServer(PORT 80);

/* Metron agent */
metron :: Metron(
	NIC               fd0,
	NIC               fd1,                 // NICs to use
	RX_MODE           $metronRxMode,       // Rx filter mode
	AGENT_IP          $agentIp,            // Own IP address
	AGENT_PORT        $agentPort,          // Own HTTP port
	DISCOVER_IP       $discoverIp,         // ONOS's REST API address
	DISCOVER_PORT     $discoverPort,       // ONOS's REST API port
	DISCOVER_USER     $discoverUser,       // ONOS's REST API username
	DISCOVER_PASSWORD $discoverPass,       // ONOS's REST API password
	SLAVE_DPDK_ARGS   "-w$ifacePCI0",      // Arguments passed to secondary DPDK processes
	SLAVE_DPDK_ARGS   "-w$ifacePCI1",
	SLAVE_DPDK_ARGS   "--log-level=eal,8",
	MONITORING        $monitoringMode,     // Defines the monitoring behaviour of the agent
	VERBOSE           $agentVerbosity      // Defines the verbosity level of the agent
);

/**
 * NICs are inactive at boot time.
 * Slave processes will be dynamically created by the Metron controller.
 */
fd0 :: FromDPDKDevice(PORT $ifacePCI0, MODE $dpdkRxMode, VF_POOLS $vfPools, N_QUEUES $queues, ACTIVE false)
	-> Idle
	-> ToDPDKDevice(PORT $ifacePCI0, N_QUEUES $queues);

fd1 :: FromDPDKDevice(PORT $ifacePCI1, MODE $dpdkRxMode, VF_POOLS $vfPools, N_QUEUES $queues, ACTIVE false)
	-> Idle
	-> ToDPDKDevice(PORT $ifacePCI1, N_QUEUES $queues);
