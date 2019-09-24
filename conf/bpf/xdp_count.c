#include <linux/bpf.h>
#include <linux/if_ether.h>
#include <linux/if_packet.h>
#include <linux/if_vlan.h>
#include <linux/types.h>
#include <asm-generic/types.h>
#include "bpf_helpers.h"
//#include <linux/ip.h>
//#include <linux/in.h>
//#include <linux/tcp.h>


struct bpf_map_def SEC("maps") count_map = {
	.type = BPF_MAP_TYPE_PERCPU_ARRAY,
	.key_size = 4, /* u8 does not work?! */
	.value_size = 8,
	.max_entries = 512,
};

#define bpf_debug(fmt, ...)                          \
    ({                                               \
        char ____fmt[] = fmt;                        \
        bpf_trace_printk(____fmt, sizeof(____fmt),   \
            ##__VA_ARGS__);                          \
    })

SEC("prog")
int  xdp_count_program(struct xdp_md *ctx)
{
/*	void *data_end = (void *)(long)ctx->data_end;
	void *data     = (void *)(long)ctx->data;
	struct ethhdr *eth = data;*/
    __u64 *value;
/*	u16 eth_proto = 0;
	u64 l3_offset = 0;
	u32 action;

	if (!(parse_eth(eth, data_end, &eth_proto, &l3_offset))) {
		bpf_debug("Cannot parse L2: L3off:%llu proto:0x%x\n",
			  l3_offset, eth_proto);
		return XDP_PASS;
	}

	bpf_debug("Reached L3: L3off:%llu proto:0x%x\n", l3_offset, eth_proto);

	action = handle_eth_protocol(ctx, eth_proto, l3_offset);*/
    __u32 id = ctx->hash % 128;
/*    if (id % 16 != bpf_get_smp_processor_id()) {
	    bpf_debug("Received idx %u on core %d instead of %d\n", ctx->hash, bpf_get_smp_processor_id(), id % 16 );
        bpf_debug("queue %d", ctx->rx_queue_index);
    }*/

    value = bpf_map_lookup_elem(&count_map, &id);
    if (value)
        *value += 1;
//    bpf_map_update_elem(&count_map, &ctx->rx_queue_index, n);
	return XDP_PASS;
}

char _license[] SEC("license") = "GPL";
