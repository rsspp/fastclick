struct bpf_map_def SEC("maps") count_map = {
	.type = BPF_MAP_TYPE_PERCPU_ARRAY,
	.key_size = 4, /* u8 does not work?! */
	.value_size = 8,
	.max_entries = 512,
};


