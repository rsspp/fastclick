#ifndef CLICK_ETHDEVICE_HH
#define CLICK_ETHDEVICE_HH

struct EthernetDevice;

typedef Vector<unsigned> (*eth_get_rss_reta)(struct EthernetDevice* eth);
typedef int (*eth_get_rss_reta_size)(struct EthernetDevice* eth);
typedef int (*eth_set_rss_reta)(struct EthernetDevice* eth, const Vector<unsigned>&);

struct EthernetDevice {
	EthernetDevice() : get_rss_reta(0) {

	}

	eth_get_rss_reta_size get_rss_reta_size;
	eth_set_rss_reta set_rss_reta;
	eth_get_rss_reta get_rss_reta;
};

#endif
