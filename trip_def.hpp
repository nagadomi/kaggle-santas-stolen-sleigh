#pragma once
#include <vector>
#include <cmath>

typedef struct gift {
	int id;
	double lat;
	double lng;
	double x;
	double y;
	double z;
	double weight;
	int trip_id;

	gift(): id(0), lat(0), lng(0), x(0), y(0), z(0), weight(0), trip_id(0) {}
	gift(double lat, double lng, double weight)
		: id(0), lat(lat), lng(lng), weight(weight), trip_id(0) 
	{
		double lat_rads = (M_PI / 180.0) * lat;
		double lng_rads = (M_PI / 180.0) * lng;
		x = 0.5 * std::cos(lat_rads) * std::sin(lng_rads);
		y = 0.5 * std::cos(lat_rads) * std::cos(lng_rads);
		z = 0.5 * std::sin(lat_rads);
	}
} gift_t;
typedef std::vector<gift_t> gifts_t;

typedef struct trip {
	static int g_node_id;
	int id;
	gifts_t gifts;
	double weight;

	trip() : weight(0) {
		id = new_id();
	}
	static int new_id()
	{
		int id = g_node_id;
		g_node_id++;
		return id;
	}
} trip_t;
typedef std::vector<trip_t> trips_t;

typedef struct path_cache {
	double dist;
	double prev_weight;
	double sum_dist;
	double wtl;
	path_cache(double d, double p, double s, double w) 
		: dist(d), prev_weight(p), sum_dist(s), wtl(w) {}
} path_cache_t;

