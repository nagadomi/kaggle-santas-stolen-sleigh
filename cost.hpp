#include "trip_def.hpp"
#include <cmath>
#include <cfloat>

// ref: https://www.kaggle.com/c/santas-stolen-sleigh/forums/t/18049/simpler-faster-haversine-distance
static inline double 
haversine_distance(double x1, double y1, double z1,
				   double x2, double y2, double z2)
{
	const static double R = 6371.0 * 2.0;
	const double dx = x1 - x2;
	const double dy = y1 - y2;
	const double dz = z1 - z2;
	return std::asin(std::sqrt(dx * dx + dy * dy + dz * dz)) * R;
}

static inline double 
haversine_distance(const gift_t &gift1, const gift_t &gift2)
{
	return haversine_distance(gift1.x, gift1.y, gift1.z,
							  gift2.x, gift2.y, gift2.z);
}

const static gift_t g_pole(90.0, 0.0, 10.0);

// weighted trip length
static inline double
wtl(const gifts_t &gifts)
{
	const size_t n = gifts.size();

	double dist = 0.0;
	double prev_weight = g_pole.weight;
	const gift_t *prev_stop;

	prev_stop = &g_pole;
	for (auto i = gifts.begin(); i != gifts.end(); ++i) {
		prev_weight += i->weight;
	}
	for (size_t i = 0; i < n; ++i) {
		dist += haversine_distance(gifts[i], *prev_stop) * prev_weight;
		prev_stop = &gifts[i];
		prev_weight -= prev_stop->weight;
	}
	dist += haversine_distance(*prev_stop, g_pole) * prev_weight;

	return dist;
}

static inline double
wtl(const gifts_t &gifts, double weight_sum)
{
	const size_t n = gifts.size();

	double dist = 0.0;
	double prev_weight = g_pole.weight + weight_sum;
	const gift_t *prev_stop = &g_pole;
	for (size_t i = 0; i < n; ++i) {
		dist += haversine_distance(gifts[i], *prev_stop) * prev_weight;
		prev_stop = &gifts[i];
		prev_weight -= prev_stop->weight;
	}
	dist += haversine_distance(*prev_stop, g_pole) * prev_weight;

	return dist;
}

// this cost function can be computed in constant time.
static inline double
wtl_test_inserting(const gifts_t &gifts,
				   double weight_sum,
				   const std::vector<path_cache_t> &path_cache,
				   size_t pos, const gift_t &gift)
{
	double dist = 0.0;
	double prev_weight = g_pole.weight + weight_sum + gift.weight;
	const gift_t *prev_stop;
	const size_t n = gifts.size();

	prev_stop = &g_pole;
	if (pos > 0) {
		const path_cache_t *p = &path_cache[pos-1];
		dist = p->wtl + p->sum_dist * gift.weight;
		prev_stop = &gifts[pos-1];
		prev_weight = p->prev_weight + gift.weight;
	}
	dist += haversine_distance(*prev_stop, gift) * prev_weight;
	prev_stop = &gift;
	prev_weight -= prev_stop->weight;
	if (pos < n) {
		dist += haversine_distance(*prev_stop, gifts[pos]) * prev_weight + 
			(path_cache[n].wtl - path_cache[pos].wtl);
	} else {
		dist += haversine_distance(*prev_stop, g_pole) * prev_weight;
	}

	return dist;
}
