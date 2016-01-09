#pragma once
#include <vector>
#include <map>
#include <set>
#include <string>
#include <iostream>
#include <cmath>
#include <cfloat>
#include <algorithm>
#include <iterator>
#include <unordered_map>
#include <climits>
#include <list>
#include <omp.h>
#include "trip_def.hpp"
#include "cost.hpp"
#include "load.hpp"
#include "tick.hpp"

int trip::g_node_id = 1;


static void 
trip_make_cluster_from_trip_id(std::map<int, gifts_t > &clusters,
							   const gifts_t &gifts)
{
	clusters.clear();
	for (auto i = gifts.begin(); i != gifts.end(); ++i) {
		auto node = clusters.find(i->trip_id);
		if (node != clusters.end()) {
			node->second.push_back(*i);
		} else {
			gifts_t val;
			val.push_back(*i);
			clusters.insert(std::make_pair(i->trip_id, val));
		}
	}
}

static void 
make_gift_id_map(std::map<int, gift_t *> &id_map, gifts_t &gifts)
{
	id_map.clear();
	for (auto i = gifts.begin(); i != gifts.end(); ++i) {
		id_map.insert(std::make_pair(i->id, &*i));
	}
}

void
trip_build(trips_t &trips, const gifts_t &gifts)
{
	std::map<int, gifts_t > clusters;
	trip_make_cluster_from_trip_id(clusters, gifts);
	trips.clear();
	for (auto c = clusters.begin(); c != clusters.end(); ++c) {
		trips.push_back(trip_t());
		trips.back().gifts = c->second;
	}
	for (auto trip= trips.begin(); trip != trips.end(); ++trip) {
		double w = 0.0;
		for (auto gift = trip->gifts.begin(); gift != trip->gifts.end(); ++gift) {
			w += gift->weight;
		}
		trip->weight = w;
	}
}

double trip_lb_score(const trips_t &trips)
{
	double score = 0.0;
	for (auto trip = trips.begin(); trip != trips.end(); ++trip) {
		score += wtl(trip->gifts);
	}
	return score;
}

void
trip_update(trip_t &trip)
{
	double w = 0.0;
	std::sort(trip.gifts.begin(), trip.gifts.end(), 
			  [](const gift_t &a, const gift_t &b) {
				  return a.lat > b.lat;
			  });
	for (auto gift = trip.gifts.begin(); gift != trip.gifts.end(); ++gift) {
		w += gift->weight;
	}
	trip.weight = w;
}

void 
trip_split(trips_t &trips, gifts_t &gifts, std::mt19937 &mt, 
		   double th1, double th2, double max_w1, double max_w2)
{
	gifts_t group;
	double w = 0.0;
	std::uniform_real_distribution<double> genw(max_w1, max_w2);
	double max_w = genw(mt);
	//group.reserve(100);
	trips.clear();
	//trips.reserve(2000);
	size_t c = 0;
	std::sort(gifts.begin(), gifts.end(),
			  [](const gift_t &a, const gift_t &b) {
				  if  (a.lng > b.lng) {
					  return true;
				  } else if (a.lng < b.lng) {
					  return false;
				  } else {
					  return a.lat > b.lat;
				  }
			  });
	for (auto gift = gifts.begin(); gift != gifts.end(); ++gift) {
		if (!(th1 >= gift->lat && th2 < gift->lat)) {
			continue;
		}
		if (w + gift->weight > max_w) {
			if (group.size() > 0) {
				trips.push_back(trip_t());
				trips.back().gifts = group;
			}
			group.clear();
			w = gift->weight;
			group.push_back(*gift);
			if (++c % 2 == 0) {
				max_w = max_w1;
			} else {
				max_w = max_w2;
			}
		} else {
			w += gift->weight;
			group.push_back(*gift);
		}
	}
	if (w != 0.0) {
		trips.push_back(trip_t());
		trips.back().gifts = group;
	}
	printf("trips: %zd\n", trips.size());
#pragma omp parallel for
	for (size_t i = 0; i < trips.size(); ++i) {
		trip_update(trips[i]);
	}
}

void
trip_id_set(trips_t &trips, gifts_t &gifts)
{
	std::map<int, gift_t *> id_map;
	make_gift_id_map(id_map, gifts);
	for (auto trip = trips.begin(); trip != trips.end(); ++trip) {
		for (auto gift = trip->gifts.begin(); gift != trip->gifts.end(); ++gift) {
			auto p = id_map.find(gift->id);
			if (p != id_map.end()) {
				p->second->trip_id = trip->id;
			} else {
				throw "Error unknown gift id";
			}
		}
	}
}

void 
trip_merge(trips_t &a, trips_t &b)
{
	std::map<std::pair<int, int>, double> cache;
	std::copy(b.begin(), b.end(), std::back_inserter(a));
}

void 
trip_to_gifts(gifts_t &gifts, trips_t &trips)
{
	gifts.clear();
	for (auto c = trips.begin(); c != trips.end(); ++c) {
		for (auto g = c->gifts.begin(); g != c->gifts.end(); ++g) {
			g->trip_id = c->id;
			gifts.push_back(*g);
		}
	}
}

void 
trip_save_gifts(trips_t &trips, const char *file)
{
	gifts_t gifts;
	trip_to_gifts(gifts, trips);
	save_gifts(file, gifts);
}

double
trip_update_path_cache(const trip_t &trip, std::vector<path_cache_t> &path_cache)
{
	const static gift_t pole(90.0, 0.0, 10.0);
	const gift_t *prev_stop = &pole;
	double sum_dist = 0.0;
	double seg_wtl = 0.0;
	double prev_weight = pole.weight + trip.weight;
	double dist;
	path_cache.clear();
	for (size_t i = 0; i < trip.gifts.size(); ++i) {
		dist = haversine_distance(*prev_stop, trip.gifts[i]);
		sum_dist += dist;
		seg_wtl += dist * prev_weight;
		prev_stop = &trip.gifts[i];
		prev_weight -= prev_stop->weight;
		path_cache.push_back(path_cache_t(dist, prev_weight, sum_dist, seg_wtl));
	}
	dist = haversine_distance(*prev_stop, pole);
	sum_dist += dist;
	seg_wtl += dist * prev_weight;
	path_cache.push_back(path_cache_t(dist, 0, sum_dist, seg_wtl));
	return seg_wtl;
}

void
trips_update_path_cache(const trips_t &trips, std::vector<path_cache_t> *path_cache)
{
	for (size_t i = 0; i < trips.size(); ++i) {
		trip_update_path_cache(trips[i], path_cache[i]);
	}
}

int trip_collect_garbage(const trips_t &trips, 
						 bool *zero_trips,
						 std::vector<size_t> *links,
						 std::set<size_t> *link_keys)
{
	int removed = 0;
	for (size_t i = 0; i < trips.size(); ++i) {
		if (trips[i].gifts.size() == 0) {
			if (!zero_trips[i]) {
				zero_trips[i] = true;
				++removed;
				for (size_t j = 0; j < trips.size(); ++j) {
					if (i == j) continue;
					for (size_t k = 0; k < links[j].size(); ++k) {
						if (links[j][k] == i) {
							links[j].erase(links[j].begin() + k);
							link_keys[j].erase(i);
							link_keys[i].erase(j);
							break;
						}
					}
				}
			}
		}
	}
	return removed;
}

size_t
trips_size(const trips_t &trips, bool *zero_trips)
{
	size_t s = 0;
	for (size_t i = 0; i < trips.size(); ++i) {
		if (trips[i].gifts.size() > 0) {
			++s;
			zero_trips[i] = false;
		}
	}
	return s;
}

void 
trip_optimize(trips_t &trips, const char *filename)
{
	std::mt19937 mt(71);
	std::uniform_int_distribution<int> geni(0, trips.size() - 1);
	std::uniform_real_distribution<double> genf(0, 1.0 - 1.0e-7);
	std::normal_distribution<> gennorm(0, 1.0);
	const static size_t PROCS = omp_get_num_procs();
	double wtl_cache[trips.size()];
	double remove_gift_scale = 4.5;
	double remove_trip_scale = 16;
	double global_search_rate = 0.025;
	double last_score = 0.0;
	size_t update_count = 0;
	size_t local_update_count = 0;
	size_t local_update_count2 = 0;
	size_t epoch = 0;
	size_t batch_size = 500;
	size_t destroy_trip_count = 0;
	size_t stucked_count = 0;
	size_t loop = 0;
	bool link_initialized = false;
	bool trip_destroy = true;
	size_t trip_destroy_enable_epoch = 0;
	std::vector<path_cache_t> path_cache[trips.size()];
	bool zero_trips[trips.size()];
	std::vector<size_t> links[trips.size()];
	std::set<size_t> link_keys[trips.size()];
	std::vector<size_t> trip_index;
	unsigned long t = tick();
	unsigned long start_time = tick();

	typedef struct backup {
		trip_t trip;
		std::vector<size_t> link;
		std::set<size_t> link_key;
		double wtl_score;
		std::vector<path_cache_t> path_cache;
		backup(const trip_t &t,
			   double w,
			   const std::vector<path_cache_t> &p,
			   const std::vector<size_t> &l,
			   const std::set<size_t> &k)
		{
			trip = t;
			wtl_score = w;
			path_cache = p;
			link = l;
			link_key = k;
		}
	} backup_t;

	for (size_t i = 0; i < trips.size(); ++i) {
		wtl_cache[i] = wtl(trips[i].gifts, trips[i].weight);
		links[i].push_back(i);
		link_keys[i].insert(i);
		trip_index.push_back(i);
		zero_trips[i] = false;
	}
	trips_update_path_cache(trips, path_cache);

	while (true) {
		++loop;
		if (loop > batch_size) {
			int removed = trip_collect_garbage(trips, zero_trips, links, link_keys);
			double current_score = trip_lb_score(trips);
			double accept_rate = update_count / (double)loop;
			size_t link_sum = 0;
			unsigned long now = tick();

			for (size_t i = 0; i < trips.size(); ++i) {
				link_sum += links[i].size();
			}
			std::printf("%4zd: %f, %4zd(local opt: %4zd, global opt: %4zd) trips: %zd, links: %.2f, accept: %f, destroy_trip: %zd, removed: %d, stucked: %zd, time: %.2fs, %02lu:%02lu\n",
						epoch++, current_score, update_count,
						local_update_count, local_update_count2,
						trips_size(trips, zero_trips),
						link_sum / (double)trips.size(),
						accept_rate,
						destroy_trip_count,
						removed,
						stucked_count,
						(now - t) / 1000.0,
						(now - start_time) / 3600000,
						((now - start_time) % 3600000) / 1000 / 60
				);
			if (!link_initialized && (link_sum / (double)trips.size()) > 10.0) {
				link_initialized = true;
				batch_size = 20000;
				std::printf("link_initialized\n");
			}
			trip_save_gifts(trips, filename);
			t = tick();
			loop = update_count = local_update_count = local_update_count2 = destroy_trip_count = 0;

			if (!trip_destroy &&
				epoch > trip_destroy_enable_epoch)
			{
				trip_destroy = true;
			}
			if (last_score > 0.0 && std::abs(last_score - current_score) < 1.0) {
				if (stucked_count < 150) {
					if (++stucked_count % 25 == 0) {
						remove_gift_scale += 1.0;
						if (stucked_count % 50 == 0) {
							remove_trip_scale += 1.0;
						}
						std::printf("stucked: remove_trip_scale -> %f,  remove_gift_scale -> %f\n",
									remove_trip_scale, remove_gift_scale);
					}
				}
			}
			last_score = current_score;
#if TICK
			std::printf("tick: %f\n", (tick() - t) / 1000.0);
			t = tick();
#endif
		}
		if (trip_destroy) {
			// destroy trip
			double best_imporoved_trip = -FLT_MAX;
			std::map<size_t, backup_t> best_state;
			std::shuffle(trip_index.begin(), trip_index.end(), mt);
			const size_t n = stucked_count == 0 ? (size_t)(trips.size() * 0.6) : trips.size();
			for (size_t i = 0; i < n; ++i) {
				size_t r1 = trip_index[i];
				std::vector<std::pair<size_t, gift_t> > removed_gifts;
				std::map<size_t, backup_t> backup;

				if (trips[r1].gifts.size() == 0) {
					continue;
				}
				backup.insert(std::make_pair(r1, backup_t(trips[r1], wtl_cache[r1], path_cache[r1], links[r1], link_keys[r1])));

				for (auto g = trips[r1].gifts.begin(); g != trips[r1].gifts.end(); ++g) {
					removed_gifts.push_back(std::make_pair(r1, *g));
				}
				trips[r1].gifts.clear();
				trips[r1].weight = 0;
				wtl_cache[r1] = trip_update_path_cache(trips[r1], path_cache[r1]);

				std::sort(removed_gifts.begin(), removed_gifts.end(),
						  [](const std::pair<size_t, gift_t> &a, const std::pair<size_t, gift_t> &b) {
							  return a.second.lat < b.second.lat;
						  });
				bool failback = false;
				for (auto rg = removed_gifts.begin(); rg != removed_gifts.end(); ++rg) {
					size_t best_i[PROCS];
					size_t best_j[PROCS];
					double best_improved[PROCS];
					const double w = rg->second.weight;
					const auto &index = trip_index;

					std::fill_n(best_improved, PROCS, FLT_MAX);
#pragma omp parallel for schedule(dynamic, 1)
					for (size_t l = 0; l < index.size(); ++l) {
						const size_t i = index[l];
						const trip_t *trip = &trips[i];
						const std::vector<path_cache_t> *pc = &path_cache[i];
						if (trip->weight + w < 1000.0) {
							const int tid = omp_get_thread_num();
							const size_t n = trip->gifts.size();
							const double base_score = wtl_cache[i];
							for (size_t j = 0; j <= n; ++j) {
								double inc_score = wtl_test_inserting(trip->gifts, trip->weight, *pc, j, rg->second) - base_score;
								if (best_improved[tid] > inc_score) {
									best_improved[tid] = inc_score;
									best_i[tid] = i;
									best_j[tid] = j;
								}
							}
						}
					}
					for (size_t i = 1; i < PROCS; ++i) {
						if (best_improved[0] > best_improved[i]) {
							best_improved[0] = best_improved[i];
							best_i[0] = best_i[i];
							best_j[0] = best_j[i];
						}
					}
					if (best_improved[0] < FLT_MAX) {
						size_t r2 = best_i[0];
						trip_t *trip = &trips[r2];
						backup.insert(std::make_pair(r2, backup_t(*trip, wtl_cache[r2], path_cache[r2], links[r2], link_keys[r2])));
						trip->gifts.insert(trip->gifts.begin() + best_j[0], rg->second);
						trip->weight += w;
						wtl_cache[r2] = trip_update_path_cache(*trip, path_cache[r2]);
						if (link_keys[rg->first].find(r2) == link_keys[rg->first].end()) {
							link_keys[rg->first].insert(r2);
							links[rg->first].push_back(r2);
							link_keys[r2].insert(rg->first);
							links[r2].push_back(rg->first);
						}
					} else {
						failback = true;
						break;
					}
				}
				if (failback) {
					for (auto b = backup.begin(); b != backup.end(); ++b) {
						trips[b->first] = b->second.trip;
						wtl_cache[b->first] = b->second.wtl_score;
						path_cache[b->first] = b->second.path_cache;
						links[b->first] = b->second.link;
						link_keys[b->first] = b->second.link_key;
					}
				} else {
					double score_diff = 0.0;
					for (auto b = backup.begin(); b != backup.end(); ++b) {
						score_diff += b->second.wtl_score - wtl_cache[b->first];
					}
					if (best_imporoved_trip < score_diff) {
						best_state.clear();
						best_imporoved_trip = score_diff;
						for (auto b = backup.begin(); b != backup.end(); ++b) {
							best_state.insert(std::make_pair(b->first, backup_t(trips[b->first], wtl_cache[b->first], path_cache[b->first], links[b->first], link_keys[b->first])));
						}
					}
					// restore
					for (auto b = backup.begin(); b != backup.end(); ++b) {
						trips[b->first] = b->second.trip;
						wtl_cache[b->first] = b->second.wtl_score;
						path_cache[b->first] = b->second.path_cache;
						links[b->first] = b->second.link;
						link_keys[b->first] = b->second.link_key;
					}
				}
			}
			if (best_state.size() > 0 &&
				best_imporoved_trip +
				 std::min<double>((stucked_count + 1) * 10000.0, 150000.0) > 0.0)
			{
				for (auto b = best_state.begin(); b != best_state.end(); ++b) {
					trips[b->first] = b->second.trip;
					wtl_cache[b->first] = b->second.wtl_score;
					path_cache[b->first] = b->second.path_cache;
					links[b->first] = b->second.link;
					link_keys[b->first] = b->second.link_key;
				}
				++destroy_trip_count;
			}
			trip_destroy = false; // disable
			if (stucked_count < 2) {
				trip_destroy_enable_epoch = epoch + 4;
			} else {
				trip_destroy_enable_epoch = epoch + 25;
			}
		} else {
			// remove & re-insert
			std::map<size_t, backup_t> backup;
			size_t num_trips = std::min<size_t>(std::max<size_t>(std::abs(gennorm(mt)) * remove_trip_scale, 1), remove_trip_scale * 2);
			std::set<size_t> trip_mutex;
			std::vector<std::pair<size_t, gift_t> > removed_gifts;
			double lp = genf(mt);
			bool neighbor_break = false;
			size_t base_trip = 0;

			if (epoch % 4 != 0 && link_initialized) {
				neighbor_break = true;
				do {
					base_trip = geni(mt);
				} while (trips[base_trip].gifts.size() == 0);
				std::shuffle(links[base_trip].begin(), links[base_trip].end(), mt);
			}

			// removing randomly selected gifts

			for (size_t k = 0; k < num_trips; ++k) {
				size_t r1;
				gifts_t r1_removed_gifts;

				if (neighbor_break) {
					if (k < links[base_trip].size()) {
						r1 = links[base_trip][k];
						if (trips[r1].gifts.size() == 0) {
							continue;
						}
					} else {
						break;
					}
				} else {
					do {
						r1 = geni(mt);
					} while (trips[r1].gifts.size() == 0 ||
							 trip_mutex.find(r1) != trip_mutex.end());
				}
				trip_mutex.insert(r1);
				backup.insert(std::make_pair(r1, backup_t(trips[r1], wtl_cache[r1], path_cache[r1], links[r1], link_keys[r1])));

				size_t num_gifts = std::min<size_t>(std::max<size_t>(std::abs(gennorm(mt)) * remove_gift_scale, 1),
													remove_gift_scale * 2);
				if (trips[r1].weight < 400.0 && stucked_count < 2) {
					num_gifts = (size_t)(num_gifts * 0.333);
				}
				lp = genf(mt);
				if (lp < 0.25) {
					std::set<int> gift_mutex;
					for (size_t i = 0; i < num_gifts; ++i) {
						size_t pos = (size_t)(genf(mt) * trips[r1].gifts.size());
						if (pos < trips[r1].gifts.size()) {
							if (gift_mutex.find(trips[r1].gifts[pos].id) == gift_mutex.end()) {
								r1_removed_gifts.push_back(trips[r1].gifts[pos]);
								trips[r1].weight -= trips[r1].gifts[pos].weight;
								gift_mutex.insert(trips[r1].gifts[pos].id);
								trips[r1].gifts.erase(trips[r1].gifts.begin() + pos);
							}
						}
					}
				} else {
					size_t pos = (size_t)(genf(mt) * trips[r1].gifts.size());
					for (size_t i = 0; i < num_gifts; ++i) {
						if (pos >= trips[r1].gifts.size()) {
							break;
						}
						r1_removed_gifts.push_back(trips[r1].gifts[pos]);
						trips[r1].weight -= trips[r1].gifts[pos].weight;
						trips[r1].gifts.erase(trips[r1].gifts.begin() + pos);
					}
				}
				for (auto gift = r1_removed_gifts.begin(); gift != r1_removed_gifts.end(); ++gift) {
					removed_gifts.push_back(std::make_pair(r1, *gift));
				}
				wtl_cache[r1] = trip_update_path_cache(trips[r1], path_cache[r1]);
			}
			if (removed_gifts.size() == 0) {
				continue;
			}
			if (stucked_count > 1 && genf(mt) < 0.5) {
				// worse but more variation
				std::shuffle(removed_gifts.begin(), removed_gifts.end(), mt);
			} else {
				// better
				std::sort(removed_gifts.begin(), removed_gifts.end(),
						  [](const std::pair<size_t, gift_t> &a, const std::pair<size_t, gift_t> &b) {
							  return a.second.lat < b.second.lat;
						  });
			}

			// re-inserting removed gifts with greedy algorithm

			bool neighbor_search;
			bool failback = false;

			lp = genf(mt);
			if (link_initialized) {
				neighbor_search = lp < (1.0 - (global_search_rate));
			} else {
				neighbor_search = false;
			}
			for (auto rg = removed_gifts.begin(); rg != removed_gifts.end(); ++rg) {
				size_t best_i[PROCS];
				size_t best_j[PROCS];
				double best_improved[PROCS];
				const double w = rg->second.weight;
				const auto &index = neighbor_search ? links[rg->first] : trip_index;

				std::fill_n(best_improved, PROCS, FLT_MAX);
#pragma omp parallel for schedule(dynamic, 1)
				for (size_t l = 0; l < index.size(); ++l) {
					const size_t i = index[l];
					const trip_t *trip = &trips[i];
					const std::vector<path_cache_t> *pc = &path_cache[i];
					if (trip->weight + w < 1000.0) {
						const int tid = omp_get_thread_num();
						const size_t n = trip->gifts.size();
						const double base_score = wtl_cache[i];
						for (size_t j = 0; j <= n; ++j) {
							double inc_score = wtl_test_inserting(trip->gifts, trip->weight, *pc, j, rg->second) - base_score;
							if (best_improved[tid] > inc_score) {
								best_improved[tid] = inc_score;
								best_i[tid] = i;
								best_j[tid] = j;
							}
						}
					}
				}
				for (size_t i = 1; i < PROCS; ++i) {
					if (best_improved[0] > best_improved[i]) {
						best_improved[0] = best_improved[i];
						best_i[0] = best_i[i];
						best_j[0] = best_j[i];
					}
				}
				if (best_improved[0] < FLT_MAX) {
					size_t r2 = best_i[0];
					trip_t *trip = &trips[r2];
					backup.insert(std::make_pair(r2, backup_t(*trip, wtl_cache[r2], path_cache[r2], links[r2], link_keys[r2])));
					trip->gifts.insert(trip->gifts.begin() + best_j[0], rg->second);
					trip->weight += w;
					wtl_cache[r2] = trip_update_path_cache(*trip, path_cache[r2]);
					if (!neighbor_search) {
						if (link_keys[rg->first].find(r2) == link_keys[rg->first].end()) {
							link_keys[rg->first].insert(r2);
							links[rg->first].push_back(r2);
							link_keys[r2].insert(rg->first);
							links[r2].push_back(rg->first);
						}
					}
				} else {
					failback = true;
					break;
				}
			}
			if (failback) {
				for (auto b = backup.begin(); b != backup.end(); ++b) {
					trips[b->first] = b->second.trip;
					wtl_cache[b->first] = b->second.wtl_score;
					path_cache[b->first] = b->second.path_cache;
					links[b->first] = b->second.link;
					link_keys[b->first] = b->second.link_key;
				}
			} else {
				double score_diff = 0.0;
				for (auto b = backup.begin(); b != backup.end(); ++b) {
					score_diff += b->second.wtl_score - wtl_cache[b->first];
				}
				if (score_diff > 1.0e-5) {
					// improved
					++update_count;
					if (neighbor_search) {
						++local_update_count;
					} else {
						++local_update_count2;
					}
				} else {
					// restore backup
					for (auto b = backup.begin(); b != backup.end(); ++b) {
						trips[b->first] = b->second.trip;
						wtl_cache[b->first] = b->second.wtl_score;
						path_cache[b->first] = b->second.path_cache;
						links[b->first] = b->second.link;
						link_keys[b->first] = b->second.link_key;
					}
				}
			}
		}
	}
}
