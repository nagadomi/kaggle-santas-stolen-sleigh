#pragma once
#include <random>
#include "trip_def.hpp"

static void
tsp_hc(gifts_t &gifts, std::mt19937 &mt, size_t tol)
{
	std::uniform_int_distribution<size_t> genn(0, gifts.size() - 1);
	std::uniform_int_distribution<size_t> genm(1, 2);
	size_t no_updated = 0;

	double best_score = wtl(gifts);
	while (no_updated < tol) {
		++no_updated;
		gifts_t tmp = gifts;
		for (int j = 0; j < 2; ++j) {
			size_t r1 = genn(mt);
			size_t r2 = genn(mt);
			gift_t a = tmp[r1];
			tmp.erase(tmp.begin() + r1);
			tmp.insert(tmp.begin() + r2, a);
		}
		double score = wtl(tmp);
		if (score < best_score) {
			best_score = score;
			gifts = tmp;
			no_updated = 0;
		}
	}
}
