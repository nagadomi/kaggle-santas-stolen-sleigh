#include "load.hpp"
#include "trip.hpp"
#include "tsp.hpp"
#include <cstdio>
using namespace std;

int main(void)
{
	gifts_t gifts, bottom_gifts;
	load_gifts("data/gifts.csv", gifts);
	printf("%lu\n", gifts.size());
	std::mt19937 mt(71);
	trips_t trips[4];

	// initial solution
	printf("make initial solution..\n");
	trip_split(trips[0], gifts, mt, 100, 6.0, 975.0, 500.0);
	trip_split(trips[1], gifts, mt, 6.0, -48.0, 975.0, 500.0);
	trip_split(trips[2], gifts, mt, -48, -100.0, 975.0, 500.0);
	trip_merge(trips[0], trips[1]);
	trip_merge(trips[0], trips[2]);

#pragma omp parallel for schedule (dynamic, 1)
	for (size_t i = 0; i < trips[0].size(); ++i) {
		std::mt19937 lmt(71);
		tsp_hc(trips[0][i].gifts, lmt, 200000);
	}
	printf("init: %f\n", trip_lb_score(trips[0]));

	// iterative improvement
	trip_optimize(trips[0], "out.csv"); 	// never return

	return 0;
}
