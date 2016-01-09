GXX=g++
CXXFLAGS=-std=c++11 -fopenmp -march=native -Wno-unused-function -Ofast -g0 -Wall -DNDEBUG
#CXXFLAGS=-std=c++11 -fopenmp -Og -g -Wall -D_GLIBCXX_DEBUG

all: santa

clean:
	rm -fr santa

santa: santa.cpp tick.hpp trip.hpp tsp.hpp load.hpp cost.hpp trip_def.hpp
	$(GXX) santa.cpp -o santa $(CXXFLAGS)
