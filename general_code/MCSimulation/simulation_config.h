#ifndef SIMULATION_CONFIG_INCLUDED
#define SIMULATION_CONFIG_INCLUDED

#include <string>
#include "global_constants.h"

using namespace std;

constexpr int TOTAL_NUMBER_OF_PARTICLES = 1000;

constexpr double DENSITY = 0.75;

const string OUTPUT_DIRECTORY = "../../data/sgcmc_shear/fresh_data/N="+to_string(TOTAL_NUMBER_OF_PARTICLES);;

constexpr int NUMBER_OF_RUNS = 1;
constexpr int NUMBER_OF_THREADS = 1;

constexpr int UPDATE_TIME_INTERVAL = 60;
constexpr int POT_ENERGY_UPDATE_INTERVAL = 200;

constexpr int MAX_RUNTIME_IN_MINUTES = 2850;

constexpr int NUMBER_OF_SAVED_STATES_PER_TEMPERATURE = 1;

constexpr double MAXIMUM_DISPLACEMENT = 0.04;
constexpr double DISPLACEMENT_PROBABILITY = 0.9;

constexpr double MAXIMUM_VOLUME_CHANGE = 0.6;
constexpr double VOLUME_CHANGE_PROBABILITY = 0.01;

constexpr double MAX_VERLET_DIST = 1.1*CUTOFF;
constexpr double MAX_VERLET_DIST_SQUARED = MAX_VERLET_DIST * MAX_VERLET_DIST;
constexpr double SKINDISTANCE = MAX_VERLET_DIST - CUTOFF;

#endif // SIMULATION_CONFIG_INCLUDED
