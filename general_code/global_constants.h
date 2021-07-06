#ifndef GLOBAL_VARIABLES_INCLUDED
#define GLOBAL_VARIABLES_INCLUDED

#include <cmath>

using namespace std;

constexpr int DIMENSION = 2;

constexpr double CUTOFF = 2.5;
constexpr double CUTOFF_SQUARED = CUTOFF * CUTOFF;
constexpr double INVERSE_CUTOFF = 1.0/CUTOFF;
constexpr double POTENTIAL_CONSTANT_1 = pow(INVERSE_CUTOFF, 6.0) - pow(INVERSE_CUTOFF, 12.0);
constexpr double POTENTIAL_CONSTANT_2 = 6.0 * pow(INVERSE_CUTOFF, 7.0) - 12.0 * pow(INVERSE_CUTOFF, 13.0);

constexpr double AA_INTERACTION_STRENGTH = 1.0;
constexpr double AB_INTERACTION_STRENGTH = 0.1;

#endif //GLOBAL_VARIABLES_INCLUDED
