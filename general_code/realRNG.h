#ifndef REAL_RNG_INCLUDED
#define REAL_RNG_INCLUDED

#include <random>

using namespace std;

class realRNG{
	private:
		std::mt19937 rng;
		std::uniform_real_distribution<double> UnifRealDist;
	public:
		realRNG():
			UnifRealDist(0.0,1.0){
				random_device rd;
				seed_seq sd{rd(),rd()};
				rng = mt19937(sd);
		}

		double drawRandomNumber(){
			return UnifRealDist(rng);
		}

		double drawRandomNumber(double Min, double Max){
			return (Max - Min) * UnifRealDist(rng) + Min;
		}
};
#endif /* REAL_RNG_INCLUDED */

