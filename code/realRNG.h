#ifndef REAL_RNG_INCLUDED
#define REAL_RNG_INCLUDED

#include <random>

using namespace std;

class realUniformRNG{
	private:
		mt19937 rng;
		uniform_real_distribution<double> UnifRealDist;
	public:
		realUniformRNG():
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

class realGaussianRNG{
	private:
		mt19937 rng;
		normal_distribution<double> GaussianRealDist;
	public:
		realGaussianRNG(double Mean, double StandardDeviation):
			GaussianRealDist(Mean, StandardDeviation){
				random_device rd;
				seed_seq sd{rd(),rd()};
				rng = mt19937(sd);
		}

		double drawRandomNumber(){
			return GaussianRealDist(rng);
		}
};
#endif /* REAL_RNG_INCLUDED */

