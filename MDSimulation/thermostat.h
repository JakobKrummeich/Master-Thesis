#ifndef THERMOSTAT_INCLUDED
#define THERMOSTAT_INCLUDED

#include <cmath>
#include "particles.h"
#include "../realRNG.h"

using namespace std;

class BussiThermostat {
	private:
		double Temperature;
		double InverseTimescale;
		int DegreesOfFreedom;
		realGaussianRNG GaussianRNG;

	public:
		BussiThermostat(double Temperature, double Timescale, int DegreesOfFreedom):
			Temperature(Temperature),
			InverseTimescale(1.0/Timescale),
			DegreesOfFreedom(DegreesOfFreedom),
			GaussianRNG(0.0,1.0)
		{
		}

		double computeRescalingFactor(const Particles& P, double Stepsize) {			
			double R1 = GaussianRNG.drawRandomNumber();
			double RandomSum = R1 * R1;
			for (int i = 2; i <= DegreesOfFreedom; i++){
				double NewR = GaussianRNG.drawRandomNumber();
				RandomSum += NewR * NewR;
			}
			double CurrentKineticEnergy = P.computeKineticEnergy();
			double ExponentialTerm = exp(-Stepsize*InverseTimescale);
			double TermInRoot = Temperature/(2.0*CurrentKineticEnergy) * (1.0 - ExponentialTerm);
			double AlphaSquared = ExponentialTerm + TermInRoot * RandomSum + 2.0 * R1 * exp(-0.5*Stepsize*InverseTimescale) * sqrt(TermInRoot);
			return sqrt(AlphaSquared);
		}
};

#endif //THERMOSTAT_INCLUDED
