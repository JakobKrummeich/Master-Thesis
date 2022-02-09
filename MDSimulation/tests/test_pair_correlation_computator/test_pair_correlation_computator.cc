#include <iostream>
#include <stdlib.h>
#include <string>
#include <chrono>
#include "../../../realRNG.h"
#include "../../particles.h"
#include "../../pair_correlation_computator.h"
#include "../../thermostat.h"

using namespace std;

int main(int argc, char* argv[]){
	int MaxNumberOfEquilibrationSweeps = 10000;
	int MaxNumberOfShearSweeps = 15000;

	int NumberOfyValues = 15;
	int NumberOfStressSubdivisions = 14;
	const double Temperature = 2.0;
	const double Beta = 1.0/Temperature;
	double ShearRate = 0.04;
	double ThermostatTime = 0.001;
	const double Stepsize = 0.001;
	const int numberOfPositions = 1000;

	realUniformRNG RNG;

	PairCorrelationComputator PCC(100,5.0,sqrt(static_cast<double>(TOTAL_NUMBER_OF_PARTICLES) / DENSITY));
	BussiThermostat BT(Temperature, ThermostatTime, DIMENSION*TOTAL_NUMBER_OF_PARTICLES);

	for (int positionsIndex = 0; positionsIndex < numberOfPositions; positionsIndex++){

		cerr << "position " << positionsIndex << endl;

		const auto StartTime = chrono::steady_clock::now();

		Particles P(500, DENSITY, Temperature, 0.0, 0.0);

		int NextUpdateTime = UPDATE_TIME_INTERVAL;

		cerr << "Equilibration started.\n";

		for (int StepNumber = 0; StepNumber < MaxNumberOfEquilibrationSweeps; StepNumber++){
			P.applyVelocityVerlet(Stepsize);
			double Alpha = BT.computeRescalingFactor(P, Stepsize);
			P.rescaleVelocities(Alpha);

			if (chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() >= NextUpdateTime){
				int Progress = MaxNumberOfEquilibrationSweeps >= 100 ? (StepNumber / (MaxNumberOfEquilibrationSweeps/100)) : StepNumber;
				cerr << "Progress: " << Progress << "%| StepNumber: " << StepNumber << "| Time passed: " <<  chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() << " s" << endl;
				NextUpdateTime += UPDATE_TIME_INTERVAL;
			}
		}

		cerr << "Equilibration finished. Shearing with rate= " << ShearRate << " enabled.\n";
		P.setShearRate(ShearRate);
		for (int StepNumber = 0; StepNumber < MaxNumberOfShearSweeps; StepNumber++){
			P.applyVelocityVerlet(Stepsize);
			double Alpha = BT.computeRescalingFactor(P, Stepsize);
			P.rescaleVelocities(Alpha);

			P.moveImageBoxes(Stepsize);

			if (chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() >= NextUpdateTime){
				int Progress = MaxNumberOfShearSweeps >= 100 ? (StepNumber / (MaxNumberOfShearSweeps/100)) : StepNumber;
				cerr << "Progress: " << Progress << "%| StepNumber: " << StepNumber << "| Time passed: " <<  chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() << " s" << endl;
				NextUpdateTime += UPDATE_TIME_INTERVAL;
			}
		}
		PCC.computeImg22(P);
		cerr << "#VerletListBuilds: " << P.getNumberOfVerletListBuilds() << endl;
		cerr << "Computation time: " << chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() << " s for " << MaxNumberOfShearSweeps+MaxNumberOfEquilibrationSweeps << " MCSweeps." <<  endl << endl;
	}

	PCC.writeResults("./");
}

