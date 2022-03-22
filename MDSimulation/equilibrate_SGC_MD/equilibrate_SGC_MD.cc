#include <iostream>
#include <stdlib.h>
#include <string>
#include <chrono>
#include <thread>
#include "../../realRNG.h"
#include "../particles.h"
#include "../stress_computator.h"
#include "../thermostat.h"

using namespace std;

int main(int argc, char* argv[]){
	const double Temperature = stod(argv[1]);
	string outputDirectory = argv[2];

	const int numberOfEquilibrationSweeps = stoi(argv[3]);
	const int numberOfDataTakingSweeps = stoi(argv[4]);
	const int numberOfAttemptedTypeSwitches = ceil(static_cast<double>(TOTAL_NUMBER_OF_PARTICLES)*0.1);

	const double Beta = 1.0/Temperature;

	const double ThermostatTime = 0.001;
	const double Stepsize = 0.0005;

	realUniformRNG RNG;

	Particles P(static_cast<int>(static_cast<double>(TOTAL_NUMBER_OF_PARTICLES)*0.5), DENSITY, Temperature, 0.0, 0.0);

	BussiThermostat BT(Temperature, ThermostatTime, TOTAL_NUMBER_OF_PARTICLES);

	vector<double> energySeries;

	const auto StartTime = chrono::steady_clock::now();
	int nextEnergyComputation = ENERGY_UPDATE_INTERVAL;

	cerr << "Equilibration started.\n";

	for (int StepNumber = 0; StepNumber < numberOfEquilibrationSweeps; StepNumber++){
		P.applyVelocityVerlet(Stepsize);
		double Alpha = BT.computeRescalingFactor(P, Stepsize);
		P.rescaleVelocities(Alpha);

		for (int typeChangeCounter = 0; typeChangeCounter < numberOfAttemptedTypeSwitches; typeChangeCounter++){
			P.runTypeChange(RNG, 0, TOTAL_NUMBER_OF_PARTICLES, Beta);
		}

		if (StepNumber == nextEnergyComputation){
			nextEnergyComputation += ENERGY_UPDATE_INTERVAL;

			energySeries.push_back(P.computePotentialEnergy());
			energySeries.push_back(P.computeKineticEnergy());
			energySeries.push_back(P.computeEnergy());
		}
	}

	cerr << "Data taking started.\n";

	nextEnergyComputation = ENERGY_UPDATE_INTERVAL;

	uint32_t histogramNA [TOTAL_NUMBER_OF_PARTICLES+1]{};

	for (int StepNumber = 0; StepNumber < numberOfDataTakingSweeps; StepNumber++){
		P.applyVelocityVerlet(Stepsize);
		double Alpha = BT.computeRescalingFactor(P, Stepsize);
		P.rescaleVelocities(Alpha);

		for (int typeChangeCounter = 0; typeChangeCounter < numberOfAttemptedTypeSwitches; typeChangeCounter++){
			P.runTypeChange(RNG, 0, TOTAL_NUMBER_OF_PARTICLES, Beta);
		}

		int numberOfAParticles = P.getNumberOfAParticles();
		histogramNA[numberOfAParticles]++;

		if (StepNumber == nextEnergyComputation){
			nextEnergyComputation += ENERGY_UPDATE_INTERVAL;

			energySeries.push_back(P.computePotentialEnergy());
			energySeries.push_back(P.computeKineticEnergy());
			energySeries.push_back(P.computeEnergy());
		}
	}

	//this_thread::sleep_for(chrono::seconds(static_cast<int>(RNG.drawRandomNumber(0.0,600.0))));

	writeEnergySeriesFile(outputDirectory, energySeries);
	writeHistogramFile(outputDirectory, histogramNA, TOTAL_NUMBER_OF_PARTICLES);
	writeFinalStateFile(outputDirectory, P);

	cerr << "#VerletListBuilds: " << P.getNumberOfVerletListBuilds() << endl;
	cerr << "Computation time: " << chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() << " s for " << numberOfDataTakingSweeps+numberOfEquilibrationSweeps << " MCSweeps." <<  endl << endl;
}

