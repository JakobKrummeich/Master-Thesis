#include <iostream>
#include <stdlib.h>
#include <string>
#include <chrono>
#include "../../realRNG.h"
#include "../particles.h"
#include "../stress_computator.h"
#include "../thermostat.h"

using namespace std;

int main(int argc, char* argv[]){
	const double Temperature = stod(argv[1]);
	string outputDirectory = argv[2];

	const int numberOfEquilibrationSweeps = 10000;
	const int numberOfDataTakingSweeps = 10000;
	const int numberOfAttemptedTypeSwitches = 100;

	const double Beta = 1.0/Temperature;

	const double ThermostatTime = 0.001;
	const double Stepsize = 0.001;

	realUniformRNG RNG;

	Particles P(500, DENSITY, Temperature, 0.0, 0.0);

	BussiThermostat BT(Temperature, ThermostatTime, DIMENSION*TOTAL_NUMBER_OF_PARTICLES);

	const auto StartTime = chrono::steady_clock::now();
	int NextUpdateTime = UPDATE_TIME_INTERVAL;

	cerr << "Equilibration started.\n";
	for (int StepNumber = 0; StepNumber < numberOfEquilibrationSweeps; StepNumber++){
		P.applyVelocityVerlet(Stepsize);
		double Alpha = BT.computeRescalingFactor(P, Stepsize);
		P.rescaleVelocities(Alpha);

		for (int typeChangeCounter = 0; typeChangeCounter < numberOfAttemptedTypeSwitches; typeChangeCounter++){
			P.runTypeChange(RNG, 0, TOTAL_NUMBER_OF_PARTICLES, Beta);
		}

		if (chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() >= NextUpdateTime){
			int Progress = numberOfEquilibrationSweeps >= 100 ? (StepNumber / (numberOfEquilibrationSweeps/100)) : StepNumber;
			cerr << "Progress: " << Progress << "%| StepNumber: " << StepNumber << "| Time passed: " <<  chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() << " s" << endl;
			NextUpdateTime += UPDATE_TIME_INTERVAL;
		}
	}

	ofstream FileStreamToWrite;
	string FileName = outputDirectory + "energySeries.dat";
	FileStreamToWrite.open(FileName);
	FileStreamToWrite << "U\t" << "T\t" << "H\n";
	FileStreamToWrite.close();

	int nextEnergyComputation = ENERGY_UPDATE_INTERVAL;

	cerr << "Data taking started.\n";

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
			string FileName = outputDirectory + "energySeries.dat";
			FileStreamToWrite.open(FileName, ios_base::app);
			FileStreamToWrite <<  fixed << setprecision(numeric_limits<long double>::digits10+1) << P.computePotentialEnergy() << '\t' << P.computeKineticEnergy() << '\t' <<  P.computeEnergy() << endl;
			FileStreamToWrite.close();
		}

		if (chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() >= NextUpdateTime){
			int Progress = numberOfDataTakingSweeps >= 100 ? (StepNumber / (numberOfDataTakingSweeps/100)) : StepNumber;
			cerr << "Progress: " << Progress << "%| StepNumber: " << StepNumber << "| Time passed: " <<  chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() << " s" << endl;
			NextUpdateTime += UPDATE_TIME_INTERVAL;
		}
	}

	FileName = outputDirectory + "histogram_NA.dat";
	FileStreamToWrite.open(FileName);
	FileStreamToWrite << "NA\t" << "#of appearances\n";
	for (int i = 0; i <= TOTAL_NUMBER_OF_PARTICLES; i++) {
		FileStreamToWrite << i << '\t' << histogramNA[i] << '\n';
	}
	FileStreamToWrite.close();

	FileName = outputDirectory + "final_state.dat";
	FileStreamToWrite.open(FileName);
	FileStreamToWrite << P;
	FileStreamToWrite.close();

	cerr << "#VerletListBuilds: " << P.getNumberOfVerletListBuilds() << endl;
	cerr << "Computation time: " << chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() << " s for " << numberOfDataTakingSweeps+numberOfEquilibrationSweeps << " MCSweeps." <<  endl << endl;
}

