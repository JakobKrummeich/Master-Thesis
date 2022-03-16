#include <iostream>
#include <stdlib.h>
#include <string>
#include <chrono>
#include <thread>
#include "../../realRNG.h"
#include "../particles.h"
#include "../stress_computator.h"
#include "../pair_correlation_computator.h"
#include "../thermostat.h"

void writeAvgVelocityFile(string outputDirectory, const vector<double>& avgVelocities, int numberOfDataTakingSweeps, double boxLength, int numberOfyValuesForVelocities){
	ofstream ofs(outputDirectory + "avgVelocities.dat");

	double yDelta = 1.0/static_cast<double>(numberOfyValuesForVelocities);
	double Currenty = yDelta*0.5;
	for (int i = 0; i < numberOfyValuesForVelocities-1; i++){
		ofs << Currenty*boxLength << '\t';
		Currenty += yDelta;
	}
	ofs << Currenty*boxLength << '\n';
	for (int i = 0; i < DIMENSION; i++){
		for (int j = 0; j < numberOfyValuesForVelocities-1; j++){
			ofs << avgVelocities[numberOfyValuesForVelocities*i+j]*boxLength/(static_cast<double>(numberOfDataTakingSweeps)) << '\t';
		}
		ofs << avgVelocities[numberOfyValuesForVelocities*i+numberOfyValuesForVelocities-1]*boxLength/(static_cast<double>(numberOfDataTakingSweeps)) << '\n';
	}
}

using namespace std;

int main(int argc, char* argv[]){
	const double Temperature = stod(argv[1]);
	const double shearRate = stod(argv[2]);
	string initialStateFile = argv[3];
	string outputDirectory = argv[4];

	const int numberOfEquilibrationSweeps = stoi(argv[5]);
	const int numberOfDataTakingSweeps = stoi(argv[6]);
	const double stepsize = stod(argv[7]);


	const int numberOfAttemptedTypeSwitches = 100;

	const int numberOfyValuesForVelocities = 15;

	const double Beta = 1.0/Temperature;

	const double ThermostatTime = 1.0;


	realUniformRNG RNG;

	Particles P(shearRate);
	P.readInParticleState(initialStateFile, 0, 0, DENSITY);

	StressComputator SC(P,14,(numberOfDataTakingSweeps+numberOfEquilibrationSweeps)/ENERGY_UPDATE_INTERVAL,1);
	PairCorrelationComputator PCC(200,5.0,sqrt(static_cast<double>(TOTAL_NUMBER_OF_PARTICLES) / DENSITY));	

	BussiThermostat BT(Temperature, ThermostatTime, DIMENSION*TOTAL_NUMBER_OF_PARTICLES);

	vector<double> avgVelocities(DIMENSION*numberOfyValuesForVelocities, 0.0);
	vector<double> energySeries;

	const auto StartTime = chrono::steady_clock::now();
	int nextEnergyComputation = ENERGY_UPDATE_INTERVAL;

	SC.computeStresses(0);

	cerr << "Equilibration started.\n";
	for (int StepNumber = 0; StepNumber < numberOfEquilibrationSweeps; StepNumber++){
		P.applyVelocityVerlet(stepsize);
		double Alpha = BT.computeRescalingFactor(P, stepsize);
		P.rescaleVelocities(Alpha);

		for (int typeChangeCounter = 0; typeChangeCounter < numberOfAttemptedTypeSwitches; typeChangeCounter++){
			P.runTypeChange(RNG, 0, TOTAL_NUMBER_OF_PARTICLES, Beta);
		}

		if (StepNumber == nextEnergyComputation){
			SC.computeStresses(StepNumber/ENERGY_UPDATE_INTERVAL);
			nextEnergyComputation += ENERGY_UPDATE_INTERVAL;

			energySeries.push_back(P.computePotentialEnergy());
			energySeries.push_back(P.computeKineticEnergy());
			energySeries.push_back(P.computeEnergy());
		}
	}

	SC.computeStresses(numberOfEquilibrationSweeps/ENERGY_UPDATE_INTERVAL);

	cerr << "Data taking started.\n";

	nextEnergyComputation = ENERGY_UPDATE_INTERVAL;

	uint32_t histogramNA [TOTAL_NUMBER_OF_PARTICLES+1]{};

	for (int StepNumber = 0; StepNumber < numberOfDataTakingSweeps; StepNumber++){
		P.applyVelocityVerlet(stepsize);
		double Alpha = BT.computeRescalingFactor(P, stepsize);
		P.rescaleVelocities(Alpha);

		for (int typeChangeCounter = 0; typeChangeCounter < numberOfAttemptedTypeSwitches; typeChangeCounter++){
			P.runTypeChange(RNG, 0, TOTAL_NUMBER_OF_PARTICLES, Beta);
		}
		int numberOfAParticles = P.getNumberOfAParticles();
		histogramNA[numberOfAParticles]++;

		P.updateAverageVelocities(avgVelocities, numberOfyValuesForVelocities);

		if (StepNumber == nextEnergyComputation){
			SC.computeStresses((StepNumber+numberOfEquilibrationSweeps)/ENERGY_UPDATE_INTERVAL);
			nextEnergyComputation += ENERGY_UPDATE_INTERVAL;

			energySeries.push_back(P.computePotentialEnergy());
			energySeries.push_back(P.computeKineticEnergy());
			energySeries.push_back(P.computeEnergy());
		}
	}
	PCC.computeImg22(P);

	this_thread::sleep_for(chrono::seconds(static_cast<int>(RNG.drawRandomNumber(0.0,600.0))));

	PCC.writeResults(outputDirectory);
	writeAvgVelocityFile(outputDirectory, avgVelocities, numberOfDataTakingSweeps, P.getBoxLength(), numberOfyValuesForVelocities);
	writeEnergySeriesFile(outputDirectory, energySeries);
	writeHistogramFile(outputDirectory, histogramNA, TOTAL_NUMBER_OF_PARTICLES);
	writeFinalStateFile(outputDirectory, P);
	SC.writeAverageStresses(outputDirectory + "avgStresses");


	cerr << "#VerletListBuilds: " << P.getNumberOfVerletListBuilds() << endl;
	cerr << "Computation time: " << chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() << " s for " << numberOfDataTakingSweeps+numberOfEquilibrationSweeps << " MCSweeps." <<  endl << endl;
}

