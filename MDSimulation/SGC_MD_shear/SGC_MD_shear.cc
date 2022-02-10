#include <iostream>
#include <stdlib.h>
#include <string>
#include <chrono>
#include "../../realRNG.h"
#include "../particles.h"
#include "../stress_computator.h"
#include "../pair_correlation_computator.h"
#include "../thermostat.h"

using namespace std;

int main(int argc, char* argv[]){
	const double Temperature = stod(argv[1]);
	string initialStateFile = argv[2];
	string outputDirectory = argv[3];

	const int numberOfEquilibrationSweeps = 1000000;
	const int numberOfDataTakingSweeps = 1000000;
	const int numberOfAttemptedTypeSwitches = 100;

	const int numberOfyValuesForVelocities = 15;

	const double Beta = 1.0/Temperature;

	const double ThermostatTime = 0.001;
	const double Stepsize = 0.001;
	const double shearRate = 0.04;

	realUniformRNG RNG;

	Particles P(shearRate);
	P.readInParticleState(initialStateFile, 0, 0, DENSITY);

	StressComputator SC(P,14,numberOfDataTakingSweeps+numberOfEquilibrationSweeps,1);
	PairCorrelationComputator PCC(200,5.0,sqrt(static_cast<double>(TOTAL_NUMBER_OF_PARTICLES) / DENSITY));	

	BussiThermostat BT(Temperature, ThermostatTime, DIMENSION*TOTAL_NUMBER_OF_PARTICLES);

	vector<double> AvgVelocities(DIMENSION*numberOfyValuesForVelocities, 0.0);

	ofstream FileStreamToWrite;
	string FileName = outputDirectory + "avgVelocities.dat";
	FileStreamToWrite.open(FileName);
	double yDelta = 1.0/static_cast<double>(numberOfyValuesForVelocities);
	double Currenty = yDelta*0.5;
	for (int i = 0; i < numberOfyValuesForVelocities-1; i++){
		FileStreamToWrite << Currenty*P.getBoxLength() << '\t';
		Currenty += yDelta;
	}
	FileStreamToWrite << Currenty*P.getBoxLength() << '\n';
	FileStreamToWrite.close();

	const auto StartTime = chrono::steady_clock::now();
	int NextUpdateTime = UPDATE_TIME_INTERVAL;

	FileName = outputDirectory + "energySeries.dat";
	FileStreamToWrite.open(FileName);
	FileStreamToWrite << "U\t" << "T\t" << "H\n";
	FileStreamToWrite.close();

	int nextEnergyComputation = ENERGY_UPDATE_INTERVAL;

	cerr << "Equilibration started.\n";
	for (int StepNumber = 0; StepNumber < numberOfEquilibrationSweeps; StepNumber++){
		P.applyVelocityVerlet(Stepsize);
		double Alpha = BT.computeRescalingFactor(P, Stepsize);
		P.rescaleVelocities(Alpha);

		for (int typeChangeCounter = 0; typeChangeCounter < numberOfAttemptedTypeSwitches; typeChangeCounter++){
			P.runTypeChange(RNG, 0, TOTAL_NUMBER_OF_PARTICLES, Beta);
		}

		SC.computeStresses(StepNumber);

		if (StepNumber == nextEnergyComputation){
			nextEnergyComputation += ENERGY_UPDATE_INTERVAL;
			string FileName = outputDirectory + "energySeries.dat";
			FileStreamToWrite.open(FileName, ios_base::app);
			FileStreamToWrite <<  fixed << setprecision(numeric_limits<long double>::digits10+1) << P.computePotentialEnergy() << '\t' << P.computeKineticEnergy() << '\t' <<  P.computeEnergy() << endl;
			FileStreamToWrite.close();
		}

		if (chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() >= NextUpdateTime){
			int Progress = numberOfEquilibrationSweeps >= 100 ? (StepNumber / (numberOfEquilibrationSweeps/100)) : StepNumber;
			cerr << "Progress: " << Progress << "%| StepNumber: " << StepNumber << "| Time passed: " <<  chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() << " s" << endl;
			NextUpdateTime += UPDATE_TIME_INTERVAL;
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

		P.updateAverageVelocities(AvgVelocities, numberOfyValuesForVelocities);
		
		SC.computeStresses(StepNumber+numberOfEquilibrationSweeps);

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

	PCC.computeImg22(P);
	PCC.writeResults(outputDirectory);

	FileName = outputDirectory + "histogram_NA.dat";
	FileStreamToWrite.open(FileName);
	FileStreamToWrite << "NA\t" << "#of appearances\n";
	for (int i = 0; i <= TOTAL_NUMBER_OF_PARTICLES; i++) {
		FileStreamToWrite << i << '\t' << histogramNA[i] << '\n';
	}
	FileStreamToWrite.close();

	FileName = outputDirectory + "N=1000_final_state.dat";
	FileStreamToWrite.open(FileName);
	FileStreamToWrite << P;
	FileStreamToWrite.close();

	SC.writeAverageStresses(outputDirectory + "avgStresses");

	FileName = outputDirectory + "avgVelocities.dat";
	FileStreamToWrite.open(FileName, ios_base::app);
	for (int i = 0; i < DIMENSION; i++){
		for (int j = 0; j < numberOfyValuesForVelocities-1; j++){
			FileStreamToWrite << AvgVelocities[numberOfyValuesForVelocities*i+j]*P.getBoxLength()/(static_cast<double>(numberOfDataTakingSweeps)) << '\t';
		}
		FileStreamToWrite << AvgVelocities[numberOfyValuesForVelocities*i+numberOfyValuesForVelocities-1]*P.getBoxLength()/(static_cast<double>(numberOfDataTakingSweeps)) << '\n';
	}
	FileStreamToWrite.close();

	cerr << "#VerletListBuilds: " << P.getNumberOfVerletListBuilds() << endl;
	cerr << "Computation time: " << chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() << " s for " << numberOfDataTakingSweeps+numberOfEquilibrationSweeps << " MCSweeps." <<  endl << endl;
}

