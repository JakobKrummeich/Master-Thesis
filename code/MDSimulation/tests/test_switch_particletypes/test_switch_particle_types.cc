#include <iostream>
#include <stdlib.h>
#include <string>
#include <chrono>
#include "../../../realRNG.h"
#include "../../particles.h"
#include "../../stress_computator.h"
#include "../../thermostat.h"

using namespace std;

int main(int argc, char* argv[]){
	int MaxNumberOfEquilibrationSweeps = 10000;
	int MaxNumberOfShearSweeps = 300000;

	int NumberOfyValues = 15;
	int NumberOfStressSubdivisions = 14;
	const double Temperature = 0.7;
	const double Beta = 1.0/Temperature;
	double ShearRate = 0.15;
	double ThermostatTime = 0.001;
	const double Stepsize = 0.001;

	realUniformRNG RNG;

	Particles P(500, DENSITY, Temperature, 0.0, 0.0);
	StressComputator SC(P,NumberOfStressSubdivisions, MaxNumberOfShearSweeps, 1);
	BussiThermostat BT(Temperature, ThermostatTime, DIMENSION*TOTAL_NUMBER_OF_PARTICLES);

	const auto StartTime = chrono::steady_clock::now();
	int NextUpdateTime = UPDATE_TIME_INTERVAL;

	ofstream FileStreamToWrite;
	cerr << "Equilibration started.\n";
	string FileName = OUTPUT_DIRECTORY+"_EnergySeries.dat";
	FileStreamToWrite.open(FileName);
	FileStreamToWrite << "U\t" << "T\t" << "H\n";
	for (int StepNumber = 0; StepNumber < MaxNumberOfEquilibrationSweeps; StepNumber++){
		P.applyVelocityVerlet(Stepsize);
		double Alpha = BT.computeRescalingFactor(P, Stepsize);
		P.rescaleVelocities(Alpha);
		//AvgMSD.push_back(P.computeAverageMSD());

		FileStreamToWrite <<  fixed << setprecision(numeric_limits<long double>::digits10+1) << P.computePotentialEnergy() << '\t' << P.computeKineticEnergy() << '\t' <<  P.computeEnergy() << endl;
		if (chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() >= NextUpdateTime){
			int Progress = MaxNumberOfEquilibrationSweeps >= 100 ? (StepNumber / (MaxNumberOfEquilibrationSweeps/100)) : StepNumber;
			cerr << "Progress: " << Progress << "%| StepNumber: " << StepNumber << "| Time passed: " <<  chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() << " s" << endl;
			NextUpdateTime += UPDATE_TIME_INTERVAL;
		}
	}

	vector<int> NumberOfASeries;
	cerr << "Equilibration finished. Shearing with rate= " << ShearRate << " enabled.\n";
	P.setShearRate(ShearRate);
	for (int StepNumber = 0; StepNumber < MaxNumberOfShearSweeps; StepNumber++){
		P.applyVelocityVerlet(Stepsize);
		double Alpha = BT.computeRescalingFactor(P, Stepsize);
		P.rescaleVelocities(Alpha);

		SC.computeStresses(StepNumber);
		//AvgMSD.push_back(P.computeAverageMSD());

		P.moveImageBoxes(Stepsize);

		P.runTypeChange(RNG, 0, TOTAL_NUMBER_OF_PARTICLES, Beta);
		NumberOfASeries.push_back(P.getNumberOfAParticles());

		FileStreamToWrite <<  fixed << setprecision(numeric_limits<long double>::digits10+1) << P.computePotentialEnergy() << '\t' << P.computeKineticEnergy() << '\t' <<  P.computeEnergy() << endl;
		if (chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() >= NextUpdateTime){
			int Progress = MaxNumberOfShearSweeps >= 100 ? (StepNumber / (MaxNumberOfShearSweeps/100)) : StepNumber;
			cerr << "Progress: " << Progress << "%| StepNumber: " << StepNumber << "| Time passed: " <<  chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() << " s" << endl;
			NextUpdateTime += UPDATE_TIME_INTERVAL;
		}
	}
	FileStreamToWrite.close();

	FileStreamToWrite.open("NANumberSeries.dat");
	for (int i = 0; i < NumberOfASeries.size(); i++){
		FileStreamToWrite << NumberOfASeries[i] << '\n';
	}
	FileStreamToWrite.close();

	cerr << "#VerletListBuilds: " << P.getNumberOfVerletListBuilds() << endl;
	cerr << "Computation time: " << chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() << " s for " << MaxNumberOfShearSweeps+MaxNumberOfEquilibrationSweeps << " MCSweeps." <<  endl << endl;

	SC.writeAverageStresses("AvgStresses");
}

