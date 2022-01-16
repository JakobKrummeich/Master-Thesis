#include <iostream>
#include <stdlib.h>
#include <string>
#include <chrono>
#include "../particles.h"
#include "../stress_computator.h"
#include "../thermostat.h"

using namespace std;

int main(int argc, char* argv[]){
	const int MaxNumberOfEquilibrationSweeps = 1000;
	const int NumberOfDataTakingSweeps = 1000;

	const double ShearRate = atof(argv[1]);
	const int NumberOfShearEquilibrationSweeps = atoi(argv[2]);
	const string DataPath = argv[3];

	const int NumberOfStressSubdivisions = 14;
	const double Temperature = 2.0;
	const double ThermostatTime = 0.05;
	const double Stepsize = 0.001;

	const int NumberOfPositions = 100;

	Particles P(500, DENSITY, Temperature, 0.0, 0.0);
	StressComputator SC(P,NumberOfStressSubdivisions, NumberOfDataTakingSweeps, NumberOfPositions);
	BussiThermostat BT(Temperature, ThermostatTime, DIMENSION*TOTAL_NUMBER_OF_PARTICLES);

	for (int PositionsIndex = 0; PositionsIndex < NumberOfPositions; PositionsIndex++){

		cerr << "Position " << PositionsIndex << endl;

		P = Particles(500, DENSITY, Temperature, 0.0, 0.0);

		const auto StartTime = chrono::steady_clock::now();
		int NextUpdateTime = UPDATE_TIME_INTERVAL;

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

		cerr << "Initial equilibration finished. Shearing with rate= " << ShearRate << " enabled.\n";
		P.setShearRate(ShearRate);
		for (int StepNumber = 0; StepNumber < NumberOfShearEquilibrationSweeps; StepNumber++){
			P.applyVelocityVerlet(Stepsize);
			double Alpha = BT.computeRescalingFactor(P, Stepsize);
			P.rescaleVelocities(Alpha);
			P.moveImageBoxes(Stepsize);

			if (chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() >= NextUpdateTime){
				int Progress = NumberOfShearEquilibrationSweeps >= 100 ? (StepNumber / (NumberOfShearEquilibrationSweeps/100)) : StepNumber;
				cerr << "Progress: " << Progress << "%| StepNumber: " << StepNumber << "| Time passed: " <<  chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() << " s" << endl;
				NextUpdateTime += UPDATE_TIME_INTERVAL;
			}
		}
		cerr << "Shearing equilibration finished. Data taking enabled.\n";
		for (int StepNumber = 0; StepNumber < NumberOfDataTakingSweeps; StepNumber++){
			P.applyVelocityVerlet(Stepsize);
			double Alpha = BT.computeRescalingFactor(P, Stepsize);
			P.rescaleVelocities(Alpha);

			SC.computeStresses(StepNumber);

			P.moveImageBoxes(Stepsize);

			if (chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() >= NextUpdateTime){
				int Progress = NumberOfDataTakingSweeps >= 100 ? (StepNumber / (NumberOfDataTakingSweeps/100)) : StepNumber;
				cerr << "Progress: " << Progress << "%| StepNumber: " << StepNumber << "| Time passed: " <<  chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() << " s" << endl;
				NextUpdateTime += UPDATE_TIME_INTERVAL;
			}
		}
		cerr << "#VerletListBuilds: " << P.getNumberOfVerletListBuilds() << endl;
		cerr << "Computation time: " << chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() << " s for " << NumberOfShearEquilibrationSweeps+MaxNumberOfEquilibrationSweeps+NumberOfDataTakingSweeps << " MCSweeps." <<  endl << endl;
	}

	SC.writeAverageStresses(DataPath+"AvgStresses_T="+to_string(Temperature)+"_shear_rate="+to_string(ShearRate));
}

