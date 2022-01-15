#include <iostream>
#include <stdlib.h>
#include <string>
#include <chrono>
#include "../particles.h"
#include "../stress_computator.h"
#include "../thermostat.h"

using namespace std;

int main(int argc, char* argv[]){
	int MaxNumberOfEquilibrationSweeps = 100000;
	int MaxNumberOfShearSweeps = 10000;

	int NumberOfyValues = 15;
	int NumberOfStressSubdivisions = 14;
	double Temperature = 2.0;
	double ShearRate = 0.1;
	double ThermostatTime = 0.1;
	const double Stepsize = 0.0004;

	Particles P(500, DENSITY, Temperature, 0.0, 0.0);
	//Particles P(ShearRate);
	//P.readInParticleState("FinalState_long_equilibration_N=1000.dat", 0, 0, DENSITY);

	StressComputator SC(P,NumberOfStressSubdivisions);
	BussiThermostat BT(Temperature, ThermostatTime, DIMENSION*TOTAL_NUMBER_OF_PARTICLES);

	const auto StartTime = chrono::steady_clock::now();
	int NextUpdateTime = UPDATE_TIME_INTERVAL;

	vector<double> AvgVelocities;
	vector<double> AvgMSD;

	ofstream FileStreamToWrite;
	string FileName = "AvgVelocities.dat";
	FileStreamToWrite.open(FileName);
	FileStreamToWrite << endl;
	double yDelta = 1.0/static_cast<double>(NumberOfyValues);
	double Currenty = yDelta*0.5;
	for (int i = 0; i < NumberOfyValues-1; i++){
		FileStreamToWrite << Currenty*P.getBoxLength() << '\t';
		Currenty += yDelta;
	}
	FileStreamToWrite << Currenty*P.getBoxLength() << '\n';
	FileStreamToWrite.close();

	cerr << "Equilibration started.\n";
	FileName = OUTPUT_DIRECTORY+"_EnergySeries.dat";
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

	cerr << "Equilibration finished. Shearing with rate= " << ShearRate << " enabled.\n";
	P.setShearRate(ShearRate);
	for (int StepNumber = 0; StepNumber < MaxNumberOfShearSweeps; StepNumber++){
		P.applyVelocityVerlet(Stepsize);
		double Alpha = BT.computeRescalingFactor(P, Stepsize);
		P.rescaleVelocities(Alpha);

		P.updateAverageVelocities(AvgVelocities, NumberOfyValues);
		SC.computeStresses();
		//AvgMSD.push_back(P.computeAverageMSD());

		P.moveImageBoxes(Stepsize);

		FileStreamToWrite <<  fixed << setprecision(numeric_limits<long double>::digits10+1) << P.computePotentialEnergy() << '\t' << P.computeKineticEnergy() << '\t' <<  P.computeEnergy() << endl;
		if (chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() >= NextUpdateTime){
			int Progress = MaxNumberOfShearSweeps >= 100 ? (StepNumber / (MaxNumberOfShearSweeps/100)) : StepNumber;
			cerr << "Progress: " << Progress << "%| StepNumber: " << StepNumber << "| Time passed: " <<  chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() << " s" << endl;
			NextUpdateTime += UPDATE_TIME_INTERVAL;
		}
	}
	FileStreamToWrite.close();

	SC.writeAverageStresses("AvgStresses");

	const int NumberOfRows = AvgVelocities.size()/NumberOfyValues;

	FileName = "AvgVelocities.dat";
	FileStreamToWrite.open(FileName, ios_base::app);
	for (int i = 0; i < NumberOfRows; i++){
		for (int j = 0; j < NumberOfyValues-1; j++){
			FileStreamToWrite << AvgVelocities[NumberOfyValues*i+j]*P.getBoxLength() << '\t';
		}
		FileStreamToWrite << AvgVelocities[NumberOfyValues*i+NumberOfyValues-1]*P.getBoxLength() << '\n';
	}
	FileStreamToWrite.close();

	/*FileStreamToWrite.open("AvgMSD.dat");
	for (int i = 0; i < AvgMSD.size(); i++){
		FileStreamToWrite << AvgMSD[i] << endl;
	}
	FileStreamToWrite.close();*/

	FileStreamToWrite.open("FinalState_long_equilibration_N=1000.dat");
	FileStreamToWrite << P;
	FileStreamToWrite.close();
	cerr << "#VerletListBuilds: " << P.getNumberOfVerletListBuilds() << endl;
	cerr << "Computation time: " << chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() << " s for " << MaxNumberOfShearSweeps+MaxNumberOfEquilibrationSweeps << " MCSweeps." <<  endl << endl;
}

