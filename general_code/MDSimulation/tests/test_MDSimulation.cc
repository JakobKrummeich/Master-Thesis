#include <iostream>
#include <stdlib.h>
#include <string>
#include <chrono>
#include "../particles.h"

using namespace std;

int main(int argc, char* argv[]){
	int MaxNumberOfSweeps = 1000000;

	int NumberOfyValues = 15;
	int NumberOfStressSubdivisions = 14;

	Particles P;
	P.initialize(0,DENSITY,1.0,0.0,NumberOfStressSubdivisions);

	const auto StartTime = chrono::steady_clock::now();
	int NextUpdateTime = UPDATE_TIME_INTERVAL;

	double Stepsize = 0.0001;
	string FileName = OUTPUT_DIRECTORY+"_EnergySeries.dat";
	ofstream FileStreamToWrite;
	FileStreamToWrite.open(FileName);
	FileStreamToWrite << "U\t" << "T\t" << "H\n";
	for (int StepNumber = 0; StepNumber < MaxNumberOfSweeps; StepNumber++){
		P.applyTimeStep(Stepsize);
		FileStreamToWrite <<  fixed << setprecision(numeric_limits<long double>::digits10+1) << P.computePotentialEnergy() << '\t' << P.computeKineticEnergy() << '\t' <<  P.computeEnergy() << endl;
		if (chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() >= NextUpdateTime){
			int Progress = MaxNumberOfSweeps >= 100 ? (StepNumber / (MaxNumberOfSweeps/100)) : StepNumber;
			cerr << "Progress: " << Progress << "%| StepNumber: " << StepNumber << "| Time passed: " <<  chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() << " s" << endl;
			NextUpdateTime += UPDATE_TIME_INTERVAL;
		}
	}
	FileStreamToWrite.close();
	cerr << "#VerletListBuilds: " << P.getNumberOfVerletListBuilds() << endl;
	cerr << "Computation time: " << chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() << " s for " << MaxNumberOfSweeps << " MCSweeps." <<  endl << endl;
}

