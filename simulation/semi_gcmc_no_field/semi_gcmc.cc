#include <iostream>
#include <stdlib.h>
#include <string>
#include "../../general_code/simulation_manager.h"

using namespace std;

int main(int argc, char* argv[]){
	double Density = 0.6;
	AB_INTERACTION_STRENGTH = 0.1;
	double MaxTemperature = 1.0;
	double TemperatureStep = 0.1;
	double MinTemperature = 0.85;
	int NumberOfSweeps = 100;
	int NumberOfRuns = 2;
	int RunNumberOffset = 0;
	string Directory = ".";
	if (argc != 10){
		cerr << "WARNING: " << argc-1 <<  " arguments were given, but exactly 9 arguments are needed: Average density, MaxTemperature, Temperature Stepsize, MinTemperature (not included), NumberOfMCSweeps, AB_INTERACTION_STRENGTH, NumberOfRuns, RunNumberOffset, Directory for fresh data. Running with default parameters: Average density = 0.6, MaxTemperature = 1.0, Temperature Stepsize = 0.1, MinTemperature = 0.85, NumberOfMCSweeps = 100, AB_INTERACTION_STRENGTH = 0.1, NumberOfRuns = 2, RunNumberOffset = 0, Directory = ." << endl;
	}
	else {
		Density = atof(argv[1]);
		MaxTemperature = atof(argv[2]);
		TemperatureStep = atof(argv[3]);
		MinTemperature = atof(argv[4]);
		NumberOfSweeps = atoi(argv[5]);
		AB_INTERACTION_STRENGTH = atof(argv[6]);
		NumberOfRuns = atoi(argv[7]);
		RunNumberOffset = atoi(argv[8]);
		Directory = argv[9];
	}

	runSGCMCSimulationForMultipleStartStates(MaxTemperature, MinTemperature, TemperatureStep, NumberOfRuns, NumberOfSweeps, RunNumberOffset, Density, Directory);
}

