#include <iostream>
#include <stdlib.h>
#include "../../general_code/simulation_manager.h"

int main(int argc, char* argv[]){
	double InitialDensity = 0.6;
	AB_INTERACTION_STRENGTH = 0.1;
	double MaxTemperature = 1.0;
	double TemperatureStep = 0.1;
	double MinTemperature = 0.85;
	int NumberOfSweeps = 100;
	int NumberOfRuns = 2;
	int RunNumberOffset = 0;
	if (argc != 9){
		cerr << "WARNING: " << argc-1 <<  " arguments were given, but exactly 8 arguments are needed: Average density, MaxTemperature, Temperature Stepsize, MinTemperature (not included), NumberOfMCSweeps, AB_INTERACTION_STRENGTH, NumberOfRuns, RunNumberOffset. Running with default parameters: Average density = 0.6, MaxTemperature = 1.0, Temperature Stepsize = 0.1, MinTemperature = 0.85, NumberOfMCSweeps = 100, AB_INTERACTION_STRENGTH = 0.1, NumberOfRuns = 2, RunNumberOffset = 0" << endl;
	}
	else {
		InitialDensity = atof(argv[1]);
		AB_INTERACTION_STRENGTH = atof(argv[6]);
		MaxTemperature = atof(argv[2]);
		TemperatureStep = atof(argv[3]);
		MinTemperature = atof(argv[4]);
		NumberOfSweeps = atoi(argv[5]);
		NumberOfRuns = atoi(argv[7]);
		RunNumberOffset = atoi(argv[8]);
	}

	runSGCMCSimulationForMultipleStartStates(MaxTemperature, MinTemperature, TemperatureStep, NumberOfRuns, NumberOfSweeps, RunNumberOffset, InitialDensity);
}

