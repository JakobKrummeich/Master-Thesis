#include <iostream>
#include <stdlib.h>
#include "../../general_code/simulation_manager.h"

int main(int argc, char* argv[]){
	double InitialDensity = 0.6;
	double MaxPressure = 1.0;
	double PressureStep = 0.1;
	double MinPressure = 0.85;
	int NumberOfSweeps = 100;
	AB_INTERACTION_STRENGTH = 0.1;
	int NumberOfRuns = 2;
	int RunNumberOffset = 0;
	double Temperature = 1.0;
	if (argc != 10){
		cerr << "WARNING: " << argc-1 <<  " arguments were given, but exactly 9 arguments are needed: Initial density, MaxPressure, Pressure Stepsize, MinPressure (not included), NumberOfMCSweeps, AB_INTERACTION_STRENGTH, NumberOfRuns, RunNumberOffset, Temperature. Running with default parameters: Initial density = 0.6, MaxPressure = 1.0, Pressure Stepsize = 0.1, MinPressure = 0.85, NumberOfMCSweeps = 100, AB_INTERACTION_STRENGTH = 0.1, NumberOfRuns = 2, RunNumberOffset = 0, Temperature = 1.0" << endl;
	}
	else {
		InitialDensity = atof(argv[1]);
		MaxPressure = atof(argv[2]);
		PressureStep = atof(argv[3]);
		MinPressure = atof(argv[4]);
		NumberOfSweeps = atoi(argv[5]);
		AB_INTERACTION_STRENGTH = atof(argv[6]);
		NumberOfRuns = atoi(argv[7]);
		RunNumberOffset = atoi(argv[8]);
		Temperature = atof(argv[9]);
	}

	runPressureMCForMultipleStartStates(MaxPressure, MinPressure, PressureStep, NumberOfRuns, NumberOfSweeps, RunNumberOffset, InitialDensity, Temperature);
}

