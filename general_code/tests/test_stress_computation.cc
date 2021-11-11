#include <iostream>
#include <stdlib.h>
#include <string>
#include "../simulation_config.h"
#include "../simulation_manager.h"

using namespace std;

int main(int argc, char* argv[]){
	double Temperature = 1.0;
	int MaxNumberOfSweeps = 1000000;
	double ShearRate = 0.0;
	int NumberOfyValues = 15;
	int NumberOfStressSubdivisions = 8;

	SimulationManager S(0, TOTAL_NUMBER_OF_PARTICLES, OUTPUT_DIRECTORY);
	S.setShearRate(ShearRate);
	S.setMaxFieldStrength(0.0);
	S.setTemperature(Temperature);

	const auto StartTime = chrono::steady_clock::now();
	S.initializeParticles(DENSITY,0.2,NumberOfStressSubdivisions);

	S.setFileNameString(0, MCModus::SGCMC, MaxNumberOfSweeps);
	S.resetCountersAndBuffers();
	S.runSGCMCSimulationForSingleTemperatureWithShear(0, MAX_RUNTIME_IN_MINUTES, StartTime, MaxNumberOfSweeps, NumberOfyValues);
}

