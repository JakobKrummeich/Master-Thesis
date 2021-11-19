#include <iostream>
#include <stdlib.h>
#include <string>
#include "../simulation_config.h"
#include "../simulation_manager.h"

using namespace std;

int main(int argc, char* argv[]){
	double Temperature = 1.0;
	int MaxNumberOfSweeps = 100000;
	double LeesEdwardsVelocity = 0.0000;
	int NumberOfyValues = 15;
	int NumberOfStressSubdivisions = 14;

	SimulationManager S(0, TOTAL_NUMBER_OF_PARTICLES, OUTPUT_DIRECTORY);
	S.setLeesEdwardsVelocity(LeesEdwardsVelocity);
	S.setMaxFieldStrength(1.0);
	S.setTemperature(Temperature);

	const auto StartTime = chrono::steady_clock::now();
	S.initializeParticles(DENSITY,0.0,NumberOfStressSubdivisions);

	S.setFileNameString(0, MCModus::SGCMC, MaxNumberOfSweeps);
	S.resetCountersAndBuffers();
	S.runSGCMCSimulationForSingleTemperatureWithShear(0, MAX_RUNTIME_IN_MINUTES, StartTime, MaxNumberOfSweeps, NumberOfyValues);
}

