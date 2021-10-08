#include <iostream>
#include <stdlib.h>
#include <string>
#include "../../general_code/simulation_config.h"
#include "../../general_code/simulation_manager.h"

using namespace std;

int main(int argc, char* argv[]){
	double Temperature = atof(argv[1]);
	string InitialStateFile = argv[2];
	int StatesToSkipPerRun = atoi(argv[3]);
	int RunNumberOffset = atoi(argv[4]);
	int RunNumberOffsetReadIn = atoi(argv[5]);
	int MaxNumberOfSweeps = atoi(argv[6]);
	double ShearRate = atof(argv[7]);
	int NumberOfyValues = atoi(argv[8]);

	#pragma omp parallel num_threads(NUMBER_OF_THREADS)
	{
		SimulationManager S(0, TOTAL_NUMBER_OF_PARTICLES, OUTPUT_DIRECTORY);
		S.setShearRate(ShearRate);
		S.setTemperature(Temperature);

		#pragma omp for
		for (int RunCount = 0; RunCount < NUMBER_OF_RUNS; RunCount++){
			const auto StartTime = chrono::steady_clock::now();
			S.readInParticleState(InitialStateFile, RunCount*StatesToSkipPerRun+RunNumberOffsetReadIn, DENSITY);

			S.setFileNameString(RunCount+RunNumberOffset, MCModus::SGCMC, MaxNumberOfSweeps);
			S.resetCountersAndBuffers();
			S.runSGCMCSimulationForSingleTemperatureWithShear(RunCount+RunNumberOffset, MAX_RUNTIME_IN_MINUTES, StartTime, MaxNumberOfSweeps, NumberOfyValues);
		}
	}
}

