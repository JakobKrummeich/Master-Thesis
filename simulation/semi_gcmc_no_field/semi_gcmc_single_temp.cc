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

	#pragma omp parallel num_threads(NUMBER_OF_THREADS)
	{
		SimulationManager S(0, TOTAL_NUMBER_OF_PARTICLES, NUMBER_OF_SWEEPS, OUTPUT_DIRECTORY);

		#pragma omp for
		for (int RunCount = RUN_NUMBER_OFFSET; RunCount < NUMBER_OF_RUNS+RUN_NUMBER_OFFSET; RunCount++){
			const auto StartTime = chrono::steady_clock::now();
			S.readInParticleState(InitialStateFile, RunCount*StatesToSkipPerRun, DENSITY);
			S.setTemperature(Temperature);
			#pragma omp critical(WRITE_TO_ERROR_STREAM)
			{
				cerr << "Run " << RunCount << ": Time for initialization and equilibration: " << chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() << " s" << endl;
			}
			S.setFileNameString(RunCount, MCModus::SGCMC);
			S.resetCountersAndBuffers();
			S.runSGCMCSimulationForSingleTemperatureTimeControlled(RunCount, MAX_RUNTIME_IN_MINUTES);
		}
	}
}

