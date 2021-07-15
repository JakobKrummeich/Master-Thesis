#include <iostream>
#include <stdlib.h>
#include <string>
#include "../../general_code/simulation_config.h"
#include "../../general_code/simulation_manager.h"

int main(int argc, char* argv[]){
	double Temperature = atof(argv[1]);
	double xA = atof(argv[2]);

	const int NumberOfAParticles = round(xA*static_cast<double>(TOTAL_NUMBER_OF_PARTICLES));
	const int NumberOfSavedStatesPerRun = NUMBER_OF_SAVED_STATES_PER_TEMPERATURE >= NUMBER_OF_RUNS ?  NUMBER_OF_SAVED_STATES_PER_TEMPERATURE / NUMBER_OF_RUNS : 1;

	#pragma omp parallel num_threads(NUMBER_OF_THREADS)
	{
		SimulationManager S(NumberOfAParticles, NumberOfAParticles, NUMBER_OF_SWEEPS, OUTPUT_DIRECTORY);

		#pragma omp for
		for (int RunCount = RUN_NUMBER_OFFSET; RunCount < NUMBER_OF_RUNS+RUN_NUMBER_OFFSET; RunCount++){
			const auto StartTime = chrono::steady_clock::now();
			S.initializeParticles(DENSITY);
			S.randomizeInitialPosition(RunCount);
			S.setTemperature(Temperature);
			#pragma omp critical(WRITE_TO_ERROR_STREAM)
			{
				cerr << "Run " << RunCount << ": Time for initialization: " << chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() << " s" << endl;
			}
			S.setFileNameString(RunCount, MCModus::CMC);
			S.resetCountersAndBuffers();
			S.runCMCSimulationForSingleTemperature(RunCount, MAX_RUNTIME_IN_MINUTES, NumberOfSavedStatesPerRun);
		}
	}
}

