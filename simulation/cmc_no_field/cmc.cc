#include <iostream>
#include <stdlib.h>
#include <string>
#include "../../general_code/simulation_config.h"
#include "../../general_code/simulation_manager.h"

using namespace std;

int main(int argc, char* argv[]){
	double Temperature = atof(argv[1]);
	double xA = atof(argv[2]);
	double RunNumberOffset = atoi(argv[3]);
	int NumberOfInitialRandomizationSweeps = atoi(argv[4]);
	int NumberOfEquilibrationSweeps = atoi(argv[5]);
	int MaxNumberOfSweeps = atoi(argv[6]);

	const int NumberOfAParticles = round(xA*static_cast<double>(TOTAL_NUMBER_OF_PARTICLES));
	const int NumberOfSavedStatesPerRun = NUMBER_OF_SAVED_STATES_PER_TEMPERATURE >= NUMBER_OF_RUNS ? NUMBER_OF_SAVED_STATES_PER_TEMPERATURE / NUMBER_OF_RUNS : 1;

	#pragma omp parallel num_threads(NUMBER_OF_THREADS)
	{
		SimulationManager S(NumberOfAParticles, NumberOfAParticles, OUTPUT_DIRECTORY);

		#pragma omp for
		for (int RunCount = 0; RunCount < NUMBER_OF_RUNS; RunCount++){
			auto StartTime = chrono::steady_clock::now();
			S.initializeParticles(DENSITY);
			S.randomizeInitialPosition(RunCount+RunNumberOffset, NumberOfInitialRandomizationSweeps);
			S.setTemperature(Temperature);
			S.equilibrate(NumberOfEquilibrationSweeps);
			#pragma omp critical(WRITE_TO_ERROR_STREAM)
			{
				cerr << "Run " << RunCount+RunNumberOffset << ": Time for " << NumberOfInitialRandomizationSweeps << " randomization sweeps and " << NumberOfEquilibrationSweeps << " equilibration sweeps: " << chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() << " s" << endl;
			}
			S.setFileNameString(RunCount+RunNumberOffset, MCModus::CMC, MaxNumberOfSweeps);
			S.resetCountersAndBuffers();
			S.runCMCSimulationForSingleTemperature(RunCount+RunNumberOffset, MAX_RUNTIME_IN_MINUTES, NumberOfSavedStatesPerRun, MaxNumberOfSweeps);
		}
	}
}

