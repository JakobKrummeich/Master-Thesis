#include <iostream>
#include <stdlib.h>
#include <string>
#include "../../general_code/simulation_config.h"
#include "../../general_code/simulation_manager.h"

using namespace std;

int main(int argc, char* argv[]){

	int NumberOfStartingStates = atoi(argv[1]);
	int NumberOfRandomizationSteps = atoi(argv[2]);

	#pragma omp parallel num_threads(NUMBER_OF_THREADS)
	{
		SimulationManager S(0, TOTAL_NUMBER_OF_PARTICLES, OUTPUT_DIRECTORY);

		#pragma omp for
		for (int RunCount = 0; RunCount < NumberOfStartingStates; RunCount++){
			const auto StartTime = chrono::steady_clock::now();
			S.initializeParticles(DENSITY);
			S.setTemperature(1.0);
			S.randomizeInitialPosition(RunCount, NumberOfRandomizationSteps);

			#pragma omp critical(WRITE_TO_ERROR_STREAM)
			{
				cerr << "Run " << RunCount << ": Time for " << NumberOfRandomizationSteps << " randomization sweeps: " << chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() << " s" << endl;
			}
			S.writeParticleStateToFile("test"+to_string(RunCount)+".dat");
		}
	}
}

