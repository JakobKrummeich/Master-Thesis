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

	#pragma omp parallel num_threads(NUMBER_OF_THREADS)
	{
		SimulationManager S(0, TOTAL_NUMBER_OF_PARTICLES, OUTPUT_DIRECTORY);

		#pragma omp for
		for (int RunCount = 0; RunCount < NUMBER_OF_RUNS; RunCount++){
			const auto StartTime = chrono::steady_clock::now();
			S.resetCountersAndBuffers();
			S.readInParticleState(InitialStateFile, RunCount*StatesToSkipPerRun, DENSITY);
			S.setTemperature(Temperature);
			#pragma omp critical(WRITE_TO_ERROR_STREAM)
			{
				cerr << "Run " << RunCount+RunNumberOffset << ": equilibration starting." << endl;
			}

			int NumberOfEquilibrationSweeps = S.equilibrateFixedTime(MAX_RUNTIME_IN_MINUTES, StartTime);

			#pragma omp critical(WRITE_TO_ERROR_STREAM)
			{
				cerr << "Run " << RunCount+RunNumberOffset << ": Time for " << NumberOfEquilibrationSweeps << " equilibration sweeps: " << chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() << " s" << endl;
			}
			S.writeParticleStateToFile(OUTPUT_DIRECTORY+"/State_"+"N="+to_string(TOTAL_NUMBER_OF_PARTICLES)+"_T="+to_string(Temperature)+"_AvgDens="+to_string(DENSITY)+"_EquilibrationTime="+to_string(MAX_RUNTIME_IN_MINUTES)+"_epsAB="+to_string(AB_INTERACTION_STRENGTH)+"_"+to_string(RunCount)+"_"+to_string(0)+".dat");
			S.writeSimulationMetaDataToErrorStream(RunCount, chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count(), NumberOfEquilibrationSweeps);
		}
	}
}

