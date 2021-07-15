#include <iostream>
#include <stdlib.h>
#include <string>
#include "../../general_code/simulation_config.h"
#include "../../general_code/simulation_manager.h"

int main(int argc, char* argv[]){
	runSGCMCSimulationForMultipleStartStates(MAX_TEMPERATURE, MIN_TEMPERATURE, TEMPERATURE_STEP, NUMBER_OF_RUNS, NUMBER_OF_SWEEPS, RUN_NUMBER_OFFSET, DENSITY, OUTPUT_DIRECTORY);
}

