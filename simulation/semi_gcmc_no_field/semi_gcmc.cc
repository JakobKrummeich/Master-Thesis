#include <iostream>
#include <stdlib.h>
#include <string>
#include "../../general_code/simulation_config.h"
#include "../../general_code/simulation_manager.h"

int main(int argc, char* argv[]){
	runSGCMCSimulationForMultipleStartStates(MaxTemperature, MinTemperature, TemperatureStep, NumberOfRuns, NumberOfSweeps, RunNumberOffset, Density, Directory);
}

