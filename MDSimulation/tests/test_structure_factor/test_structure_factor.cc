#include <iostream>
#include <stdlib.h>
#include <string>
#include <chrono>
#include "../../structure_factor.h"

using namespace std;

int main(int argc, char* argv[]){

	StructureFactorComputator SC(1000, 36.5148371670110734044, 5.0, 14.0, 100, 20);

	SC.computeNewStructureFactorValues("N=1000_final_state.dat");
	SC.writeResultsToFile("./structure_factors.dat");
}

