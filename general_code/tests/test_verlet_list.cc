#include <iostream>
#include <stdlib.h>
#include <string>
#include "../simulation_config.h"
#include "../simulation_manager.h"

using namespace std;

int main(int argc, char* argv[]){

	Particles P;
	P.initialize(0,0.1);
	cerr << P;
	P.buildVerletList(stod(argv[1]));
}

