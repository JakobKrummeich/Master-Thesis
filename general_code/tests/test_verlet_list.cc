#include <iostream>
#include <stdlib.h>
#include <string>
#include "../simulation_config.h"
#include "../simulation_manager.h"

using namespace std;

int main(int argc, char* argv[]){

	Particles P;
	P.initialize(0,0.1,stod(argv[1]));
	cerr << P;
	cerr << endl;
	double Delta [DIMENSION]{0.0,-0.05};
	for (int i = 0; i < 10; i++){
		cerr << "Change of Pot energy: " << P.computeChangeInPotentialEnergyByMoving(0,Delta) << endl << endl;
		P.updatePosition(0,Delta);
		cerr << P << endl;
	}
}

