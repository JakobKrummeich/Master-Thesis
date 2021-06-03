#include "../General_Code/structure_factor.h"
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>

using namespace std;

double extractBoxLength(string InputFileName){
	ifstream InputFileStream;
	string CurrentString;
	InputFileStream.open(InputFileName);
	getline(InputFileStream, CurrentString, ':');
	getline(InputFileStream, CurrentString, ':');
	getline(InputFileStream, CurrentString, ':');
	getline(InputFileStream, CurrentString, '\n');
	return stod(CurrentString);
}

int main(int argc, char* argv[]){
	double kMax = 25.0;
	if (argc == 2){
		kMax = atof(argv[1]);
	}
	StructureFactorComputator SFComputator(kMax);
	string InputFileName;
	string OutputFileName = "StructureFactors";
	cin >> InputFileName;
	SFComputator.initialize(extractBoxLength(InputFileName));
	SFComputator.readInParticleState(InputFileName);
	SFComputator.computeNewStructureFactorValues();
	
	while (cin >> InputFileName){
		SFComputator.readInParticleState(InputFileName);
		SFComputator.computeNewStructureFactorValues();
	}
	SFComputator.writeResultsToFile(OutputFileName);
}
