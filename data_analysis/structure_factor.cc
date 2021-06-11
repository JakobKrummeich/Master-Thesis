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
	double kMax;
	if (argc >= 3){
		kMax = atof(argv[2]);
	}
	else if (argc < 3){
		cerr << "StructureFactorComputation failed as we need a target FileName and kMax!" << endl;
		return 0;
	}
	StructureFactorComputator SFComputator(kMax);
	string InputFileName;
	cin >> InputFileName;
	string OutputFileName = "StructureFactors_";
	OutputFileName += argv[1];
	OutputFileName += ".dat";

	SFComputator.initialize(extractBoxLength(InputFileName));
	SFComputator.readInParticleState(InputFileName);
	SFComputator.computeNewStructureFactorValues();

	while (cin >> InputFileName){
		SFComputator.readInParticleState(InputFileName);
		SFComputator.computeNewStructureFactorValues();
	}
	SFComputator.writeResultsToFile(OutputFileName);
}
