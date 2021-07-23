#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include "../general_code/structure_factor.h"

using namespace std;

double extractBoxLength(string InputFileName){
	ifstream InputFileStream;
	string CurrentString;
	InputFileStream.open(InputFileName);
	getline(InputFileStream, CurrentString, '\n');
	getline(InputFileStream, CurrentString, ':');
	getline(InputFileStream, CurrentString, ':');
	getline(InputFileStream, CurrentString, ':');
	getline(InputFileStream, CurrentString, '\n');
	return stod(CurrentString);
}

int main(int argc, char* argv[]){
	double kMax;
	int TotalNumberOfParticles;
	int NumberOfStates;
	string InputFileName;
	string OutputFileName;
	if (argc == 6){
		InputFileName = argv[1];
		OutputFileName = "StructureFactors_";
		OutputFileName += argv[2];
		OutputFileName += ".dat";
		kMax = atof(argv[3]);
		TotalNumberOfParticles = atoi(argv[4]);
		NumberOfStates = atoi(argv[5]);
	}
	else if (argc < 6){
		cerr << "StructureFactorComputation failed as we need a target File, extracted FileInfo, kMax, TotalNumberOfParticles and NumberOfStates!" << endl;
		return 0;
	}

	StructureFactorComputator SFComputator(kMax, extractBoxLength(InputFileName), TotalNumberOfParticles);

	SFComputator.computeStructureFactors(InputFileName, NumberOfStates);

	SFComputator.writeResultsToFile(OutputFileName);
}
