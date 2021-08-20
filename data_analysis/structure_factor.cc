#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include "../general_code/structure_factor_alternative.h"

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
	string InputFileName;
	string OutputFileName;
	double kMax;
	double kWidth;
	int TotalNumberOfParticles;
	int NumberOfStates;

	if (argc == 7){
		InputFileName = argv[1];
		OutputFileName = "StructureFactors_";
		OutputFileName += argv[2];
		OutputFileName += ".dat";
		kMax = atof(argv[3]);
		kWidth = atof(argv[4]);
		TotalNumberOfParticles = atoi(argv[5]);
		NumberOfStates = atoi(argv[6]);
	}
	else if (argc < 7){
		cerr << "Structure factor computation failed as we need a target File, extracted FileInfo, kMax, kWidth, TotalNumberOfParticles and NumberOfStates!" << endl;
		return 0;
	}

	StructureFactorComputator SFComputator(kMax, kWidth, extractBoxLength(InputFileName), TotalNumberOfParticles);

	SFComputator.computeStructureFactors(InputFileName, NumberOfStates);

	SFComputator.writeResultsToFile(OutputFileName);
}
