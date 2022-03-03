#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include "../../MDSimulation/structure_factor.h"

using namespace std;

double extractBoxLength(string fileName){
	ifstream ifs(fileName);
	string CurrentString;
	getline(ifs, CurrentString, '|');
	getline(ifs, CurrentString, '|');
	getline(ifs, CurrentString, '|');
	getline(ifs, CurrentString, ':');
	getline(ifs, CurrentString, '|');
	return stod(CurrentString);
}

int main(int argc, char* argv[]){
	if (argc < 7){
		cerr << "StructureFactorComputation failed as we need an outputFileName, outputDirectory, kAll, kMax, numberOfkIntervals, numberOfAngleIntervals, totalNumberOfParticles!" << endl;
		return 0;
	}

	string outputFileName = argv[2];
	outputFileName += "StructureFactors_";
	outputFileName += argv[1];
	outputFileName += ".dat";

	double kAll = atof(argv[3]);
	double kMax = atof(argv[4]);
	int numberOfkIntervals = atoi(argv[5]);
	int numberOfAngleIntervals = atoi(argv[6]);
	int totalNumberOfParticles = atoi(argv[7]);

	string inputFileName;

	double boxLength;

	if (cin >> inputFileName){
		cerr << "new state: " << inputFileName << endl;
		boxLength = extractBoxLength(inputFileName);
	}

	StructureFactorComputator SC(totalNumberOfParticles, boxLength, kAll, kMax, numberOfkIntervals, numberOfAngleIntervals);
	SC.computeNewStructureFactorValues(inputFileName);

	while (cin >> inputFileName){
		cerr << "new state: " << inputFileName << endl;
		SC.computeNewStructureFactorValues(inputFileName);
	}

	SC.writeResultsToFile(outputFileName);
}
