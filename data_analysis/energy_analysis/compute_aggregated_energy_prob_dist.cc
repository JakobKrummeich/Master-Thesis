#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <string>

#include "../energy_analysis/energy_analysis.h"

using namespace std;

int main(int argc, char* argv[]){
	if (argc < 4){
		cerr << "energyAnalyzer failed as we need numberOfBins, Output-directory, numberOfEquilibrations!" << endl;
		return 0;
	}

	int numberOfBins = atoi(argv[1]);
	int numberOfEquilibrations = atoi(argv[3]);

	vector<EnergyAnalyzer> analyzers;
	for (int i = 0; i < 4; i++){
		analyzers.emplace_back(numberOfBins);
	}

	string inputFile;
	while (cin >> inputFile){
		cerr << "new energy series: " << inputFile << endl;
		for (int i = 0; i < 4; i++){
			analyzers[i].addNewSeries(inputFile, numberOfEquilibrations, i);
		}
	}
	vector<string> labels{"U","T","H","T_y"};
	for (int i = 0; i < 4; i++){
		analyzers[i].generateDistribution();
		string fileName = argv[2];
		fileName += "/aggregatedProbDist_"+labels[i]+".dat";
		analyzers[i].writeProbabilityDistributionToFile(fileName, labels[i]);
		double firstMoment = analyzers[i].computeFirstMomentOfDistribution();
		double standardDeviation = sqrt(analyzers[i].computeSecondMoment());
		cout << labels[i] << ": " << firstMoment << '\t' << standardDeviation << endl;
	}
}

