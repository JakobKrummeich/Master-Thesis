#include "../General_Code/series_analysis.h"
#include <iostream>
#include <stdlib.h>

using namespace std;

int main(int argc, char* argv[]){
	int TotalNumberOfParticles;
	if (argc >= 3){
		TotalNumberOfParticles = atoi(argv[2]);
	}
	else if (argc < 3){
		cerr << "SeriesAnalyzer failed as we need a target FileName and TotalNumberOfParticles!" << endl;
		return 0;
	}
	SeriesAnalyzer Analyzer(TotalNumberOfParticles);

	string OutputFileName = "NAProbDist_";
	OutputFileName += argv[1];
	OutputFileName += ".dat";

	string InputFileNameNASeries;
	string InputFileNamePotEnergySeries;
	while (cin >> InputFileNameNASeries){
		cin >> InputFileNamePotEnergySeries;
		Analyzer.addNewSeries(InputFileNameNASeries, InputFileNamePotEnergySeries);
	}
	Analyzer.normalizeNADistribution();
	Analyzer.writeNAProbabilityDistributionToFile(OutputFileName);
}

