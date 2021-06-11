#include "../general_code/series_analysis.h"
#include <iostream>
#include <stdlib.h>

using namespace std;

int main(int argc, char* argv[]){
	int TotalNumberOfParticles;
	if (argc >= 5){
		TotalNumberOfParticles = atoi(argv[2]);
	}
	else if (argc < 5){
		cerr << "SeriesAnalyzer failed as we need a NADistribution-FileName, TotalNumberOfParticles, Output-directory and MinNumberOfEquilibrationSweeps!" << endl;
		return 0;
	}
	SeriesAnalyzer Analyzer(TotalNumberOfParticles);

	string OutputFileName = argv[3];
	OutputFileName += "NAProbDist_";
	OutputFileName += argv[1];
	OutputFileName += ".dat";
	int MinNumberOfEquilibrationSweeps = stoi(argv[4]);

	string InputFileNameNASeries;
	string InputFileNamePotEnergySeries;
	while (cin >> InputFileNameNASeries){
		cin >> InputFileNamePotEnergySeries;
		Analyzer.addNewSeries(InputFileNameNASeries, InputFileNamePotEnergySeries, MinNumberOfEquilibrationSweeps);
	}
	Analyzer.normalizeNADistribution();
	Analyzer.writeNAProbabilityDistributionToFile(OutputFileName);
	double FirstMoment = Analyzer.computeFirstMomentOfHalfDistribution();
	double BinderCumulant = Analyzer.computeBinderCumulant();
	cout << FirstMoment << '\t' << BinderCumulant;
}

