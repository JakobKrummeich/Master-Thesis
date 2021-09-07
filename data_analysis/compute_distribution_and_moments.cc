#include <iostream>
#include <stdlib.h>

#include "../general_code/series_analysis.h"

using namespace std;

int main(int argc, char* argv[]){
	if (argc < 5){
		cerr << "SeriesAnalyzer failed as we need a NADistribution-FileName, TotalNumberOfParticles, Output-directory and MinNumberOfEquilibrationSweeps!" << endl;
		return 0;
	}

	string OutputFileName = argv[3];
	OutputFileName += "NAProbDist_";
	OutputFileName += argv[1];
	OutputFileName += ".dat";
	int TotalNumberOfParticles = atoi(argv[2]);
	int MinNumberOfEquilibrationSweeps = stoi(argv[4]);

	SeriesAnalyzer Analyzer(TotalNumberOfParticles);

	string InputFileNameNASeries;
	while (cin >> InputFileNameNASeries){
		cerr << "NewNASeriesFile: " << InputFileNameNASeries << endl;
		Analyzer.addNewSeries(InputFileNameNASeries, MinNumberOfEquilibrationSweeps);
	}
	Analyzer.normalizeNADistribution();
	Analyzer.writeNAProbabilityDistributionToFile(OutputFileName);
	double FirstMoment = Analyzer.computeFirstMomentOfHalfDistribution();
	double BinderCumulant = Analyzer.computeBinderCumulant();
	cout << FirstMoment << '\t' << BinderCumulant;
}

