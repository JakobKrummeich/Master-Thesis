#include <iostream>
#include <stdlib.h>

#include "../general_code/series_analysis.h"

using namespace std;

int main(int argc, char* argv[]){
	if (argc < 5){
		cerr << "SeriesAnalyzer failed as we need a NADistribution-FileName, TotalNumberOfParticles, Output-directory and MinNumberOfEquilibrationSweeps!" << endl;
		return 0;
	}
	SeriesAnalyzer Analyzer(TotalNumberOfParticles);

	OutputFileName += "NAProbDist_";
	OutputFileName += argv[1];
	OutputFileName += ".dat";
	int TotalNumberOfParticles = atoi(argv[2]);
	string OutputFileName = argv[3];
	int MinNumberOfEquilibrationSweeps = stoi(argv[4]);

	string InputFileNameNASeries;
	string InputFileNamePotEnergySeries;
	while (cin >> InputFileNameNASeries){
		cerr << "NewNASeriesFile: " << InputFileNameNASeries << endl;
		cin >> InputFileNamePotEnergySeries;
		cerr << "NewPotSeriesFile: " << InputFileNamePotEnergySeries << endl;
		Analyzer.addNewSeries(InputFileNameNASeries, InputFileNamePotEnergySeries, MinNumberOfEquilibrationSweeps);
	}
	Analyzer.normalizeNADistribution();
	double Variance = Analyzer.computeStructureFactor();
	cout << Variance;
}

