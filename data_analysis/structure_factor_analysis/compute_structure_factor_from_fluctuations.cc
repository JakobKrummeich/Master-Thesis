#include <iostream>
#include <stdlib.h>

#include "../../general_code/series_analysis.h"

using namespace std;

int main(int argc, char* argv[]){
	if (argc < 6){
		cerr << "SeriesAnalyzer failed as we need a NADistribution-FileName, TotalNumberOfParticles, Output-directory, MinNumberOfEquilibrationSweeps and T<Tc!" << endl;
		return 0;
	}

	string OutputFileName = argv[3];
	OutputFileName += "NAProbDist_";
	OutputFileName += argv[1];
	OutputFileName += ".dat";
	int TotalNumberOfParticles = atoi(argv[2]);
	int MinNumberOfEquilibrationSweeps = stoi(argv[4]);
	bool RestrictToHalfDistribution = stoi(argv[5]); // below Tc we need to compute the concentration fluctuation only for [0,0.5] or [0.5,1]

	SeriesAnalyzer Analyzer(TotalNumberOfParticles);

	string InputFileNameNASeries;
	while (cin >> InputFileNameNASeries){
		cerr << "NewNASeriesFile: " << InputFileNameNASeries << endl;
		Analyzer.addNewSeries(InputFileNameNASeries, MinNumberOfEquilibrationSweeps);
	}
	Analyzer.normalizeNADistribution();
	Analyzer.writeNAProbabilityDistributionToFile(OutputFileName);
	double Variance = Analyzer.computeStructureFactor(RestrictToHalfDistribution);
	cout << Variance;
}

