#include <iostream>
#include <stdlib.h>

#include "NA_analysis.h"

using namespace std;

int main(int argc, char* argv[]){
	if (argc < 4){
		cerr << "NAAnalyzer failed as we need a NADistribution-FileName, TotalNumberOfParticles, Output-directory!" << endl;
		return 0;
	}

	string OutputFileName = argv[3];
	OutputFileName += "NAProbDist_";
	OutputFileName += argv[1];
	OutputFileName += ".dat";

	string symmetrizedDistFilename = argv[3];
	symmetrizedDistFilename += "symmetrizedNAProbDist_";
	symmetrizedDistFilename += argv[1];
	symmetrizedDistFilename += ".dat";

	int TotalNumberOfParticles = atoi(argv[2]);

	NAAnalyzer Analyzer(TotalNumberOfParticles);

	string InputFileNameNASeries;
	while (cin >> InputFileNameNASeries){
		cerr << "NewNASeriesFile: " << InputFileNameNASeries << endl;
		Analyzer.addNewHistogram(InputFileNameNASeries);
	}
	Analyzer.normalizeNADistribution();
	Analyzer.writeNAProbabilityDistributionToFile(OutputFileName);
	Analyzer.computeAndWriteSymmetrizedNADistribution(symmetrizedDistFilename);

	double FirstMoment = Analyzer.computeFirstMomentOfHalfDistribution();
	double BinderCumulant = Analyzer.computeBinderCumulant();
	cout << FirstMoment << '\t' << BinderCumulant;
}

