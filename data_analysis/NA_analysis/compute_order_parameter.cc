#include "NA_analysis.h"


int main(int argc, char* argv[]){
	string file = argv[1];
	int numberOfParticles = stoi(argv[2]);

	NAAnalyzer analyzer(numberOfParticles);
	analyzer.readDistribution(file);
	cout << computeOrderParameter(analyzer.NADistribution);
}
