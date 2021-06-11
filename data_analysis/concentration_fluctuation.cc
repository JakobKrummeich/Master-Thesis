#include <random>
#include <fstream>
#include <vector>
#include <iostream>

using namespace std;

vector<int> readInParticleNumberSeries(string FileNameToReadIn) {
	ifstream FileStreamToReadIn;
	FileStreamToReadIn.open(FileNameToReadIn);

	string CurrentString;

	getline(FileStreamToReadIn, CurrentString, '\n');

	vector<int> NumberOfASeries;

	int CurrentNumber;

	while (FileStreamToReadIn >> CurrentNumber){
		NumberOfASeries.push_back(CurrentNumber);
	}
	FileStreamToReadIn.close();
	return NumberOfASeries;
}

int main(int argc, char* argv[]){
	int TotalNumberOfParticles = atoi(argv[1]);
	vector<int> NumberOfASeries = readInParticleNumberSeries("data/NA_Series_N=1000_T=1.000000_AvgDens=0.600000_MCRuns=500000_epsAB=0.100000.dat");
	double xASquaredAverage = 0.0;
	int NumberOfValues = 0;
	for (int i = 0; i < NumberOfASeries.size(); i++){
		double NewxA = static_cast<double>(NumberOfASeries[i])/static_cast<double>(TotalNumberOfParticles);
		xASquaredAverage += NewxA*NewxA;
		xASquaredAverage += (1.0-NewxA)*(1.0-NewxA);
		NumberOfValues += 2;
	}
	double xAAverage = 0.5;
	xASquaredAverage /= static_cast<double>(NumberOfValues);
	cerr << xAAverage << endl;
	cerr << xASquaredAverage << endl;
	double Result = static_cast<double>(TotalNumberOfParticles) * (xASquaredAverage - xAAverage*xAAverage);
	cout << "cc-structure-factor by concentration fluctuation: " << Result << endl;
}
