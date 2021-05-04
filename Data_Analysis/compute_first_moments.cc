#include <iostream>
#include <cmath>
#include <chrono>
#include <random>
#include <fstream>
#include <vector>
#include <iomanip>
#include <stdlib.h>

using namespace std;

vector<int> readInSeries(string FilePath){
	
	vector<int> Series;
	
	ifstream FileStreamToReadIn;
	FileStreamToReadIn.open(FilePath);

	string CurrentString;
	getline(FileStreamToReadIn, CurrentString);
	
	while (getline(FileStreamToReadIn, CurrentString)){
		Series.push_back(stoi(CurrentString));
	}
	return Series;
}

struct FunctionValues {
	vector<double> xValues;
	vector<double> yValues;
};


FunctionValues computeNAHist(const vector<int>& NASeries, int EquilibrationIndex, int TotalNumberOfParticles){
	vector<double> Hist;
	Hist.resize(TotalNumberOfParticles + 1);
	for (int i = EquilibrationIndex; i < NASeries.size(); i++){
		Hist[NASeries[i]]++;
	}
	FunctionValues Histogram ;
	Histogram.yValues = Hist;
	for (int i = 0; i < TotalNumberOfParticles + 1; i++){
		Histogram.xValues.push_back(static_cast<double>(i)/static_cast<double>(TotalNumberOfParticles));
	}
	return Histogram;
}

int computeEquilibrationIndex(const vector<int>& NASeries){
	for (int CurrentIndex = 0; CurrentIndex < NASeries.size(); CurrentIndex++){
		double AverageValue = 0.0;
		int NumberOfPoints = 0;
		for (int j = CurrentIndex; j < NASeries.size(); j++){
			AverageValue += static_cast<double>(NASeries[j]);
			NumberOfPoints++;
		}
		AverageValue /= static_cast<double>(NumberOfPoints);
		double AverageDeviation = 0.0;
		for (int j = CurrentIndex; j < NASeries.size(); j++){
			AverageDeviation += abs(AverageValue - static_cast<double>(NASeries[j]));
		}
		AverageDeviation /= static_cast<double>(NumberOfPoints);
		if (abs(AverageValue - NASeries[CurrentIndex]) < AverageDeviation){
			return CurrentIndex;
		}
	}
	return 0;
}


double computeMomentOfDistribution(const FunctionValues& Values, double Exponent){
	double h = Values.xValues[1] - Values.xValues[0];
	double Moment = 0.5 * (pow(Values.xValues[0],Exponent) * Values.yValues[0] + pow(Values.xValues.back(), Exponent) * Values.yValues.back());
	for (int i = 1; i < Values.xValues.size(); i++){
		Moment += pow(Values.xValues[i], Exponent)*Values.yValues[i];
	}
	return h*Moment;
}

FunctionValues computeNormalizedDistribution(const FunctionValues& Function){
	double TotalIntegral = computeMomentOfDistribution(Function, 0.0);
	FunctionValues NormalizedFunction;
	NormalizedFunction.xValues = Function.xValues;
	for (int i = 0; i < NormalizedFunction.xValues.size(); i++){
		NormalizedFunction.yValues.push_back(Function.yValues[i]/TotalIntegral);
	}
	return NormalizedFunction;
}

double computeFirstMomentInSubInterval(const FunctionValues& Values, int LeftIntervalBoundary, int RightIntervalBoundary){
	FunctionValues SubIntervalValues;
	for (int i = LeftIntervalBoundary; i <= RightIntervalBoundary; i++){
		SubIntervalValues.xValues.push_back(Values.xValues[i]);
		SubIntervalValues.yValues.push_back(Values.yValues[i]);
	}
	FunctionValues NormalizedDistribution = computeNormalizedDistribution(SubIntervalValues);
	return computeMomentOfDistribution(NormalizedDistribution, 1.0);
}


int main(int argc, char* argv[]){
	string FilePath = argv[1];
	vector<int> NASeries = readInSeries(FilePath);
	int TotalNumberOfParticles = stoi(argv[2]);
	int EquilibrationIndex = computeEquilibrationIndex(NASeries);
	cerr << EquilibrationIndex << endl;
	FunctionValues NAHist = computeNAHist(NASeries, EquilibrationIndex, TotalNumberOfParticles);
	cerr << computeFirstMomentInSubInterval(NAHist, 0, static_cast<int>(TotalNumberOfParticles/2)) << endl;
}
