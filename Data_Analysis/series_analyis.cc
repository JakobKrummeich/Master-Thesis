#include <iostream>
#include <cmath>
#include <chrono>
#include <random>
#include <fstream>
#include <vector>
#include <iomanip>
#include <stdlib.h>
#include <regex>

using namespace std;

int TotalNumberOfParticles;

vector<double> readInSeries(string FilePath){
	
	vector<double> Series;
	
	ifstream FileStreamToReadIn;
	FileStreamToReadIn.open(FilePath);

	string CurrentString;
	getline(FileStreamToReadIn, CurrentString);
	
	while (getline(FileStreamToReadIn, CurrentString)){
		Series.push_back(stod(CurrentString));
	}
	return Series;
}

vector<int> computeNAHist(const vector<double>& NASeries){
	vector<int> Hist;
	return Hist;
}

int main(int argc, char* argv[]){
	string FilePath = "NA_Series_N=1000_T=0.500000_AvgDens=0.560000_MCRuns=600000.dat";
	vector<double> NASeries = readInSeries(FilePath);
	smatch sm;
	regex_search (FilePath, sm, regex (R"(N=\d*)"));
	string matchedSubstring = sm[0];
	regex_search (matchedSubstring, sm, regex (R"(\d+)"));
	TotalNumberOfParticles = stoi(sm[0]);
}
