#ifndef NA_ANALYSIS_INCLUDED
#define NA_ANALYSIS_INCLUDED

#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <string>
#include <fstream>

#include "../value_pair.h"
#include "../integration.h"

using namespace std;

class NAAnalyzer{

	private:

		unsigned long TotalNumberOfParticles;

		vector<unsigned long> readInSeries(string FilePath){
			vector<unsigned long> Series;

			ifstream FileStreamToReadIn;
			FileStreamToReadIn.open(FilePath);

			string CurrentString;
			getline(FileStreamToReadIn, CurrentString);

			while (getline(FileStreamToReadIn, CurrentString)){
				Series.push_back(stoul(CurrentString));
			}
			return Series;
		}

		vector<unsigned long> readHistogram(string filePath){
			vector<unsigned long> histogram(TotalNumberOfParticles,0);

			ifstream fS(filePath);
			if (!fS.is_open()){
				return histogram;
			}

			string currentString;
			getline(fS, currentString);

			for (int i = 0; i <= TotalNumberOfParticles; i++){
				fS >> currentString;
				fS >> currentString;
				histogram[i] = stoul(currentString);
			}
			return histogram;
		}

		void updateHistogramWithNewSeries(const vector<unsigned long>& NASeries){
			for (unsigned long i = 0; i < NASeries.size(); i++){
				NADistribution[NASeries[i]].yValue++;
			}
		}

		void updateHistogramWithNewHistogram(const vector<unsigned long>& NAHistogram){
			for (int i = 0; i <= TotalNumberOfParticles; i++){
				NADistribution[i].yValue += static_cast<double>(NAHistogram[i]);
			}
		}

		static void writeDistributionToFile(string FileName, const vector<ValuePair>& Distribution) {
			ofstream FileStreamToWriteTo;
			FileStreamToWriteTo.open(FileName);
			FileStreamToWriteTo << fixed << setprecision(numeric_limits<long double>::digits10+1) << "xA\tprobability_density\n";
			for (unsigned long i = 0; i < Distribution.size(); i++){
				FileStreamToWriteTo << Distribution[i].xValue << '\t' << Distribution[i].yValue << '\n';
			}
			FileStreamToWriteTo.close();
		}

		static vector<ValuePair> symmetrizeDistribution(const vector<ValuePair>& Distribution) {
			vector<ValuePair> SymmetrizedDistribution(Distribution);
			for (unsigned long i = 0; i < SymmetrizedDistribution.size()/2; i++){
				SymmetrizedDistribution[i].yValue += SymmetrizedDistribution[SymmetrizedDistribution.size()-1-i].yValue;
				SymmetrizedDistribution[SymmetrizedDistribution.size()-1-i].yValue = SymmetrizedDistribution[i].yValue;
			}
			if (SymmetrizedDistribution.size() % 2 != 0){
				SymmetrizedDistribution[SymmetrizedDistribution.size()/2].yValue *= 2.0;
			}
			normalizeDistribution(SymmetrizedDistribution);
			return SymmetrizedDistribution;
		}

	public:

		vector<ValuePair> NADistribution;

		NAAnalyzer(unsigned long TotalNumberOfParticles):
			TotalNumberOfParticles(TotalNumberOfParticles){
			for (unsigned long i = 0; i <= TotalNumberOfParticles; i++){
				NADistribution.emplace_back(static_cast<double>(i)/static_cast<double>(TotalNumberOfParticles), 0.0);
			}
		}

		void readDistribution(string path){
			NADistribution.clear();

			ifstream fs(path);

			string s;
			getline(fs, s);

			ValuePair vp;

			while(fs >> vp.xValue){
				fs >> vp.yValue;
				NADistribution.push_back(vp);
			}
		}

		void addNewSeries(string FileNameNASeries){
			vector<unsigned long> NewNASeries = readInSeries(FileNameNASeries);
			updateHistogramWithNewSeries(NewNASeries);
		}

		void addNewHistogram(string fileName){
			vector<unsigned long> newNAHistogram = readHistogram(fileName);
			updateHistogramWithNewHistogram(newNAHistogram);
		}

		void normalizeNADistribution(){
			normalizeDistribution(NADistribution);
		}

		void writeNAProbabilityDistributionToFile(string FileName) const {
			writeDistributionToFile(FileName, NADistribution);
		}

		void computeAndWriteSymmetrizedNADistribution(string filename) const {
			vector<ValuePair> symmetrizedDistribution = symmetrizeDistribution(NADistribution);
			writeDistributionToFile(filename, symmetrizedDistribution);
		}

		double computeFirstMomentOfHalfDistribution() const {

			vector<ValuePair> SymmetrizedDistribution = symmetrizeDistribution(NADistribution);
			vector<ValuePair> HalfDistribution;

			for (int i = 0; i < SymmetrizedDistribution.size()/2; i++){
				HalfDistribution.push_back(SymmetrizedDistribution[i]);
			}
			if (SymmetrizedDistribution.size() % 2 != 0){
				HalfDistribution.push_back(SymmetrizedDistribution[SymmetrizedDistribution.size()/2]);
			}
			normalizeDistribution(HalfDistribution);
			return computeFirstMoment(HalfDistribution);
		}

		double computeBinderCumulant() const {

			vector<ValuePair> SymmetrizedDistribution = symmetrizeDistribution(NADistribution);

			double SecondCentralMoment = computeSecondCentralMoment(SymmetrizedDistribution, 0.5);
			double FourthCentralMoment = computeFourthCentralMoment(SymmetrizedDistribution, 0.5);
			return (1.0 - FourthCentralMoment/(SecondCentralMoment*SecondCentralMoment*3.0));
		}

		double computeStructureFactor(bool RestrictToHalfDistribution) const {
			vector<ValuePair> SymmetrizedDistribution = symmetrizeDistribution(NADistribution);
			if (!RestrictToHalfDistribution){
				return static_cast<double>(TotalNumberOfParticles)*(computeSecondCentralMoment(SymmetrizedDistribution, 0.5));
			}

			vector<ValuePair> HalfDistribution;
			for (int i = 0; i < SymmetrizedDistribution.size()/2; i++){
				HalfDistribution.push_back(SymmetrizedDistribution[i]);
			}
			if (SymmetrizedDistribution.size() % 2 != 0){
				HalfDistribution.push_back(SymmetrizedDistribution[SymmetrizedDistribution.size()/2]);
			}
			normalizeDistribution(HalfDistribution);
			double Mean = computeFirstMoment(HalfDistribution);
			return static_cast<double>(TotalNumberOfParticles)*(computeSecondCentralMoment(HalfDistribution, Mean));
		}
};

#endif /* NA_ANALYSIS_INCLUDED */
