#ifndef NA_ANALYSIS_INCLUDED
#define NA_ANALYSIS_INCLUDED

#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <limits>
#include <string>
#include <fstream>

using namespace std;

class NAAnalyzer{

	private:

		struct ValuePair {
			double xValue;
			double yValue;

			ValuePair(double xValue, double yValue):
				xValue(xValue),
				yValue(yValue){
			}
		};

		unsigned long TotalNumberOfParticles;
		vector<ValuePair> NADistribution;

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
			vector<unsigned long> histogram(TotalNumberOfParticles);

			ifstream fS;
			fS.open(filePath);

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

		static double computeIntegral(const vector<ValuePair>& Distribution) {
			double h = Distribution[1].xValue - Distribution[0].xValue;
			double Integral = 0.5 * (Distribution[0].yValue + Distribution.back().yValue);
			for (unsigned long i = 1; i < Distribution.size() - 1; i++){
				Integral += Distribution[i].yValue;
			}
			return h*Integral;
		}

		static double computeFirstMoment(const vector<ValuePair>& Distribution) {
			vector<ValuePair> FirstMomentDistribution(Distribution);
			for (unsigned long i = 0; i < Distribution.size(); i++){
				FirstMomentDistribution[i].yValue = Distribution[i].yValue * Distribution[i].xValue;
			}
			return computeIntegral(FirstMomentDistribution);
		}

		static double computeSecondCentralMoment(const vector<ValuePair>& Distribution, double Mean) {
			vector<ValuePair> SecondCentralMomentDistribution(Distribution);
			for (unsigned long i = 0; i < Distribution.size(); i++){
				double Difference = (Distribution[i].xValue - Mean);
				SecondCentralMomentDistribution[i].yValue = Distribution[i].yValue * Difference * Difference;
			}
			return computeIntegral(SecondCentralMomentDistribution);
		}

		static double computeFourthCentralMoment(const vector<ValuePair>& Distribution, double Mean) {
			vector<ValuePair> FourthCentralMomentDistribution(Distribution);
			for (unsigned long i = 0; i < Distribution.size(); i++){
				double Difference = (Distribution[i].xValue - Mean);
				FourthCentralMomentDistribution[i].yValue = Distribution[i].yValue * Difference * Difference * Difference * Difference;
			}
			return computeIntegral(FourthCentralMomentDistribution);
		}

		static void normalizeDistribution(vector<ValuePair>& Distribution) {
			double TotalIntegral = computeIntegral(Distribution);
			for (unsigned long i = 0; i < Distribution.size(); i++){
				Distribution[i].yValue /= TotalIntegral;
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

		NAAnalyzer(unsigned long TotalNumberOfParticles):
			TotalNumberOfParticles(TotalNumberOfParticles){
			for (unsigned long i = 0; i <= TotalNumberOfParticles; i++){
				NADistribution.emplace_back(static_cast<double>(i)/static_cast<double>(TotalNumberOfParticles), 0.0);
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
