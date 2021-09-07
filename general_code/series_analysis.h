#ifndef SERIES_ANALYSIS_INCLUDED
#define SERIES_ANALYSIS_INCLUDED

#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <fstream>

using namespace std;

class SeriesAnalyzer{

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

		void updateHistogramWithNewSeries(const vector<unsigned long>& NASeries, unsigned long MinNumberOfSweeps ){
			for (unsigned long i = MinNumberOfSweeps; i < NASeries.size(); i++){
				NADistribution[NASeries[i]].yValue++;
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

		static double computeFirstMomentOfDistribution(const vector<ValuePair>& Distribution) {
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
			FileStreamToWriteTo << "xA\tprobability_density\n";
			for (unsigned long i = 0; i < Distribution.size(); i++){
				FileStreamToWriteTo << Distribution[i].xValue << '\t' << Distribution[i].yValue << '\n';
			}
			FileStreamToWriteTo.close();
		}

		static vector<ValuePair> symmetrizeDistribution(const vector<ValuePair>& Distribution) {
			vector<ValuePair> LeftDistribution;
			for (unsigned long i = 0; i <= (Distribution.size()-1)/2; i++){
				LeftDistribution.push_back(Distribution[i]);
			}
			double LeftHalfIntegral = computeIntegral(LeftDistribution);

			vector<ValuePair> RightDistribution;
			for (unsigned long i = Distribution.size()/2; i < Distribution.size(); i++){
				RightDistribution.push_back(Distribution[i]);
			}
			double RightHalfIntegral = computeIntegral(RightDistribution);

			vector<ValuePair> SymmetrizedDistribution(Distribution);
			if (LeftHalfIntegral > RightHalfIntegral){
				for (unsigned long i = 0; i < SymmetrizedDistribution.size()/2; i++){
					SymmetrizedDistribution[SymmetrizedDistribution.size()-1-i].yValue = SymmetrizedDistribution[i].yValue;
				}
			}
			else {
				for (unsigned long i = 0; i < SymmetrizedDistribution.size()/2; i++){
					SymmetrizedDistribution[i].yValue = SymmetrizedDistribution[SymmetrizedDistribution.size()-1-i].yValue;
				}
			}
			normalizeDistribution(SymmetrizedDistribution);
			return SymmetrizedDistribution;
		}

	public:

		SeriesAnalyzer(unsigned long TotalNumberOfParticles):
			TotalNumberOfParticles(TotalNumberOfParticles){
			for (unsigned long i = 0; i <= TotalNumberOfParticles; i++){
				NADistribution.emplace_back(static_cast<double>(i)/static_cast<double>(TotalNumberOfParticles), 0.0);
			}
		}

		void addNewSeries(string FileNameNASeries, unsigned long MinNumberOfEquilibrationSweeps){
			vector<unsigned long> NewNASeries = readInSeries(FileNameNASeries);
			updateHistogramWithNewSeries(NewNASeries, MinNumberOfEquilibrationSweeps);
		}

		void normalizeNADistribution(){
			normalizeDistribution(NADistribution);
		}

		void writeNAProbabilityDistributionToFile(string FileName) const {
			writeDistributionToFile(FileName, NADistribution);
		}

		double computeFirstMomentOfHalfDistribution() const {
			vector<ValuePair> LeftDistribution;
			for (unsigned long i = 0; i <= TotalNumberOfParticles/2; i++){
				LeftDistribution.push_back(NADistribution[i]);
			}
			double LeftHalfIntegral = computeIntegral(LeftDistribution);

			vector<ValuePair> RightDistribution;
			for (unsigned long i = (TotalNumberOfParticles-1)/2+1; i < NADistribution.size(); i++){
				RightDistribution.push_back(NADistribution[i]);
			}
			double RightHalfIntegral = computeIntegral(RightDistribution);

			if (LeftHalfIntegral > RightHalfIntegral){
				normalizeDistribution(LeftDistribution);
				return computeFirstMomentOfDistribution(LeftDistribution);
			}
			normalizeDistribution(RightDistribution);
			return (1.0 - computeFirstMomentOfDistribution(RightDistribution));
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
			for (unsigned long i = 0; i <= (SymmetrizedDistribution.size()-1)/2; i++){
				HalfDistribution.push_back(SymmetrizedDistribution[i]);
			}
			normalizeDistribution(HalfDistribution);
			double Mean = computeFirstMomentOfDistribution(HalfDistribution);
			return static_cast<double>(TotalNumberOfParticles)*(computeSecondCentralMoment(SymmetrizedDistribution, Mean));
		}
};

#endif /* SERIES_ANALYSIS_INCLUDED */
