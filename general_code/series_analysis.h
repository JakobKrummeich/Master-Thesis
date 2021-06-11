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

		int TotalNumberOfParticles;
		vector<ValuePair> NADistribution;
		
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
		
		void updateHistogramWithNewSeries(const vector<double>& NASeries, int EquilibrationIndex){
			for (int i = EquilibrationIndex; i < NASeries.size(); i++){
				NADistribution[static_cast<int>(NASeries[i])].yValue++;
			}
		}

		static int computeEquilibrationIndex(const vector<double>& Series, int MinNumberOfThrowAwayPoints) {
			for (int CurrentIndex = (MinNumberOfThrowAwayPoints < Series.size() ? MinNumberOfThrowAwayPoints : 0); CurrentIndex < Series.size(); CurrentIndex++){
				double AverageValue = 0.0;
				for (int j = CurrentIndex; j < Series.size(); j++){
					AverageValue += static_cast<double>(Series[j]);
				}
				AverageValue /= static_cast<double>(Series.size()-CurrentIndex);
				double AverageDeviation = 0.0;
				for (int j = CurrentIndex; j < Series.size(); j++){
					AverageDeviation += abs(AverageValue - static_cast<double>(Series[j]));
				}
				AverageDeviation /= static_cast<double>(Series.size()-CurrentIndex);
				if (abs(AverageValue - Series[CurrentIndex]) < AverageDeviation){
					return CurrentIndex;
				}
			}
			return 0;
		}
		
		static double computeCentralMomentOfDistribution(const vector<ValuePair>& Distribution, double Mean, double Exponent){
			double h = Distribution[1].xValue - Distribution[0].xValue;
			double Moment = 0.5 * (pow(Distribution[0].xValue - Mean,Exponent) * Distribution[0].yValue + pow(Distribution.back().xValue - Mean, Exponent) * Distribution.back().yValue);
			for (int i = 1; i < Distribution.size(); i++){
				Moment += pow(Distribution[i].xValue - Mean, Exponent)*Distribution[i].yValue;
			}
			return h*Moment;
		}
		
		static double computeMomentOfDistribution(const vector<ValuePair>& Distribution, double Exponent) {
			return computeCentralMomentOfDistribution(Distribution, 0.0, Exponent);
		}

		static void normalizeDistribution(vector<ValuePair>& Distribution) {
			double TotalIntegral = computeMomentOfDistribution(Distribution, 0.0);
			for (int i = 0; i < Distribution.size(); i++){
				Distribution[i].yValue /= TotalIntegral;
			}
		}

		static double computeFirstMomentInSubInterval(const vector<ValuePair>& Distribution, int LeftIntervalBoundary, int RightIntervalBoundary) {
			vector<ValuePair> SubIntervalDistribution;
			for (int i = LeftIntervalBoundary; i <= RightIntervalBoundary; i++){
				SubIntervalDistribution.emplace_back(Distribution[i].xValue,Distribution[i].yValue);
			}
			normalizeDistribution(SubIntervalDistribution);
			return computeMomentOfDistribution(SubIntervalDistribution, 1.0);
		}
		
		static void writeDistributionToFile(string FileName, const vector<ValuePair> Distribution) {
			ofstream FileStreamToWriteTo;
			FileStreamToWriteTo.open(FileName);
			FileStreamToWriteTo << "xA\tprobability_density\n";
			for (int i = 0; i < Distribution.size(); i++){
				FileStreamToWriteTo << Distribution[i].xValue << '\t' << Distribution[i].yValue << '\n';
			}
			FileStreamToWriteTo.close();
		}

	public:

		SeriesAnalyzer(int TotalNumberOfParticles):
			TotalNumberOfParticles(TotalNumberOfParticles){
			for (int i = 0; i < TotalNumberOfParticles + 1; i++){
				NADistribution.emplace_back(static_cast<double>(i)/static_cast<double>(TotalNumberOfParticles), 0.0);
			}
		}
		
		void addNewSeries(string FileNameNASeries, string FileNamePotEnergySeries, int MinNumberOfEquilibrationSweeps){
			vector<double> NewNASeries = readInSeries(FileNameNASeries);
			vector<double> NewPotEnergySeries = readInSeries(FileNamePotEnergySeries);
			int IndexConversionFactor = ceil(static_cast<double>(NewNASeries.size())/(static_cast<double>(NewPotEnergySeries.size())));
			int EquilibriumIndexNASeries = computeEquilibrationIndex(NewNASeries, MinNumberOfEquilibrationSweeps);
			int EquilibriumIndexPotEnergySeries = computeEquilibrationIndex(NewPotEnergySeries, MinNumberOfEquilibrationSweeps/IndexConversionFactor);
			int EquilibriumIndex = max(EquilibriumIndexNASeries, IndexConversionFactor*EquilibriumIndexPotEnergySeries);
			updateHistogramWithNewSeries(NewNASeries, EquilibriumIndex);
		}
		
		void normalizeNADistribution(){
			normalizeDistribution(NADistribution);
		}
		
		void writeNAProbabilityDistributionToFile(string FileName) const {
			writeDistributionToFile(FileName, NADistribution);
		}

		double computeFirstMomentOfHalfDistribution() const {
			vector<ValuePair> LeftDistribution;
			for (int i = 0; i < TotalNumberOfParticles/2; i++){
				LeftDistribution.push_back(NADistribution[i]);
			}
			double LeftHalfIntegral = computeMomentOfDistribution(LeftDistribution, 0.0);
			vector<ValuePair> RightDistribution;
			for (int i = TotalNumberOfParticles/2; i < NADistribution.size(); i++){
				RightDistribution.push_back(NADistribution[i]);
			}
			double RightHalfIntegral = computeMomentOfDistribution(RightDistribution, 0.0);
			if (LeftHalfIntegral > RightHalfIntegral){
				normalizeDistribution(LeftDistribution);
				return computeMomentOfDistribution(LeftDistribution, 1.0);
			}
			normalizeDistribution(RightDistribution);
			return (1.0 - computeMomentOfDistribution(RightDistribution, 1.0));
		}

		double computeBinderCumulant() const {
			vector<ValuePair> LeftDistribution;
			for (int i = 0; i < TotalNumberOfParticles/2; i++){
				LeftDistribution.push_back(NADistribution[i]);
			}
			double LeftHalfIntegral = computeMomentOfDistribution(LeftDistribution, 0.0);
			vector<ValuePair> RightDistribution;
			for (int i = TotalNumberOfParticles/2; i < NADistribution.size(); i++){
				RightDistribution.push_back(NADistribution[i]);
			}
			double RightHalfIntegral = computeMomentOfDistribution(RightDistribution, 0.0);
			vector<ValuePair> SymmetrizedDistribution(NADistribution);
			if (LeftHalfIntegral > RightHalfIntegral){
				for (int i = 0; i < SymmetrizedDistribution.size()/2; i++){
					SymmetrizedDistribution[SymmetrizedDistribution.size()-1-i].yValue = SymmetrizedDistribution[i].yValue;
				}
			}
			else {
				for (int i = 0; i < SymmetrizedDistribution.size()/2; i++){
					SymmetrizedDistribution[i].yValue = SymmetrizedDistribution[SymmetrizedDistribution.size()-1-i].yValue;
				}
			}
			normalizeDistribution(SymmetrizedDistribution);
			double SecondCentralMoment = computeCentralMomentOfDistribution(SymmetrizedDistribution, 0.5, 2.0);
			double FourthCentralMoment = computeCentralMomentOfDistribution(SymmetrizedDistribution, 0.5, 4.0);
			return (1.0 - FourthCentralMoment/(SecondCentralMoment*SecondCentralMoment*3.0));
		}
};

#endif /* SERIES_ANALYSIS_INCLUDED */
