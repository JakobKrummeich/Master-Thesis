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

		static int computeEquilibrationIndex(const vector<double>& Series){
			for (int CurrentIndex = 0; CurrentIndex < Series.size(); CurrentIndex++){
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

		double computeMomentOfDistribution(const vector<ValuePair>& Distribution, double Exponent){
			double h = Distribution[1].xValue - Distribution[0].xValue;
			double Moment = 0.5 * (pow(Distribution[0].xValue,Exponent) * Distribution[0].yValue + pow(Distribution.back().xValue, Exponent) * Distribution.back().yValue);
			for (int i = 1; i < Distribution.size(); i++){
				Moment += pow(Distribution[i].xValue, Exponent)*Distribution[i].yValue;
			}
			return h*Moment;
		}

		void normalizeDistribution(vector<ValuePair>& Distribution){
			double TotalIntegral = computeMomentOfDistribution(Distribution, 0.0);
			for (int i = 0; i < Distribution.size(); i++){
				Distribution[i].yValue /= TotalIntegral;
			}
		}

		double computeFirstMomentInSubInterval(const vector<ValuePair>& Distribution, int LeftIntervalBoundary, int RightIntervalBoundary){
			vector<ValuePair> SubIntervalDistribution;
			for (int i = LeftIntervalBoundary; i <= RightIntervalBoundary; i++){
				SubIntervalDistribution.emplace_back(Distribution[i].xValue,Distribution[i].yValue);
			}
			normalizeDistribution(SubIntervalDistribution);
			return computeMomentOfDistribution(SubIntervalDistribution, 1.0);
		}

	public:

		SeriesAnalyzer(int TotalNumberOfParticles):
			TotalNumberOfParticles(TotalNumberOfParticles){
			for (int i = 0; i < TotalNumberOfParticles + 1; i++){
				NADistribution.emplace_back(static_cast<double>(i)/static_cast<double>(TotalNumberOfParticles), 0.0);
			}
		}
		
		void addNewSeries(string FileNameNASeries, string FileNamePotEnergySeries){
			vector<double> NewNASeries = readInSeries(FileNameNASeries);
			vector<double> NewPotEnergySeries = readInSeries(FileNamePotEnergySeries);
			int IndexConversionFactor = ceil(static_cast<double>(NewNASeries.size())/(static_cast<double>(NewPotEnergySeries.size())));
			int EquilibriumIndexNASeries = computeEquilibrationIndex(NewNASeries);
			int EquilibriumIndexPotEnergySeries = computeEquilibrationIndex(NewPotEnergySeries);
			int EquilibriumIndex = max(EquilibriumIndexNASeries, IndexConversionFactor*EquilibriumIndexPotEnergySeries);
			updateHistogramWithNewSeries(NewNASeries, EquilibriumIndex);
		}
		
		void normalizeNADistribution(){
			normalizeDistribution(NADistribution);
		}
		
		void writeNAProbabilityDistributionToFile(string FileName){
			ofstream FileStreamToWriteTo;
			FileStreamToWriteTo.open(FileName);
			FileStreamToWriteTo << "xA\tprobability_density\n";
			for (int i = 0; i < NADistribution.size(); i++){
				FileStreamToWriteTo << NADistribution[i].xValue << '\t' << NADistribution[i].yValue << '\n';
			}
			FileStreamToWriteTo.close();
		}

		double computeFirstMomentOfDistribution(){
			return computeMomentOfDistribution(NADistribution, 1.0);
		}
};

#endif /* SERIES_ANALYSIS_INCLUDED */
