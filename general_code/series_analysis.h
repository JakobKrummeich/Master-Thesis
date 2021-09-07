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

		void updateHistogramWithNewSeries(const vector<double>& NASeries, unsigned long MinNumberOfSweeps ){
			for (unsigned long i = MinNumberOfSweeps; i < NASeries.size(); i++){
				NADistribution[static_cast<unsigned long>(NASeries[i])].yValue++;
			}
		}

		static double computeCentralMomentOfDistribution(const vector<ValuePair>& Distribution, double Mean, double Exponent){
			double h = Distribution[1].xValue - Distribution[0].xValue;
			double Moment = 0.5 * (pow(Distribution[0].xValue - Mean, Exponent) * Distribution[0].yValue + pow(Distribution.back().xValue - Mean, Exponent) * Distribution.back().yValue);
			for (unsigned long i = 1; i < Distribution.size() - 1; i++){
				Moment += pow(Distribution[i].xValue - Mean, Exponent)*Distribution[i].yValue;
			}
			return h*Moment;
		}

		static double computeMomentOfDistribution(const vector<ValuePair>& Distribution, double Exponent) {
			return computeCentralMomentOfDistribution(Distribution, 0.0, Exponent);
		}

		static void normalizeDistribution(vector<ValuePair>& Distribution) {
			double TotalIntegral = computeMomentOfDistribution(Distribution, 0.0);
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
			double LeftHalfIntegral = computeMomentOfDistribution(LeftDistribution, 0.0);

			vector<ValuePair> RightDistribution;
			for (unsigned long i = Distribution.size()/2; i < Distribution.size(); i++){
				RightDistribution.push_back(Distribution[i]);
			}
			double RightHalfIntegral = computeMomentOfDistribution(RightDistribution, 0.0);

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
			vector<double> NewNASeries = readInSeries(FileNameNASeries);
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
			double LeftHalfIntegral = computeMomentOfDistribution(LeftDistribution, 0.0);

			vector<ValuePair> RightDistribution;
			for (unsigned long i = (TotalNumberOfParticles-1)/2+1; i < NADistribution.size(); i++){
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

			vector<ValuePair> SymmetrizedDistribution = symmetrizeDistribution(NADistribution);

			double SecondCentralMoment = computeCentralMomentOfDistribution(SymmetrizedDistribution, 0.5, 2.0);
			double FourthCentralMoment = computeCentralMomentOfDistribution(SymmetrizedDistribution, 0.5, 4.0);
			return (1.0 - FourthCentralMoment/(SecondCentralMoment*SecondCentralMoment*3.0));
		}

		double computeStructureFactor(bool RestrictToHalfDistribution) const {
			vector<ValuePair> SymmetrizedDistribution = symmetrizeDistribution(NADistribution);
			if (!RestrictToHalfDistribution){
				return static_cast<double>(TotalNumberOfParticles)*(computeCentralMomentOfDistribution(SymmetrizedDistribution, 0.5, 2.0));
			}

			vector<ValuePair> HalfDistribution;
			for (unsigned long i = 0; i <= (SymmetrizedDistribution.size()-1)/2; i++){
				HalfDistribution.push_back(SymmetrizedDistribution[i]);
			}
			normalizeDistribution(HalfDistribution);
			double Mean = computeMomentOfDistribution(HalfDistribution, 1.0);
			return static_cast<double>(TotalNumberOfParticles)*(computeCentralMomentOfDistribution(SymmetrizedDistribution, Mean, 2.0));
		}
};

#endif /* SERIES_ANALYSIS_INCLUDED */
