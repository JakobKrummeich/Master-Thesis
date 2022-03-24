#ifndef ENERGY_ANALYSIS_INCLUDED
#define ENERGY_ANALYSIS_INCLUDED

#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <fstream>

#include "../value_pair.h"

using namespace std;

class EnergyAnalyzer{

	private:

		vector<double> values;

		vector<ValuePair> energyDistribution;

		vector<unsigned long> readInSeries(string FilePath, int numberOfEquilibrationValues){
			ifstream FileStreamToReadIn(FilePath);

			string CurrentString;
			getline(FileStreamToReadIn, CurrentString);

			for (int i = 0; i < numberOfEquilibrationValues; i++){
				getline(FileStreamToReadIn, CurrentString);
			}

			while (getline(FileStreamToReadIn, CurrentString)){
				values.push_back(stod(CurrentString));
			}
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

#endif /* ENERGY_ANALYSIS_INCLUDED */
