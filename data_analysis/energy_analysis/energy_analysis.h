#ifndef ENERGY_ANALYSIS_INCLUDED
#define ENERGY_ANALYSIS_INCLUDED

#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>

#include "../value_pair.h"
#include "../integration.h"

using namespace std;

class EnergyAnalyzer{

	private:

		vector<double> values;

		int numberOfBins;

		vector<ValuePair> energyDistribution;

		void readInSeries(string FilePath, int numberOfEquilibrationValues, int columnIndex){
			ifstream FileStreamToReadIn(FilePath);

			string CurrentString;
			getline(FileStreamToReadIn, CurrentString);

			for (int i = 0; i < numberOfEquilibrationValues; i++){
				getline(FileStreamToReadIn, CurrentString);
			}

			while (getline(FileStreamToReadIn, CurrentString)){
				istringstream iss(CurrentString);
				double newValue;
				for (int i = 0; i <= columnIndex; i++){
					iss >> newValue;
				}
				values.push_back(newValue);
			}
		}

		static void writeDistributionToFile(string FileName, const vector<ValuePair>& Distribution, string firstColumnName) {
			ofstream FileStreamToWriteTo(FileName);
			FileStreamToWriteTo << fixed << setprecision(numeric_limits<long double>::digits10+1) << firstColumnName << "\tprobability_density\n";
			for (unsigned long i = 0; i < Distribution.size(); i++){
				FileStreamToWriteTo << Distribution[i].xValue << '\t' << Distribution[i].yValue << '\n';
			}
		}

	public:

		EnergyAnalyzer(int numberOfBins):
			numberOfBins(numberOfBins)
		{
		}

		void addNewSeries(string FileNameNASeries, int numberOfEquilibrationValues, int columnIndex){
			readInSeries(FileNameNASeries, numberOfEquilibrationValues, columnIndex);
		}

		void writeProbabilityDistributionToFile(string FileName, string firstColumnName) const {
			writeDistributionToFile(FileName, energyDistribution, firstColumnName);
		}

		void generateDistribution() {
			double maxValue = values[0];
			double minValue = values[0];
			double averageValue = 0.0;
			for (int i = 0; i < values.size(); i++){
				if (values[i] > maxValue){
					maxValue = values[i];
				}
				else if (values[i] < minValue){
					minValue = values[i];
				}
			}
			double intervalWidth = (maxValue - minValue)/static_cast<double>(numberOfBins);
			double currentx = 0.5*intervalWidth+minValue;
			for (int i = 0; i < numberOfBins; i++, currentx+=intervalWidth){
				energyDistribution.emplace_back(currentx,0.0);
			}
			for (int i = 0; i < values.size(); i++){
				int index = static_cast<int>((values[i]-minValue)/intervalWidth);
				if (index >= energyDistribution.size()){
					index = energyDistribution.size() - 1;
				}
				energyDistribution[index].yValue++;
			}
			normalizeDistribution(energyDistribution);
		}

		double computeFirstMomentOfDistribution() const {
			return computeFirstMoment(energyDistribution);
		}

		double computeSecondMoment() const {
			double firstMoment = computeFirstMomentOfDistribution();
			return computeSecondCentralMoment(energyDistribution, firstMoment);
		}
};

#endif /* ENERGY_ANALYSIS_INCLUDED */
