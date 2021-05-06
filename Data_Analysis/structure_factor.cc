#include <iostream>
#include <cmath>
#include <chrono>
#include <random>
#include <fstream>
#include <vector>
#include <iomanip>
#include <stdlib.h>

using namespace std;

static const int DIMENSION = 2;

class realRNG{
	private:
		mt19937 rng;
		uniform_real_distribution<double> UnifRealDist;
	public:
		realRNG():
			UnifRealDist(0.0,1.0){
			random_device rd;
			seed_seq sd{rd(),rd()};
			rng = mt19937(sd);
		}

		double drawRandomNumber(){
			return UnifRealDist(rng);
		}

		double drawRandomNumber(double Min, double Max){
			return (Max - Min) * UnifRealDist(rng) + Min;
		}
};

realRNG RNG;
vector<double> Positions;
vector<int> ParticleTypes;
int NumberOfAParticles;
int NumberOfBParticles;
int TotalNumberOfParticles;
double BoxLength;
vector<double> APositions;
vector<double> BPositions;


void readInParticleState(string FileNameToReadIn) {
	ifstream FileStreamToReadIn;
	FileStreamToReadIn.open(FileNameToReadIn);

	string CurrentString;

	getline(FileStreamToReadIn, CurrentString, ':');
	getline(FileStreamToReadIn, CurrentString, '|');
	NumberOfAParticles =  stod(CurrentString);
	getline(FileStreamToReadIn, CurrentString, ':');
	getline(FileStreamToReadIn, CurrentString, '\n');
	NumberOfBParticles = stod(CurrentString);
	TotalNumberOfParticles = NumberOfAParticles + NumberOfBParticles;

	for (int ParticleIndex = 0; ParticleIndex < TotalNumberOfParticles; ParticleIndex++){
		getline(FileStreamToReadIn, CurrentString, '\t');
		double CurrentPosition [DIMENSION];
		for (int j = 0; j < DIMENSION; j++){
			getline(FileStreamToReadIn, CurrentString, '\t');
			CurrentPosition[j] = stod(CurrentString)*BoxLength;
			Positions.push_back(CurrentPosition[j]);
		}
		getline(FileStreamToReadIn, CurrentString, '\n');
		if (CurrentString == "A"){
			ParticleTypes.push_back(0);
			for (int j = 0; j < DIMENSION; j++){
				APositions.push_back(CurrentPosition[j]);
			}
		}
		else {
			ParticleTypes.push_back(1);
			for (int j = 0; j < DIMENSION; j++){
				BPositions.push_back(CurrentPosition[j]);
			}
		}
	}
	FileStreamToReadIn.close();
}

double computeCosSum(const vector<double>& Positions, double kx, double ky){
	double Sum = 0.0;
	for (int i = 0; i < Positions.size(); i += DIMENSION){
		Sum += cos(kx*Positions[i]+ky*Positions[i+1]);
	}
	return Sum;
}

double computeSinSum(const vector<double>& Positions, double kx, double ky){
	double Sum = 0.0;
	for (int i = 0; i < Positions.size(); i += DIMENSION){
		Sum += sin(kx*Positions[i]+ky*Positions[i+1]);
	}
	return Sum;
}

double computeConcentrationFactorValue(double AAStructureFactorValue, double BBStructureFactorValue, double ABStructureFactorValue){
	double xA = static_cast<double>(NumberOfAParticles)/static_cast<double>(TotalNumberOfParticles);
	double xB = static_cast<double>(NumberOfBParticles)/static_cast<double>(TotalNumberOfParticles);
	return xB*xB*AAStructureFactorValue + xA*xA*BBStructureFactorValue - 2.0 * xA * xB * ABStructureFactorValue;
}

int main(int argc, char* argv[]){
	BoxLength = atof(argv[2]);
	BoxLength = 40.8248905;
	readInParticleState("FinalParticleConfig_N=1000_T=0.335000_AvgDens=0.600000_MCRuns=500000_epsAB=0.100000.dat");

	double kMin = 2.0*M_PI/BoxLength;
	double kMax = 9.0;
	int NumberOfkValues = 400;
	double kDelta = (kMax - kMin)/static_cast<double>(NumberOfkValues);
	double kWidth = kDelta*0.5;
	double CurrentkMag = kMin;
	int NumberOfAveragesPerk = 1000;
	
	double AAStructureFactorValues [NumberOfkValues];
	double BBStructureFactorValues [NumberOfkValues];
	double ABStructureFactorValues [NumberOfkValues];
	double ConcentrationStructureValues [NumberOfkValues];
	double kValues [NumberOfkValues];
	bool kValuesOnGridFound [NumberOfkValues];

	for (int i = 0; i < NumberOfkValues; i++){
		kValues[i] = CurrentkMag;
		CurrentkMag += kDelta;
	}

	for (int i = 0; i < NumberOfkValues; i++){
		double AverageStructureFactorAA = 0.0;
		double AverageStructureFactorBB = 0.0;
		double AverageStructureFactorAB = 0.0;
		int NumberOfSuccessfulAverages = 0;
		for (int j = 0; j < NumberOfAveragesPerk; j++){
			double RandomAngle = RNG.drawRandomNumber(0.0, 2.0 * M_PI);
			double kMagnitude = RNG.drawRandomNumber(-kWidth,kWidth) + kValues[i];
			int nx = round(BoxLength*kMagnitude*cos(RandomAngle)/(2.0*M_PI));
			int ny = round(BoxLength*kMagnitude*sin(RandomAngle)/(2.0*M_PI));
			double GridkMagnitude = 2.0*M_PI/BoxLength*sqrt(static_cast<double>(nx*nx+ny*ny));
			if (GridkMagnitude <= (kValues[i]+kWidth) && GridkMagnitude >= (kValues[i]-kWidth)){
				NumberOfSuccessfulAverages++;
				double kx = nx * 2.0*M_PI/BoxLength;
				double ky = ny * 2.0*M_PI/BoxLength;
				double cosSumA;
				double sinSumA;
				double cosSumB;
				double sinSumB;
				if (NumberOfAParticles > 0){
					cosSumA = computeCosSum(APositions, kx, ky);
					sinSumA = computeSinSum(APositions, kx, ky);
					AverageStructureFactorAA += (cosSumA*cosSumA+sinSumA*sinSumA);
				}
				if (NumberOfBParticles > 0){
					cosSumB = computeCosSum(BPositions, kx, ky);
					sinSumB = computeSinSum(BPositions, kx, ky);
					AverageStructureFactorBB += (cosSumB*cosSumB+sinSumB*sinSumB);
				}
				if (NumberOfAParticles > 0 && NumberOfBParticles > 0){
					AverageStructureFactorAB += (cosSumA*cosSumB + sinSumA*sinSumB);
				}
			}
		}
		if (NumberOfSuccessfulAverages > 0){
			kValuesOnGridFound[i] = true;
			AAStructureFactorValues[i] = AverageStructureFactorAA / (static_cast<double>(NumberOfSuccessfulAverages)*static_cast<double>(TotalNumberOfParticles));
			BBStructureFactorValues[i] = AverageStructureFactorBB / (static_cast<double>(NumberOfSuccessfulAverages)*static_cast<double>(TotalNumberOfParticles));
			ABStructureFactorValues[i] = AverageStructureFactorAB / (static_cast<double>(NumberOfSuccessfulAverages)*static_cast<double>(TotalNumberOfParticles));
			ConcentrationStructureValues [i] = computeConcentrationFactorValue(AAStructureFactorValues[i],BBStructureFactorValues[i],ABStructureFactorValues[i]);
		}
		else {
			kValuesOnGridFound[i] = false;
		}
	}

	string FileName("structure_factor_T=0.335_Roh=0.6_epsAB=0.1.dat");
	ofstream FileStreamToWrite;
	FileStreamToWrite.open(FileName);
	FileStreamToWrite << "k" << '\t' << "AAStructureFactor\t" << "BBStructureFactor\t" << "ABStructureFactor\t" << "ConcentrationStructureFactor\n";
	for (int i = 0; i < NumberOfkValues; i++){
		if (kValuesOnGridFound[i]){
			FileStreamToWrite << kValues[i] << '\t';
			if (NumberOfAParticles > 0){
				FileStreamToWrite << AAStructureFactorValues[i];
			}
			else {
				FileStreamToWrite << 0.0;
			}
			FileStreamToWrite << '\t';
			if (NumberOfBParticles > 0){
				FileStreamToWrite << BBStructureFactorValues[i];
			}
			else {
				FileStreamToWrite << 0.0;
			}
			FileStreamToWrite << '\t';
			if (NumberOfAParticles > 0 && NumberOfBParticles > 0){
				FileStreamToWrite << ABStructureFactorValues[i] << '\t' << ConcentrationStructureValues[i];
			}
			else {
				FileStreamToWrite << 0.0 << '\t' << 0.0;
			}
			FileStreamToWrite << '\n';
		}
	}
	FileStreamToWrite.close();
}
