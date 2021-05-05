#include <iostream>
#include <cmath>
#include <chrono>
#include <random>
#include <fstream>
#include <vector>
#include <iomanip>
#include <stdlib.h>
#include <omp.h>

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
vector<double> rABDifferences;


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

void computeABDistances(){
	for (int i = 0; i < NumberOfAParticles; i++){
		for (int j = 0; j < NumberOfBParticles; j++){
			for (int k = 0; k < DIMENSION; k++){
				double CoordinateDifference = APositions[DIMENSION * i + k] - BPositions[DIMENSION * j + k];
				rABDifferences.push_back(CoordinateDifference);
			}
		}
	}
}

void computeSelfDifferences(const vector<double>& Positions, vector<double>& SelfDifferences){
	SelfDifferences.clear();
	for (int i = 0; i < Positions.size(); i += DIMENSION){
		for (int j = 0; j < i; j += DIMENSION){
			for (int k = 0; k < DIMENSION; k++){
				double CoordinateDifference = Positions[i + k] - Positions[j + k];
				SelfDifferences.push_back(CoordinateDifference);
			}
		}
	}
}

double computeSelfStructureFactor(const vector<double>& Positions, double kx, double ky){
	double cosSum = 0.0;
	double sinSum = 0.0;
	for (int i = 0; i < Positions.size(); i += DIMENSION){
		cosSum += cos(kx*Positions[i]+ky*Positions[i+1]);
		sinSum += sin(kx*Positions[i]+ky*Positions[i+1]);
	}
	return (cosSum*cosSum+sinSum*sinSum)/static_cast<int>(TotalNumberOfParticles);
}

double computeAverageSelfStructureFactor(const vector<double>& Positions, double kMagnitudeCenter, double kMagnitudeWidth, int NumberOfAverageValues){
	double AverageStructureFactor = 0.0;
	for (int i = 0; i < NumberOfAverageValues; i++){
		double RandomAngle = RNG.drawRandomNumber(0.0, 2.0 * M_PI);
		double kMagnitude = RNG.drawRandomNumber(-kMagnitudeWidth,kMagnitudeWidth) + kMagnitudeCenter;
		AverageStructureFactor += computeSelfStructureFactor(Positions, kMagnitude*cos(RandomAngle), kMagnitude*sin(RandomAngle));
	}
	return AverageStructureFactor/static_cast<double>(NumberOfAverageValues);
}

double computeSelfStructureFactorAlternative(const vector<double>& SelfDifferences, int NumberOfParticles, double kx, double ky){
	double Result = static_cast<double>(NumberOfParticles)/static_cast<double>(TotalNumberOfParticles);
	for (int i = 0; i < SelfDifferences.size(); i += DIMENSION){
		Result += 2.0/static_cast<double>(TotalNumberOfParticles)*cos(kx*SelfDifferences[i]+ky*SelfDifferences[i+1]);
	}
	return Result;
}

double computeAverageSelfStructureFactorAlternative(const vector<double>& SelfDifferences, int NumberOfParticles, double kMagnitudeCenter, double kMagnitudeWidth, int NumberOfAverageValues){
	double AverageStructureFactor = 0.0;
	for (int i = 0; i < NumberOfAverageValues; i++){
		double RandomAngle = RNG.drawRandomNumber(0.0, 2.0 * M_PI);
		double kMagnitude = RNG.drawRandomNumber(-kMagnitudeWidth,kMagnitudeWidth) + kMagnitudeCenter;
		AverageStructureFactor += computeSelfStructureFactorAlternative(SelfDifferences, NumberOfParticles, kMagnitude*cos(RandomAngle), kMagnitude*sin(RandomAngle));
	}
	return AverageStructureFactor/static_cast<double>(NumberOfAverageValues);
}

double computeABStructureFactor(const vector<double>& ABDifferences, double kx, double ky){
	double Result = 0.0;
	for (int i = 0; i < ABDifferences.size(); i += DIMENSION){
		Result += cos(kx*ABDifferences[i]+ky*ABDifferences[i+1]);
	}
	return Result/static_cast<double>(TotalNumberOfParticles);
}

double computeAverageABStructureFactor(const vector<double>& ABDifferences, double kMagnitudeCenter, double kMagnitudeWidth, int NumberOfAverageValues){
	double AverageStructureFactor = 0.0;
	for (int i = 0; i < NumberOfAverageValues; i++){
		double RandomAngle = RNG.drawRandomNumber(0.0, 2.0 * M_PI);
		double kMagnitude = RNG.drawRandomNumber(-kMagnitudeWidth,kMagnitudeWidth) + kMagnitudeCenter;
		AverageStructureFactor += computeABStructureFactor(Positions, kMagnitude*cos(RandomAngle), kMagnitude*sin(RandomAngle));
	}
	return AverageStructureFactor/static_cast<double>(NumberOfAverageValues);
}

double computeConcentrationFactorValue(double AAStructureFactorValue, double BBStructureFactorValue, double ABStructureFactorValue){
	double xA = static_cast<double>(NumberOfAParticles)/static_cast<double>(TotalNumberOfParticles);
	double xB = static_cast<double>(NumberOfBParticles)/static_cast<double>(TotalNumberOfParticles);
	return xB*xB*AAStructureFactorValue + xA*xA*BBStructureFactorValue - 2.0 * xA * xB * ABStructureFactorValue;
}

int main(int argc, char* argv[]){
	BoxLength = atof(argv[2]);
	BoxLength = 35.3553;
	readInParticleState("FinalParticleConfig_N=1000_T=0.900000_AvgDens=0.600000_MCRuns=600000.dat");

	computeABDistances();

	double kMin = 2.0*M_PI/BoxLength;
	double kMax = 9.0;
	int NumberOfkValues = 400;
	double kDelta = (kMax - kMin)/static_cast<double>(NumberOfkValues);
	double kWidth = kDelta*0.5;
	double CurrentkMag = kMin;
	int NumberOfAveragesPerk = 2000;
	
	double AAStructureFactorValues [NumberOfkValues];
	double BBStructureFactorValues [NumberOfkValues];
	double ABStructureFactorValues [NumberOfkValues];
	double ConcentrationStructureValues [NumberOfkValues];
	double kValues [NumberOfkValues];

	for (int i = 0; i < NumberOfkValues; i++){
		kValues[i] = CurrentkMag;
		CurrentkMag += kDelta;
	}

	#pragma omp parallel num_threads(2)
	{
		#pragma omp critical 
		{
			cerr << omp_get_thread_num() << endl;
		}
		
		#pragma omp for
		for (int i = 0; i < NumberOfkValues; i++){
			if (NumberOfAParticles > 0){
				AAStructureFactorValues[i] = computeAverageSelfStructureFactor(APositions, kValues[i], kWidth, NumberOfAveragesPerk);
			}
			if (NumberOfBParticles > 0){
				BBStructureFactorValues[i] = computeAverageSelfStructureFactor(BPositions, kValues[i], kWidth, NumberOfAveragesPerk);
			}
			if (NumberOfAParticles > 0 && NumberOfBParticles > 0){
				ABStructureFactorValues[i] = computeAverageABStructureFactor(rABDifferences, kValues[i], kWidth, NumberOfAveragesPerk);
				ConcentrationStructureValues[i] = computeConcentrationFactorValue(AAStructureFactorValues[i], BBStructureFactorValues[i], ABStructureFactorValues[i]);
			}
		}
	}

	string FileName("structure_factor_T=0.9.dat");
	ofstream FileStreamToWrite;
	FileStreamToWrite.open(FileName);
	FileStreamToWrite << "k" << '\t' << "AAStructureFactor\t" << "BBStructureFactor\t" << "ABStructureFactor\t" << "ConcentrationStructureFactor\n";
	for (int i = 0; i < NumberOfkValues; i++){
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
	FileStreamToWrite.close();
}
