#include "../General_Code/realRNG.h"
#include <fstream>
#include <vector>
#include <algorithm>

using namespace std;

static const int DIMENSION = 2;

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

double computeConcentrationFactorValue(double AAStructureFactor, double BBStructureFactor, double ABStructureFactor){
	double xA = static_cast<double>(NumberOfAParticles)/static_cast<double>(TotalNumberOfParticles);
	double xB = static_cast<double>(NumberOfBParticles)/static_cast<double>(TotalNumberOfParticles);
	return xB*xB*AAStructureFactor + xA*xA*BBStructureFactor - 2.0 * xA * xB * ABStructureFactor;
}

struct ResultEntry{
	double kMagnitude;
	double AAStructureFactor;
	double BBStructureFactor;
	double ABStructureFactor;
	double ConcentrationStructureFactor;
};

int main(int argc, char* argv[]){
	BoxLength = atof(argv[2]);
	BoxLength = 40.8248;
	readInParticleState("unsorted_data/FinalParticleConfig_N=1000_T=1.000000_AvgDens=0.600000_MCRuns=10000_epsAB=0.100000.dat");

	double kMin = 2.0*M_PI/BoxLength;
	double kMax = 9.0;
	int NumberOfkValueSubdivisions = 400;
	double kDelta = (kMax - kMin)/static_cast<double>(NumberOfkValueSubdivisions);
	double kWidth = kDelta*0.5;
	double CurrentkMag = kMin;
	int NumberOfAttemptedAveragesPerk = 1000;

	vector<ResultEntry> Results;

	vector<int> CombinationsFound;
	vector<int> MagnitudesFound;
	vector<ResultEntry> IntermediateResults;
	vector<int> NumberOfDifferentkPerMagnitude;

	for (int i = 0; i < NumberOfkValueSubdivisions; i++){
		int NumberOfCombinationsFound = 0;
		CombinationsFound.clear();
		MagnitudesFound.clear();
		IntermediateResults.clear();
		NumberOfDifferentkPerMagnitude.clear();

		for (int j = 0; j < NumberOfAttemptedAveragesPerk; j++){
			double RandomAngle = RNG.drawRandomNumber(0.0, 0.5 * M_PI);
			double kMagnitude = RNG.drawRandomNumber(-kWidth,kWidth) + CurrentkMag;
			int nx = round(BoxLength*kMagnitude*cos(RandomAngle)/(2.0*M_PI));
			int ny = round(BoxLength*kMagnitude*sin(RandomAngle)/(2.0*M_PI));
			double GridkMagnitude = 2.0*M_PI/BoxLength*sqrt(static_cast<double>(nx*nx+ny*ny));
			if (GridkMagnitude <= (CurrentkMag+kWidth) && GridkMagnitude >= (CurrentkMag-kWidth) && GridkMagnitude >= kMin){
				bool sameCombinationAlreadyFoundBefore = false;
				for (int k = 0; k < NumberOfCombinationsFound && !sameCombinationAlreadyFoundBefore; k++){
					if ((CombinationsFound[DIMENSION*k] == nx && CombinationsFound[DIMENSION*k+1] == ny) || (CombinationsFound[DIMENSION*k] == ny && CombinationsFound[DIMENSION*k+1] == nx)){
						sameCombinationAlreadyFoundBefore = true;
					}
				}
				if (!sameCombinationAlreadyFoundBefore){
					NumberOfCombinationsFound++;
					CombinationsFound.push_back(nx);
					CombinationsFound.push_back(ny);
					int NewGridMagnitude = nx*nx+ny*ny;
					bool sameMagnitudeFoundBefore = false;
					int IndexOfMagnitude;
					for (int k = 0; k < MagnitudesFound.size(); k++){
						if (NewGridMagnitude == MagnitudesFound[k]){
							sameMagnitudeFoundBefore = true;
							IndexOfMagnitude = k;
							break;
						}
					}
					ResultEntry NewEntry{};
					NewEntry.kMagnitude = GridkMagnitude;
					int NumberOfkWithThisCombination = 0;
					for (int k = 0; k < (abs(nx) == abs(ny) ? 1 : 2); k++){
						for (int xSign = 1; xSign >= (nx == 0 ? 1 : -1); xSign-=2){
							for (int ySign = 1; ySign >= (ny == 0 ? 1 : -1); ySign-=2){
								NumberOfkWithThisCombination++;
								double kx = xSign * nx * 2.0*M_PI/BoxLength;
								double ky = ySign * ny * 2.0*M_PI/BoxLength;
								double cosSumA;
								double sinSumA;
								double cosSumB;
								double sinSumB;
								if (NumberOfAParticles > 0){
									cosSumA = computeCosSum(APositions, kx, ky);
									sinSumA = computeSinSum(APositions, kx, ky);
									NewEntry.AAStructureFactor += (cosSumA*cosSumA+sinSumA*sinSumA);
								}
								if (NumberOfBParticles > 0){
									cosSumB = computeCosSum(BPositions, kx, ky);
									sinSumB = computeSinSum(BPositions, kx, ky);
									NewEntry.BBStructureFactor += (cosSumB*cosSumB+sinSumB*sinSumB);
								}
								if (NumberOfAParticles > 0 && NumberOfBParticles > 0){
									NewEntry.ABStructureFactor += (cosSumA*cosSumB + sinSumA*sinSumB);
								}
							}
						}
						swap(nx,ny);
					}
					if (!sameMagnitudeFoundBefore){
						MagnitudesFound.push_back(NewGridMagnitude);
						NumberOfDifferentkPerMagnitude.push_back(NumberOfkWithThisCombination);
						IntermediateResults.push_back(NewEntry);
						for (int CurrentIndex = MagnitudesFound.size() - 1; CurrentIndex > 0 && MagnitudesFound[CurrentIndex-1] > MagnitudesFound[CurrentIndex]; CurrentIndex--){
							swap(MagnitudesFound[CurrentIndex], MagnitudesFound[CurrentIndex-1]);
							swap(IntermediateResults[CurrentIndex], IntermediateResults[CurrentIndex-1]);
							swap(NumberOfDifferentkPerMagnitude[CurrentIndex], NumberOfDifferentkPerMagnitude[CurrentIndex-1]);
						}
					}
					else {
						IntermediateResults[IndexOfMagnitude].AAStructureFactor += NewEntry.AAStructureFactor;
						IntermediateResults[IndexOfMagnitude].BBStructureFactor += NewEntry.BBStructureFactor;
						IntermediateResults[IndexOfMagnitude].ABStructureFactor += NewEntry.ABStructureFactor;
						NumberOfDifferentkPerMagnitude[IndexOfMagnitude] += NumberOfkWithThisCombination;
					}
				}
			}
		}
		for (int i = 0; i < IntermediateResults.size(); i++){
			IntermediateResults[i].AAStructureFactor /= (static_cast<double>(NumberOfDifferentkPerMagnitude[i])*static_cast<double>(TotalNumberOfParticles));
			IntermediateResults[i].BBStructureFactor /= (static_cast<double>(NumberOfDifferentkPerMagnitude[i])*static_cast<double>(TotalNumberOfParticles));
			IntermediateResults[i].ABStructureFactor /= (static_cast<double>(NumberOfDifferentkPerMagnitude[i])*static_cast<double>(TotalNumberOfParticles));
			IntermediateResults[i].ConcentrationStructureFactor = computeConcentrationFactorValue(IntermediateResults[i].AAStructureFactor,IntermediateResults[i].BBStructureFactor,IntermediateResults[i].ABStructureFactor);
			Results.push_back(IntermediateResults[i]);
		}
		CurrentkMag += kDelta;
	}

	vector<ResultEntry> SmoothedResults;
	double SmoothingWindowSize = 0.5;
	for (int i = 0; i < Results.size(); i++){
		ResultEntry NewEntry = Results[i];
		int AveragedValues = 1;
		int SearchIndex = i-1;
		while (SearchIndex >= 0 && Results[i].kMagnitude - Results[SearchIndex].kMagnitude < SmoothingWindowSize*0.5){
			NewEntry.AAStructureFactor += Results[SearchIndex].AAStructureFactor;
			NewEntry.BBStructureFactor += Results[SearchIndex].BBStructureFactor;
			NewEntry.ABStructureFactor += Results[SearchIndex].ABStructureFactor;
			NewEntry.ConcentrationStructureFactor += Results[SearchIndex].ConcentrationStructureFactor;
			SearchIndex--;
			AveragedValues++;
		}
		SearchIndex = i+1;
		while (SearchIndex < Results.size() && Results[SearchIndex].kMagnitude - Results[i].kMagnitude < SmoothingWindowSize*0.5){
			NewEntry.AAStructureFactor += Results[SearchIndex].AAStructureFactor;
			NewEntry.BBStructureFactor += Results[SearchIndex].BBStructureFactor;
			NewEntry.ABStructureFactor += Results[SearchIndex].ABStructureFactor;
			NewEntry.ConcentrationStructureFactor += Results[SearchIndex].ConcentrationStructureFactor;
			AveragedValues++;
			SearchIndex++;
		}
		NewEntry.AAStructureFactor /= static_cast<double>(AveragedValues);
		NewEntry.BBStructureFactor /= static_cast<double>(AveragedValues);
		NewEntry.ABStructureFactor /= static_cast<double>(AveragedValues);
		NewEntry.ConcentrationStructureFactor /= static_cast<double>(AveragedValues);
		SmoothedResults.push_back(NewEntry);
	}
	

	string FileName("structure_factor_T=1.0_Roh=0.6_epsAB=0.1.dat");
	ofstream FileStreamToWrite;
	FileStreamToWrite.open(FileName);
	FileStreamToWrite << "k" << '\t' << "AAStructureFactor\tBBStructureFactor\tABStructureFactor\tConcentrationStructureFactor\tSmoothedAAStructureFactor\tSmoothedBBStructureFactor\tSmoothedABStructureFactor\tSmoothedConcentrationStructureFactor\n";
	for (int i = 0; i < Results.size(); i++){
		FileStreamToWrite << Results[i].kMagnitude << '\t' << Results[i].AAStructureFactor << '\t' << Results[i].BBStructureFactor << '\t' << Results[i].ABStructureFactor << '\t' << Results[i].ConcentrationStructureFactor << '\t' << SmoothedResults[i].AAStructureFactor << '\t' << SmoothedResults[i].BBStructureFactor << '\t' << SmoothedResults[i].ABStructureFactor << '\t' << SmoothedResults[i].ConcentrationStructureFactor << '\n';
	}
	FileStreamToWrite.close();
}
