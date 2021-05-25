#include <random>
#include <fstream>
#include <vector>

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
	BoxLength = 40.8248905;
	readInParticleState("data/FinalParticleConfig_N=1000_T=0.335000_AvgDens=0.600000_MCRuns=500000_epsAB=0.100000.dat");

	double kMin = 2.0*M_PI/BoxLength;
	double kMax = 9.0;
	int NumberOfkValueSubdivisions = 200;
	double kDelta = (kMax - kMin)/static_cast<double>(NumberOfkValueSubdivisions);
	double kWidth = kDelta*0.5;
	double CurrentkMag = kMin;
	int NumberOfAttemptedAveragesPerk = 1000;

	vector<ResultEntry> Results;

	vector<int> CombinationsFound;

	for (int i = 0; i < NumberOfkValueSubdivisions; i++){
		int NumberOfCombinationsFound = 0;
		CombinationsFound.clear();
		int NumberOfSuccessfulAverages = 0;
		ResultEntry NewEntry{};
		NewEntry.kMagnitude = CurrentkMag;
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
					for (int k = 0; k < (abs(nx) == abs(ny) ? 1 : 2); k++){
						for (int xSign = 1; xSign >= (nx == 0 ? 1 : -1); xSign-=2){
							for (int ySign = 1; ySign >= (ny == 0 ? 1 : -1); ySign-=2){
								NumberOfSuccessfulAverages++;
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
						int temp = nx;
						nx = ny;
						ny = temp;
					}
				}
			}
		}
		if (NumberOfSuccessfulAverages > 0){
			NewEntry.AAStructureFactor /= (static_cast<double>(NumberOfSuccessfulAverages)*static_cast<double>(TotalNumberOfParticles));
			NewEntry.BBStructureFactor /= (static_cast<double>(NumberOfSuccessfulAverages)*static_cast<double>(TotalNumberOfParticles));
			NewEntry.ABStructureFactor /= (static_cast<double>(NumberOfSuccessfulAverages)*static_cast<double>(TotalNumberOfParticles));
			NewEntry.ConcentrationStructureFactor = computeConcentrationFactorValue(NewEntry.AAStructureFactor,NewEntry.BBStructureFactor,NewEntry.ABStructureFactor);
			Results.push_back(NewEntry);
		}
		CurrentkMag += kDelta;
	}

	string FileName("structure_factor_T=0.335_Roh=0.6_epsAB=0.1.dat");
	ofstream FileStreamToWrite;
	FileStreamToWrite.open(FileName);
	FileStreamToWrite << "k" << '\t' << "AAStructureFactor\t" << "BBStructureFactor\t" << "ABStructureFactor\t" << "ConcentrationStructureFactor\n";
	for (int i = 0; i < Results.size(); i++){
		FileStreamToWrite << Results[i].kMagnitude << '\t';
		if (NumberOfAParticles > 0){
			FileStreamToWrite << Results[i].AAStructureFactor;
		}
		else {
			FileStreamToWrite << 0.0;
		}
		FileStreamToWrite << '\t';
		if (NumberOfBParticles > 0){
			FileStreamToWrite << Results[i].BBStructureFactor;
		}
		else {
			FileStreamToWrite << 0.0;
		}
		FileStreamToWrite << '\t';
		if (NumberOfAParticles > 0 && NumberOfBParticles > 0){
			FileStreamToWrite << Results[i].ABStructureFactor << '\t' << Results[i].ConcentrationStructureFactor;
		}
		else {
			FileStreamToWrite << 0.0 << '\t' << 0.0;
		}
		FileStreamToWrite << '\n';
	}
	FileStreamToWrite.close();
}
