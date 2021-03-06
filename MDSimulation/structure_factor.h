#ifndef STRUCTURE_FACTOR_INCLUDED
#define STRUCTURE_FACTOR_INCLUDED

#include <vector>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <iostream>
#include <cmath>

#include "../realRNG.h"
#include "../utility_functions.h"

using namespace std;

static const int DIMENSION = 2;

static const int NUMBER_OF_THREADS = 2;

class StructureFactorComputator{
	private:

		realUniformRNG RNG;

		double BoxLength;
		double kAll;
		int numberOfAngleIntervals;
		double angleSpacing;
		int numberOfkIntervals;
		double kSpacing;

		vector<double> APositions;
		vector<double> BPositions;
		int TotalNumberOfParticles;
		int NumberOfAParticles;
		int NumberOfBParticles;

		struct kVector {
			double kx;
			double ky;
			double kMagnitude;
			double angle;
			int z1;
			int z2;
		};

		struct Entry {
			double kMag;
			double angle;
			double SAA;
			double SBB;
			double SAB;
			double Scc;
			int NumberOfDataPoints;
			bool kVectorsFound;
			vector<kVector> kVectors;

			Entry(double kMag, double angle):
				kMag(kMag),
				angle(angle),
				SAA(0.0),
				SBB(0.0),
				SAB(0.0),
				Scc(0.0),
				kVectorsFound(false),
				NumberOfDataPoints(0){
			}
		};

		vector<Entry> Results;

		bool kVectorFoundAlready(int z1, int z2, int indexInResults) const {
			for (int i = 0; i < Results[indexInResults].kVectors.size(); i++){
				if (z1 == Results[indexInResults].kVectors[i].z1 && z2 == Results[indexInResults].kVectors[i].z2){
					return true;
				}
			}
			return false;
		}

		double computeCosSum(const vector<double>& Positions, double kx, double ky) const{
			double Sum = 0.0;
			for (int i = 0; i < Positions.size(); i += DIMENSION){
				Sum += cos( kx * Positions[i] + ky * Positions[i+1]);
			}
			return Sum;
		}

		double computeSinSum(const vector<double>& Positions, double kx, double ky) const{
			double Sum = 0.0;
			for (int i = 0; i < Positions.size(); i += DIMENSION){
				Sum += sin( kx * Positions[i] + ky * Positions[i+1]);
			}
			return Sum;
		}

		void findkValuesOnGrid(double xDisplacement){

			const double kMin = 4.0*M_PI/BoxLength; // minimum k Value given by half box length due to periodic boundary conditions
			const double k1 [DIMENSION]{2.0*M_PI/BoxLength, -2.0*M_PI*xDisplacement/(BoxLength*BoxLength)}; //lattice basis in reciprocal lattice
			const double k2 [DIMENSION]{0.0, 2.0*M_PI/BoxLength};

			int Maxn = ceil(BoxLength*kAll*0.5/(M_PI));
			for (int z2 = -Maxn; z2 <= Maxn; z2++){
				for (int z1 = 0; z1 <= Maxn; z1++){
					if (z1 != 0 || z2 < 0){
						double kx = static_cast<double>(z1)*k1[0]+static_cast<double>(z2)*k2[0];
						double ky = static_cast<double>(z1)*k1[1]+static_cast<double>(z2)*k2[1];
						double magnitude = sqrt(kx*kx+ky*ky);
						if (magnitude > kAll){
						 break;
						}
						else if (magnitude >= kMin){
							kVector newkVector;
							newkVector.kx = kx;
							newkVector.ky = ky;
							newkVector.kMagnitude = magnitude;
							newkVector.z1 = z1;
							newkVector.z2 = z2;
							if (kx == 0.0){
								newkVector.angle = -M_PI*0.5;
							}
							else {
								newkVector.angle = atan2(ky,kx);
							}
							double shiftedAngle = newkVector.angle;
							shiftedAngle += 0.5*angleSpacing+0.5*M_PI;
							while (shiftedAngle < 0.0){
								shiftedAngle += M_PI;
							}
							while (shiftedAngle >= M_PI){
								shiftedAngle -= M_PI;
							}
							int angleIndex = static_cast<int>(shiftedAngle/angleSpacing);
							if (angleIndex < 0 || angleIndex >= numberOfAngleIntervals){
								cerr << "ERROR! angleIndex is wrong" << endl;
							}
							int kMagIndex = static_cast<int>( (newkVector.kMagnitude - kMin + 0.5*kSpacing)/(kSpacing));
							if (kMagIndex < 0 || kMagIndex >= numberOfkIntervals){
								cerr << "ERROR! kMagIndex is wrong" << endl;
							}
							int indexInResults = numberOfkIntervals * angleIndex + kMagIndex;
							Results[indexInResults].kVectors.push_back(newkVector);
							Results[indexInResults].kVectorsFound = true;
						}
					}
				}
			}

			const int numberOfTriesPerBin = 10*static_cast<int>(ceil(angleSpacing*kSpacing*BoxLength*BoxLength*kAll*0.25/(M_PI*M_PI)));

			for (int i = 0; i < Results.size(); i++){
				double kMag = Results[i].kMag;
				if (kMag >= kAll){
					double angle = Results[i].angle;
					double kMin = kMag - kSpacing*0.5;
					double kMax = kMag + kSpacing*0.5;
					double angleMin = angle - angleSpacing*0.5;
					double angleMax = angle + angleSpacing*0.5;
					for (int i = 0; i < numberOfTriesPerBin; i++){
						double randomkMag = RNG.drawRandomNumber(kMin, kMax);
						double randomAngle = RNG.drawRandomNumber(angleMin, angleMax);
						double randomkx = randomkMag * cos(randomAngle);
						double randomky = randomkMag * sin(randomAngle);
						int z1 = round(randomkx*BoxLength*0.5 / (M_PI));
						int z2 = round(randomky*BoxLength*0.5 / (M_PI) + static_cast<double>(z1)*xDisplacement/(BoxLength));
						double kxLattice = static_cast<double>(z1)*2.0*M_PI/BoxLength;
						double kyLattice = static_cast<double>(z1)*2.0*M_PI*(-xDisplacement)/(BoxLength*BoxLength) + static_cast<double>(z2)*2.0*M_PI/BoxLength;

						kVector newkVector;
						newkVector.kx = kxLattice;
						newkVector.ky = kyLattice;
						newkVector.kMagnitude = sqrt(kxLattice*kxLattice + kyLattice*kyLattice);
						newkVector.z1 = z1;
						newkVector.z2 = z2;
						if (newkVector.kMagnitude < kMax && newkVector.kMagnitude >= kMin){
							if (kxLattice == 0.0){
								newkVector.angle = -M_PI*0.5;
							}
							else {
								newkVector.angle = atan2(kyLattice,kxLattice);
							}
							double shiftedAngle = newkVector.angle;
							shiftedAngle += 0.5*angleSpacing+0.5*M_PI;
							while (shiftedAngle < 0.0){
								shiftedAngle += M_PI;
							}
							while (shiftedAngle >= M_PI){
								shiftedAngle -= M_PI;
							}
							int angleIndex = static_cast<int>(shiftedAngle/angleSpacing);
							if (angleIndex < 0 || angleIndex >= numberOfAngleIntervals){
								cerr << "ERROR! angleIndex is wrong" << endl;
							}
							int kMagIndex = static_cast<int>( (newkVector.kMagnitude - 4.0*M_PI/BoxLength + 0.5*kSpacing)/(kSpacing));
							if (kMagIndex < 0 || kMagIndex >= numberOfkIntervals){
								cerr << "ERROR! kMagIndex is wrong" << endl;
							}
							int indexInResults = numberOfkIntervals * angleIndex + kMagIndex;
							if (!kVectorFoundAlready(z1,z2,indexInResults)){
								Results[indexInResults].kVectors.push_back(newkVector);
								Results[indexInResults].kVectorsFound = true;
							}
						}
					}
				}
			}
		}

		bool checkIfAllKVectorsAreInCorrectBin() const {
			for (int i = 0; i < Results.size(); i++){
				double angle = Results[i].angle;
				double kMag = Results[i].kMag;
				double kMin = kMag - kSpacing*0.5;
				double kMax = kMag + kSpacing*0.5;
				double angleMin = angle - angleSpacing*0.5;
				double angleMax = angle + angleSpacing*0.5;
				for (int j = 0; j < Results[i].kVectors.size(); j++){
					double kx = Results[i].kVectors[j].kx;
					double ky = Results[i].kVectors[j].ky;
					double mag = sqrt(kx*kx + ky*ky);
					if (mag < kMin || mag >= kMax){
						return false;
					}
					double actualAngle;
					if (kx == 0.0){
						actualAngle = -0.5*M_PI;
					}
					else {
						actualAngle = atan2(ky,kx);
					}
					while (actualAngle < -M_PI*0.5-angleSpacing*0.5){
						actualAngle += M_PI;
					}
					while (actualAngle >= M_PI*0.5-angleSpacing*0.5){
						actualAngle -= M_PI;
					}
					if (actualAngle < angleMin || actualAngle >= angleMax){
						return false;
					}
				}
			}
			return true;
		}

	public:

		StructureFactorComputator(int TotalNumberOfParticles, double BoxLength, double kAll, double kMax, int numberOfkValues, int numberOfAngleBins):
			TotalNumberOfParticles(TotalNumberOfParticles),
			kAll(kAll),
			BoxLength(BoxLength),
			numberOfkIntervals(numberOfkValues),
			kSpacing((kMax - 4.0*M_PI/(BoxLength))/static_cast<double>(numberOfkValues)), // start at kMin = 4*M_PI/(L)
			numberOfAngleIntervals(numberOfAngleBins),
			angleSpacing(M_PI/static_cast<double>(numberOfAngleBins)){
				double currentAngle = -M_PI*0.5;
				for (int angleIndex = 0; angleIndex < numberOfAngleBins; angleIndex++, currentAngle += angleSpacing){
					for (double currentkMag = 4.0*M_PI/BoxLength, kIndex = 0; kIndex < numberOfkValues; currentkMag += kSpacing, kIndex++){
						Results.emplace_back(currentkMag, currentAngle);
					}
				}
			}

		void readInParticleState(string FileNameToReadIn, double& xDisplacement) {
			ifstream FileStreamToReadIn;
			FileStreamToReadIn.open(FileNameToReadIn);

			string CurrentString;

			APositions.clear();
			BPositions.clear();

			getline(FileStreamToReadIn, CurrentString, ':');
			getline(FileStreamToReadIn, CurrentString, '|');
			NumberOfAParticles = stoi(CurrentString);

			getline(FileStreamToReadIn, CurrentString, ':');
			getline(FileStreamToReadIn, CurrentString, '|');
			NumberOfBParticles = stoi(CurrentString);

			getline(FileStreamToReadIn, CurrentString, ':');
			getline(FileStreamToReadIn, CurrentString, '|');
			BoxLength = stod(CurrentString);

			getline(FileStreamToReadIn, CurrentString, ':');
			getline(FileStreamToReadIn, CurrentString, '\n');
			xDisplacement = stod(CurrentString);

			if (xDisplacement > 0.5){
				xDisplacement -= 1.0;
			}
			xDisplacement *= BoxLength;

			for (int ParticleIndex = 0; ParticleIndex < TotalNumberOfParticles; ParticleIndex++){
				double CurrentPosition [DIMENSION];
				for (int j = 0; j < DIMENSION; j++){
					getline(FileStreamToReadIn, CurrentString, '\t');
					CurrentPosition[j] = stod(CurrentString)*BoxLength;
				}
				getline(FileStreamToReadIn, CurrentString, '\t');
				getline(FileStreamToReadIn, CurrentString, '\t');
				getline(FileStreamToReadIn, CurrentString, '\n');
				if (CurrentString == "A"){
					for (int j = 0; j < DIMENSION; j++){
						APositions.push_back(CurrentPosition[j]);
					}
				}
				else {
					for (int j = 0; j < DIMENSION; j++){
						BPositions.push_back(CurrentPosition[j]);
					}
				}
			}
			FileStreamToReadIn.close();
		}

		void computeNewStructureFactorValues(string stateFile){
			double xDisplacement;
			readInParticleState(stateFile, xDisplacement);
			const auto StartTime = chrono::steady_clock::now();
			const double xA = static_cast<double>(NumberOfAParticles)/static_cast<double>(TotalNumberOfParticles);
			const double xB = static_cast<double>(NumberOfBParticles)/static_cast<double>(TotalNumberOfParticles);

			findkValuesOnGrid(xDisplacement);

			#pragma omp parallel num_threads(NUMBER_OF_THREADS)
			{

				#pragma omp for
				for (int i = 0; i < Results.size(); i++){
					for (int j = 0; j < Results[i].kVectors.size(); j++){
						double kx = Results[i].kVectors[j].kx;
						double ky = Results[i].kVectors[j].ky;

						double cosSumA = computeCosSum(APositions, kx, ky);
						double sinSumA = computeSinSum(APositions, kx, ky);
						double cosSumB = computeCosSum(BPositions, kx, ky);
						double sinSumB = computeSinSum(BPositions, kx, ky);

						double SAA = cosSumA * cosSumA + sinSumA * sinSumA;
						double SBB = cosSumB * cosSumB + sinSumB * sinSumB;
						double SAB = cosSumA * cosSumB + sinSumA * sinSumB;

						#pragma omp critical(UPDATE_RESULTS)
						{
							Results[i].SAA += SAA;
							Results[i].SBB += SBB;
							Results[i].SAB += SAB;
							Results[i].Scc += xB * xB * SAA + xA * xA * SBB - 2.0 * xA * xB * SAB;
							Results[i].NumberOfDataPoints++;
						}
					}
					Results[i].kVectors.clear();
				}
			}
			cerr << "Computation time for structure factors: " << chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() << " s" << endl;
		}

		void computeStructureFactors(string inputFileName, int NumberOfStates){
			for (int StateIndex = 0; StateIndex < NumberOfStates; StateIndex++){
				computeNewStructureFactorValues(inputFileName);
			}
		}

		void writeResultsToFile(string fileName) {
			for (int i = 0; i < Results.size(); i++){
				if (Results[i].kVectorsFound){
					Results[i].SAA /= (static_cast<double>(Results[i].NumberOfDataPoints) * static_cast<double>(TotalNumberOfParticles));
					Results[i].SBB /= (static_cast<double>(Results[i].NumberOfDataPoints) * static_cast<double>(TotalNumberOfParticles));
					Results[i].SAB /= (static_cast<double>(Results[i].NumberOfDataPoints) * static_cast<double>(TotalNumberOfParticles));
					Results[i].Scc /= (static_cast<double>(Results[i].NumberOfDataPoints) * static_cast<double>(TotalNumberOfParticles));
				}
			}
			ofstream FileStreamToWriteTo(fileName);
			FileStreamToWriteTo << "angle\tk\tAAStructureFactor\tBBStructureFactor\tABStructureFactor\tConcentrationStructureFactor\n";
			for (int i = 0; i < Results.size(); i++){
				if (Results[i].kVectorsFound){
					FileStreamToWriteTo << setprecision(6) << Results[i].angle << '\t' << Results[i].kMag << '\t' << Results[i].SAA << '\t' << Results[i].SBB << '\t' << Results[i].SAB << '\t' << Results[i].Scc << '\n';
				}
			}
		}
};

#endif /* STRUCTURE_FACTOR_INCLUDED */
