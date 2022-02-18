#ifndef STRUCTURE_FACTOR_INCLUDED
#define STRUCTURE_FACTOR_INCLUDED

#include <vector>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <iostream>
#include <cmath>
#include <algorithm>

#include "../realRNG.h"
#include "../utility_functions.h"

using namespace std;

static const int DIMENSION = 2;

static const int NUMBER_OF_THREADS = 2;

class StructureFactorComputator{
	private:

		realRNG RNG;

		double BoxLength;
		double angleSpacing;
		double kSpacing;

		vector<double> APositions;
		vector<double> BPositions;
		int TotalNumberOfParticles;
		int NumberOfAParticles;
		int NumberOfBParticles;

		struct Entry {
			double kMag;
			double angle;
			double SAA;
			double SBB;
			double SAB;
			double Scc;
			int NumberOfDataPoints;

			Entry(double kMag, double angle):
				kMag(kMag),
				angle(angle),
				SAA(0.0),
				SBB(0.0),
				SAB(0.0),
				Scc(0.0),
				NumberOfDataPoints(0){
			}
		};

		struct kVector {
			double kx;
			double ky;
			double kMagnitude;
			double angle;
		};

		vector<Entry> Results;

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

		vector<kVector> findkValuesOnGrid(double kAll, double kMax, double xDisplacement){

			const double kMin = 4.0*M_PI/BoxLength; // minimum k Value given by half box length due to periodic boundary conditions
			const double k1 [DIMENSION]{2.0*M_PI/BoxLength, -2.0*M_PI*xDisplacement/(BoxLength*BoxLength)}; //lattice basis in reciprocal lattice
			const double k2 [DIMENSION]{0.0, 2.0*M_PI/BoxLength};

			const double k1Mag = sqrt(k1[0]*k1[0]+k1[1]*k1[1]);
			const double k2Mag = sqrt(k2[0]*k2[0]+k2[1]*k2[1]);

			vector<kVector> kValuesOnGrid;

			int Maxn = ceil(BoxLength*kAll*0.5/(M_PI));
			for (int n1 = -Maxn; n1 <= Maxn; n1++){
				for (int n2 = 0; n2 <= Maxn; n2++){
					double kx = static_cast<double>(n1)*k1[0]+static_cast<double>(n2)*k2[0];
					double ky = static_cast<double>(n1)*k1[1]+static_cast<double>(n2)*k2[1];
					double magnitude = sqrt(kx*kx+ky*ky);
					if (magnitude > kAll){
					 break;
					}
					else if (magnitude > kMin){
						kVector newkVector;
						newkVector.kx = kx;
						newkVector.ky = ky;
						newkVector.kMagnitude = magnitude;
						if (kx == 0.0){
							newkVector.angle = M_PI*0.5;
						}
						else {
							newkVector.angle = atan2(ky,kx);
						}
						kValuesOnGrid.push_back(newkVector);
					}
				}
			}

			for (int i = 0; i < Results.size(); i++){
				double kMag = Results[i].kMag;
				
				double angle = Results[i].angle;
				
			}
		}

	public:

		StructureFactorComputator(double BoxLength, int TotalNumberOfParticles, double kMax, int kSpacing, int numberOfAngleBins):
			BoxLength(BoxLength),
			TotalNumberOfParticles(TotalNumberOfParticles),
			kSpacing(kSpacing),
			angleSpacing(M_PI/static_cast<double>(numberOfAngleBins)){
				double currentAngle = -M_PI*0.5;
				for (int angleIndex = 0; angleIndex < numberOfAngleBins; angleIndex++, currentAngle += angleSpacing){
					for (double currentkMag = 4.0*_M_PI/BoxLength; currentkMag < kMax; currentkMag += kSpacing){ // start at kMin = 4*M_PI/(L)
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

			vector<kVector> chosenkValues = findkValuesOnGrid(kAll, kMax, xDisplacement);

			#pragma omp parallel num_threads(NUMBER_OF_THREADS)
			{

				#pragma omp for
				for (int i = 0; i < kCombinationMapping.size(); i++){
					for (int j = 0; j < kCombinationMapping[i].Combinations.size(); j++){
						int nx = kCombinationMapping[i].Combinations[j].nx;
						int ny = kCombinationMapping[i].Combinations[j].ny;
						for (int k = 0; k < (abs(nx) == abs(ny) ? 1 : 2); k++){
							for (int Sign = 1; Sign >= ((nx == 0 || ny == 0) ? 1 : -1); Sign -= 2){
								double kx = Sign * nx * 2.0*M_PI/BoxLength;
								double ky = ny * 2.0*M_PI/BoxLength;
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
							swap(nx,ny);
						}
					}
				}
			}
			cerr << "Computation time for structure factors: " << chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() << " s" << endl;
		}

		void computeStructureFactors(string InputFileName, int NumberOfStates){
			for (int StateIndex = 0; StateIndex < NumberOfStates; StateIndex++){
				readInParticleState(InputFileName, StateIndex);
				computeNewStructureFactorValues();
			}
		}

		void writeResultsToFile(string FileName) {
			for (int kIndex = 0; kIndex < kCombinationMapping.size(); kIndex++){
				Results[kIndex].SAA /= (static_cast<double>(Results[kIndex].NumberOfDataPoints) * static_cast<double>(TotalNumberOfParticles));
				Results[kIndex].SBB /= (static_cast<double>(Results[kIndex].NumberOfDataPoints) * static_cast<double>(TotalNumberOfParticles));
				Results[kIndex].SAB /= (static_cast<double>(Results[kIndex].NumberOfDataPoints) * static_cast<double>(TotalNumberOfParticles));
				Results[kIndex].Scc /= (static_cast<double>(Results[kIndex].NumberOfDataPoints) * static_cast<double>(TotalNumberOfParticles));
			}
			ofstream FileStreamToWriteTo;
			FileStreamToWriteTo.open(FileName);
			FileStreamToWriteTo << "k\tAAStructureFactor\tBBStructureFactor\tABStructureFactor\tConcentrationStructureFactor\n";
			for (int i = 0; i < kCombinationMapping.size(); i++){
				FileStreamToWriteTo << setprecision(6) << kCombinationMapping[i].kMagnitude << '\t' << Results[i].SAA << '\t' << Results[i].SBB << '\t' << Results[i].SAB << '\t' << Results[i].Scc << '\n';
			}
			FileStreamToWriteTo.close();
		}
};

#endif /* STRUCTURE_FACTOR_INCLUDED */
