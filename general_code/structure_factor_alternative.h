#ifndef STRUCTURE_FACTOR_INCLUDED
#define STRUCTURE_FACTOR_INCLUDED

#include <vector>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <iostream>
#include <cmath>

#include "realRNG.h"
#include "utility_functions.h"

using namespace std;

static const int DIMENSION = 2;

static const int NUMBER_OF_THREADS = 2;

class StructureFactorComputator{
	private:
		double kMin;
		const double kMax;

		realRNG RNG;

		double BoxLength;

		vector<double> APositions;
		vector<double> BPositions;
		int TotalNumberOfParticles;
		int NumberOfAParticles;
		int NumberOfBParticles;


		struct RowEntry {
			double kMagnitude;
			double SAA;
			double SBB;
			double SAB;
			double Scc;
			int NumberOfDataPoints;

			RowEntry():
				kMagnitude(0.0),
				SAA(0.0),
				SBB(0.0),
				SAB(0.0),
				Scc(0.0),
				NumberOfDataPoints(0){
			}
		};

		vector<RowEntry> Results;

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

	public:

		StructureFactorComputator(double kMax, double kWidth, double BoxLength, int TotalNumberOfParticles):
			kMax(kMax),
			BoxLength(BoxLength),
			TotalNumberOfParticles(TotalNumberOfParticles){
				kMin = 2.0*M_PI/BoxLength;
				Results = vector<RowEntry>(int((kMax - kMin)/(2.0*kWidth))+1, RowEntry());
				double CurrentkMag = kMin+kWidth;
				for (int i = 0; i < Results.size(); i++){
					Results[i].kMagnitude = CurrentkMag;
					CurrentkMag += 2.0*kWidth;
				}
		}

		void readInParticleState(string FileNameToReadIn, int StateNumber) {
			ifstream FileStreamToReadIn;
			FileStreamToReadIn.open(FileNameToReadIn);

			string CurrentString;

			APositions.clear();
			BPositions.clear();

			skipLines(FileStreamToReadIn, 1+StateNumber*(TotalNumberOfParticles+2));

			getline(FileStreamToReadIn, CurrentString, ':');
			getline(FileStreamToReadIn, CurrentString, '|');
			NumberOfAParticles =  stoi(CurrentString);
			getline(FileStreamToReadIn, CurrentString, ':');
			getline(FileStreamToReadIn, CurrentString, '|');
			NumberOfBParticles = stoi(CurrentString);
			getline(FileStreamToReadIn, CurrentString, '\n');

			for (int ParticleIndex = 0; ParticleIndex < TotalNumberOfParticles; ParticleIndex++){
				double CurrentPosition [DIMENSION];
				for (int j = 0; j < DIMENSION; j++){
					getline(FileStreamToReadIn, CurrentString, '\t');
					CurrentPosition[j] = stod(CurrentString)*BoxLength;
				}
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

		void computeNewStructureFactorValues(){
			const auto StartTime = chrono::steady_clock::now();
			double xA = static_cast<double>(NumberOfAParticles)/static_cast<double>(TotalNumberOfParticles);
			double xB = static_cast<double>(NumberOfBParticles)/static_cast<double>(TotalNumberOfParticles);

			#pragma omp parallel num_threads(NUMBER_OF_THREADS)
			{

				static const double kWidth = 0.01;
				static const int NumberOfAveragesPerk = 1000;

				#pragma omp for
				for (int i = 0; i < Results.size(); i++){
					double CurrentkMag = Results[i].kMagnitude;

					for (int j = 0; j < NumberOfAveragesPerk; j++){

						double RandomAngle = RNG.drawRandomNumber(0.0, M_PI);
						double RandomkMagnitude = RNG.drawRandomNumber(-kWidth,kWidth) + CurrentkMag;
						double kx = cos(RandomAngle)*RandomkMagnitude;
						double ky = sin(RandomAngle)*RandomkMagnitude;

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
				}
			}
			cerr << "Computation time for structure factors: " << chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() << " s" << endl;
		}

		void computeStructureFactors(string InputFileName, int NumberOfStates){
			for (int StateIndex = 0; StateIndex < NumberOfStates; StateIndex++){
				readInParticleState(InputFileName, StateIndex);
				cerr << "State " << StateIndex << " structure factor computation starting. ";
				computeNewStructureFactorValues();
			}
		}

		void writeResultsToFile(string FileName) {
			for (int kIndex = 0; kIndex < Results.size(); kIndex++){
				Results[kIndex].SAA /= (static_cast<double>(Results[kIndex].NumberOfDataPoints) * static_cast<double>(TotalNumberOfParticles));
				Results[kIndex].SBB /= (static_cast<double>(Results[kIndex].NumberOfDataPoints) * static_cast<double>(TotalNumberOfParticles));
				Results[kIndex].SAB /= (static_cast<double>(Results[kIndex].NumberOfDataPoints) * static_cast<double>(TotalNumberOfParticles));
				Results[kIndex].Scc /= (static_cast<double>(Results[kIndex].NumberOfDataPoints) * static_cast<double>(TotalNumberOfParticles));
			}
			ofstream FileStreamToWriteTo;
			FileStreamToWriteTo.open(FileName);
			FileStreamToWriteTo << "k\tAAStructureFactor\tBBStructureFactor\tABStructureFactor\tConcentrationStructureFactor\n";
			for (int i = 0; i < Results.size(); i++){
				FileStreamToWriteTo << setprecision(6) << Results[i].kMagnitude << '\t' << Results[i].SAA << '\t' << Results[i].SBB << '\t' << Results[i].SAB << '\t' << Results[i].Scc << '\n';
			}
			FileStreamToWriteTo.close();
		}
};

#endif /* STRUCTURE_FACTOR_INCLUDED */
