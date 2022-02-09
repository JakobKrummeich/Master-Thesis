#ifndef STRUCTURE_FACTOR_INCLUDED
#define STRUCTURE_FACTOR_INCLUDED

#include <vector>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <iostream>
#include <cmath>
#include <algorithm>

#include "realRNG.h"
#include "utility_functions.h"

using namespace std;

static const int DIMENSION = 2;

static const int NUMBER_OF_THREADS = 2;

class StructureFactorComputator{
	private:

		realRNG RNG;

		double BoxLength;

		vector<double> APositions;
		vector<double> BPositions;
		int TotalNumberOfParticles;
		int NumberOfAParticles;
		int NumberOfBParticles;

		struct Combination{
			int nx;
			int ny;

			Combination(){
			}

			Combination(int nx, int ny):
				nx(nx),
				ny(ny){
			}

			bool operator==(const Combination& RHS) const {
				return ((nx == RHS.nx && ny == RHS.ny) || (nx == RHS.ny && ny == RHS.nx));
			}
		};

		struct SquaredCombinationMappingEntry{
			long MagnitudeSquared;
			Combination AssociatedCombination;

			static bool isSmaller(const SquaredCombinationMappingEntry& LHS, const SquaredCombinationMappingEntry& RHS){
				return (LHS.MagnitudeSquared < RHS.MagnitudeSquared);
			}
		};

		struct kCombinationMappingEntry{
			double kMagnitude;
			vector<Combination> Combinations;
		};

		struct RowEntry {
			double SAA;
			double SBB;
			double SAB;
			double Scc;
			int NumberOfDataPoints;

			RowEntry():
				SAA(0.0),
				SBB(0.0),
				SAB(0.0),
				Scc(0.0),
				NumberOfDataPoints(0){
			}
		};

		vector<kCombinationMappingEntry> kCombinationMapping;
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

		void findkValuesOnGrid(double kAll, double kMax){

			vector<SquaredCombinationMappingEntry> SquareCombinationMapping;
			int Maxn = ceil(BoxLength*kAll*0.5/(M_PI));
			for (int nx = 1; nx <= Maxn; nx++){
				for (int ny = 0; ny <= nx; ny++){
					long Magnitude = nx*nx+ny*ny;
					if (Magnitude < Maxn*Maxn){
						SquareCombinationMapping.push_back(SquaredCombinationMappingEntry{Magnitude,Combination{nx,ny}});
					}
					else {
						break;
					}
				}
			}
			sort(SquareCombinationMapping.begin(), SquareCombinationMapping.end(), SquaredCombinationMappingEntry::isSmaller);

			for (int i = 0; i < SquareCombinationMapping.size(); i++){
				kCombinationMappingEntry NewMapping{};
				NewMapping.kMagnitude = 2.0*M_PI/(BoxLength)*sqrt(static_cast<double>(SquareCombinationMapping[i].MagnitudeSquared));
				NewMapping.Combinations.push_back(SquareCombinationMapping[i].AssociatedCombination);
				while (i+1 < SquareCombinationMapping.size() && SquareCombinationMapping[i].MagnitudeSquared == SquareCombinationMapping[i+1].MagnitudeSquared){
					NewMapping.Combinations.push_back(SquareCombinationMapping[i+1].AssociatedCombination);
					i++;
				}
				kCombinationMapping.push_back(NewMapping);
			}

			double kWidth = 0.01;
			double kAllChosen = 2.0*M_PI/(BoxLength)*static_cast<double>(Maxn);
			double CurrentkMag = kAllChosen;
			int NumberOfAttemptedAveragesPerk = 1000;
			int MaxNumberOfCombinationsPerInterval = 20;
			vector<Combination> CombinationsFound;

			while (CurrentkMag < kMax){
				CombinationsFound.clear();
				kCombinationMappingEntry NewMapping{};
				bool MagnitudeFound = false;
				for (int i = 0; i < NumberOfAttemptedAveragesPerk && CombinationsFound.size() < MaxNumberOfCombinationsPerInterval; i++){
					double RandomAngle = RNG.drawRandomNumber(0.0, 0.5 * M_PI);
					double RandomkMagnitude = RNG.drawRandomNumber(-kWidth,kWidth) + CurrentkMag;
					int nx = round(BoxLength*RandomkMagnitude*cos(RandomAngle)/(2.0*M_PI));
					int ny = round(BoxLength*RandomkMagnitude*sin(RandomAngle)/(2.0*M_PI));
					double GridkMagnitude = 2.0*M_PI/BoxLength*sqrt(static_cast<double>(nx*nx+ny*ny));
					if (GridkMagnitude <= (CurrentkMag + kWidth) && GridkMagnitude >= (CurrentkMag - kWidth) && GridkMagnitude >= kAllChosen){
						bool sameCombinationAlreadyFoundBefore = false;
						Combination NewCombination(nx,ny);
						for (int j = 0; j < CombinationsFound.size() && !sameCombinationAlreadyFoundBefore; j++){
							if (CombinationsFound[j] == NewCombination){
								sameCombinationAlreadyFoundBefore = true;
							}
						}
						if (!sameCombinationAlreadyFoundBefore){
							CombinationsFound.push_back(NewCombination);
							if (!MagnitudeFound){
								NewMapping.kMagnitude = GridkMagnitude;
								NewMapping.Combinations.push_back(NewCombination);
								MagnitudeFound = true;
							}
							else {
								NewMapping.kMagnitude = CurrentkMag;
								NewMapping.Combinations.push_back(NewCombination);
							}
						}
					}
				}
				if (MagnitudeFound){
					kCombinationMapping.push_back(NewMapping);
				}
				CurrentkMag += 2.0*kWidth;
			}

			for (int i = 0; i < kCombinationMapping.size(); i++){
				cerr << kCombinationMapping[i].kMagnitude << ": ";
				for (int j = 0; j < kCombinationMapping[i].Combinations.size(); j++){
					cerr << "(" << kCombinationMapping[i].Combinations[j].nx << "," << kCombinationMapping[i].Combinations[j].ny << "),";
				}
				cerr << endl;
			}
		}

	public:

		StructureFactorComputator(double kAll, double kMax, double BoxLength, int TotalNumberOfParticles):
			BoxLength(BoxLength),
			TotalNumberOfParticles(TotalNumberOfParticles){
				if (kAll <= 2.0*M_PI/BoxLength){
					cerr << "We need to have kAll >= 2.0*PI/L!" << endl;
					exit;
				}
				else if (kMax < kAll){
					cerr << "We need to have kAll <= kMax!" << endl;
					exit;
				}
				findkValuesOnGrid(kAll,kMax);
				Results = vector<RowEntry>(kCombinationMapping.size(), RowEntry());
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
