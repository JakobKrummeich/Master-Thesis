#include <iostream>
#include <cmath>
#include <chrono>
#include <fstream>
#include <vector>
#include <iomanip>
#include <stdlib.h>
#include "../../General_Code/realRNG.h"
#include "../../General_Code/fvec.h"

using namespace std;

const int DIMENSION = 2;
const int TOTAL_NUMBER_OF_PARTICLES = 1000;

const double AA_INTERACTION_STRENGTH = 1.0;
double AB_INTERACTION_STRENGTH;
const double CUTOFF = 2.5;
const double CUTOFF_SQUARED = CUTOFF * CUTOFF;
const double INVERSE_CUTOFF = 1.0/CUTOFF;
const double POTENTIAL_CONSTANT_1 = pow(INVERSE_CUTOFF, 6.0) - pow(INVERSE_CUTOFF, 12.0);
const double POTENTIAL_CONSTANT_2 = 6.0 * pow(INVERSE_CUTOFF, 7.0) - 12.0 * pow(INVERSE_CUTOFF, 13.0);

double DENSITY;
double BOX_LENGTH;
double BOX_LENGTH_SQUARED;
double INVERSE_BOX_LENGTH;

const double MAXIMUM_DISPLACEMENT = 0.1;
const double MAX_VERLET_DIST = 1.3*CUTOFF;
const double MAX_VERLET_DIST_SQUARED = MAX_VERLET_DIST * MAX_VERLET_DIST;
const double SKINDISTANCE = MAX_VERLET_DIST - CUTOFF;
int NUMBER_OF_SUBDIVISIONS;

const double DISPLACEMENT_PROBABILITY = 0.9;
const int UPDATE_TIME_INTERVAL = 60;
const int POT_ENERGY_UPDATE_INTERVAL = 200;
const int NUMBER_OF_STRUCTURE_FACTOR_AVERAGES = 1000;

void initializeBox(){
	BOX_LENGTH = sqrt(static_cast<double>(TOTAL_NUMBER_OF_PARTICLES) / DENSITY);
	BOX_LENGTH_SQUARED = BOX_LENGTH * BOX_LENGTH;
	INVERSE_BOX_LENGTH = 1.0/BOX_LENGTH;
	NUMBER_OF_SUBDIVISIONS = static_cast<int>(BOX_LENGTH/MAX_VERLET_DIST) > 3 ? static_cast<int>(BOX_LENGTH/MAX_VERLET_DIST) : 1;
}

enum class ParticleType{
	A = 0,
	B = 1
};

ostream& operator<<(ostream& OStream, ParticleType T){
	if (T == ParticleType::A){
		OStream << "A";
	}
	else {
		OStream << "B";
	}
	return OStream;
}

struct Particles {
	private:
	double Positions [DIMENSION*TOTAL_NUMBER_OF_PARTICLES];

	ParticleType ParticleTypes [TOTAL_NUMBER_OF_PARTICLES];

	fvec<int, TOTAL_NUMBER_OF_PARTICLES> TypeAParticleIndices;
	fvec<int, TOTAL_NUMBER_OF_PARTICLES> TypeBParticleIndices;

	vector<int> CellListHead;
	int CellListIndices [TOTAL_NUMBER_OF_PARTICLES];
	int VerletListHead [2*TOTAL_NUMBER_OF_PARTICLES];
	vector<int> VerletIndicesOfNeighbors;

	double ChangeInCoordinates [DIMENSION*TOTAL_NUMBER_OF_PARTICLES];
	double MostTraveledDistances [2];
	int MostTraveledParticleIndex;

	friend class StructureFactorComputator;
	friend class SimulationManager;

	public:

	int NumberOfVerletListBuilds;

	Particles(){
		for (int i = 0; i < TOTAL_NUMBER_OF_PARTICLES; i++){
			for (int j = 0; j < DIMENSION; j++){
				ChangeInCoordinates[DIMENSION*i+j] = 0.0;
			}
		}
		for (int i = 0; i < 2; i++){
			MostTraveledDistances[i] = 0.0;
		}
		NumberOfVerletListBuilds = 0;
	}

	void initialize(int InitialNumberOfAParticles){
		int InitialNumberOfBParticles = TOTAL_NUMBER_OF_PARTICLES - InitialNumberOfAParticles;
		double FractionOfAParticles = static_cast<double>(InitialNumberOfAParticles)/static_cast<double>(TOTAL_NUMBER_OF_PARTICLES);
		realRNG RNG;

		int NumberOfParticlesInARow(ceil(pow(static_cast<double>(TOTAL_NUMBER_OF_PARTICLES),1.0/static_cast<double>(DIMENSION))));
		double Distance(1.0/static_cast<double>(NumberOfParticlesInARow));
		cerr << "NumberOfParticlesInRow: " << NumberOfParticlesInARow << endl;
		cerr << "Distance: " << Distance*BOX_LENGTH << endl;
		double CurrentPosition [DIMENSION];
		for (int i = 0; i < DIMENSION; i++){
			CurrentPosition[i] = Distance*0.5;
		}
		int NumberOfAParticlesInitialized = 0;
		for (int ParticlesInitialized = 0; ParticlesInitialized < TOTAL_NUMBER_OF_PARTICLES; ParticlesInitialized++){
			for (int i = 0; i < DIMENSION; i++){
				Positions[ParticlesInitialized*DIMENSION + i] = CurrentPosition[i];
			}
			CurrentPosition[0] += Distance;
			for (int i = 1; i < DIMENSION; i++){
				if (CurrentPosition[i-1] >= 1.0){
					CurrentPosition[i-1] = Distance*0.5;
					CurrentPosition[i] += Distance;
				}
			}
			if ((RNG.drawRandomNumber() <= FractionOfAParticles && NumberOfAParticlesInitialized < InitialNumberOfAParticles) || (TOTAL_NUMBER_OF_PARTICLES - ParticlesInitialized <= InitialNumberOfAParticles - NumberOfAParticlesInitialized)){
				ParticleTypes[ParticlesInitialized] = ParticleType::A;
				TypeAParticleIndices.push_back(ParticlesInitialized);
				NumberOfAParticlesInitialized++;
			}
			else {
				ParticleTypes[ParticlesInitialized] = ParticleType::B;
				TypeBParticleIndices.push_back(ParticlesInitialized);
			}
		}
	}

	void readInParticleState(string FileNameToReadIn) {
		ifstream FileStreamToReadIn;
		FileStreamToReadIn.open(FileNameToReadIn);

		string CurrentString;

		getline(FileStreamToReadIn,CurrentString, '\n');

		for (int ParticleIndex = 0; ParticleIndex < TOTAL_NUMBER_OF_PARTICLES; ParticleIndex++){
			getline(FileStreamToReadIn, CurrentString, '\t');
			for (int j = 0; j < DIMENSION; j++){
				getline(FileStreamToReadIn, CurrentString, '\t');
				Positions[DIMENSION*ParticleIndex+j] = stod(CurrentString);
			}
			getline(FileStreamToReadIn, CurrentString, '\n');
			if (CurrentString == "A"){
				ParticleTypes[ParticleIndex] = ParticleType::A;
				TypeAParticleIndices.push_back(ParticleIndex);
			}
			else {
				ParticleTypes[ParticleIndex] = ParticleType::B;
				TypeBParticleIndices.push_back(ParticleIndex);
			}
		}
		FileStreamToReadIn.close();
	}

	double getPosition(int ParticleIndex, int Coordinate) const {
		return Positions[DIMENSION*ParticleIndex+Coordinate];
	}

	ParticleType getParticleType(int ParticleIndex) const {
		return ParticleTypes[ParticleIndex];
	}

	void switchParticleType(int ParticleID, int IndexInTypeArray) {
		if (ParticleTypes[ParticleID] == ParticleType::A){
			TypeAParticleIndices.erase(IndexInTypeArray);
			TypeBParticleIndices.push_back(ParticleID);
		}
		else {
			TypeBParticleIndices.erase(IndexInTypeArray);
			TypeAParticleIndices.push_back(ParticleID);
		}
		ParticleTypes[ParticleID] = ParticleTypes[ParticleID] == ParticleType::A ? ParticleType::B : ParticleType::A;
	}

	int getNumberOfAParticles() const {
		return TypeAParticleIndices.size();
	}

	int getNumberOfBParticles() const {
		return TypeBParticleIndices.size();
	}

	static double computePairwiseParticlePotentialEnergy(double DimensionlessDistanceSquared) {
		double InverseDistanceSquared = 1.0/(DimensionlessDistanceSquared*BOX_LENGTH_SQUARED);
		double InverseDistanceToThePowerOfSix = InverseDistanceSquared * InverseDistanceSquared * InverseDistanceSquared;
		return 4.0*(InverseDistanceToThePowerOfSix * InverseDistanceToThePowerOfSix - InverseDistanceToThePowerOfSix + POTENTIAL_CONSTANT_1 - (sqrt(DimensionlessDistanceSquared) * BOX_LENGTH - CUTOFF) * POTENTIAL_CONSTANT_2);
	}

	static double computePairwiseParticlePotentialEnergy(const double* Position0, const double* Position1) {
		double DistanceSquared = 0.0;
		for (int i = 0; i < DIMENSION; i++){
			double CoordinateDifference = *(Position0 + i) - *(Position1 + i);
			if (CoordinateDifference > 0.5){
				CoordinateDifference -= 1.0;
			}
			else if (CoordinateDifference <= -0.5){
				CoordinateDifference += 1.0;
			}
			DistanceSquared += CoordinateDifference*CoordinateDifference;
		}
		if (DistanceSquared*BOX_LENGTH_SQUARED >= CUTOFF_SQUARED){
			return 0.0;
		}
		return computePairwiseParticlePotentialEnergy(DistanceSquared);
	}

	double computeChangeInPotentialEnergyByMoving(int ParticleIndex, const double* Delta) const {
		double PotEnergyChange = 0.0;
		double UpdatedCoordinates [DIMENSION];
		for (int i = 0; i < DIMENSION; i++){
			UpdatedCoordinates[i] = Positions[DIMENSION*ParticleIndex + i] + *(Delta + i);
			if (UpdatedCoordinates[i] < 0.0){
				UpdatedCoordinates[i] += 1.0;
			}
			else if (UpdatedCoordinates[i] > 1.0){
				UpdatedCoordinates[i] -= 1.0;
			}
		}
		for (int i = 0; i < VerletListHead[2*ParticleIndex+1]; i++){
			int OtherParticleIndex = VerletIndicesOfNeighbors[VerletListHead[2*ParticleIndex]+i];
			double InteractionStrength = AB_INTERACTION_STRENGTH;
			if (ParticleTypes[ParticleIndex] == ParticleTypes[OtherParticleIndex]) {
				InteractionStrength = AA_INTERACTION_STRENGTH;
			}
			PotEnergyChange += InteractionStrength * (computePairwiseParticlePotentialEnergy(&Positions[DIMENSION*OtherParticleIndex], UpdatedCoordinates) 
																								- computePairwiseParticlePotentialEnergy(&Positions[DIMENSION*OtherParticleIndex], &Positions[DIMENSION*ParticleIndex]));
		}
		return PotEnergyChange;
	}

	double computeChangeInPotentialEnergyBySwitching(int ParticleIndex) const {
		double PotEnergyChange = 0.0;
		for (int i = 0; i < VerletListHead[2*ParticleIndex+1]; i++){
			int OtherParticleIndex = VerletIndicesOfNeighbors[VerletListHead[2*ParticleIndex]+i];
			double PrefactorDifference = AA_INTERACTION_STRENGTH - AB_INTERACTION_STRENGTH;
			if (ParticleTypes[ParticleIndex] == ParticleTypes[OtherParticleIndex]){
				PrefactorDifference *= -1.0;
			}
			PotEnergyChange += PrefactorDifference * computePairwiseParticlePotentialEnergy(&Positions[DIMENSION*OtherParticleIndex], &Positions[DIMENSION*ParticleIndex]);
		}
		return PotEnergyChange;
	}

	double computePotentialEnergy() const {
		double PotEnergy = 0.0;
		for (int ParticleIndex = 0; ParticleIndex < TOTAL_NUMBER_OF_PARTICLES; ParticleIndex++){
			for (int i = 0; i < VerletListHead[2*ParticleIndex+1]; i++){
				int OtherParticleIndex = VerletIndicesOfNeighbors[VerletListHead[2*ParticleIndex]+i];
				if (OtherParticleIndex < ParticleIndex){
					double InteractionStrength = AB_INTERACTION_STRENGTH;
					if (ParticleTypes[ParticleIndex] == ParticleTypes[OtherParticleIndex]) {
						InteractionStrength = AA_INTERACTION_STRENGTH;
					}
					PotEnergy += InteractionStrength * computePairwiseParticlePotentialEnergy(&Positions[DIMENSION*OtherParticleIndex], &Positions[DIMENSION*ParticleIndex]);
				}
			}
		}
		return PotEnergy;
	}

	void buildCellList(){
		int NumberOfSubcells = NUMBER_OF_SUBDIVISIONS;
		for (int i = 0; i < DIMENSION-1; i++){
			NumberOfSubcells *= NUMBER_OF_SUBDIVISIONS;
		}
		CellListHead.clear();
		CellListHead.resize(NumberOfSubcells,-1);
		int CurrentCellIndex;
		int IndexFactor;
		for (int ParticleIndex = 0; ParticleIndex < TOTAL_NUMBER_OF_PARTICLES; ParticleIndex++){
			CurrentCellIndex = 0;
			IndexFactor = 1;
			for (int j = 0; j < DIMENSION; j++){
				CurrentCellIndex += static_cast<int>(static_cast<double>(NUMBER_OF_SUBDIVISIONS)*Positions[DIMENSION*ParticleIndex+j])*IndexFactor;
				IndexFactor *= NUMBER_OF_SUBDIVISIONS;
			}
			CellListIndices[ParticleIndex] = CellListHead[CurrentCellIndex];
			CellListHead[CurrentCellIndex] = ParticleIndex;
		}
	}

	void printCellList() const {
		for (int i = 0; i < CellListHead.size(); i++)	{
			cerr << "Cell " << i << ": Head: " << CellListHead[i] << endl;
			int CurrentParticleIndex(CellListHead[i]);
			while (CurrentParticleIndex >= 0){
				cerr << CurrentParticleIndex << ',';
				CurrentParticleIndex = CellListIndices[CurrentParticleIndex];
			}
			cerr << endl;
		}
	}

	void buildVerletList() {
		buildCellList();
		VerletIndicesOfNeighbors.clear();
		VerletIndicesOfNeighbors.reserve(40*TOTAL_NUMBER_OF_PARTICLES);
		int CurrentIndexInVerletIndices = 0;

		int Indices [DIMENSION];
		for (int i = 0; i < DIMENSION; i++){
			Indices[i] = 0;
		}
		while (Indices[DIMENSION - 1] < NUMBER_OF_SUBDIVISIONS){
			int Cell = 0;
			int IndexFactor = 1;
			for (int i = 0; i < DIMENSION; i++){
				Cell += Indices[i]*IndexFactor;
				IndexFactor *= NUMBER_OF_SUBDIVISIONS;

			}
			int CurrentParticleIndex = CellListHead[Cell];
			while (CurrentParticleIndex >= 0){
				int NumberOfNeighbors = 0;
				VerletListHead[2*CurrentParticleIndex] = CurrentIndexInVerletIndices;
				int IndicesOffsets [DIMENSION];
				for (int i = 0; i < DIMENSION; i++){
					IndicesOffsets[i] = -1;
				}
				while (IndicesOffsets[DIMENSION - 1] < 2){
					int NeighborCell = 0;
					int IndexFactor = 1;
					for (int i = 0; i < DIMENSION; i++){
						int OtherCellIndex = Indices[i]+IndicesOffsets[i];
						if (OtherCellIndex < 0){
							OtherCellIndex += NUMBER_OF_SUBDIVISIONS;
						}
						else if (OtherCellIndex >= NUMBER_OF_SUBDIVISIONS){
							OtherCellIndex -= NUMBER_OF_SUBDIVISIONS;
						}
						NeighborCell += IndexFactor*OtherCellIndex;
						IndexFactor *= NUMBER_OF_SUBDIVISIONS;
					}
					int OtherParticleIndex = CellListHead[NeighborCell];
					while (OtherParticleIndex >= 0){
						if (OtherParticleIndex != CurrentParticleIndex){
							double DistanceSquared = 0.0;
							for (int k = 0; k < DIMENSION; k++){
								double CoordinateDifference = Positions[DIMENSION * CurrentParticleIndex + k] - Positions[DIMENSION * OtherParticleIndex + k];
								if (CoordinateDifference > 0.5){
									CoordinateDifference -= 1.0;
								}
								else if (CoordinateDifference <= -0.5){
									CoordinateDifference += 1.0;
								}
								DistanceSquared += CoordinateDifference * CoordinateDifference;
							}
							if (DistanceSquared * BOX_LENGTH_SQUARED <= MAX_VERLET_DIST_SQUARED){
								NumberOfNeighbors++;
								VerletIndicesOfNeighbors.push_back(OtherParticleIndex);
								CurrentIndexInVerletIndices++;
							}
						}
						OtherParticleIndex = CellListIndices[OtherParticleIndex];
					}
					IndicesOffsets[0]++;
					for (int i = 0; i < DIMENSION - 1; i++){
						if (IndicesOffsets[i] >= 2){
							IndicesOffsets[i] = -1;
							IndicesOffsets[i+1]++;
						}
					}
				}
				VerletListHead[2*CurrentParticleIndex+1] = NumberOfNeighbors;
				CurrentParticleIndex = CellListIndices[CurrentParticleIndex];
			}
			Indices[0]++;
			for (int i = 0; i < DIMENSION - 1; i++){
				if (Indices[i] >= NUMBER_OF_SUBDIVISIONS){
					Indices[i] = 0;
					Indices[i+1]++;
				}
			}
		}
		NumberOfVerletListBuilds++;
	}

	void printVerletList() const{
		for (int i = 0; i < TOTAL_NUMBER_OF_PARTICLES; i++){
			cerr << VerletListHead[2*i] << "," << VerletListHead[2*i+1] << '|';
		}
		cerr << endl;
		for (int i = 0; i < VerletIndicesOfNeighbors.size(); i++){
			cerr << VerletIndicesOfNeighbors[i] << ',';
		}
			cerr << endl;
	}

	void updatePosition(int ParticleIndex, const double* Deltas){
		double CurrentTraveledDistance = 0.0;
		for (int k = 0; k < DIMENSION; k++){
			Positions[DIMENSION * ParticleIndex + k] += Deltas[k];
			if (Positions[DIMENSION * ParticleIndex + k] < 0.0){
				Positions[DIMENSION * ParticleIndex + k] += 1.0;
			}
			else if (Positions[DIMENSION * ParticleIndex + k] > 1.0){
				Positions[DIMENSION * ParticleIndex + k] -= 1.0;
			}
			ChangeInCoordinates[DIMENSION * ParticleIndex + k] += Deltas[k];
			CurrentTraveledDistance += ChangeInCoordinates[DIMENSION * ParticleIndex + k] * ChangeInCoordinates[DIMENSION * ParticleIndex + k];
		}
		CurrentTraveledDistance = sqrt(CurrentTraveledDistance);
		if (CurrentTraveledDistance >= MostTraveledDistances[0]) {
			if (ParticleIndex == MostTraveledParticleIndex) {
				MostTraveledDistances[0] = CurrentTraveledDistance;
			}
			else {
				MostTraveledDistances[1] = MostTraveledDistances[0];
				MostTraveledDistances[0] = CurrentTraveledDistance;
				MostTraveledParticleIndex = ParticleIndex;
			}
		}
		else if (CurrentTraveledDistance > MostTraveledDistances[1]) {
			MostTraveledDistances[1] = CurrentTraveledDistance;
		}
		if (BOX_LENGTH*(MostTraveledDistances[0]+MostTraveledDistances[1]) > SKINDISTANCE){
			buildVerletList();
			for (int k = 0; k < DIMENSION * TOTAL_NUMBER_OF_PARTICLES; k++){
				ChangeInCoordinates[k] = 0.0;
			}
			for (int i = 0; i < 2; i++){
				MostTraveledDistances[i] = 0.0;
			}
		}
	}
};

ostream& operator<<(ostream& OStream, const Particles& State){
	OStream << "#ID\tX       Y       Type | #AParticles:  " << State.getNumberOfAParticles() << "| #BParticles: " << State.getNumberOfBParticles() << "| BOX_LENGTH: " << BOX_LENGTH << endl;
	for (int i = 0; i < TOTAL_NUMBER_OF_PARTICLES; i++){
		OStream << i << "\t";
		OStream << fixed << setprecision(5) << State.getPosition(i,0) << "\t" << State.getPosition(i,1) << "\t" << State.getParticleType(i) <<  endl;
	}
	return OStream;
}

class StructureFactorComputator{
	private:
		const double kMin;
		const double kMax;

		realRNG RNG;

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

		struct RowEntry {
			double kMagnitude;
			vector<Combination> Combinations;
			double SAA;
			double SBB;
			double SAB;
			double Scc;
			int NumberOfDataPoints;
		};

		vector<RowEntry> Results;

		double computeCosSum(const double* const Positions, const fvec<int, TOTAL_NUMBER_OF_PARTICLES>& ParticleIndices, double kx, double ky) const{
			double Sum = 0.0;
			for (int i = 0; i < ParticleIndices.size(); i++){
				Sum += cos(kx* (*(Positions+DIMENSION*ParticleIndices[i])) + ky * (*(Positions+DIMENSION*ParticleIndices[i]+1)));
			}
			return Sum;
		}

		double computeSinSum(const double* const Positions, const fvec<int, TOTAL_NUMBER_OF_PARTICLES>& ParticleIndices, double kx, double ky) const{
			double Sum = 0.0;
			for (int i = 0; i < ParticleIndices.size(); i++){
				Sum += sin(kx* (*(Positions+DIMENSION*ParticleIndices[i])) + ky * (*(Positions+DIMENSION*ParticleIndices[i]+1)));
			}
			return Sum;
		}

	public:
		StructureFactorComputator(double kMax):
			kMin(2.0*M_PI/BOX_LENGTH),
			kMax(kMax){
		}

		void findkValuesOnGrid(){
			double kWidth = 0.01;
			double CurrentkMag = kMin;
			int NumberOfAttemptedAveragesPerk = 1000;
			vector<Combination> CombinationsFound;
			vector<int> MagnitudesFound;
			vector<RowEntry> IntermediateResults;

			while (CurrentkMag < kMax){
				CombinationsFound.clear();
				MagnitudesFound.clear();
				IntermediateResults.clear();
				for (int i = 0; i < NumberOfAttemptedAveragesPerk; i++){
					double RandomAngle = RNG.drawRandomNumber(0.0, 0.5 * M_PI);
					double RandomkMagnitude = RNG.drawRandomNumber(-kWidth,kWidth) + CurrentkMag;
					int nx = round(BOX_LENGTH*RandomkMagnitude*cos(RandomAngle)/(2.0*M_PI));
					int ny = round(BOX_LENGTH*RandomkMagnitude*sin(RandomAngle)/(2.0*M_PI));
					double GridkMagnitude = 2.0*M_PI/BOX_LENGTH*sqrt(static_cast<double>(nx*nx+ny*ny));
					if (GridkMagnitude <= (CurrentkMag+kWidth) && GridkMagnitude >= (CurrentkMag-kWidth) && GridkMagnitude >= kMin){
						bool sameCombinationAlreadyFoundBefore = false;
						Combination NewCombination(nx,ny);
						for (int j = 0; j < CombinationsFound.size() && !sameCombinationAlreadyFoundBefore; j++){
							if (CombinationsFound[j] == NewCombination){
								sameCombinationAlreadyFoundBefore = true;
							}
						}
						if (!sameCombinationAlreadyFoundBefore){
							CombinationsFound.push_back(NewCombination);
							int NewGridMagnitude = nx*nx+ny*ny;
							bool sameMagnitudeFoundBefore = false;
							for (int j = 0; j < MagnitudesFound.size(); j++){
								if (NewGridMagnitude == MagnitudesFound[j]){
									sameMagnitudeFoundBefore = true;
									IntermediateResults[j].Combinations.push_back(NewCombination);
									break;
								}
							}
							if (!sameMagnitudeFoundBefore){
								RowEntry NewEntry{};
								NewEntry.kMagnitude = GridkMagnitude;
								NewEntry.Combinations.push_back(NewCombination);
								IntermediateResults.push_back(NewEntry);
								MagnitudesFound.push_back(NewGridMagnitude);
								for (int CurrentIndex = MagnitudesFound.size() - 1; CurrentIndex > 0 && MagnitudesFound[CurrentIndex-1] > MagnitudesFound[CurrentIndex]; CurrentIndex--){
									swap(MagnitudesFound[CurrentIndex], MagnitudesFound[CurrentIndex-1]);
									swap(IntermediateResults[CurrentIndex], IntermediateResults[CurrentIndex-1]);
								}
							}
						}
					}
				}
				for (int i = 0; i < IntermediateResults.size(); i++){
					Results.push_back(IntermediateResults[i]);
				}
				CurrentkMag += 2.0*kWidth;
			}
		}

		void computeNewStructureFactorValues(const double* const Positions, const fvec<int, TOTAL_NUMBER_OF_PARTICLES>& TypeAParticleIndices, const fvec<int, TOTAL_NUMBER_OF_PARTICLES>& TypeBParticleIndices){
			double xA = static_cast<double>(TypeAParticleIndices.size())/static_cast<double>(TOTAL_NUMBER_OF_PARTICLES);
			double xB = static_cast<double>(TypeBParticleIndices.size())/static_cast<double>(TOTAL_NUMBER_OF_PARTICLES);
			for (int i = 0; i < Results.size(); i++){
				for (int j = 0; j < Results[i].Combinations.size(); j++){
					int nx = Results[i].Combinations[j].nx;
					int ny = Results[i].Combinations[j].ny;
					for (int k = 0; k < (abs(nx) == abs(ny) ? 1 : 2); k++){
						for (int Sign = 1; Sign >= ((nx == 0 || ny == 0)? 1: -1); Sign -= 2){
							double kx = Sign * nx * 2.0*M_PI/BOX_LENGTH;
							double ky = ny * 2.0*M_PI/BOX_LENGTH;
							double cosSumA = computeCosSum(Positions, TypeAParticleIndices, kx, ky);
							double sinSumA = computeSinSum(Positions, TypeAParticleIndices, kx, ky);
							double cosSumB = computeCosSum(Positions, TypeBParticleIndices, kx, ky);
							double sinSumB = computeSinSum(Positions, TypeBParticleIndices, kx, ky);

							double SAA = cosSumA*cosSumA + sinSumA * sinSumA;
							double SBB = cosSumB*cosSumB + sinSumB * sinSumB;
							double SAB = cosSumA*cosSumB + sinSumA * sinSumB;

							Results[i].SAA += SAA;
							Results[i].SBB += SBB;
							Results[i].SAB += SAB;
							Results[i].Scc += xB*xB*SAA + xA*xA*SBB - 2.0 * xA * xB * SAB;
							Results[i].NumberOfDataPoints++;
						}
						swap(nx,ny);
					}
				}
			}
		}

		void reset(){
			for (int i = 0; i < Results.size(); i++){
				Results[i].SAA = 0.0;
				Results[i].SBB = 0.0;
				Results[i].SAB = 0.0;
				Results[i].Scc = 0.0;
				Results[i].NumberOfDataPoints = 0;
			}
		}

		void writeResultsToFile(string FileName) {
			for (int i = 0; i < Results.size(); i++){
				Results[i].SAA /= (static_cast<double>(Results[i].NumberOfDataPoints) * static_cast<double>(TOTAL_NUMBER_OF_PARTICLES));
				Results[i].SBB /= (static_cast<double>(Results[i].NumberOfDataPoints) * static_cast<double>(TOTAL_NUMBER_OF_PARTICLES));
				Results[i].SAB /= (static_cast<double>(Results[i].NumberOfDataPoints) * static_cast<double>(TOTAL_NUMBER_OF_PARTICLES));
				Results[i].Scc /= (static_cast<double>(Results[i].NumberOfDataPoints) * static_cast<double>(TOTAL_NUMBER_OF_PARTICLES));
			}
			ofstream FileStreamToWriteTo;
			FileStreamToWriteTo.open(FileName);
			FileStreamToWriteTo << "k\tAAStructureFactor\tBBStructureFactor\tABStructureFactor\tConcentrationStructureFactor\n";
			for (int i = 0; i < Results.size(); i++){
				FileStreamToWriteTo << Results[i].kMagnitude << '\t' << Results[i].SAA << '\t' << Results[i].SBB << '\t' << Results[i].SAB << '\t' << Results[i].Scc << '\n';
			}
			FileStreamToWriteTo.close();
		}
};

struct SimulationManager {
	Particles P;
	double Temperature;
	double Beta;
	realRNG RNG;
	double ChemicalPotentialDiff;

	int MinNumberOfA;
	int MaxNumberOfA;
	int NumberOfMCSweeps;

	vector<int> NumberOfABuffer;
	vector<double> PotEnergyBuffer;

	int NumberOfTriedDisplacements;
	int NumberOfAcceptedDisplacements;

	int NumberOfTriedTypeChanges;
	int NumberOfAcceptedTypeChanges;

	string FileNameString;

	StructureFactorComputator SFComputator;

	SimulationManager(double _Temperature, double _ChemicalPotentialDiff, int _MinNumberOfA, int _MaxNumberOfA, int _NumberOfMCSweeps, double kMax):
		SFComputator(kMax),
		Temperature(_Temperature),
		Beta(1.0/Temperature),
		ChemicalPotentialDiff(_ChemicalPotentialDiff),
		MinNumberOfA(_MinNumberOfA),
		MaxNumberOfA(_MaxNumberOfA),
		NumberOfMCSweeps(_NumberOfMCSweeps),
		NumberOfTriedDisplacements(0),
		NumberOfAcceptedDisplacements(0),
		NumberOfTriedTypeChanges(0),
		NumberOfAcceptedTypeChanges(0),
		FileNameString("N="+to_string(TOTAL_NUMBER_OF_PARTICLES)+"_T="+to_string(Temperature)+"_AvgDens="+to_string(DENSITY)+"_MCRuns="+to_string(NumberOfMCSweeps)+"_epsAB="+to_string(AB_INTERACTION_STRENGTH)+".dat"){
	}

	void initialize() {
		P.initialize(MinNumberOfA);
		P.buildVerletList();
		SFComputator.findkValuesOnGrid();
	}

	void initialize(string FileNameInitialConfiguration) {
		P.readInParticleState(FileNameInitialConfiguration);
		P.buildVerletList();
	}

	void changeTemperature(double NewTemperature){
		Temperature = NewTemperature;
		Beta = 1.0/Temperature;
		FileNameString = "N="+to_string(TOTAL_NUMBER_OF_PARTICLES)+"_T="+to_string(Temperature)+"_AvgDens="+to_string(DENSITY)+"_MCRuns="+to_string(NumberOfMCSweeps)+"_epsAB="+to_string(AB_INTERACTION_STRENGTH)+".dat";
	}

	void reset(){
		NumberOfABuffer.clear();
		PotEnergyBuffer.clear();
		NumberOfTriedDisplacements = 0;
		NumberOfAcceptedDisplacements = 0;
		NumberOfTriedTypeChanges = 0;
		NumberOfAcceptedTypeChanges = 0;
		P.NumberOfVerletListBuilds = 0;
	}

	void runDisplacementStep() {
			NumberOfTriedDisplacements++;
			int RandomParticleID = static_cast<int>(RNG.drawRandomNumber()*static_cast<double>(TOTAL_NUMBER_OF_PARTICLES));
			double Deltas [DIMENSION];
			for (int i = 0; i < DIMENSION; i++){
				Deltas[i] = RNG.drawRandomNumber(-MAXIMUM_DISPLACEMENT, MAXIMUM_DISPLACEMENT)*INVERSE_BOX_LENGTH;
			}
			double AcceptanceProbability = exp(-P.computeChangeInPotentialEnergyByMoving(RandomParticleID, Deltas)*Beta);
			if (AcceptanceProbability >= 1.0 || (RNG.drawRandomNumber() < AcceptanceProbability)){
				P.updatePosition(RandomParticleID, Deltas);
				NumberOfAcceptedDisplacements++;
			}
	}

	void runTypeChange() {
		NumberOfTriedTypeChanges++;
		int NumberOfAParticles = P.getNumberOfAParticles();
		int NumberOfBParticles = P.getNumberOfBParticles();
		int RandomParticleIndexInTypeArray;
		int RandomParticleID;
		double ParticleNumbersPrefactor;
		double ChemicalPotentialSign;
		bool ParticleSwitchAllowed = false;

		if (RNG.drawRandomNumber() <= 0.5){
			if (NumberOfAParticles > MinNumberOfA){
				ParticleSwitchAllowed = true;
				RandomParticleIndexInTypeArray = static_cast<int>(RNG.drawRandomNumber()*static_cast<double>(NumberOfAParticles));
				RandomParticleID = P.TypeAParticleIndices[RandomParticleIndexInTypeArray];
				ParticleNumbersPrefactor = static_cast<double>(NumberOfAParticles)/static_cast<double>(NumberOfBParticles+1);
				ChemicalPotentialSign = 1.0;
			}
		}
		else {
			if (NumberOfAParticles < MaxNumberOfA){
				ParticleSwitchAllowed = true;
				RandomParticleIndexInTypeArray = static_cast<int>(RNG.drawRandomNumber()*static_cast<double>(NumberOfBParticles));
				RandomParticleID = P.TypeBParticleIndices[RandomParticleIndexInTypeArray];
				ParticleNumbersPrefactor = static_cast<double>(NumberOfBParticles)/static_cast<double>(NumberOfAParticles+1);
				ChemicalPotentialSign = -1.0;
			}
		}
		if (ParticleSwitchAllowed){
			double AcceptanceProbability = ParticleNumbersPrefactor*exp(-Beta*(P.computeChangeInPotentialEnergyBySwitching(RandomParticleID)  + ChemicalPotentialSign * ChemicalPotentialDiff));
			if (AcceptanceProbability >= 1.0 || (RNG.drawRandomNumber() < AcceptanceProbability)){
				P.switchParticleType(RandomParticleID, RandomParticleIndexInTypeArray);
				NumberOfAcceptedTypeChanges++;
			}
		}
	}

	void runSimulation() {
		const auto StartTime = chrono::steady_clock::now();
		int NextUpdateTime = UPDATE_TIME_INTERVAL;
		int NextPotEnergyComputation = POT_ENERGY_UPDATE_INTERVAL;
		int StructureFactorComputationInterval = (NumberOfMCSweeps > 2*NUMBER_OF_STRUCTURE_FACTOR_AVERAGES ? static_cast<int>(0.5 * NumberOfMCSweeps/NUMBER_OF_STRUCTURE_FACTOR_AVERAGES) : 1);
		int NextStructureFactorComputation = (NumberOfMCSweeps > 2*NUMBER_OF_STRUCTURE_FACTOR_AVERAGES ? static_cast<int>(0.5 * NumberOfMCSweeps) : 1);
		writeMetaData();
		cerr << "T = " << Temperature << ". Simulation running. Progress: ";
		for (int i = 0; i < NumberOfMCSweeps; i++){
			for (int j = 0; j < TOTAL_NUMBER_OF_PARTICLES; j++){
				if (RNG.drawRandomNumber() <= DISPLACEMENT_PROBABILITY){
					runDisplacementStep();
				}
				else {
					runTypeChange();
				}
			}
			NumberOfABuffer.push_back(P.getNumberOfAParticles());
			if (i == NextPotEnergyComputation){
				PotEnergyBuffer.push_back(P.computePotentialEnergy());
				NextPotEnergyComputation += POT_ENERGY_UPDATE_INTERVAL;
			}
			if (chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() == NextUpdateTime){
				cerr << i / (NumberOfMCSweeps/100) << "%|";
				cerr.flush();
				writeResults();
				NumberOfABuffer.clear();
				PotEnergyBuffer.clear();
				NextUpdateTime += UPDATE_TIME_INTERVAL;
			}
			if (i >= static_cast<int>(NumberOfMCSweeps*0.5) && i == NextStructureFactorComputation){
				NextStructureFactorComputation += StructureFactorComputationInterval;
				SFComputator.computeNewStructureFactorValues(P.Positions, P.TypeAParticleIndices, P.TypeBParticleIndices);
			}
		}
		writeResults();
		NumberOfABuffer.clear();
		PotEnergyBuffer.clear();
		cerr << endl << "Tried displacements: " << NumberOfTriedDisplacements << "| Accepted displacements:" << NumberOfAcceptedDisplacements << "| Ratio of accepted displacements: " << static_cast<double>(NumberOfAcceptedDisplacements)/static_cast<double>(NumberOfTriedDisplacements) << endl;
		cerr << "Tried type changes: " << NumberOfTriedTypeChanges << "| Accepted type changes:" << NumberOfAcceptedTypeChanges << "| Ratio of accepted type changes: " << static_cast<double>(NumberOfAcceptedTypeChanges)/static_cast<double>(NumberOfTriedTypeChanges) << endl;
		cerr << "#VerletListBuilds: " << P.NumberOfVerletListBuilds << endl;
		cerr << "Computation time: " << chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() << " s for " << NumberOfMCSweeps << " runs" <<  endl;
	}

	void writeMetaData() const {
		string FileName("unsorted_data/NA_Series_"+FileNameString);
		string MetaDataString("BOX_LENGTH = "+to_string(BOX_LENGTH)+"AA_INTERACTION_STRENGTH = "+to_string(AA_INTERACTION_STRENGTH)+" | AB_INTERACTION_STRENGTH = "+to_string(AB_INTERACTION_STRENGTH)+" | MAXIMUM_DISPLACEMENT = "+to_string(MAXIMUM_DISPLACEMENT)+" | DISPLACEMENT_PROBABILITY = "+to_string(DISPLACEMENT_PROBABILITY)+" | MinNA = "+to_string(MinNumberOfA)+" | MaxNA "+to_string(MaxNumberOfA)+'\n');
		ofstream FileStreamToWrite;
		FileStreamToWrite.open(FileName);
		FileStreamToWrite << MetaDataString;
		FileStreamToWrite.close();
		FileName = "unsorted_data/PotEnergySeries_"+FileNameString;
		FileStreamToWrite.open(FileName);
		FileStreamToWrite << MetaDataString;
		FileStreamToWrite.close();
	}

	void writeResults() const {
		string FileName("unsorted_data/NA_Series_"+FileNameString);
		ofstream FileStreamToWrite;
		FileStreamToWrite.open(FileName, ios_base::app);
		for (int i = 0; i < NumberOfABuffer.size(); i++){
			FileStreamToWrite << NumberOfABuffer[i] << '\n';
		}
		FileStreamToWrite.close();
		FileName = "unsorted_data/PotEnergySeries_"+FileNameString;
		FileStreamToWrite.open(FileName, ios_base::app);
		for (int i = 0; i < PotEnergyBuffer.size(); i++){
			FileStreamToWrite << PotEnergyBuffer[i] << '\n';
		}
		FileStreamToWrite.close();
	}

	void writeParticleConfigurationToFile(string FileName) const {
		ofstream FileStreamToWrite;
		FileStreamToWrite.open(FileName);
		FileStreamToWrite << P;
		FileStreamToWrite.close();
	}
};

int main(int argc, char* argv[]){
	if (argc != 7){
		cerr << argc-1 <<  " arguments were given, but exactly 6 arguments are needed: Average density, MaxTemperature, Temperature Stepsize, MinTemperature (not included), NumberOfMCSweeps, AB_INTERACTION_STRENGTH. Stopping now." << endl;
		return 0;
	}
	DENSITY = atof(argv[1]);
	AB_INTERACTION_STRENGTH = atof(argv[6]);
	double CurrentTemperature = atof(argv[2]);
	double TemperatureStep = atof(argv[3]);
	double MinTemperature = atof(argv[4]);
	int NumberOfSweeps = atof(argv[5]);

	initializeBox();
	SimulationManager S(CurrentTemperature, 0.0, 0, TOTAL_NUMBER_OF_PARTICLES, NumberOfSweeps, 8.0);
	S.initialize();

	while (CurrentTemperature > MinTemperature){
		S.reset();
		S.runSimulation();
		S.writeParticleConfigurationToFile("unsorted_data/FinalParticleConfig_"+S.FileNameString);
		S.SFComputator.writeResultsToFile("unsorted_data/StructureFactors_"+S.FileNameString);
		S.SFComputator.reset();
		CurrentTemperature -= TemperatureStep;
		S.changeTemperature(CurrentTemperature);
	}
}

