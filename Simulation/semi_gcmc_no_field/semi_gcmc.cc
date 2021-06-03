#include <iostream>
#include <cmath>
#include <chrono>
#include <fstream>
#include <vector>
#include <iomanip>
#include <stdlib.h>
#include <algorithm>
#include "../../General_Code/realRNG.h"
#include "../../General_Code/fvec.h"

using namespace std;

const int DIMENSION = 2;
const int TOTAL_NUMBER_OF_PARTICLES = 16000;

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
const int NUMBER_OF_INITIAL_RANDOMIZATION_SWEEPS = 100;

const int UPDATE_TIME_INTERVAL = 60;
const int POT_ENERGY_UPDATE_INTERVAL = 200;

const int NUMBER_OF_SAVED_STATES_PER_TEMPERATURE = 10;

const int NUMBER_OF_THREADS = 20;

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
	double MostTraveledDistancesSquared [2];
	int MostTraveledParticleIndices [2];

	friend class SimulationManager;

	public:

	int NumberOfVerletListBuilds;

	void initialize(int InitialNumberOfAParticles){
		int InitialNumberOfBParticles = TOTAL_NUMBER_OF_PARTICLES - InitialNumberOfAParticles;
		double FractionOfAParticles = static_cast<double>(InitialNumberOfAParticles)/static_cast<double>(TOTAL_NUMBER_OF_PARTICLES);
		realRNG RNG;
		TypeAParticleIndices.clear();
		TypeBParticleIndices.clear();
		for (int i = 0; i < TOTAL_NUMBER_OF_PARTICLES; i++){
			for (int j = 0; j < DIMENSION; j++){
				ChangeInCoordinates[DIMENSION*i+j] = 0.0;
			}
		}
		for (int i = 0; i < 2; i++){
			MostTraveledDistancesSquared[i] = 0.0;
			MostTraveledParticleIndices[i] = i;
		}

		NumberOfVerletListBuilds = 0;

		int NumberOfParticlesInARow(ceil(pow(static_cast<double>(TOTAL_NUMBER_OF_PARTICLES),1.0/static_cast<double>(DIMENSION))));
		double Distance(1.0/static_cast<double>(NumberOfParticlesInARow));
		#pragma omp critical(WRITE_TO_ERROR_STREAM)
		{
			cerr << "NumberOfParticlesInRow: " << NumberOfParticlesInARow << endl;
			cerr << "Distance: " << Distance*BOX_LENGTH << endl;
		}
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

	double computePairwiseParticlePotentialEnergy(double DimensionlessDistanceSquared) const {
		double InverseDistanceSquared = 1.0/(DimensionlessDistanceSquared*BOX_LENGTH_SQUARED);
		double InverseDistanceToThePowerOfSix = InverseDistanceSquared * InverseDistanceSquared * InverseDistanceSquared;
		return 4.0*(InverseDistanceToThePowerOfSix * InverseDistanceToThePowerOfSix - InverseDistanceToThePowerOfSix + POTENTIAL_CONSTANT_1 - (sqrt(DimensionlessDistanceSquared) * BOX_LENGTH - CUTOFF) * POTENTIAL_CONSTANT_2);
	}

	double computePairwiseParticlePotentialEnergy(const double* Position0, const double* Position1) const {
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
			double InteractionStrength = (ParticleTypes[ParticleIndex] == ParticleTypes[OtherParticleIndex]) ? AA_INTERACTION_STRENGTH : AB_INTERACTION_STRENGTH;
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
					double InteractionStrength = (ParticleTypes[ParticleIndex] == ParticleTypes[OtherParticleIndex] ? AA_INTERACTION_STRENGTH : AB_INTERACTION_STRENGTH);
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
		double CurrentTraveledDistanceSquared = 0.0;
		for (int k = 0; k < DIMENSION; k++){
			Positions[DIMENSION * ParticleIndex + k] += Deltas[k];
			if (Positions[DIMENSION * ParticleIndex + k] < 0.0){
				Positions[DIMENSION * ParticleIndex + k] += 1.0;
			}
			else if (Positions[DIMENSION * ParticleIndex + k] > 1.0){
				Positions[DIMENSION * ParticleIndex + k] -= 1.0;
			}
			ChangeInCoordinates[DIMENSION * ParticleIndex + k] += Deltas[k];
			CurrentTraveledDistanceSquared += ChangeInCoordinates[DIMENSION * ParticleIndex + k] * ChangeInCoordinates[DIMENSION * ParticleIndex + k];
		}
		bool TraveledDistanceIncreased = false;
		if (MostTraveledParticleIndices[0] == ParticleIndex){
			if (CurrentTraveledDistanceSquared > MostTraveledDistancesSquared[0]){
				MostTraveledDistancesSquared[0] = CurrentTraveledDistanceSquared;
				TraveledDistanceIncreased = true;
			}
			else if (CurrentTraveledDistanceSquared > MostTraveledDistancesSquared[1]){
				MostTraveledDistancesSquared[0] = CurrentTraveledDistanceSquared;
			}
			else {
				MostTraveledDistancesSquared[0] = MostTraveledDistancesSquared[1];
				MostTraveledParticleIndices[0] = MostTraveledParticleIndices[1];
				MostTraveledDistancesSquared[1] = 0.0;
				for (int i = 0; i < TOTAL_NUMBER_OF_PARTICLES; i++){
					double NewTraveledDistanceSquared = 0.0;
					for (int j = 0; j < DIMENSION; j++){
						NewTraveledDistanceSquared += ChangeInCoordinates[DIMENSION * i + j] * ChangeInCoordinates[DIMENSION * i + j];
					}
					if (NewTraveledDistanceSquared > MostTraveledDistancesSquared[1] && i != MostTraveledParticleIndices[0]){
						MostTraveledParticleIndices[1] = i;
						MostTraveledDistancesSquared[1] = NewTraveledDistanceSquared;
					}
				}
			}
		}
		else if (MostTraveledParticleIndices[1] == ParticleIndex){
			if (CurrentTraveledDistanceSquared > MostTraveledDistancesSquared[0]){
				swap(MostTraveledParticleIndices[0], MostTraveledParticleIndices[1]);
				MostTraveledDistancesSquared[1] = MostTraveledDistancesSquared[0];
				MostTraveledDistancesSquared[0] = CurrentTraveledDistanceSquared;
				TraveledDistanceIncreased = true;
			}
			else if (CurrentTraveledDistanceSquared > MostTraveledDistancesSquared[1]){
				MostTraveledDistancesSquared[1] = CurrentTraveledDistanceSquared;
				TraveledDistanceIncreased = true;
			}
			else {
				MostTraveledDistancesSquared[1] = 0.0;
				for (int i = 0; i < TOTAL_NUMBER_OF_PARTICLES; i++){
					double NewTraveledDistanceSquared = 0.0;
					for (int j = 0; j < DIMENSION; j++){
						NewTraveledDistanceSquared += ChangeInCoordinates[DIMENSION * i + j] * ChangeInCoordinates[DIMENSION * i + j];
					}
					if (NewTraveledDistanceSquared > MostTraveledDistancesSquared[1] && i != MostTraveledParticleIndices[0]){
						MostTraveledParticleIndices[1] = i;
						MostTraveledDistancesSquared[1] = NewTraveledDistanceSquared;
					}
				}
			}
		}
		else if (CurrentTraveledDistanceSquared > MostTraveledDistancesSquared[0]){
			MostTraveledDistancesSquared[1] = MostTraveledDistancesSquared[0];
			MostTraveledParticleIndices[1] = MostTraveledParticleIndices[0];
			MostTraveledDistancesSquared[0] = CurrentTraveledDistanceSquared;
			MostTraveledParticleIndices[0] = ParticleIndex;
			TraveledDistanceIncreased = true;
		}
		else if (CurrentTraveledDistanceSquared > MostTraveledDistancesSquared[1]){
			MostTraveledDistancesSquared[1] = CurrentTraveledDistanceSquared;
			MostTraveledParticleIndices[1] = ParticleIndex;
			TraveledDistanceIncreased = true;
		}
		if (TraveledDistanceIncreased && BOX_LENGTH*(sqrt(MostTraveledDistancesSquared[0])+sqrt(MostTraveledDistancesSquared[1])) > SKINDISTANCE){
			buildVerletList();
			for (int i = 0; i < TOTAL_NUMBER_OF_PARTICLES; i++){
				for (int j = 0; j < DIMENSION; j++){
					ChangeInCoordinates[DIMENSION*i+j] = 0.0;
				}
			}
			for (int i = 0; i < 2; i++){
				MostTraveledDistancesSquared[i] = 0.0;
				MostTraveledParticleIndices[i] = i;
			}
		}
	}
};

ostream& operator<<(ostream& OStream, const Particles& State){
	OStream << "#ID\tX       Y       Type | #AParticles:  " << State.getNumberOfAParticles() << "| #BParticles: " << State.getNumberOfBParticles() << "| BOX_LENGTH: " << BOX_LENGTH << endl;
	for (int i = 0; i < TOTAL_NUMBER_OF_PARTICLES; i++){
		OStream << i << "\t";
		OStream << fixed << setprecision(6) << State.getPosition(i,0) << "\t" << State.getPosition(i,1) << "\t" << State.getParticleType(i) <<  endl;
	}
	return OStream;
}

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
	string DirectoryString;

	SimulationManager(double ChemicalPotentialDiff, int MinNumberOfA, int MaxNumberOfA, int NumberOfMCSweeps):
		ChemicalPotentialDiff(ChemicalPotentialDiff),
		MinNumberOfA(MinNumberOfA),
		MaxNumberOfA(MaxNumberOfA),
		NumberOfMCSweeps(NumberOfMCSweeps),
		NumberOfTriedDisplacements(0),
		NumberOfAcceptedDisplacements(0),
		NumberOfTriedTypeChanges(0),
		NumberOfAcceptedTypeChanges(0),
		DirectoryString("fresh_data/N="+to_string(TOTAL_NUMBER_OF_PARTICLES)+"/"){
	}

	void initializeParticles() {
		P.initialize(MinNumberOfA);
		P.buildVerletList();
	}

	void readInParticleState(string FileNameInitialState) {
		P.readInParticleState(FileNameInitialState);
		P.buildVerletList();
	}

	void setTemperature(double NewTemperature){
		Temperature = NewTemperature;
		Beta = 1.0/Temperature;
	}

	void setFileNameString(int RunNumber){
		FileNameString = "N="+to_string(TOTAL_NUMBER_OF_PARTICLES)+"_T="+to_string(Temperature)+"_AvgDens="+to_string(DENSITY)+"_MCRuns="+to_string(NumberOfMCSweeps)+"_epsAB="+to_string(AB_INTERACTION_STRENGTH)+"_"+to_string(RunNumber);
	}

	void resetCountersAndBuffers(){
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

	void runSimulationForSingleTemperature(int RunCount, int NumberOfSavedStatesPerRun) {
		const auto StartTime = chrono::steady_clock::now();
		int NextUpdateTime = UPDATE_TIME_INTERVAL;
		int NextPotEnergyComputation = POT_ENERGY_UPDATE_INTERVAL;
		int NextStateSave = NumberOfMCSweeps / 2;
		int StateSaveInterval = (NumberOfMCSweeps / 2) / NumberOfSavedStatesPerRun;
		int SavedStateCount = 0;
		writeMetaData();
		#pragma omp critical(WRITE_TO_ERROR_STREAM)
		{
			cerr << "Run " << RunCount <<  ": T=" << Temperature << ". Simulation started." << endl << endl;
		}
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
			if (i == NextStateSave){
				writeParticleStateToFile(DirectoryString+"State_"+FileNameString+"_"+to_string(SavedStateCount)+".dat");
				NextStateSave += StateSaveInterval;
				SavedStateCount++;
			}
			if (chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() >= NextUpdateTime){
				#pragma omp critical(WRITE_TO_ERROR_STREAM)
				{
					cerr << "Run " << RunCount << ": T = " << Temperature <<  ". Progress: " << i / (NumberOfMCSweeps/100) << "%| Time passed: " << chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() << " s" << endl;
				}
				writeResults();
				NumberOfABuffer.clear();
				PotEnergyBuffer.clear();
				NextUpdateTime += UPDATE_TIME_INTERVAL;
			}
		}
		writeResults();
		NumberOfABuffer.clear();
		PotEnergyBuffer.clear();
		#pragma omp critical(WRITE_TO_ERROR_STREAM)
		{
			cerr << "Run " << RunCount << ": Simulation of T=" << Temperature << " finished. Simulation-metadata: " << endl;
			cerr << "Tried displacements: " << NumberOfTriedDisplacements << "| Accepted displacements: " << NumberOfAcceptedDisplacements << "| Ratio of accepted displacements: " << static_cast<double>(NumberOfAcceptedDisplacements)/static_cast<double>(NumberOfTriedDisplacements) << endl;
			cerr << "Tried type changes: " << NumberOfTriedTypeChanges << "| Accepted type changes: " << NumberOfAcceptedTypeChanges << "| Ratio of accepted type changes: " << static_cast<double>(NumberOfAcceptedTypeChanges)/static_cast<double>(NumberOfTriedTypeChanges) << endl;
			cerr << "#VerletListBuilds: " << P.NumberOfVerletListBuilds << endl;
			cerr << "Computation time: " << chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() << " s for " << NumberOfMCSweeps << " runs." <<  endl << endl;
		}
	}

	void randomizeInitialPosition(int RunCount){
		setTemperature(AA_INTERACTION_STRENGTH*10.0);
		for (int i = 0; i < NUMBER_OF_INITIAL_RANDOMIZATION_SWEEPS; i++){
			for (int j = 0; j < TOTAL_NUMBER_OF_PARTICLES; j++){
				if (RNG.drawRandomNumber() <= DISPLACEMENT_PROBABILITY){
					runDisplacementStep();
				}
				else {
					runTypeChange();
				}
			}
		}
		#pragma omp critical(WRITE_TO_ERROR_STREAM)
		{
			cerr << "Run " << RunCount << ": Initial randomization finished. Tried displacements: " << NumberOfTriedDisplacements << "| Accepted displacements:" << NumberOfAcceptedDisplacements << "| Ratio of accepted displacements: " << static_cast<double>(NumberOfAcceptedDisplacements)/static_cast<double>(NumberOfTriedDisplacements) << endl;
			cerr << "Tried type changes: " << NumberOfTriedTypeChanges << "| Accepted type changes:" << NumberOfAcceptedTypeChanges << "| Ratio of accepted type changes: " << static_cast<double>(NumberOfAcceptedTypeChanges)/static_cast<double>(NumberOfTriedTypeChanges) << endl;
			cerr << "#VerletListBuilds: " << P.NumberOfVerletListBuilds << endl << endl;
		}
	}

	void writeMetaData() const {
		string FileName(DirectoryString+"/NA_Series_"+FileNameString+".dat");
		string MetaDataString("BOX_LENGTH = "+to_string(BOX_LENGTH)+"AA_INTERACTION_STRENGTH = "+to_string(AA_INTERACTION_STRENGTH)+" | AB_INTERACTION_STRENGTH = "+to_string(AB_INTERACTION_STRENGTH)+" | MAXIMUM_DISPLACEMENT = "+to_string(MAXIMUM_DISPLACEMENT)+" | DISPLACEMENT_PROBABILITY = "+to_string(DISPLACEMENT_PROBABILITY)+" | MinNA = "+to_string(MinNumberOfA)+" | MaxNA "+to_string(MaxNumberOfA)+'\n');
		ofstream FileStreamToWrite;
		FileStreamToWrite.open(FileName);
		FileStreamToWrite << MetaDataString;
		FileStreamToWrite.close();
		FileName = DirectoryString+"/PotEnergySeries_"+FileNameString+".dat";
		FileStreamToWrite.open(FileName);
		FileStreamToWrite << MetaDataString;
		FileStreamToWrite.close();
	}

	void writeResults() const {
		string FileName(DirectoryString+"/NA_Series_"+FileNameString+".dat");
		ofstream FileStreamToWrite;
		FileStreamToWrite.open(FileName, ios_base::app);
		for (int i = 0; i < NumberOfABuffer.size(); i++){
			FileStreamToWrite << NumberOfABuffer[i] << '\n';
		}
		FileStreamToWrite.close();
		FileName = DirectoryString+"/PotEnergySeries_"+FileNameString+".dat";
		FileStreamToWrite.open(FileName, ios_base::app);
		for (int i = 0; i < PotEnergyBuffer.size(); i++){
			FileStreamToWrite << PotEnergyBuffer[i] << '\n';
		}
		FileStreamToWrite.close();
	}

	void writeParticleStateToFile(string FileName) const {
		ofstream FileStreamToWrite;
		FileStreamToWrite.open(FileName);
		FileStreamToWrite << P;
		FileStreamToWrite.close();
	}
};

void runSimulationForMultipleStartStates(double MaxTemperature, double MinTemperature, double TemperatureStep, int NumberOfRuns, int NumberOfMCSweeps) {

	const int NumberOfSavedStatesPerRun = NUMBER_OF_SAVED_STATES_PER_TEMPERATURE / NumberOfRuns;

	#pragma omp parallel num_threads(NUMBER_OF_THREADS)
	{
		SimulationManager S(0.0, 0, TOTAL_NUMBER_OF_PARTICLES, NumberOfMCSweeps);

		#pragma omp for
		for (int RunCount = 0; RunCount < NumberOfRuns; RunCount++){
			S.initializeParticles();
			S.randomizeInitialPosition(RunCount);
			for (double CurrentTemperature = MaxTemperature; CurrentTemperature > MinTemperature; CurrentTemperature -= TemperatureStep){
				S.setTemperature(CurrentTemperature);
				S.setFileNameString(RunCount);
				S.resetCountersAndBuffers();
				S.runSimulationForSingleTemperature(RunCount, NumberOfSavedStatesPerRun);
			}
		}
	}
}

int main(int argc, char* argv[]){
	double MaxTemperature;
	double TemperatureStep;
	double MinTemperature;
	int NumberOfSweeps;
	int NumberOfRuns;
	if (argc != 8){
		cerr << "WARNING: " << argc-1 <<  " arguments were given, but exactly 7 arguments are needed: Average density, MaxTemperature, Temperature Stepsize, MinTemperature (not included), NumberOfMCSweeps, AB_INTERACTION_STRENGTH, NumberOfRuns. Running with default parameters: Average density = 0.6, MaxTemperature = 1.0, Temperature Stepsize = 0.1, MinTemperature = 0.85, NumberOfMCSweeps = 100, AB_INTERACTION_STRENGTH = 0.1, NumberOfRuns = 2" << endl;
		DENSITY = 0.6;
		AB_INTERACTION_STRENGTH = 0.1;
		MaxTemperature = 1.0;
		TemperatureStep = 0.1;
		MinTemperature = 0.85;
		NumberOfSweeps = 100;
		NumberOfRuns = 2;
	}
	else {
		DENSITY = atof(argv[1]);
		AB_INTERACTION_STRENGTH = atof(argv[6]);
		MaxTemperature = atof(argv[2]);
		TemperatureStep = atof(argv[3]);
		MinTemperature = atof(argv[4]);
		NumberOfSweeps = atoi(argv[5]);
		NumberOfRuns = atoi(argv[7]);
	}

	initializeBox();

	runSimulationForMultipleStartStates(MaxTemperature, MinTemperature, TemperatureStep, NumberOfRuns, NumberOfSweeps);
}

