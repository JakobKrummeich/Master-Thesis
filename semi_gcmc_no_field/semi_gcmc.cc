#include <iostream>
#include <cmath>
#include <chrono>
#include <random>
#include <fstream>
#include <vector>
#include <array>
#include <iomanip>

using namespace std;

static const int DIMENSION = 2;
static const int TOTAL_NUMBER_OF_PARTICLES = 33;
static const double AA_INTERACTION_STRENGTH = 1.0;
static const double AB_INTERACTION_STRENGTH = 0.5;
static const double CUTOFF = 2.5;
static const double CUTOFF_SQUARED = CUTOFF * CUTOFF;
static const double INVERSE_CUTOFF = 1.0/CUTOFF;
static const double BOX_LENGTH = 13.0;
static const double BOX_LENGTH_SQUARED = BOX_LENGTH * BOX_LENGTH;
static const double INVERSE_BOX_LENGTH = 1.0/BOX_LENGTH;
static const double MAXIMUM_DISPLACEMENT = 0.1;
static const double MAX_VERLET_DIST = 1.3*CUTOFF;
static const double MAX_VERLET_DIST_SQUARED = MAX_VERLET_DIST * MAX_VERLET_DIST;
static const int NUMBER_OF_SUBDIVISIONS = static_cast<int>(BOX_LENGTH/MAX_VERLET_DIST);
static const int MIN_NUMBER_OF_SUBDIVISIONS = 4;


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
		
		double drawRandomNumber(double Min, double Max){
			return (Max - Min) * UnifRealDist(rng) + Min;
		}
};

struct Particles {
	double Positions [DIMENSION*TOTAL_NUMBER_OF_PARTICLES];
	int ParticleTypeBoundaryIndex;
	
	vector<int> CellListHead;
	int CellListIndices [TOTAL_NUMBER_OF_PARTICLES];
	int VerletListHead [2*TOTAL_NUMBER_OF_PARTICLES];
	vector<int> VerletIndicesOfNeighbors;

	void initialize(){
		ParticleTypeBoundaryIndex = static_cast<int>(round(0.5*static_cast<double>(TOTAL_NUMBER_OF_PARTICLES)) - 1.0);

		int NumberOfParticlesInARow(ceil(pow(static_cast<double>(TOTAL_NUMBER_OF_PARTICLES),1.0/static_cast<double>(DIMENSION))));
		double Distance(1.0/static_cast<double>(NumberOfParticlesInARow));
		cerr << "NumberOfParticlesInRow: " << NumberOfParticlesInARow << endl;
		cerr << "Distance: " << Distance*BOX_LENGTH << endl;
		array<double,DIMENSION> CurrentPosition;
		CurrentPosition.fill(Distance*0.5);
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
		}
	}
	
	double getPosition(int ParticleIndex, int Coordinate) const {
		return Positions[DIMENSION*ParticleIndex+Coordinate];
	}
	
	void switchParticleType(int ParticleIndex) {
		int SwapParticleIndex;
		if (ParticleIndex <= ParticleTypeBoundaryIndex){
			SwapParticleIndex = ParticleTypeBoundaryIndex;
			ParticleTypeBoundaryIndex--;
		}
		else {
			SwapParticleIndex = ParticleTypeBoundaryIndex + 1;
			ParticleTypeBoundaryIndex++;
		}
		double Temp [DIMENSION];
		for (int i = 0; i < DIMENSION; i++){
			Temp[i] = Positions[DIMENSION*ParticleIndex + i];
			Positions[ParticleIndex*DIMENSION + i] = Positions[SwapParticleIndex*DIMENSION + i];
			Positions[SwapParticleIndex*DIMENSION + i] = Temp[i];
		}
	}
	
	int getNumberOfAParticles() const {
		return ParticleTypeBoundaryIndex + 1;
	}
	
	int getNumberOfBParticles() const {
		return TOTAL_NUMBER_OF_PARTICLES - getNumberOfAParticles();
	}
	
	static double computePairwiseParticlePotentialEnergy(double DimensionlessDistance) {
		double InverseDistance(1.0/(DimensionlessDistance*BOX_LENGTH));
		return (pow(InverseDistance, 12.0) - pow(InverseDistance, 6.0) - pow(INVERSE_CUTOFF, 12.0) + pow(INVERSE_CUTOFF, 6.0) - (DimensionlessDistance * BOX_LENGTH - CUTOFF) * ((-12.0) * pow(INVERSE_CUTOFF, 13.0) + 6.0 * pow(INVERSE_CUTOFF, 7.0)));
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
		return computePairwiseParticlePotentialEnergy(sqrt(DistanceSquared));
	}
	
	double computeChangeInPotentialEnergyByMoving(int ParticleIndex, const double* Delta) const {
		double PotEnergyDiff = 0.0;
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
		for (int OtherParticleIndex = 0; OtherParticleIndex < TOTAL_NUMBER_OF_PARTICLES; OtherParticleIndex++){
			if (OtherParticleIndex != ParticleIndex){
				double InteractionStrength = AB_INTERACTION_STRENGTH;
				if (((OtherParticleIndex <= ParticleTypeBoundaryIndex) && (ParticleIndex <= ParticleTypeBoundaryIndex)) || ((OtherParticleIndex > ParticleTypeBoundaryIndex) && (ParticleIndex > ParticleTypeBoundaryIndex)) ) {
					InteractionStrength = AA_INTERACTION_STRENGTH;
				}
				PotEnergyDiff += InteractionStrength * (computePairwiseParticlePotentialEnergy(&Positions[DIMENSION*OtherParticleIndex], &Positions[DIMENSION*ParticleIndex]) - computePairwiseParticlePotentialEnergy(&Positions[DIMENSION*OtherParticleIndex], UpdatedCoordinates));
			}
		}
		return PotEnergyDiff;
	}
	
	double computeChangeInPotentialEnergyBySwitching(int ParticleIndex) const {
		double PotEnergyChange = 0.0;
		for (int OtherParticleIndex = 0; OtherParticleIndex < TOTAL_NUMBER_OF_PARTICLES; OtherParticleIndex++){
			if (OtherParticleIndex != ParticleIndex){
				double PrefactorDifference = AA_INTERACTION_STRENGTH - AB_INTERACTION_STRENGTH;
				if ((OtherParticleIndex <= ParticleTypeBoundaryIndex && ParticleIndex > ParticleIndex) || (OtherParticleIndex > ParticleTypeBoundaryIndex && ParticleIndex <= ParticleIndex)){
					PrefactorDifference *= -1.0;
				}
				PotEnergyChange += PrefactorDifference * computePairwiseParticlePotentialEnergy(&Positions[DIMENSION*OtherParticleIndex], &Positions[DIMENSION*ParticleIndex]);
			}
		}
		return PotEnergyChange;
	}
	
	void buildCellList(){
		int NumberOfSubcells = NUMBER_OF_SUBDIVISIONS > 3 ? NUMBER_OF_SUBDIVISIONS : 1;
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
		VerletIndicesOfNeighbors.reserve(30*TOTAL_NUMBER_OF_PARTICLES);
		int CurrentIndexInVerletIndices = 0;
		
		int Indices [DIMENSION];
		for (int i = 0; i < DIMENSION; i++){
			Indices[i] = 0;
		}
		while (Indices[DIMENSION - 1] < NUMBER_OF_SUBDIVISIONS){
			int Cell = 0;
			int IndexFactor = 1;
			int IndicesOffsets [DIMENSION];
			for (int i = 0; i < DIMENSION; i++){
				Cell += Indices[i]*IndexFactor;
				IndexFactor *= NUMBER_OF_SUBDIVISIONS;
				IndicesOffsets[i] = -1;
			}
			int CurrentParticleIndex = CellListHead[Cell];
			while (CurrentParticleIndex >= 0){
				int NumberOfNeighbors = 0;
				VerletListHead[2*CurrentParticleIndex] = CurrentIndexInVerletIndices;
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
							if (DistanceSquared * BOX_LENGTH * BOX_LENGTH <= MAX_VERLET_DIST_SQUARED){
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
	
};

ostream& operator<<(ostream& OStream, const Particles& State){
	OStream << "#ID\tX       Y       ParticleTypeBoundaryIndex: " << State.ParticleTypeBoundaryIndex << endl;
	for (int i = 0; i < TOTAL_NUMBER_OF_PARTICLES; i++){
		OStream << i << "\t";
		OStream << fixed << setprecision(5) << State.getPosition(i,0) << "\t" << State.getPosition(i,1) << endl;
	}
	return OStream;
}

struct SimulationManager {
	Particles P;
	double Temperature;
	double Beta;
	realRNG RNG;
	
	void runCanonicalSteps(int NumberOfSteps) {
		for (int i = 0; i < NumberOfSteps; i++){
			for (int j = 0; j < TOTAL_NUMBER_OF_PARTICLES; j++){
				int RandomParticleID = static_cast<int>(RNG.drawRandomNumber(0.0,1.0)*static_cast<double>(TOTAL_NUMBER_OF_PARTICLES));
				double Deltas [DIMENSION];
				for (int i = 0; i < DIMENSION; i++){
				 Deltas[i] = RNG.drawRandomNumber(-MAXIMUM_DISPLACEMENT, MAXIMUM_DISPLACEMENT)*INVERSE_BOX_LENGTH;
				}
				double PotentialEnergyChange = P.computeChangeInPotentialEnergyByMoving(RandomParticleID, Deltas);
				double AcceptanceProbability = exp(-PotentialEnergyChange*Beta);
				if (AcceptanceProbability >= 1.0 || (RNG.drawRandomNumber(0.0, 1.0) < AcceptanceProbability)){
					for (int k = 0; k < DIMENSION; k++){
						P.Positions[DIMENSION * RandomParticleID + k] += Deltas[k];
						if (P.Positions[DIMENSION * RandomParticleID + k] < 0.0){
							P.Positions[DIMENSION * RandomParticleID + k] += 1.0;
						}
						else if (P.Positions[DIMENSION * RandomParticleID + k] > 1.0){
							P.Positions[DIMENSION * RandomParticleID + k] -= 1.0;
						}
					}
				}
			}
		}
	}

};

int main(){
	SimulationManager S;
	S.P.initialize();
	cerr << S.P;
	cerr << "Number Of A particles: " << S.P.getNumberOfAParticles() << endl;
	
	S.P.buildCellList();
	S.P.printCellList();
	S.P.buildVerletList();
	S.P.printVerletList();
}

