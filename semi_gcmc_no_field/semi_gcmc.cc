#include <iostream>
#include <cmath>
#include <chrono>
#include <random>
#include <fstream>
#include <vector>

using namespace std;

static const int DIMENSION = 2;
static const int TOTAL_NUMBER_OF_PARTICLES = 1000;
static const double AA_INTERACTION_STRENGTH = 1.0;
static const double AB_INTERACTION_STRENGTH = 0.5;
static const double CUTOFF = 2.5;
static const double CUTOFF_SQUARED = CUTOFF * CUTOFF;
static const double INVERSE_CUTOFF = 1.0/CUTOFF;
static const double BOX_LENGTH = 34.0;
static const double INVERSE_BOX_LENGTH = 1.0/BOX_LENGTH;
static const double MAXIMUM_DISPLACEMENT = 0.1;
static const double MAX_VERLET_DIST = 1.3*CUTOFF;
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
	

	void initialize(){
		ParticleTypeBoundaryIndex = static_cast<int>(round(0.5*static_cast<double>(TOTAL_NUMBER_OF_PARTICLES)) - 1.0);

		int NumberOfParticlesInARow(ceil(sqrt(static_cast<double>(TOTAL_NUMBER_OF_PARTICLES))));
		double Distance(BOX_LENGTH/static_cast<double>(NumberOfParticlesInARow));
		double Currentx(Distance*0.5);
		double Currenty(Distance*0.5);
		int ParticlesInitialized(0);
		while (ParticlesInitialized < TOTAL_NUMBER_OF_PARTICLES){
			if ((Currentx > BOX_LENGTH) || (Currenty > BOX_LENGTH)){
				cerr << Currentx << "," << Currenty << endl;
			}
			Positions[ParticlesInitialized*DIMENSION] = Currentx;
			Positions[ParticlesInitialized*DIMENSION+1] = Currenty;
			Currentx += Distance;
			if (Currentx >= BOX_LENGTH){
				Currentx = Distance*0.5;
				Currenty += Distance;
			}
			ParticlesInitialized++;
		}
	}
	
	double getPosition(int ParticleIndex, int Coordinate) const {
		return Positions[DIMENSION*ParticleIndex+Coordinate];
	}
	
	void switchParticleType(int ParticleIndex) {
		double Temp [DIMENSION] = {Positions[DIMENSION*ParticleIndex], Positions[DIMENSION*ParticleIndex + 1]};
		if (ParticleIndex <= ParticleTypeBoundaryIndex){
			Positions[ParticleIndex*DIMENSION] = Positions[ParticleTypeBoundaryIndex*DIMENSION];
			Positions[ParticleIndex*DIMENSION + 1] = Positions[ParticleTypeBoundaryIndex*DIMENSION + 1];
			Positions[ParticleTypeBoundaryIndex*DIMENSION] = Temp[0];
			Positions[ParticleTypeBoundaryIndex*DIMENSION + 1] = Temp[1];
			ParticleTypeBoundaryIndex--;
		}
		else {
			Positions[ParticleIndex*DIMENSION] = Positions[(ParticleTypeBoundaryIndex + 1)*DIMENSION];
			Positions[ParticleIndex*DIMENSION + 1] = Positions[(ParticleTypeBoundaryIndex + 1)*DIMENSION + 1];
			Positions[(ParticleTypeBoundaryIndex + 1)*DIMENSION] = Temp[0];
			Positions[(ParticleTypeBoundaryIndex + 1)*DIMENSION + 1] = Temp[1];
			ParticleTypeBoundaryIndex++;
		}
	}
	
	int getNumberOfAParticles() const {
		return ParticleTypeBoundaryIndex + 1;
	}
	
	int getNumberOfBParticles() const {
		return TOTAL_NUMBER_OF_PARTICLES - getNumberOfAParticles();
	}
	
	static double computePairwiseParticlePotentialEnergy(double Distance) {
		double InverseDistance(1.0/Distance);
		return (pow(InverseDistance, 12.0) - pow(InverseDistance, 6.0) - pow(INVERSE_CUTOFF, 12.0) + pow(INVERSE_CUTOFF, 6.0) - (Distance - CUTOFF) * ((-12.0) * pow(InverseDistance, 13.0) + 6.0 * pow(InverseDistance, 7.0)));
	}
	
	double computePairwiseParticlePotentialEnergy(int ParticleIndex, double xPos, double yPos) const {
		double CoordinateDifferences [DIMENSION] = {Positions[DIMENSION*ParticleIndex] - xPos, Positions[DIMENSION*ParticleIndex + 1] - yPos};
		for (int i = 0; i < DIMENSION; i++){
			if (CoordinateDifferences[i] > 0.5 * BOX_LENGTH){
				CoordinateDifferences[i] -= BOX_LENGTH;
			}
			else if (CoordinateDifferences[i] <= - 0.5 * BOX_LENGTH ){
				CoordinateDifferences[i] += BOX_LENGTH;
			}
		}
		double DistanceSquared  = CoordinateDifferences[0] * CoordinateDifferences[0] + CoordinateDifferences[1] * CoordinateDifferences[1];
		if (DistanceSquared >= CUTOFF_SQUARED){
			return 0.0;
		}
		return computePairwiseParticlePotentialEnergy(sqrt(DistanceSquared));
	}
	
	double computeChangeInPotentialEnergyByMoving(int ParticleIndex, double Deltax, double Deltay) const {
		double PotEnergyBefore = 0.0;
		double PotEnergyAfter = 0.0;
		double UpdatedCoordinates [DIMENSION] = {Positions[DIMENSION*ParticleIndex] + Deltax, Positions[DIMENSION*ParticleIndex + 1] + Deltay};
		for (int i = 0; i < DIMENSION; i++){
			if (UpdatedCoordinates[i] < 0.0){
				UpdatedCoordinates[i] += BOX_LENGTH;
			}
			else if (UpdatedCoordinates[i] > BOX_LENGTH){
				UpdatedCoordinates[i] -= BOX_LENGTH;
			}
		}
		for (int i = 0; i < TOTAL_NUMBER_OF_PARTICLES; i++){
			if (i != ParticleIndex){
				double InteractionStrength = AB_INTERACTION_STRENGTH;
				if (((i <= ParticleTypeBoundaryIndex) && (ParticleIndex <= ParticleTypeBoundaryIndex)) || ((i > ParticleTypeBoundaryIndex) && (ParticleIndex > ParticleTypeBoundaryIndex)) ) {
					InteractionStrength = AA_INTERACTION_STRENGTH;
				}
				PotEnergyBefore += InteractionStrength * computePairwiseParticlePotentialEnergy(i, Positions[DIMENSION*ParticleIndex], Positions[DIMENSION*ParticleIndex + 1]);
				PotEnergyAfter += InteractionStrength * computePairwiseParticlePotentialEnergy(i, UpdatedCoordinates[0], UpdatedCoordinates[1]);
			}
		}
		return PotEnergyAfter - PotEnergyBefore;
	}
	
	double computeChangeInPotentialEnergyBySwitching(int ParticleIndex) const {
		double PotEnergyChange = 0.0;
		for (int i = 0; i < TOTAL_NUMBER_OF_PARTICLES; i++){
			if (i != ParticleIndex){
				double PrefactorDifference = AA_INTERACTION_STRENGTH - AB_INTERACTION_STRENGTH;
				if ((i <= ParticleTypeBoundaryIndex && ParticleIndex > ParticleIndex) || (i > ParticleTypeBoundaryIndex && ParticleIndex <= ParticleIndex)){
					PrefactorDifference *= -1.0;
				}
				PotEnergyChange += PrefactorDifference * computePairwiseParticlePotentialEnergy(i, Positions[DIMENSION*ParticleIndex], Positions[DIMENSION*ParticleIndex + 1]);
			}
		}
		return PotEnergyChange;
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
		for (int i = 0; i < TOTAL_NUMBER_OF_PARTICLES; i++){
			CurrentCellIndex = 0;
			IndexFactor = 1;
			for (int j = 0; j < DIMENSION; j++){
				CurrentCellIndex += static_cast<int>(static_cast<double>(NUMBER_OF_SUBDIVISIONS)*Positions[DIMENSION*i+j]*INVERSE_BOX_LENGTH)*IndexFactor;
				IndexFactor *= NUMBER_OF_SUBDIVISIONS;
			}
				CellListIndices[i] = CellListHead[CurrentCellIndex];
				CellListHead[CurrentCellIndex] = i;
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
			cerr << '\n';
		}
	}
};

ostream& operator<<(ostream& OStream, const Particles& State){
	OStream << "#ID     X       Y       ParticleTypeBoundaryIndex: " << State.ParticleTypeBoundaryIndex << endl;
	for (int i = 0; i < TOTAL_NUMBER_OF_PARTICLES; i++){
		OStream << i << "\t";
		OStream << State.getPosition(i,0) << "\t" << State.getPosition(i,1) << endl;
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
				double Deltas [DIMENSION] = {RNG.drawRandomNumber(-MAXIMUM_DISPLACEMENT, MAXIMUM_DISPLACEMENT), RNG.drawRandomNumber(-MAXIMUM_DISPLACEMENT, MAXIMUM_DISPLACEMENT)};
				double PotentialEnergyChange = P.computeChangeInPotentialEnergyByMoving(RandomParticleID, Deltas[0], Deltas[1]);
				double AcceptanceProbability = exp(-PotentialEnergyChange*Beta);
				if (AcceptanceProbability >= 1.0 || (RNG.drawRandomNumber(0.0, 1.0) < AcceptanceProbability)){
					P.Positions[DIMENSION * RandomParticleID] += Deltas[0];
					P.Positions[DIMENSION * RandomParticleID + 1] += Deltas[1];
					for (int k = 0; k < DIMENSION; k++){
						if (P.Positions[DIMENSION * RandomParticleID + k] < 0.0){
							P.Positions[DIMENSION * RandomParticleID + k] += BOX_LENGTH;
						}
						else if (P.Positions[DIMENSION * RandomParticleID + k] > BOX_LENGTH){
							P.Positions[DIMENSION * RandomParticleID + k] -= BOX_LENGTH;
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
}

