#include <iostream>
#include <cmath>
#include <chrono>
#include <random>
#include <fstream>

using namespace std;

static const int DIMENSION = 2;
static const int TOTAL_NUMBER_OF_PARTICLES = 1000;
static const double AA_INTERACTION_STRENGTH = 1.0;
static const double AB_INTERACTION_STRENGTH = 0.5;
static const double CUTOFF = 2.5;
static const double CUTOFF_SQUARED = CUTOFF * CUTOFF;
static const double INVERSE_CUTOFF = 1.0/CUTOFF;
static const double BOX_LENGTH = 34.0;
static const double MAXIMUM_DISPLACEMENT = 0.1;

class realRNG{
	private:
		mt19937 rng;
		uniform_real_distribution<double> UnifRealDist;
	public:
		realRNG(double Minimum, double Maximum):
			UnifRealDist(Minimum,Maximum){
			random_device rd;
			seed_seq sd{rd(),rd()};
			rng = mt19937(sd);
		}
		
		double drawRandomNumber(){
			return UnifRealDist(rng);
		}
};

struct Particles {
	double Positions [DIMENSION*TOTAL_NUMBER_OF_PARTICLES];
	int ParticleTypeBoundaryIndex;
	
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
};

struct SimulationManager{
	Particles P;
	double Temperature;
	
};

int main(){
}

