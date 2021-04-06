#include <iostream>
#include <cmath>
#include <chrono>
#include <random>
#include <fstream>

static const int DIMENSION(2);
static const int TOTAL_NUMBER_OF_PARTICLES;


struct Particles {
	double Positions [DIMENSION*TOTAL_NUMBER_OF_PARTICLES];
	bool Types [TOTAL_NUMBER_OF_PARTICLES];
	int NumberOfAParticles;
	
	double getPosition(int ParticleIndex, int Coordinate) const {
		return Positions[DIMENSION*ParticleIndex+Coordinate];
	}
	
	void switchParticleType(int ParticleIndex) {
		if (Types[ParticleIndex]){
			NumberOfAParticles--;
		}
		else {
			NumberOfAParticles++;
		}
		Types[ParticleIndex] = !Types[ParticleIndex];
	}
	
	double computeChangeInPotentialEnergyByMoving(int ParticleIndex, double Deltax, double Deltay) const {
		
	}
};


