#ifndef PARTICLES_INCLUDED
#define PARTICLES_INCLUDED

#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <stdlib.h>
#include <algorithm>

#include "../global_constants.h"
#include "simulation_config.h"
#include "../particle_type.h"
#include "../fvec.h"
#include "../realRNG.h"
#include "../utility_functions.h"

using namespace std;


class Particles {
	private:
		double Positions [DIMENSION*TOTAL_NUMBER_OF_PARTICLES];
		double Velocities [DIMENSION*TOTAL_NUMBER_OF_PARTICLES];

		ParticleType ParticleTypes [TOTAL_NUMBER_OF_PARTICLES];

		fvec<int, TOTAL_NUMBER_OF_PARTICLES> TypeAParticleIndices;
		fvec<int, TOTAL_NUMBER_OF_PARTICLES> TypeBParticleIndices;

		double Density;
		double BoxLength;
		double BoxLengthSquared;
		double InverseBoxLength;
		double Volume;
		int NumberOfSubdivisions;

		double xDisplacement;
		double ShearRate;

		vector<int> CellListHead;
		int CellListIndices [TOTAL_NUMBER_OF_PARTICLES];
		int VerletListHead [2*TOTAL_NUMBER_OF_PARTICLES];
		vector<int> VerletIndicesOfNeighbors;

		vector<int> CellNeighborsHead;
		vector<int> CellNeighborsIndices;

		double ChangeInCoordinates [DIMENSION*TOTAL_NUMBER_OF_PARTICLES];
		double MostTraveledDistancesSquared [2];
		int MostTraveledParticleIndices [2];

		double TotalChangeInCoordinates [DIMENSION*TOTAL_NUMBER_OF_PARTICLES];

		double NumberOfVerletListBuilds;

		void buildCellList(){
			CellListHead.clear();
			CellListHead.resize(NumberOfSubdivisions*NumberOfSubdivisions,-1);
			int CurrentCellIndex;
			int IndexFactor;
			for (int ParticleIndex = 0; ParticleIndex < TOTAL_NUMBER_OF_PARTICLES; ParticleIndex++){
				CurrentCellIndex = 0;
				IndexFactor = 1;
				for (int j = 0; j < DIMENSION; j++){
					CurrentCellIndex += static_cast<int>(static_cast<double>(NumberOfSubdivisions)*Positions[DIMENSION*ParticleIndex+j])*IndexFactor;
					IndexFactor *= NumberOfSubdivisions;
				}
				CellListIndices[ParticleIndex] = CellListHead[CurrentCellIndex];
				CellListHead[CurrentCellIndex] = ParticleIndex;
			}
		}

		void buildVerletList() {
			buildCellList();

			CellNeighborsHead.clear();
			CellNeighborsHead.reserve(2*NumberOfSubdivisions*NumberOfSubdivisions);

			CellNeighborsIndices.clear();
			CellNeighborsIndices.reserve(10*NumberOfSubdivisions*NumberOfSubdivisions);
			if (NumberOfSubdivisions == 1){
				CellNeighborsHead.push_back(0);
				CellNeighborsHead.push_back(1);
				CellNeighborsIndices.push_back(0);
			}
			else {
				int CurrentIndexInCellNeighbors = 0;

				// cells at the bottom of the box
				int Indices [DIMENSION]{0,0};
				for (; Indices[0] < NumberOfSubdivisions; Indices[0]++){
					CellNeighborsHead.push_back(CurrentIndexInCellNeighbors);
					int NumberOfNeighbors = 0;

					int IndicesOffsets [DIMENSION]{-1,0};
					while (IndicesOffsets[DIMENSION - 1] < 2){
						int NeighborCell = 0;
						int IndexFactor = 1;
						for (int i = 0; i < DIMENSION; i++){
							int OtherCellIndex = Indices[i]+IndicesOffsets[i];
							if (OtherCellIndex < 0){
								OtherCellIndex += NumberOfSubdivisions;
							}
							else if (OtherCellIndex >= NumberOfSubdivisions){
								OtherCellIndex -= NumberOfSubdivisions;
							}
							NeighborCell += IndexFactor*OtherCellIndex;
							IndexFactor *= NumberOfSubdivisions;
						}
						CellNeighborsIndices.push_back(NeighborCell);
						NumberOfNeighbors++;
						CurrentIndexInCellNeighbors++;

						IndicesOffsets[0]++;
						for (int i = 0; i < DIMENSION - 1; i++){
							if (IndicesOffsets[i] >= 2){
								IndicesOffsets[i] = -1;
								IndicesOffsets[i+1]++;
							}
						}
					}
					int xOffsetShear = static_cast<int>(xDisplacement * static_cast<double>(NumberOfSubdivisions));
					for (int xOffset = -1; xOffset < 3; xOffset++){
						int NeighborCell = NumberOfSubdivisions*(NumberOfSubdivisions-1); // shift to top row of the cells
						int xOtherCellIndex = xOffset+xOffsetShear+Indices[0];
						if (xOtherCellIndex < 0){
							xOtherCellIndex += NumberOfSubdivisions;
						}
						else if (xOtherCellIndex >= NumberOfSubdivisions){
							xOtherCellIndex = xOtherCellIndex % NumberOfSubdivisions;
						}
						NeighborCell += xOtherCellIndex;

						CellNeighborsIndices.push_back(NeighborCell);
						NumberOfNeighbors++;
						CurrentIndexInCellNeighbors++;
					}

					CellNeighborsHead.push_back(NumberOfNeighbors);
				}

				// cells in the middle of the box
				Indices[0] = 0;
				Indices[1] = 1;
				while (Indices[1] < NumberOfSubdivisions - 1){
					CellNeighborsHead.push_back(CurrentIndexInCellNeighbors);
					int NumberOfNeighbors = 0;

					int IndicesOffsets [DIMENSION]{-1,-1};
					while (IndicesOffsets[DIMENSION - 1] < 2){
						int NeighborCell = 0;
						int IndexFactor = 1;
						for (int i = 0; i < DIMENSION; i++){
							int OtherCellIndex = Indices[i]+IndicesOffsets[i];
							if (OtherCellIndex < 0){
								OtherCellIndex += NumberOfSubdivisions;
							}
							else if (OtherCellIndex >= NumberOfSubdivisions){
								OtherCellIndex -= NumberOfSubdivisions;
							}
							NeighborCell += IndexFactor*OtherCellIndex;
							IndexFactor *= NumberOfSubdivisions;
						}
						CellNeighborsIndices.push_back(NeighborCell);
						NumberOfNeighbors++;
						CurrentIndexInCellNeighbors++;
						IndicesOffsets[0]++;
						for (int i = 0; i < DIMENSION - 1; i++){
							if (IndicesOffsets[i] >= 2){
								IndicesOffsets[i] = -1;
								IndicesOffsets[i+1]++;
							}
						}
					}
					CellNeighborsHead.push_back(NumberOfNeighbors);

					Indices[0]++;
					if (Indices[0] >= NumberOfSubdivisions){
						Indices[0] = 0;
						Indices[1]++;
					}
				}

				// cells at the top of the box
				Indices[0] = 0;
				Indices[1] = NumberOfSubdivisions-1;
				for (; Indices[0] < NumberOfSubdivisions; Indices[0]++){
					CellNeighborsHead.push_back(CurrentIndexInCellNeighbors);
					int NumberOfNeighbors = 0;

					int IndicesOffsets [DIMENSION]{-1,-1};
					while (IndicesOffsets[DIMENSION - 1] < 1){
						int NeighborCell = 0;
						int IndexFactor = 1;
						for (int i = 0; i < DIMENSION; i++){
							int OtherCellIndex = Indices[i]+IndicesOffsets[i];
							if (OtherCellIndex < 0){
								OtherCellIndex += NumberOfSubdivisions;
							}
							else if (OtherCellIndex >= NumberOfSubdivisions){
								OtherCellIndex -= NumberOfSubdivisions;
							}
							NeighborCell += IndexFactor*OtherCellIndex;
							IndexFactor *= NumberOfSubdivisions;
						}
						CellNeighborsIndices.push_back(NeighborCell);
						NumberOfNeighbors++;
						CurrentIndexInCellNeighbors++;

						IndicesOffsets[0]++;
						for (int i = 0; i < DIMENSION - 1; i++){
							if (IndicesOffsets[i] >= 2){
								IndicesOffsets[i] = -1;
								IndicesOffsets[i+1]++;
							}
						}
					}
					int xOffsetShear = NumberOfSubdivisions-static_cast<int>(xDisplacement * static_cast<double>(NumberOfSubdivisions));
					for (int xOffset = -2; xOffset < 2; xOffset++){
						int NeighborCell = 0;
						int xOtherCellIndex = xOffset+xOffsetShear+Indices[0];
						if (xOtherCellIndex < 0){
							xOtherCellIndex += NumberOfSubdivisions;
						}
						else if (xOtherCellIndex >= NumberOfSubdivisions){
							xOtherCellIndex = xOtherCellIndex % NumberOfSubdivisions;
						}
						NeighborCell += xOtherCellIndex;

						CellNeighborsIndices.push_back(NeighborCell);
						NumberOfNeighbors++;
						CurrentIndexInCellNeighbors++;
					}

					CellNeighborsHead.push_back(NumberOfNeighbors);
				}
			}

			VerletIndicesOfNeighbors.clear();
			VerletIndicesOfNeighbors.reserve(40*TOTAL_NUMBER_OF_PARTICLES);
			int CurrentIndexInVerletIndices = 0;

			for (int i = 0; i < CellNeighborsHead.size(); i+=2){
				int CurrentParticleIndex = CellListHead[i/2];
				while (CurrentParticleIndex >= 0){
					int NumberOfNeighbors = 0;
					VerletListHead[2*CurrentParticleIndex] = CurrentIndexInVerletIndices;
					for (int j = 0; j < CellNeighborsHead[i+1]; j++){
						int OtherParticleIndex = CellListHead[CellNeighborsIndices[CellNeighborsHead[i] + j]];
						while (OtherParticleIndex >= 0){
							if (OtherParticleIndex != CurrentParticleIndex){
								double xCoordinateDifference = Positions[DIMENSION * CurrentParticleIndex] - Positions[DIMENSION * OtherParticleIndex];
								double yCoordinateDifference = Positions[DIMENSION * CurrentParticleIndex + 1] - Positions[DIMENSION * OtherParticleIndex + 1];
								if (yCoordinateDifference > 0.5){
									yCoordinateDifference -= 1.0;
									xCoordinateDifference -= xDisplacement;
								}
								else if (yCoordinateDifference <= -0.5){
									yCoordinateDifference += 1.0;
									xCoordinateDifference += xDisplacement;
								}
								while (xCoordinateDifference > 0.5){
									xCoordinateDifference -= 1.0;
								}
								while (xCoordinateDifference <= -0.5){
									xCoordinateDifference += 1.0;
								}
								double DistanceSquared = xCoordinateDifference * xCoordinateDifference + yCoordinateDifference * yCoordinateDifference;
								if (DistanceSquared * BoxLengthSquared <= MAX_VERLET_DIST_SQUARED){
									NumberOfNeighbors++;
									VerletIndicesOfNeighbors.push_back(OtherParticleIndex);
									CurrentIndexInVerletIndices++;
								}
							}
							OtherParticleIndex = CellListIndices[OtherParticleIndex];
						}
					}
					VerletListHead[2*CurrentParticleIndex+1] = NumberOfNeighbors;
					CurrentParticleIndex = CellListIndices[CurrentParticleIndex];
				}
			}

			NumberOfVerletListBuilds++;
			resetTraveledDistances();
		}

		void resetTraveledDistances() {
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

		void updateBoxParametersWithDensity(double NewDensity){
			Density = NewDensity;
			BoxLength = sqrt(static_cast<double>(TOTAL_NUMBER_OF_PARTICLES) / Density);
			Volume = BoxLength*BoxLength;
			BoxLengthSquared = BoxLength * BoxLength;
			InverseBoxLength = 1.0/BoxLength;
			NumberOfSubdivisions = static_cast<int>(BoxLength/MAX_VERLET_DIST) > 3 ? static_cast<int>(BoxLength/MAX_VERLET_DIST) : 1;
		}

		void updateBoxParametersWithVolumeChange(double VolumeChange){
			Volume += VolumeChange;
			Density = static_cast<double>(TOTAL_NUMBER_OF_PARTICLES)/Volume;
			BoxLength = sqrt(Volume);
			BoxLengthSquared = BoxLength * BoxLength;
			InverseBoxLength = 1.0/BoxLength;
			NumberOfSubdivisions = static_cast<int>(BoxLength/MAX_VERLET_DIST) > 3 ? static_cast<int>(BoxLength/MAX_VERLET_DIST) : 1;
		}

		void computeForceOnParticle(int ParticleIndex, double* ForceOnParticle) {
			*(ForceOnParticle  ) = 0.0;
			*(ForceOnParticle+1) = 0.0;

			double x0 = Positions[DIMENSION*ParticleIndex];
			double y0 = Positions[DIMENSION*ParticleIndex+1];

			for (int i = 0; i < VerletListHead[2*ParticleIndex+1]; i++){
				int OtherParticleIndex = VerletIndicesOfNeighbors[VerletListHead[2*ParticleIndex]+i];

				double x1 = Positions[DIMENSION*OtherParticleIndex];
				double y1 = Positions[DIMENSION*OtherParticleIndex+1];
				double Deltax = x1 - x0;
				double Deltay = y1 - y0;
				if (Deltay > 0.5){
					Deltay -= 1.0;
					Deltax -= xDisplacement;
				}
				else if (Deltay <= -0.5){
					Deltay += 1.0;
					Deltax += xDisplacement;
				}
				while (Deltax > 0.5){
					Deltax -= 1.0;
				}
				while (Deltax <= -0.5){
					Deltax += 1.0;
				}
				double DistanceSquared = (Deltax*Deltax+Deltay*Deltay)*BoxLengthSquared;
				if (DistanceSquared < CUTOFF_SQUARED){
					double InteractionStrength = (ParticleTypes[ParticleIndex] == ParticleTypes[OtherParticleIndex] ? AA_INTERACTION_STRENGTH : AB_INTERACTION_STRENGTH);
					double MagnitudeOfForce = computePairwiseMagnitudeOfForce(DistanceSquared);
					*(ForceOnParticle  ) += -MagnitudeOfForce*InteractionStrength*Deltax*BoxLength;
					*(ForceOnParticle+1) += -MagnitudeOfForce*InteractionStrength*Deltay*BoxLength;
				}
			}
		}

		double computePairwiseParticlePotentialEnergy(double DimensionlessDistanceSquared) const {
			double InverseDistanceSquared = 1.0/(DimensionlessDistanceSquared*BoxLengthSquared);
			double InverseDistanceToThePowerOfSix = InverseDistanceSquared * InverseDistanceSquared * InverseDistanceSquared;
			return 4.0*(InverseDistanceToThePowerOfSix * InverseDistanceToThePowerOfSix - InverseDistanceToThePowerOfSix + POTENTIAL_CONSTANT_1 - (sqrt(DimensionlessDistanceSquared) * BoxLength - CUTOFF) * POTENTIAL_CONSTANT_2);
		}

		double computePairwiseParticlePotentialEnergy(const double* Position0, const double* Position1) const {
			double xCoordinateDifference = Position0[0] - Position1[0];
			double yCoordinateDifference = Position0[1] - Position1[1];
			if (yCoordinateDifference > 0.5){
				yCoordinateDifference -= 1.0;
				xCoordinateDifference -= xDisplacement;
			}
			else if (yCoordinateDifference <= -0.5){
				yCoordinateDifference += 1.0;
				xCoordinateDifference += xDisplacement;
			}
			while (xCoordinateDifference > 0.5){
				xCoordinateDifference -= 1.0;
			}
			while (xCoordinateDifference <= -0.5){
				xCoordinateDifference += 1.0;
			}
			double DistanceSquared = xCoordinateDifference * xCoordinateDifference + yCoordinateDifference * yCoordinateDifference;
			if (DistanceSquared*BoxLengthSquared >= CUTOFF_SQUARED){
				return 0.0;
			}
			return computePairwiseParticlePotentialEnergy(DistanceSquared);
		}

	public:

		Particles(double ShearRate):
			NumberOfVerletListBuilds(0.0),
			xDisplacement(0.0),
			ShearRate(ShearRate)
		{
			for (int i = 0; i < TOTAL_NUMBER_OF_PARTICLES; i++){
				for (int j = 0; j < DIMENSION; j++){
					ChangeInCoordinates[DIMENSION*i+j] = 0.0;
					TotalChangeInCoordinates[DIMENSION*i+j] = 0.0;
				}
			}
		}

		Particles(int NumberOfAParticles, double InitialDensity, double InitialTemperature, double InitialxDisplacement, double InitialShearRate){
			xDisplacement = InitialxDisplacement;
			ShearRate = InitialShearRate;

			NumberOfVerletListBuilds = 0.0;
			updateBoxParametersWithDensity(InitialDensity);

			int NumberOfBParticles = TOTAL_NUMBER_OF_PARTICLES - NumberOfAParticles;
			double FractionOfAParticles = static_cast<double>(NumberOfAParticles)/static_cast<double>(TOTAL_NUMBER_OF_PARTICLES);

			TypeAParticleIndices.clear();
			TypeBParticleIndices.clear();
			for (int i = 0; i < TOTAL_NUMBER_OF_PARTICLES; i++){
				for (int j = 0; j < DIMENSION; j++){
					ChangeInCoordinates[DIMENSION*i+j] = 0.0;
					TotalChangeInCoordinates[DIMENSION*i+j] = 0.0;
				}
			}
			for (int i = 0; i < 2; i++){
				MostTraveledDistancesSquared[i] = 0.0;
				MostTraveledParticleIndices[i] = i;
			}

			int NumberOfParticlesInARow(ceil(pow(static_cast<double>(TOTAL_NUMBER_OF_PARTICLES),1.0/static_cast<double>(DIMENSION))));
			double Distance(1.0/static_cast<double>(NumberOfParticlesInARow));
			#pragma omp critical(WRITE_TO_ERROR_STREAM)
			{
				cerr << "NumberOfParticlesInRow: " << NumberOfParticlesInARow << endl;
				cerr << "Distance: " << Distance*BoxLength << endl;
			}
			double CurrentPosition [DIMENSION];
			for (int i = 0; i < DIMENSION; i++){
				CurrentPosition[i] = Distance*0.5;
			}

			realUniformRNG UniformRNG;
			initializeVelocities(InitialTemperature);
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
				if ((UniformRNG.drawRandomNumber() <= FractionOfAParticles && NumberOfAParticlesInitialized < NumberOfAParticles) || (TOTAL_NUMBER_OF_PARTICLES - ParticlesInitialized <= NumberOfAParticles - NumberOfAParticlesInitialized)){
					ParticleTypes[ParticlesInitialized] = ParticleType::A;
					TypeAParticleIndices.push_back(ParticlesInitialized);
					NumberOfAParticlesInitialized++;
				}
				else {
					ParticleTypes[ParticlesInitialized] = ParticleType::B;
					TypeBParticleIndices.push_back(ParticlesInitialized);
				}
			}
			buildVerletList();
			cerr << "Initial potential energy: " << computePotentialEnergy() << endl;
		}

		void readInParticleState(string FileNameToReadIn, int Skiplines, int StateNumber, double Density) {
			updateBoxParametersWithDensity(Density);

			TypeAParticleIndices.clear();
			TypeBParticleIndices.clear();
			ifstream FileStreamToReadIn;
			FileStreamToReadIn.open(FileNameToReadIn);

			skipLines(FileStreamToReadIn, Skiplines+StateNumber*(TOTAL_NUMBER_OF_PARTICLES+2));

			string CurrentString;
			for (int i = 0; i < 4; i++){
				getline(FileStreamToReadIn, CurrentString, '|');
			}
			getline(FileStreamToReadIn, CurrentString, ':');
			getline(FileStreamToReadIn, CurrentString, '\n');
			xDisplacement = stod(CurrentString);

			for (int ParticleIndex = 0; ParticleIndex < TOTAL_NUMBER_OF_PARTICLES; ParticleIndex++){
				for (int j = 0; j < DIMENSION; j++){
					getline(FileStreamToReadIn, CurrentString, '\t');
					Positions[DIMENSION*ParticleIndex+j] = stod(CurrentString);
				}
				for (int j = 0; j < DIMENSION; j++){
					getline(FileStreamToReadIn, CurrentString, '\t');
					Velocities[DIMENSION*ParticleIndex+j] = stod(CurrentString);
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
			buildVerletList();
		}

		void readInParticlePositions(string FileNameToReadIn, int Skiplines, int StateNumber, double Density) {
			updateBoxParametersWithDensity(Density);
			xDisplacement = 0.0;

			TypeAParticleIndices.clear();
			TypeBParticleIndices.clear();
			ifstream FileStreamToReadIn;
			FileStreamToReadIn.open(FileNameToReadIn);

			skipLines(FileStreamToReadIn, Skiplines+StateNumber*(TOTAL_NUMBER_OF_PARTICLES+2));

			string CurrentString;

			for (int ParticleIndex = 0; ParticleIndex < TOTAL_NUMBER_OF_PARTICLES; ParticleIndex++){
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
			buildVerletList();
			cerr << "Initial potential energy: " << computePotentialEnergy() << endl;
		}

		void initializeVelocities(double InitialTemperature){
			realGaussianRNG GaussianRNG(0.0,sqrt(InitialTemperature));
			double TotalMomentum [DIMENSION]{};
			for (int ParticlesInitialized = 0; ParticlesInitialized < TOTAL_NUMBER_OF_PARTICLES; ParticlesInitialized++){
				for (int i = 0; i < DIMENSION; i++){
					Velocities[ParticlesInitialized*DIMENSION + i] = GaussianRNG.drawRandomNumber()*InverseBoxLength;
					TotalMomentum[i] += Velocities[ParticlesInitialized*DIMENSION + i];
				}
			}
			double KineticEnergy = 0.0;
			for (int ParticleIndex = 0; ParticleIndex < TOTAL_NUMBER_OF_PARTICLES; ParticleIndex++){
				for (int i = 0; i < DIMENSION; i++){
					Velocities[DIMENSION*ParticleIndex + i] -= TotalMomentum[i]/static_cast<double>(TOTAL_NUMBER_OF_PARTICLES);
					KineticEnergy += Velocities[DIMENSION*ParticleIndex + i] * Velocities[DIMENSION*ParticleIndex + i];
				}
			}
			for (int ParticleIndex = 0; ParticleIndex < TOTAL_NUMBER_OF_PARTICLES; ParticleIndex++){
				for (int i = 0; i < DIMENSION; i++){
					Velocities[DIMENSION*ParticleIndex + i] *= sqrt(2.0 * static_cast<double>(TOTAL_NUMBER_OF_PARTICLES) * InitialTemperature / (KineticEnergy * BoxLengthSquared));
				}
			}
			cerr << "Initial kinetic energy: " << computeKineticEnergy() << endl;
			computeTotalMomentum(TotalMomentum);
			cerr << fixed << setprecision(numeric_limits<long double>::digits10+1) << "TotalMomentum: (" << TotalMomentum[0] << ',' << TotalMomentum[1] << ')' << endl;
		}

		static double computePairwiseMagnitudeOfForce(double DistanceSquared) { //intentionally off by a factor r, so we can just multiply with r-vector to get the correct result
			double InverseDistanceSquared = 1.0/DistanceSquared;
			double InverseDistanceToThePowerOfSix = InverseDistanceSquared * InverseDistanceSquared * InverseDistanceSquared;
			return (48.0*InverseDistanceToThePowerOfSix*InverseDistanceToThePowerOfSix*InverseDistanceSquared-24.0*InverseDistanceToThePowerOfSix*InverseDistanceSquared+4.0*POTENTIAL_CONSTANT_2*sqrt(InverseDistanceSquared));
		}

		double getPosition(int ParticleIndex, int Coordinate) const {
			return Positions[DIMENSION*ParticleIndex+Coordinate];
		}

		double getVelocity(int ParticleIndex, int Coordinate) const {
			return Velocities[DIMENSION*ParticleIndex+Coordinate];
		}

		ParticleType getParticleType(int ParticleIndex) const {
			return ParticleTypes[ParticleIndex];
		}

		double getBoxLength() const {
			return BoxLength;
		}

		double getBoxLengthSquared() const {
			return BoxLengthSquared;
		}

		double getInverseBoxLength() const {
			return InverseBoxLength;
		}

		double getDensity() const {
			return Density;
		}

		double getVolume() const {
			return Volume;
		}

		double getxDisplacement() const {
			return xDisplacement;
		}

		void setShearRate(double NewShearRate) {
			ShearRate = NewShearRate;
		}

		void resetCounters() {
			NumberOfVerletListBuilds = 0.0;
		}

		int getNumberOfAParticles() const {
			return TypeAParticleIndices.size();
		}

		int getNumberOfBParticles() const {
			return TypeBParticleIndices.size();
		}

		double getNumberOfVerletListBuilds() const {
			return NumberOfVerletListBuilds;
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

		double computePotentialEnergy() const {
			double TotPotEnergy = 0.0;
			#pragma omp parallel num_threads(NUMBER_OF_THREADS)
			{
				double PotEnergy = 0.0;
				#pragma omp for
				for (int ParticleIndex = 0; ParticleIndex < TOTAL_NUMBER_OF_PARTICLES; ParticleIndex++){
					for (int i = 0; i < VerletListHead[2*ParticleIndex+1]; i++){
						int OtherParticleIndex = VerletIndicesOfNeighbors[VerletListHead[2*ParticleIndex]+i];
						if (OtherParticleIndex < ParticleIndex){
							double InteractionStrength = (ParticleTypes[ParticleIndex] == ParticleTypes[OtherParticleIndex] ? AA_INTERACTION_STRENGTH : AB_INTERACTION_STRENGTH);
							double PairwisePotEnergy = InteractionStrength * computePairwiseParticlePotentialEnergy(&Positions[DIMENSION*OtherParticleIndex], &Positions[DIMENSION*ParticleIndex]);
							PotEnergy += PairwisePotEnergy;
						}
					}
				}
				#pragma omp atomic
				TotPotEnergy += PotEnergy;
			}
			return TotPotEnergy;
		}

		double computeKineticEnergy() const {
			double KineticEnergy = 0.0;
			for (int ParticleIndex = 0; ParticleIndex < TOTAL_NUMBER_OF_PARTICLES; ParticleIndex++){
					double xVelWithoutDrift = Velocities[DIMENSION*ParticleIndex] - ShearRate*(Positions[DIMENSION*ParticleIndex+1]-0.5);
					KineticEnergy += xVelWithoutDrift * xVelWithoutDrift;
					KineticEnergy += Velocities[DIMENSION*ParticleIndex+1]*Velocities[DIMENSION*ParticleIndex+1];
			}
			return KineticEnergy*BoxLengthSquared*0.5;
		}

		double computeEnergy() const {
			return computePotentialEnergy()+computeKineticEnergy();
		}

		void computeTotalMomentum(double* Momentum) {
			*(Momentum + 0) = 0.0;
			*(Momentum + 1) = 0.0;
			for (int ParticleIndex = 0; ParticleIndex < TOTAL_NUMBER_OF_PARTICLES; ParticleIndex++) {
				for (int Coordinate = 0; Coordinate < DIMENSION; Coordinate++) {
					*(Momentum + Coordinate) += Velocities[DIMENSION*ParticleIndex+Coordinate];
				}
			}
		}

		void updatePosition(int ParticleIndex, const double* Deltas){
			Positions[DIMENSION * ParticleIndex] += Deltas[0];
			Positions[DIMENSION * ParticleIndex + 1] += Deltas[1];
			if (Positions[DIMENSION * ParticleIndex + 1] < 0.0){
				Positions[DIMENSION * ParticleIndex] += xDisplacement;
				Positions[DIMENSION * ParticleIndex + 1] += 1.0;
				Velocities[DIMENSION * ParticleIndex] += ShearRate;
			}
			else if (Positions[DIMENSION * ParticleIndex + 1] >= 1.0){
				Positions[DIMENSION * ParticleIndex] -= xDisplacement;
				Positions[DIMENSION * ParticleIndex + 1] -= 1.0;
				Velocities[DIMENSION * ParticleIndex] -= ShearRate;
			}
			while (Positions[DIMENSION * ParticleIndex] < 0.0){
				Positions[DIMENSION * ParticleIndex] += 1.0;
			}
			while (Positions[DIMENSION * ParticleIndex] >= 1.0){
				Positions[DIMENSION * ParticleIndex] -= 1.0;
			}
			if (Positions[DIMENSION * ParticleIndex] >= 1.0 || Positions[DIMENSION * ParticleIndex] < 0.0 || Positions[DIMENSION * ParticleIndex + 1] >= 1.0 || Positions[DIMENSION * ParticleIndex + 1] < 0.0)
			{
				cerr << "ERROR: Particle " << ParticleIndex << " has wrong coordinates!" << endl;
				cerr << "Position: (" << Positions[DIMENSION*ParticleIndex] << "," << Positions[DIMENSION * ParticleIndex + 1] << "), Deltas: (" << Deltas[0] << "," << Deltas[1] << ")" << endl;
				cerr << "X       Y       Type | #AParticles:  " << fixed << setprecision(numeric_limits<long double>::digits10+1) << getNumberOfAParticles() << "| #BParticles: " << getNumberOfBParticles() << "| BoxLength: " << getBoxLength() << "| xDisplacement: " << getxDisplacement() << endl;
				for (int i = 0; i < TOTAL_NUMBER_OF_PARTICLES; i++){
					cerr << getPosition(i,0) << "\t" << getPosition(i,1) << "\t" << getParticleType(i) <<  endl;
				}
			}
			ChangeInCoordinates[DIMENSION * ParticleIndex] += Deltas[0];
			ChangeInCoordinates[DIMENSION * ParticleIndex + 1] += Deltas[1];
			TotalChangeInCoordinates[DIMENSION * ParticleIndex] += Deltas[0];
			TotalChangeInCoordinates[DIMENSION * ParticleIndex + 1] += Deltas[1];
			double CurrentTraveledDistanceSquared = ChangeInCoordinates[DIMENSION * ParticleIndex] * ChangeInCoordinates[DIMENSION * ParticleIndex] + ChangeInCoordinates[DIMENSION * ParticleIndex + 1] * ChangeInCoordinates[DIMENSION * ParticleIndex + 1];
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
			if (TraveledDistanceIncreased && BoxLength*(sqrt(MostTraveledDistancesSquared[0])+sqrt(MostTraveledDistancesSquared[1])) > SKINDISTANCE){
				buildVerletList();
			}
		}

		void computeForces(double* Forces){
			#pragma omp parallel num_threads(NUMBER_OF_THREADS)
			{
				#pragma omp for
				for (int ParticleIndex = 0; ParticleIndex < TOTAL_NUMBER_OF_PARTICLES; ParticleIndex++){
					computeForceOnParticle(ParticleIndex, Forces+DIMENSION*ParticleIndex);
				}
			}
		}

		void applyVelocityVerlet(double Stepsize) {
			double CurrentForces [DIMENSION*TOTAL_NUMBER_OF_PARTICLES];
			computeForces(CurrentForces);
			for (int ParticleIndex = 0; ParticleIndex < TOTAL_NUMBER_OF_PARTICLES; ParticleIndex++) {
				double Deltas [DIMENSION];
				for (int i = 0; i < DIMENSION; i++) {
					Deltas[i] = Velocities[DIMENSION*ParticleIndex+i]*Stepsize + 0.5*Stepsize*Stepsize*InverseBoxLength*CurrentForces[DIMENSION*ParticleIndex+i];
				}
				updatePosition(ParticleIndex, Deltas);
			}
			double NextForces [DIMENSION*TOTAL_NUMBER_OF_PARTICLES];
			computeForces(NextForces);
			for (int ParticleIndex = 0; ParticleIndex < TOTAL_NUMBER_OF_PARTICLES; ParticleIndex++){
				for (int i = 0; i < DIMENSION; i++){
					Velocities[DIMENSION*ParticleIndex+i] += 0.5*Stepsize*InverseBoxLength*(CurrentForces[DIMENSION*ParticleIndex+i]+NextForces[DIMENSION*ParticleIndex+i]);
				}
			}
		}

		void moveImageBoxes(double Stepsize){
			xDisplacement += ShearRate*Stepsize;
			xDisplacement -= static_cast<double>(static_cast<int>(xDisplacement));
			buildVerletList();
		}

		void rescaleVelocities(double RescalingFactor) {
			for (int ParticleIndex = 0; ParticleIndex < TOTAL_NUMBER_OF_PARTICLES; ParticleIndex++){
					//Velocities[DIMENSION*ParticleIndex] -= ShearRate*(Positions[DIMENSION*ParticleIndex+1]-0.5);
					//Velocities[DIMENSION*ParticleIndex] *= RescalingFactor;
					//Velocities[DIMENSION*ParticleIndex] += ShearRate*(Positions[DIMENSION*ParticleIndex+1]-0.5);
					Velocities[DIMENSION*ParticleIndex + 1] *= RescalingFactor;
			}
		}

		void updateAverageVelocities(vector<double>& AvgVelocities, const int NumberOfyValues) const {
			vector<double> NewVelocityEntries(NumberOfyValues,0.0);
			vector<int> NumberOfValues(NumberOfyValues,0);
			for (int ParticleID = 0; ParticleID < TOTAL_NUMBER_OF_PARTICLES; ParticleID++){
				int Index = static_cast<int>(getPosition(ParticleID,1)*static_cast<double>(NumberOfyValues));
				NewVelocityEntries[Index] += Velocities[DIMENSION*ParticleID];
				NumberOfValues[Index]++;
			}
			for (int i = 0; i < NumberOfyValues; i++){
				AvgVelocities.push_back(NewVelocityEntries[i] / static_cast<double>(NumberOfValues[i]));
			}
		}

		double computeAverageMSD() const {
			double AverageMSD = 0.0;
			for (int ParticleIndex = 0; ParticleIndex < TOTAL_NUMBER_OF_PARTICLES; ParticleIndex++){
				AverageMSD += TotalChangeInCoordinates[DIMENSION * ParticleIndex]*TotalChangeInCoordinates[DIMENSION * ParticleIndex] + TotalChangeInCoordinates[DIMENSION * ParticleIndex + 1] * TotalChangeInCoordinates[DIMENSION * ParticleIndex + 1];
			}
			return AverageMSD*BoxLengthSquared;
		}
};

ostream& operator<<(ostream& OStream, const Particles& State){
	OStream << "x\ty\tvX\tvy\tType\t| #AParticles:  " << fixed << setprecision(numeric_limits<long double>::digits10+1) << State.getNumberOfAParticles() << "| #BParticles: " << State.getNumberOfBParticles() << "| BoxLength: " << State.getBoxLength() << "| xDisplacement: " << State.getxDisplacement() << endl;
	for (int i = 0; i < TOTAL_NUMBER_OF_PARTICLES; i++){
		OStream << State.getPosition(i,0) << '\t' << State.getPosition(i,1) << '\t' << State.getVelocity(i,0) << '\t' << State.getVelocity(i,1) << '\t' << State.getParticleType(i) <<  endl;
	}
	return OStream;
}

#endif //PARTICLES_INCLUDED
