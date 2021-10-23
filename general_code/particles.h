#ifndef PARTICLES_INCLUDED
#define PARTICLES_INCLUDED

#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <stdlib.h>
#include <algorithm>

#include "global_constants.h"
#include "simulation_config.h"
#include "particle_type.h"
#include "fvec.h"
#include "realRNG.h"
#include "utility_functions.h"

using namespace std;


class Particles {
	public:
		double Positions [DIMENSION*TOTAL_NUMBER_OF_PARTICLES];

		ParticleType ParticleTypes [TOTAL_NUMBER_OF_PARTICLES];

		fvec<int, TOTAL_NUMBER_OF_PARTICLES> TypeAParticleIndices;
		fvec<int, TOTAL_NUMBER_OF_PARTICLES> TypeBParticleIndices;

		double Density;
		double BoxLength;
		double BoxLengthSquared;
		double InverseBoxLength;
		double Volume;
		int NumberOfSubdivisions;

		vector<int> CellListHead;
		int CellListIndices [TOTAL_NUMBER_OF_PARTICLES];
		int VerletListHead [2*TOTAL_NUMBER_OF_PARTICLES];
		vector<int> VerletIndicesOfNeighbors;

		double ChangeInCoordinates [DIMENSION*TOTAL_NUMBER_OF_PARTICLES];
		double MostTraveledDistancesSquared [2];
		int MostTraveledParticleIndices [2];

		double NumberOfVerletListBuilds;

		void buildCellList(){
			int NumberOfSubcells = NumberOfSubdivisions;
			for (int i = 0; i < DIMENSION-1; i++){
				NumberOfSubcells *= NumberOfSubdivisions;
			}
			CellListHead.clear();
			CellListHead.resize(NumberOfSubcells,-1);
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

		void buildVerletList(double xDisplacement) {
			buildCellList();

			vector<int> CellNeighborsHead;
			CellNeighborsHead.reserve(2*NumberOfSubdivisions*NumberOfSubdivisions);
			vector<int> CellNeighborsIndices;
			CellNeighborsIndices.reserve(10*NumberOfSubdivisions*NumberOfSubdivisions);
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

			cerr << "xDisplacement = " << xDisplacement << endl;
			for (int i = 0; i < CellNeighborsHead.size(); i+=2){
				cerr << i/2 << ": ";
				for (int j = 0; j < CellNeighborsHead[i+1]; j++){
					cerr << CellNeighborsIndices[CellNeighborsHead[i] + j] << ",";
				}
				cerr << endl;
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
			cerr << endl << "Resulting neighbors: " << endl;
			for (int i = 0; i < TOTAL_NUMBER_OF_PARTICLES; i++){
				cerr << i << ": ";
				for (int j = 0; j < VerletListHead[2*i+1]; j++){
					cerr << VerletIndicesOfNeighbors[VerletListHead[2*i]+j] << ",";
				}
				cerr << endl;
			}
			cerr << endl;
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

		double computePairwiseParticlePotentialEnergy(double DimensionlessDistanceSquared) const {
			double InverseDistanceSquared = 1.0/(DimensionlessDistanceSquared*BoxLengthSquared);
			double InverseDistanceToThePowerOfSix = InverseDistanceSquared * InverseDistanceSquared * InverseDistanceSquared;
			return 4.0*(InverseDistanceToThePowerOfSix * InverseDistanceToThePowerOfSix - InverseDistanceToThePowerOfSix + POTENTIAL_CONSTANT_1 - (sqrt(DimensionlessDistanceSquared) * BoxLength - CUTOFF) * POTENTIAL_CONSTANT_2);
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
			if (DistanceSquared*BoxLengthSquared >= CUTOFF_SQUARED){
				return 0.0;
			}
			return computePairwiseParticlePotentialEnergy(DistanceSquared);
		}

	public:

		double getPosition(int ParticleIndex, int Coordinate) const {
			return Positions[DIMENSION*ParticleIndex+Coordinate];
		}

		ParticleType getParticleType(int ParticleIndex) const {
			return ParticleTypes[ParticleIndex];
		}

		double getBoxLength() const {
			return BoxLength;
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

		void initialize(int InitialNumberOfAParticles, double InitialDensity){
			updateBoxParametersWithDensity(InitialDensity);
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

			NumberOfVerletListBuilds = 0.0;

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
			buildVerletList(0.0);
		}

		void readInParticleState(string FileNameToReadIn, int StateNumber, double Density) {
			updateBoxParametersWithDensity(Density);
			TypeAParticleIndices.clear();
			TypeBParticleIndices.clear();
			ifstream FileStreamToReadIn;
			FileStreamToReadIn.open(FileNameToReadIn);

			skipLines(FileStreamToReadIn, 2+StateNumber*(TOTAL_NUMBER_OF_PARTICLES+2));

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
			buildVerletList(0.0);
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

		double computeChangeInPotentialEnergyByMoving(int ParticleIndex, const double* Delta) const {
			double PotEnergyChange = 0.0;
			double UpdatedCoordinates [DIMENSION];
			for (int i = 0; i < DIMENSION; i++){
				UpdatedCoordinates[i] = Positions[DIMENSION*ParticleIndex + i] + *(Delta + i);
				if (UpdatedCoordinates[i] < 0.0){
					UpdatedCoordinates[i] += static_cast<double>(static_cast<int>(abs(UpdatedCoordinates[i]))+1);
				}
				else if (UpdatedCoordinates[i] >= 1.0){
					UpdatedCoordinates[i] -= static_cast<double>(static_cast<int>(abs(UpdatedCoordinates[i])));
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

		double computeChangeInPotentialEnergyBySwitching(int ParticleIndexInTypeArray, ParticleType ParticleTypeBefore) const {
			int ParticleIndex = (ParticleTypeBefore == ParticleType::A) ? TypeAParticleIndices[ParticleIndexInTypeArray] : TypeBParticleIndices[ParticleIndexInTypeArray];
			double PotEnergyChange = 0.0;
			for (int i = 0; i < VerletListHead[2*ParticleIndex+1]; i++){
				int OtherParticleIndex = VerletIndicesOfNeighbors[VerletListHead[2*ParticleIndex]+i];
				double PrefactorDifference = (ParticleTypes[ParticleIndex] == ParticleTypes[OtherParticleIndex]) ? AB_INTERACTION_STRENGTH - AA_INTERACTION_STRENGTH : AA_INTERACTION_STRENGTH - AB_INTERACTION_STRENGTH;
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

		double computeChangeInPotentialEnergyByChangingVolume(double VolumeChange) {
			double PotEnergyBefore = computePotentialEnergy();
			updateBoxParametersWithVolumeChange(VolumeChange);
			buildVerletList(0.0);
			resetTraveledDistances();
			return computePotentialEnergy() - PotEnergyBefore;
		}

		void updatePosition(int ParticleIndex, const double* Deltas){
			double CurrentTraveledDistanceSquared = 0.0;
			for (int k = 0; k < DIMENSION; k++){
				Positions[DIMENSION * ParticleIndex + k] += Deltas[k];

				if (Positions[DIMENSION * ParticleIndex + k] < 0.0){
					Positions[DIMENSION * ParticleIndex + k] += static_cast<double>(static_cast<int>(abs(Positions[DIMENSION * ParticleIndex + k]))+1);
				}
				else if (Positions[DIMENSION * ParticleIndex + k] >= 1.0){
					Positions[DIMENSION * ParticleIndex + k] -= static_cast<double>(static_cast<int>(abs(Positions[DIMENSION * ParticleIndex + k])));
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
			if (TraveledDistanceIncreased && BoxLength*(sqrt(MostTraveledDistancesSquared[0])+sqrt(MostTraveledDistancesSquared[1])) > SKINDISTANCE){
				buildVerletList(0.0);
				resetTraveledDistances();
			}
		}

		void switchParticleType(int ParticleIndexInTypeArray, ParticleType ParticleTypeBefore) {
			if (ParticleTypeBefore == ParticleType::A){
				TypeBParticleIndices.push_back(TypeAParticleIndices[ParticleIndexInTypeArray]);
				ParticleTypes[TypeAParticleIndices[ParticleIndexInTypeArray]] = ParticleType::B;
				TypeAParticleIndices.erase(ParticleIndexInTypeArray);
			}
			else {
				TypeAParticleIndices.push_back(TypeBParticleIndices[ParticleIndexInTypeArray]);
				ParticleTypes[TypeBParticleIndices[ParticleIndexInTypeArray]] = ParticleType::A;
				TypeBParticleIndices.erase(ParticleIndexInTypeArray);
			}
		}

		void changeVolume(double VolumeChange) {
			updateBoxParametersWithVolumeChange(VolumeChange);
			buildVerletList(0.0);
		}
};

ostream& operator<<(ostream& OStream, const Particles& State){
	OStream << "X       Y       Type | #AParticles:  " << fixed << setprecision(numeric_limits<long double>::digits10+1) << State.getNumberOfAParticles() << "| #BParticles: " << State.getNumberOfBParticles() << "| BoxLength: " << State.getBoxLength() << endl;
	for (int i = 0; i < TOTAL_NUMBER_OF_PARTICLES; i++){
		OStream << State.getPosition(i,0) << "\t" << State.getPosition(i,1) << "\t" << State.getParticleType(i) <<  endl;
	}
	return OStream;
}

#endif //PARTICLES_INCLUDED
