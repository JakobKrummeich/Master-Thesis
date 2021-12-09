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
	private:
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

		double xDisplacement;

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

		class StressComputator {
			private:
				const Particles* P;
				vector<double> StressOfEdgesInyDirection;
				vector<int> NumberOfForceValuesInyDirection;
				vector<double> StressOfEdgesInxDirection;
				vector<int> NumberOfForceValuesInxDirection;
				int NumberOfAverages;

				double EdgeLength;
				int NumberOfSubdivisions;
				double DimensionlessEdgeLength;

				vector<int> CellListHead;
				int CellListIndices [TOTAL_NUMBER_OF_PARTICLES];

				vector<int> CellListHeadDisplacedCells;
				int CellListIndicesDisplacedCells [TOTAL_NUMBER_OF_PARTICLES];

				void buildCellLists(){
					CellListHead.clear();
					CellListHead.resize(NumberOfSubdivisions*NumberOfSubdivisions,-1);
					CellListHeadDisplacedCells.clear();
					CellListHeadDisplacedCells.resize(2*NumberOfSubdivisions,-1);
					int CurrentCellIndex;
					int IndexFactor;
					for (int ParticleIndex = 0; ParticleIndex < TOTAL_NUMBER_OF_PARTICLES; ParticleIndex++){
						CurrentCellIndex = 0;
						IndexFactor = 1;
						for (int j = 0; j < DIMENSION; j++){
							CurrentCellIndex += static_cast<int>(static_cast<double>(NumberOfSubdivisions)*P->Positions[DIMENSION*ParticleIndex+j])*IndexFactor;
							IndexFactor *= NumberOfSubdivisions;
						}
						if (static_cast<int>(static_cast<double>(NumberOfSubdivisions)*P->Positions[DIMENSION*ParticleIndex+1]) == 0) {
							double DisplacedxPosition = P->Positions[DIMENSION*ParticleIndex]+P->xDisplacement;
							if (DisplacedxPosition >= 1.0){
								DisplacedxPosition -= 1.0;
							}
							int DisplacedxIndex = static_cast<int>(static_cast<double>(NumberOfSubdivisions)*DisplacedxPosition);
							CellListIndicesDisplacedCells[ParticleIndex] = CellListHeadDisplacedCells[DisplacedxIndex+NumberOfSubdivisions];
							CellListHeadDisplacedCells[DisplacedxIndex+NumberOfSubdivisions] = ParticleIndex;
						}
						else if (static_cast<int>(static_cast<double>(NumberOfSubdivisions)*P->Positions[DIMENSION*ParticleIndex+1]) >= NumberOfSubdivisions-1) {
							double DisplacedxPosition = P->Positions[DIMENSION*ParticleIndex]-P->xDisplacement;
							if (DisplacedxPosition < 0.0){
								DisplacedxPosition += 1.0;
							}
							int DisplacedxIndex = static_cast<int>(static_cast<double>(NumberOfSubdivisions)*DisplacedxPosition);
							CellListIndicesDisplacedCells[ParticleIndex] = CellListHeadDisplacedCells[DisplacedxIndex];
							CellListHeadDisplacedCells[DisplacedxIndex] = ParticleIndex;
						}
						CellListIndices[ParticleIndex] = CellListHead[CurrentCellIndex];
						CellListHead[CurrentCellIndex] = ParticleIndex;
					}
				}

				void printCellLists() const {
					cerr << "normal CellList: " << endl;
					for (int i = 0; i < CellListHead.size(); i++){
						cerr << "Cell " << i << ": " << endl;
						int CurrentParticleIndex = CellListHead[i];
						while (CurrentParticleIndex >= 0){
							cerr << CurrentParticleIndex << ": (" << P->Positions[2*CurrentParticleIndex] << "," << P->Positions[2*CurrentParticleIndex+1] << ") : " << "(" << static_cast<int>(static_cast<double>(NumberOfSubdivisions)*P->Positions[2*CurrentParticleIndex]) << "," << static_cast<int>(static_cast<double>(NumberOfSubdivisions)*P->Positions[2*CurrentParticleIndex+1]) << ")" << endl;
							CurrentParticleIndex = CellListIndices[CurrentParticleIndex];
						}
					}

					cerr << "displaced CellList with xDisplacement=  " << P->xDisplacement << endl;
					for (int i = 0; i < CellListHeadDisplacedCells.size(); i++){
						cerr << "Cell " << i << ": " << endl;
						int CurrentParticleIndex = CellListHeadDisplacedCells[i];
						while (CurrentParticleIndex >= 0){
							cerr << CurrentParticleIndex << ": (" << P->Positions[2*CurrentParticleIndex] << "," << P->Positions[2*CurrentParticleIndex+1] << ") : " << "(" << static_cast<int>(static_cast<double>(NumberOfSubdivisions)*P->Positions[2*CurrentParticleIndex]) << "," << static_cast<int>(static_cast<double>(NumberOfSubdivisions)*P->Positions[2*CurrentParticleIndex+1]) << ")" << endl;
							CurrentParticleIndex = CellListIndicesDisplacedCells[CurrentParticleIndex];
						}
					}
				}

				double computePairwiseMagnitudeOfForce(double DistanceSquared) const { //intentionally off by a factor r, so we can just multiply with r-vector to get the correct result
					double InverseDistanceSquared = 1.0/DistanceSquared;
					double InverseDistanceToThePowerOfSix = InverseDistanceSquared * InverseDistanceSquared * InverseDistanceSquared;
					return (48.0*InverseDistanceToThePowerOfSix*InverseDistanceToThePowerOfSix*InverseDistanceSquared-24.0*InverseDistanceToThePowerOfSix*InverseDistanceSquared-4.0*POTENTIAL_CONSTANT_2*sqrt(InverseDistanceSquared));
				}

				void computeForceThroughyEdge(int CurrentParticleIndex, int OtherParticleIndex, double xOffset, double yOffset, int xEdgeIndex, int yEdgeIndex, vector<double>& StressOfEdgesInyDirection) {
					int EdgeIndex = xEdgeIndex+NumberOfSubdivisions*yEdgeIndex;
					double xPositionOfEdge = static_cast<double>(xEdgeIndex+1)*DimensionlessEdgeLength;
					double LoweryOfEdge = static_cast<double>(yEdgeIndex)*DimensionlessEdgeLength;
					double UpperyOfEdge = static_cast<double>(yEdgeIndex+1)*DimensionlessEdgeLength;

					double x0 = P->Positions[DIMENSION*CurrentParticleIndex];
					double y0 = P->Positions[DIMENSION*CurrentParticleIndex+1];
					double x1 = P->Positions[DIMENSION*OtherParticleIndex];
					double y1 = P->Positions[DIMENSION*OtherParticleIndex+1];
					double Deltax = x1 - x0;
					double Deltay = y1 - y0;
					if (Deltay > 0.5){
						Deltay -= 1.0;
						Deltax -= P->xDisplacement;
					}
					else if (Deltay <= -0.5){
						Deltay += 1.0;
						Deltax += P->xDisplacement;
					}
					while (Deltax > 0.5){
						Deltax -= 1.0;
					}
					while (Deltax <= -0.5){
						Deltax += 1.0;
					}
					if (Deltax > 0.0){
						double DistanceSquared = (Deltax * Deltax + Deltay * Deltay)*P->BoxLengthSquared;
						if (DistanceSquared < CUTOFF_SQUARED){
							double Shiftedx0 = x0 + xOffset;
							if (Shiftedx0 < 0.0){
								Shiftedx0 += 1.0;
							}
							else if (Shiftedx0 >= 1.0){
								Shiftedx0 -= 1.0;
							}
							double yIntersect = Deltay/Deltax*(xPositionOfEdge-Shiftedx0)+y0+yOffset;
							if (LoweryOfEdge <= yIntersect && UpperyOfEdge > yIntersect){ // the force intersects the edge in question, compute the force of the 2 particles and add it to the total force through the edge
								double InteractionStrength = (P->ParticleTypes[CurrentParticleIndex] == P->ParticleTypes[OtherParticleIndex] ? AA_INTERACTION_STRENGTH : AB_INTERACTION_STRENGTH);
								double MagnitudeOfForce = computePairwiseMagnitudeOfForce(DistanceSquared);
								StressOfEdgesInyDirection[EdgeIndex*DIMENSION] += InteractionStrength*MagnitudeOfForce*P->BoxLength*Deltax;
								StressOfEdgesInyDirection[EdgeIndex*DIMENSION+1] += InteractionStrength*MagnitudeOfForce*P->BoxLength*Deltay;
								NumberOfForceValuesInyDirection[EdgeIndex]++;
							}
						}
					}
				}

				void computeForcesThroughyEdge(int xEdgeIndex, int yEdgeIndex, int yIndexOtherCell, int CurrentParticleIndex, double xOffset, double yOffset, vector<double>& StressOfEdgesInyDirection) {
					int xIndexOtherCell = xEdgeIndex + 1 < NumberOfSubdivisions ? xEdgeIndex + 1 : 0;
					if (yIndexOtherCell >= NumberOfSubdivisions){
						int OtherParticleIndex = CellListHeadDisplacedCells[NumberOfSubdivisions+xIndexOtherCell];
						while (OtherParticleIndex >= 0){
							computeForceThroughyEdge(CurrentParticleIndex, OtherParticleIndex, xOffset, yOffset, xEdgeIndex, yEdgeIndex, StressOfEdgesInyDirection);
							OtherParticleIndex = CellListIndicesDisplacedCells[OtherParticleIndex];
						}
					}
					else if (yIndexOtherCell < 0){
						int OtherParticleIndex = CellListHeadDisplacedCells[xIndexOtherCell];
						while (OtherParticleIndex >= 0){
							computeForceThroughyEdge(CurrentParticleIndex, OtherParticleIndex, xOffset, yOffset, xEdgeIndex, yEdgeIndex, StressOfEdgesInyDirection);
							OtherParticleIndex = CellListIndicesDisplacedCells[OtherParticleIndex];
						}
					}
					else {
						int OtherCell = xIndexOtherCell+NumberOfSubdivisions*yIndexOtherCell;
						int OtherParticleIndex = CellListHead[OtherCell];
						while (OtherParticleIndex >= 0){
							computeForceThroughyEdge(CurrentParticleIndex, OtherParticleIndex, xOffset, yOffset, xEdgeIndex, yEdgeIndex, StressOfEdgesInyDirection);
							OtherParticleIndex = CellListIndices[OtherParticleIndex];
						}
					}
				}

				void computeForceThroughyEdge(int xEdgeIndex, int yEdgeIndex, vector<double>& StressOfEdgesInyDirection){
					int CurrentCell = xEdgeIndex+NumberOfSubdivisions*yEdgeIndex;
					int CurrentParticleIndex = CellListHead[CurrentCell];
					while (CurrentParticleIndex >= 0){
						for (int yOffset = -1; yOffset < 2; yOffset++){
							computeForcesThroughyEdge(xEdgeIndex, yEdgeIndex, yEdgeIndex+yOffset, CurrentParticleIndex, 0.0, 0.0, StressOfEdgesInyDirection);
						}
						CurrentParticleIndex = CellListIndices[CurrentParticleIndex];
					}

					double xOffset = 0.0;
					double yOffset = 0.0;
					int yOtherCellIndex = yEdgeIndex+1;
					int* CellListIndicesPointer;
					if (yOtherCellIndex >= NumberOfSubdivisions){
						CurrentParticleIndex = CellListHeadDisplacedCells[xEdgeIndex+NumberOfSubdivisions];
						CellListIndicesPointer = CellListIndicesDisplacedCells;
						xOffset = P->xDisplacement;
						yOffset = 1.0;
					}
					else {
						CurrentParticleIndex = CellListHead[xEdgeIndex+NumberOfSubdivisions*yOtherCellIndex];
						CellListIndicesPointer = CellListIndices;
					}
					while (CurrentParticleIndex >= 0){
						for (int yIndexOffset = -2; yIndexOffset < 0; yIndexOffset++){
							computeForcesThroughyEdge(xEdgeIndex, yEdgeIndex, yOtherCellIndex+yIndexOffset, CurrentParticleIndex, xOffset, yOffset, StressOfEdgesInyDirection);
						}
						CurrentParticleIndex = CellListIndicesPointer[CurrentParticleIndex];
					}

					xOffset = 0.0;
					yOffset = 0.0;
					yOtherCellIndex = yEdgeIndex-1;
					if (yOtherCellIndex < 0){
						CurrentParticleIndex = CellListHeadDisplacedCells[xEdgeIndex];
						CellListIndicesPointer = CellListIndicesDisplacedCells;
						xOffset = -P->xDisplacement;
						yOffset = -1.0;
					}
					else {
						CurrentParticleIndex = CellListHead[xEdgeIndex+NumberOfSubdivisions*yOtherCellIndex];
						CellListIndicesPointer = CellListIndices;
					}
					while (CurrentParticleIndex >= 0){
						for (int yIndexOffset = 1; yIndexOffset < 3; yIndexOffset++){
							computeForcesThroughyEdge(xEdgeIndex, yEdgeIndex, yOtherCellIndex+yIndexOffset, CurrentParticleIndex, xOffset, yOffset, StressOfEdgesInyDirection);
						}
						CurrentParticleIndex = CellListIndicesPointer[CurrentParticleIndex];
					}
				}

				void computeForceThroughxEdge(int CurrentParticleIndex, int OtherParticleIndex, double xOffset, int xEdgeIndex, int yEdgeIndex, vector<double>& StressOfEdgesInxDirection) {
					int EdgeIndex = xEdgeIndex+NumberOfSubdivisions*yEdgeIndex;
					double yPositionOfEdge = static_cast<double>(yEdgeIndex+1)*DimensionlessEdgeLength;
					double LowerxOfEdge = static_cast<double>(xEdgeIndex)*DimensionlessEdgeLength;
					double UpperxOfEdge = static_cast<double>(xEdgeIndex+1)*DimensionlessEdgeLength;

					double x0 = P->Positions[DIMENSION*CurrentParticleIndex];
					double y0 = P->Positions[DIMENSION*CurrentParticleIndex+1];
					double x1 = P->Positions[DIMENSION*OtherParticleIndex];
					double y1 = P->Positions[DIMENSION*OtherParticleIndex+1];
					double Deltax = x1 - x0;
					double Deltay = y1 - y0;
					if (Deltay > 0.5){
						Deltay -= 1.0;
						Deltax -= P->xDisplacement;
					}
					else if (Deltay <= -0.5){
						Deltay += 1.0;
						Deltax += P->xDisplacement;
					}
					while (Deltax > 0.5){
						Deltax -= 1.0;
					}
					while (Deltax <= -0.5){
						Deltax += 1.0;
					}
					if (Deltay > 0.0){
						double DistanceSquared = (Deltax * Deltax + Deltay * Deltay)*P->BoxLengthSquared;
						if (DistanceSquared < CUTOFF_SQUARED){
							double xIntersect = Deltax/Deltay*(yPositionOfEdge-y0)+x0+xOffset;
							if (LowerxOfEdge <= xIntersect && UpperxOfEdge > xIntersect){ // the force intersects the edge in question, compute the force of the 2 particles and add it to the total force through the edge
								double InteractionStrength = (P->ParticleTypes[CurrentParticleIndex] == P->ParticleTypes[OtherParticleIndex] ? AA_INTERACTION_STRENGTH : AB_INTERACTION_STRENGTH);
								double MagnitudeOfForce = computePairwiseMagnitudeOfForce(DistanceSquared);
								StressOfEdgesInxDirection[EdgeIndex*DIMENSION] += InteractionStrength*MagnitudeOfForce*P->BoxLength*Deltay;
								StressOfEdgesInxDirection[EdgeIndex*DIMENSION+1] += InteractionStrength*MagnitudeOfForce*P->BoxLength*Deltax;
								NumberOfForceValuesInxDirection[EdgeIndex]++;
							}
						}
					}
				}

				void computeForcesThroughxEdge(int xEdgeIndex, int yEdgeIndex, int xIndexOtherCell, int CurrentParticleIndex, double xOffset, vector<double>& StressOfEdgesInxDirection) {
					if (xIndexOtherCell < 0){
						xIndexOtherCell += NumberOfSubdivisions;
					}
					else if (xIndexOtherCell >= NumberOfSubdivisions){
						xIndexOtherCell -= NumberOfSubdivisions;
					}
					int yIndexOtherCell = yEdgeIndex + 1;
					if (yIndexOtherCell >= NumberOfSubdivisions){
						int OtherParticleIndex = CellListHeadDisplacedCells[NumberOfSubdivisions+xIndexOtherCell];
						while (OtherParticleIndex >= 0){
							computeForceThroughxEdge(CurrentParticleIndex, OtherParticleIndex, xOffset, xEdgeIndex, yEdgeIndex, StressOfEdgesInxDirection);
							OtherParticleIndex = CellListIndicesDisplacedCells[OtherParticleIndex];
						}
					}
					else {
						int OtherCell = xIndexOtherCell+NumberOfSubdivisions*yIndexOtherCell;
						int OtherParticleIndex = CellListHead[OtherCell];
						while (OtherParticleIndex >= 0){
							computeForceThroughxEdge(CurrentParticleIndex, OtherParticleIndex, xOffset, xEdgeIndex, yEdgeIndex, StressOfEdgesInxDirection);
							OtherParticleIndex = CellListIndices[OtherParticleIndex];
						}
					}
				}

				void computeForceThroughxEdge(int xEdgeIndex, int yEdgeIndex, vector<double>& StressOfEdgesInxDirection){
					int CurrentParticleIndex = CellListHead[xEdgeIndex+NumberOfSubdivisions*yEdgeIndex];
					while (CurrentParticleIndex >= 0){
						for (int xIndexOffset = -1; xIndexOffset < 2; xIndexOffset++){
							computeForcesThroughxEdge(xEdgeIndex, yEdgeIndex, xEdgeIndex+xIndexOffset, CurrentParticleIndex, 0.0, StressOfEdgesInxDirection);
						}
						CurrentParticleIndex = CellListIndices[CurrentParticleIndex];
					}

					double xOffset = 0.0;
					int xCellIndex = xEdgeIndex+1;
					if (xCellIndex >= NumberOfSubdivisions){
						xCellIndex = 0;
						xOffset = 1.0;
					}
					CurrentParticleIndex = CellListHead[xCellIndex+NumberOfSubdivisions*yEdgeIndex];
					while (CurrentParticleIndex >= 0){
						for (int xIndexOffset = -2; xIndexOffset < 0; xIndexOffset++){
							computeForcesThroughxEdge(xEdgeIndex, yEdgeIndex, xCellIndex+xIndexOffset, CurrentParticleIndex, xOffset, StressOfEdgesInxDirection);
						}
						CurrentParticleIndex = CellListIndices[CurrentParticleIndex];
					}

					xOffset = 0.0;
					xCellIndex = xEdgeIndex-1;
					if (xCellIndex < 0){
						xCellIndex = NumberOfSubdivisions - 1;
						xOffset = -1.0;
					}
					CurrentParticleIndex = CellListHead[xCellIndex+NumberOfSubdivisions*yEdgeIndex];
					while (CurrentParticleIndex >= 0){
						for (int xIndexOffset = 1; xIndexOffset < 3; xIndexOffset++){
							computeForcesThroughxEdge(xEdgeIndex, yEdgeIndex, xCellIndex+xIndexOffset, CurrentParticleIndex, xOffset, StressOfEdgesInxDirection);
						}
						CurrentParticleIndex = CellListIndices[CurrentParticleIndex];
					}
				}

			public:

				void initialize(Particles* OuterP, int IntendedNumberOfSubdivisions){
					P = OuterP;
					NumberOfAverages = 0;
					NumberOfSubdivisions = P->BoxLength/static_cast<double>(IntendedNumberOfSubdivisions) > CUTOFF ? IntendedNumberOfSubdivisions : static_cast<int>(P->BoxLength / CUTOFF);
					EdgeLength = P->BoxLength/static_cast<double>(NumberOfSubdivisions);
					DimensionlessEdgeLength = 1.0/static_cast<double>(NumberOfSubdivisions);
					StressOfEdgesInyDirection.assign(NumberOfSubdivisions*NumberOfSubdivisions*DIMENSION,0.0);
					StressOfEdgesInxDirection.assign(NumberOfSubdivisions*NumberOfSubdivisions*DIMENSION,0.0);
					NumberOfForceValuesInxDirection.assign(NumberOfSubdivisions*NumberOfSubdivisions,0);
					NumberOfForceValuesInyDirection.assign(NumberOfSubdivisions*NumberOfSubdivisions,0);
				}

				void computeStresses() {
					NumberOfAverages++;
					buildCellLists();
					for (int xEdgeIndex = 0; xEdgeIndex < NumberOfSubdivisions; xEdgeIndex++){
						for (int yEdgeIndex = 0; yEdgeIndex < NumberOfSubdivisions; yEdgeIndex++){
							computeForceThroughyEdge(xEdgeIndex,yEdgeIndex,StressOfEdgesInyDirection);
							computeForceThroughxEdge(xEdgeIndex,yEdgeIndex,StressOfEdgesInxDirection);
						}
					}
				}

				void writeAverageStresses(string FilePath) const {
					ofstream FileStreamToWrite;
					FileStreamToWrite.open(FilePath+"_yEdges.dat");
					FileStreamToWrite << "y stresses (i.e. stresses through edges in y direction)" << endl;
					FileStreamToWrite << "xPos\tlowery\tuppery\tnormal_stress\ttangential_stress" << endl;
					for (int yEdgeIndex = 0; yEdgeIndex < NumberOfSubdivisions; yEdgeIndex++){
						double LoweryOfEdge = static_cast<double>(yEdgeIndex)*DimensionlessEdgeLength;
						double UpperyOfEdge = static_cast<double>(yEdgeIndex+1)*DimensionlessEdgeLength;
						for (int xEdgeIndex = 0; xEdgeIndex < NumberOfSubdivisions; xEdgeIndex++){
							double xPositionOfEdge = static_cast<double>(xEdgeIndex+1)*DimensionlessEdgeLength;
							FileStreamToWrite << xPositionOfEdge << "\t" << LoweryOfEdge << "\t" << UpperyOfEdge << "\t" << StressOfEdgesInyDirection[DIMENSION*(xEdgeIndex+NumberOfSubdivisions*yEdgeIndex)]/(static_cast<double>(NumberOfAverages)*EdgeLength) << "\t" << StressOfEdgesInyDirection[DIMENSION*(xEdgeIndex+NumberOfSubdivisions*yEdgeIndex)+1]/(static_cast<double>(NumberOfAverages)*EdgeLength) << endl;
						}
					}
					FileStreamToWrite.close();
					FileStreamToWrite.open(FilePath+"_xEdges.dat");
					FileStreamToWrite << "x stresses (i.e. stresses through edges in x direction)" << endl;
					FileStreamToWrite << "yPos\tlowerx\tupperx\tnormal_stress\ttangential_stress" << endl;
					for (int xEdgeIndex = 0; xEdgeIndex < NumberOfSubdivisions; xEdgeIndex++){
						double LowerxOfEdge = static_cast<double>(xEdgeIndex)*DimensionlessEdgeLength;
						double UpperxOfEdge = static_cast<double>(xEdgeIndex+1)*DimensionlessEdgeLength;
						for (int yEdgeIndex = 0; yEdgeIndex < NumberOfSubdivisions; yEdgeIndex++){
							double yPositionOfEdge = static_cast<double>(yEdgeIndex+1)*DimensionlessEdgeLength;
							FileStreamToWrite << yPositionOfEdge << "\t" << LowerxOfEdge << "\t" << UpperxOfEdge << "\t" << StressOfEdgesInxDirection[DIMENSION*(xEdgeIndex+NumberOfSubdivisions*yEdgeIndex)]/(static_cast<double>(NumberOfAverages)*EdgeLength) << "\t" << StressOfEdgesInxDirection[DIMENSION*(xEdgeIndex+NumberOfSubdivisions*yEdgeIndex)+1]/(static_cast<double>(NumberOfAverages)*EdgeLength) << endl;
						}
					}
					FileStreamToWrite.close();
					cerr << "Avg number of force values per edge:\n";
					cerr << "xEdges:\n";
					double TotalAverage = 0.0;
					for (int yEdgeIndex = 0; yEdgeIndex < NumberOfSubdivisions; yEdgeIndex++){
						for (int xEdgeIndex = 0; xEdgeIndex < NumberOfSubdivisions; xEdgeIndex++){
							cerr << static_cast<double>(NumberOfForceValuesInxDirection[xEdgeIndex+NumberOfSubdivisions*yEdgeIndex])/(NumberOfAverages) << " ";
							TotalAverage += static_cast<double>(NumberOfForceValuesInxDirection[xEdgeIndex+NumberOfSubdivisions*yEdgeIndex]);
						}
						cerr << endl;
					}
					cerr << "Total avg: " << TotalAverage/(NumberOfAverages*NumberOfSubdivisions*NumberOfSubdivisions) << endl;
					cerr << "yEdges:\n";
					TotalAverage = 0.0;
					for (int yEdgeIndex = 0; yEdgeIndex < NumberOfSubdivisions; yEdgeIndex++){
						for (int xEdgeIndex = 0; xEdgeIndex < NumberOfSubdivisions; xEdgeIndex++){
							cerr << static_cast<double>(NumberOfForceValuesInyDirection[xEdgeIndex+NumberOfSubdivisions*yEdgeIndex])/(NumberOfAverages) << " ";
							TotalAverage += static_cast<double>(NumberOfForceValuesInyDirection[xEdgeIndex+NumberOfSubdivisions*yEdgeIndex]);
						}
						cerr << endl;
					}
					cerr << "Total avg: " << TotalAverage/(NumberOfAverages*NumberOfSubdivisions*NumberOfSubdivisions) << endl;
					cerr << endl;
				}

		};

		StressComputator SC;

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

			/*cerr << "xDisplacement = " << xDisplacement << endl;
			cerr << "CellNeighborsHead.size() = " << CellNeighborsHead.size() << endl;
			for (int i = 0; i < CellNeighborsHead.size(); i += 2){
				cerr << i/2 << ": ";
				for (int j = 0; j < CellNeighborsHead[i+1]; j++){
					cerr << CellNeighborsIndices[CellNeighborsHead[i] + j] << ",";
				}
				cerr << endl;
			}*/

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
			/*cerr << endl << "Size of VerletIndicesOfNeighbors: " << VerletIndicesOfNeighbors.size() << endl;
			cerr << endl << "Resulting neighbors: " << endl;
			for (int i = 0; i < TOTAL_NUMBER_OF_PARTICLES; i++){
				cerr << i << ": ";
				for (int j = 0; j < VerletListHead[2*i+1]; j++){
					cerr << VerletIndicesOfNeighbors[VerletListHead[2*i]+j] << ",";
				}
				cerr << endl;
			}
			cerr << endl;*/
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

		double getxDisplacement() const {
			return xDisplacement;
		}

		void applyShear(double ShearRate, double* ChangeInCoordinates) {
			for (int ParticleID = 0; ParticleID < TOTAL_NUMBER_OF_PARTICLES; ParticleID++) {
				double& xPos = Positions[DIMENSION*ParticleID];
				double& yPos = Positions[DIMENSION*ParticleID+1];
				xPos += ShearRate*(yPos-0.5);
				ChangeInCoordinates[DIMENSION*ParticleID] += ShearRate*(yPos-0.5);
				if (xPos < 0.0) {
					xPos += 1.0;
				}
				else if (xPos >= 1.0) {
					xPos -= 1.0;
				}
			}
			xDisplacement += ShearRate;
			xDisplacement -= static_cast<double>(static_cast<int>(xDisplacement));
			buildVerletList();
		}

		void writeAverageStresses(string FilePath) const {
			SC.writeAverageStresses(FilePath);
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

		void initialize(int InitialNumberOfAParticles, double InitialDensity, double InitialxDisplacement, int NumberOfIntendedStressSubdivisions){
			updateBoxParametersWithDensity(InitialDensity);
			SC.initialize(this, NumberOfIntendedStressSubdivisions);
			int InitialNumberOfBParticles = TOTAL_NUMBER_OF_PARTICLES - InitialNumberOfAParticles;
			double FractionOfAParticles = static_cast<double>(InitialNumberOfAParticles)/static_cast<double>(TOTAL_NUMBER_OF_PARTICLES);
			realRNG RNG;
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

			xDisplacement = InitialxDisplacement;

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
			buildVerletList();
		}

		void readInParticleState(string FileNameToReadIn, int StateNumber, double Density) {
			updateBoxParametersWithDensity(Density);
			xDisplacement = 0.0;

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
			buildVerletList();
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
			double UpdatedCoordinates [DIMENSION]{Positions[DIMENSION*ParticleIndex]+Delta[0],Positions[DIMENSION*ParticleIndex+1]+Delta[1]};
			if (UpdatedCoordinates[1] < 0.0){
				UpdatedCoordinates[0] += xDisplacement;
				UpdatedCoordinates[1] += 1.0;
			}
			else if (UpdatedCoordinates[1] >= 1.0){
				UpdatedCoordinates[0] -= xDisplacement;
				UpdatedCoordinates[1] -= 1.0;
			}
			while (UpdatedCoordinates[0] < 0.0){
				UpdatedCoordinates[0] += 1.0;
			}
			while (UpdatedCoordinates[0] >= 1.0){
				UpdatedCoordinates[0] -= 1.0;
			}
			double PotEnergyChange = 0.0;
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

		void computeStresses() {
			SC.computeStresses();
		}

		double computeChangeInPotentialEnergyByChangingVolume(double VolumeChange) {
			double PotEnergyBefore = computePotentialEnergy();
			updateBoxParametersWithVolumeChange(VolumeChange);
			buildVerletList();
			return computePotentialEnergy() - PotEnergyBefore;
		}

		void updatePosition(int ParticleIndex, const double* Deltas){
			Positions[DIMENSION * ParticleIndex] += Deltas[0];
			Positions[DIMENSION * ParticleIndex + 1] += Deltas[1];
			if (Positions[DIMENSION * ParticleIndex + 1] < 0.0){
				Positions[DIMENSION * ParticleIndex] += xDisplacement;
				Positions[DIMENSION * ParticleIndex + 1] += 1.0;
			}
			else if (Positions[DIMENSION * ParticleIndex + 1] >= 1.0){
				Positions[DIMENSION * ParticleIndex] -= xDisplacement;
				Positions[DIMENSION * ParticleIndex + 1] -= 1.0;
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
				exit(EXIT_FAILURE);
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
			buildVerletList();
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
	OStream << "X       Y       Type | #AParticles:  " << fixed << setprecision(numeric_limits<long double>::digits10+1) << State.getNumberOfAParticles() << "| #BParticles: " << State.getNumberOfBParticles() << "| BoxLength: " << State.getBoxLength() << "| xDisplacement: " << State.getxDisplacement() << endl;
	for (int i = 0; i < TOTAL_NUMBER_OF_PARTICLES; i++){
		OStream << State.getPosition(i,0) << "\t" << State.getPosition(i,1) << "\t" << State.getParticleType(i) <<  endl;
	}
	return OStream;
}

#endif //PARTICLES_INCLUDED
