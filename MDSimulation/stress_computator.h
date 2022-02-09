#ifndef STRESS_COMPUTATOR_INCLUDED
#define STRESS_COMPUTATOR_INCLUDED

#include <iostream>
#include <vector>
#include <string>

#include "particles.h"

using namespace std;

class StressComputator {
	private:
		const Particles& P;
		vector<double> StressOfEdgesInyDirection;
		vector<int> NumberOfForceValuesInyDirection;
		vector<double> StressOfEdgesInxDirection;
		vector<int> NumberOfForceValuesInxDirection;

		int NumberOfAverages;
		int NumberOfConfigurations;

		vector<double> seriesOfAvgShearStresses;

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
					CurrentCellIndex += static_cast<int>(static_cast<double>(NumberOfSubdivisions)*P.getPosition(ParticleIndex,j))*IndexFactor;
					IndexFactor *= NumberOfSubdivisions;
				}
				if (static_cast<int>(static_cast<double>(NumberOfSubdivisions)*P.getPosition(ParticleIndex,1)) == 0) {
					double DisplacedxPosition = P.getPosition(ParticleIndex,0)+P.getxDisplacement();
					if (DisplacedxPosition >= 1.0){
						DisplacedxPosition -= 1.0;
					}
					int DisplacedxIndex = static_cast<int>(static_cast<double>(NumberOfSubdivisions)*DisplacedxPosition);
					CellListIndicesDisplacedCells[ParticleIndex] = CellListHeadDisplacedCells[DisplacedxIndex+NumberOfSubdivisions];
					CellListHeadDisplacedCells[DisplacedxIndex+NumberOfSubdivisions] = ParticleIndex;
				}
				else if (static_cast<int>(static_cast<double>(NumberOfSubdivisions)*P.getPosition(ParticleIndex,1)) >= NumberOfSubdivisions-1) {
					double DisplacedxPosition = P.getPosition(ParticleIndex,0) - P.getxDisplacement();
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
					cerr << CurrentParticleIndex << ": (" << P.getPosition(CurrentParticleIndex,0) << "," << P.getPosition(CurrentParticleIndex,1) << ") : " << "(" << static_cast<int>(static_cast<double>(NumberOfSubdivisions)*P.getPosition(CurrentParticleIndex,0)) << "," << static_cast<int>(static_cast<double>(NumberOfSubdivisions)*P.getPosition(CurrentParticleIndex,1)) << ")" << endl;
					CurrentParticleIndex = CellListIndices[CurrentParticleIndex];
				}
			}

			cerr << "displaced CellList with xDisplacement=  " << P.getxDisplacement() << endl;
			for (int i = 0; i < CellListHeadDisplacedCells.size(); i++){
				cerr << "Cell " << i << ": " << endl;
				int CurrentParticleIndex = CellListHeadDisplacedCells[i];
				while (CurrentParticleIndex >= 0){
					cerr << CurrentParticleIndex << ": (" << P.getPosition(CurrentParticleIndex,0) << "," << P.getPosition(CurrentParticleIndex,1) << ") : " << "(" << static_cast<int>(static_cast<double>(NumberOfSubdivisions)*P.getPosition(CurrentParticleIndex,0)) << "," << static_cast<int>(static_cast<double>(NumberOfSubdivisions)*P.getPosition(CurrentParticleIndex,1)) << ")" << endl;
					CurrentParticleIndex = CellListIndicesDisplacedCells[CurrentParticleIndex];
				}
			}
		}

		void computeForceThroughyEdge(int CurrentParticleIndex, int OtherParticleIndex, double xOffset, double yOffset, int xEdgeIndex, int yEdgeIndex, vector<double>& StressOfEdgesInyDirection) {
			int EdgeIndex = xEdgeIndex+NumberOfSubdivisions*yEdgeIndex;
			double xPositionOfEdge = static_cast<double>(xEdgeIndex+1)*DimensionlessEdgeLength;
			double LoweryOfEdge = static_cast<double>(yEdgeIndex)*DimensionlessEdgeLength;
			double UpperyOfEdge = static_cast<double>(yEdgeIndex+1)*DimensionlessEdgeLength;

			double x0 = P.getPosition(CurrentParticleIndex,0);
			double y0 = P.getPosition(CurrentParticleIndex,1);
			double x1 = P.getPosition(OtherParticleIndex,0);
			double y1 = P.getPosition(OtherParticleIndex,1);
			double Deltax = x1 - x0;
			double Deltay = y1 - y0;
			if (Deltay > 0.5){
				Deltay -= 1.0;
				Deltax -= P.getxDisplacement();
			}
			else if (Deltay <= -0.5){
				Deltay += 1.0;
				Deltax += P.getxDisplacement();
			}
			while (Deltax > 0.5){
				Deltax -= 1.0;
			}
			while (Deltax <= -0.5){
				Deltax += 1.0;
			}
			if (Deltax > 0.0){
				double DistanceSquared = (Deltax * Deltax + Deltay * Deltay)*P.getBoxLengthSquared();
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
						double InteractionStrength = (P.getParticleType(CurrentParticleIndex) == P.getParticleType(OtherParticleIndex) ? AA_INTERACTION_STRENGTH : AB_INTERACTION_STRENGTH);
						double MagnitudeOfForce = P.computePairwiseMagnitudeOfForce(DistanceSquared);
						StressOfEdgesInyDirection[EdgeIndex*DIMENSION] += InteractionStrength*MagnitudeOfForce*P.getBoxLength()*Deltax;
						StressOfEdgesInyDirection[EdgeIndex*DIMENSION+1] += InteractionStrength*MagnitudeOfForce*P.getBoxLength()*Deltay;
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
				xOffset = P.getxDisplacement();
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
				xOffset = -P.getxDisplacement();
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

			double x0 = P.getPosition(CurrentParticleIndex,0);
			double y0 = P.getPosition(CurrentParticleIndex,1);
			double x1 = P.getPosition(OtherParticleIndex,0);
			double y1 = P.getPosition(OtherParticleIndex,1);
			double Deltax = x1 - x0;
			double Deltay = y1 - y0;
			if (Deltay > 0.5){
				Deltay -= 1.0;
				Deltax -= P.getxDisplacement();
			}
			else if (Deltay <= -0.5){
				Deltay += 1.0;
				Deltax += P.getxDisplacement();
			}
			while (Deltax > 0.5){
				Deltax -= 1.0;
			}
			while (Deltax <= -0.5){
				Deltax += 1.0;
			}
			if (Deltay > 0.0){
				double DistanceSquared = (Deltax * Deltax + Deltay * Deltay)*P.getBoxLengthSquared();
				if (DistanceSquared < CUTOFF_SQUARED){
					double xIntersect = Deltax/Deltay*(yPositionOfEdge-y0)+x0+xOffset;
					if (LowerxOfEdge <= xIntersect && UpperxOfEdge > xIntersect){ // the force intersects the edge in question, compute the force of the 2 particles and add it to the total force through the edge
						double InteractionStrength = (P.getParticleType(CurrentParticleIndex) == P.getParticleType(OtherParticleIndex) ? AA_INTERACTION_STRENGTH : AB_INTERACTION_STRENGTH);
						double MagnitudeOfForce = P.computePairwiseMagnitudeOfForce(DistanceSquared);
						StressOfEdgesInxDirection[EdgeIndex*DIMENSION] += InteractionStrength*MagnitudeOfForce*P.getBoxLength()*Deltay;
						StressOfEdgesInxDirection[EdgeIndex*DIMENSION+1] += InteractionStrength*MagnitudeOfForce*P.getBoxLength()*Deltax;

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

		double computeTotalShearStress() const {
			double TotShearStress = 0.0;
			#pragma omp parallel num_threads(NUMBER_OF_THREADS)
			{
				double ShearStress = 0.0;
				#pragma omp for
				for (int ParticleIndex = 0; ParticleIndex < TOTAL_NUMBER_OF_PARTICLES; ParticleIndex++){
					for (int i = 0; i < P.VerletListHead[2*ParticleIndex+1]; i++){
						int OtherParticleIndex = P.VerletIndicesOfNeighbors[P.VerletListHead[2*ParticleIndex]+i];
						if (OtherParticleIndex < ParticleIndex){
							double x0 = P.getPosition(ParticleIndex,0);
							double y0 = P.getPosition(ParticleIndex,1);
							double x1 = P.getPosition(OtherParticleIndex,0);
							double y1 = P.getPosition(OtherParticleIndex,1);
							double Deltax = x1 - x0;
							double Deltay = y1 - y0;
							if (Deltay > 0.5){
								Deltay -= 1.0;
								Deltax -= P.getxDisplacement();
							}
							else if (Deltay <= -0.5){
								Deltay += 1.0;
								Deltax += P.getxDisplacement();
							}
							while (Deltax > 0.5){
								Deltax -= 1.0;
							}
							while (Deltax <= -0.5){
								Deltax += 1.0;
							}
							if (Deltay < 0.0){
								Deltax *= -1.0;
								Deltay *= -1.0;
							}
							double DistanceSquared = (Deltax * Deltax + Deltay * Deltay)*P.getBoxLengthSquared();
							if (DistanceSquared < CUTOFF_SQUARED){
								double InteractionStrength = (P.getParticleType(ParticleIndex) == P.getParticleType(OtherParticleIndex) ? AA_INTERACTION_STRENGTH : AB_INTERACTION_STRENGTH);
								double MagnitudeOfForce = P.computePairwiseMagnitudeOfForce(DistanceSquared);
								ShearStress += InteractionStrength*MagnitudeOfForce*Deltax*Deltay;
							}
						}
					}
				}
				#pragma omp atomic
				TotShearStress += ShearStress;
			}
			return TotShearStress;
		}

	public:

		StressComputator(const Particles& P, int IntendedNumberOfSubdivisions, int NumberOfTimeSteps, int NumberOfAverageConfigurations):
			P(P),
			NumberOfAverages(0),
			NumberOfConfigurations(NumberOfAverageConfigurations),
			seriesOfAvgShearStresses(NumberOfTimeSteps,0.0)
		{
			NumberOfSubdivisions = P.getBoxLength()/static_cast<double>(IntendedNumberOfSubdivisions) > CUTOFF ? IntendedNumberOfSubdivisions : static_cast<int>(P.getBoxLength() / CUTOFF);
			EdgeLength = P.getBoxLength()/static_cast<double>(NumberOfSubdivisions);
			DimensionlessEdgeLength = 1.0/static_cast<double>(NumberOfSubdivisions);
			StressOfEdgesInyDirection.assign(NumberOfSubdivisions*NumberOfSubdivisions*DIMENSION,0.0);
			StressOfEdgesInxDirection.assign(NumberOfSubdivisions*NumberOfSubdivisions*DIMENSION,0.0);
			NumberOfForceValuesInxDirection.assign(NumberOfSubdivisions*NumberOfSubdivisions,0);
			NumberOfForceValuesInyDirection.assign(NumberOfSubdivisions*NumberOfSubdivisions,0);
		}

		void computeStresses(int TimeStep) {
			NumberOfAverages++;
			buildCellLists();

			for (int xEdgeIndex = 0; xEdgeIndex < NumberOfSubdivisions; xEdgeIndex++){
				for (int yEdgeIndex = 0; yEdgeIndex < NumberOfSubdivisions; yEdgeIndex++){
					computeForceThroughyEdge(xEdgeIndex, yEdgeIndex, StressOfEdgesInyDirection);
					computeForceThroughxEdge(xEdgeIndex, yEdgeIndex, StressOfEdgesInxDirection);
				}
			}
			double currentTotalShearStress = computeTotalShearStress();
			seriesOfAvgShearStresses[TimeStep] += currentTotalShearStress;
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

			double avgShearStressTotal = 0.0;
			FileStreamToWrite.open(FilePath+"_ShearStressSeries.dat");
			FileStreamToWrite << "avg shear stresses series" << endl;
			for (int i = 0; i < seriesOfAvgShearStresses.size(); i++){
				double avgShearStressThisTimeStep = seriesOfAvgShearStresses[i]/static_cast<double>(NumberOfConfigurations);
				FileStreamToWrite << avgShearStressThisTimeStep << '\n';
				avgShearStressTotal += avgShearStressThisTimeStep;
			}
			FileStreamToWrite.close();

			cerr << fixed << setprecision(5);
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

			cerr << "Avg shear stress total=" << avgShearStressTotal/(static_cast<double>(seriesOfAvgShearStresses.size())) << endl;
		}
};

#endif // STRESS_COMPUTATOR_INCLUDED
