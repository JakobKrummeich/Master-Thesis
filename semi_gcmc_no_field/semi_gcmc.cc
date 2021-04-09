#include <iostream>
#include <cmath>
#include <chrono>
#include <random>
#include <fstream>
#include <vector>
#include <iomanip>

using namespace std;

static const int DIMENSION = 2;
static const int TOTAL_NUMBER_OF_PARTICLES = 500;
static const double AA_INTERACTION_STRENGTH = 1.0;
static const double AB_INTERACTION_STRENGTH = 0.5;
static const double CUTOFF = 2.5;
static const double CUTOFF_SQUARED = CUTOFF * CUTOFF;
static const double INVERSE_CUTOFF = 1.0/CUTOFF;
static const double DENSITY = 0.9;
static const double BOX_LENGTH = sqrt(static_cast<double>(TOTAL_NUMBER_OF_PARTICLES) / DENSITY);
static const double BOX_LENGTH_SQUARED = BOX_LENGTH * BOX_LENGTH;
static const double INVERSE_BOX_LENGTH = 1.0/BOX_LENGTH;
static const double MAXIMUM_DISPLACEMENT = 0.1;
static const double MAX_VERLET_DIST = 1.3*CUTOFF;
static const double MAX_VERLET_DIST_SQUARED = MAX_VERLET_DIST * MAX_VERLET_DIST;
static const double SKINDISTANCE = MAX_VERLET_DIST - CUTOFF;
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
		
		double drawRandomNumber(){
			return UnifRealDist(rng);
		}
		
		double drawRandomNumber(double Min, double Max){
			return (Max - Min) * UnifRealDist(rng) + Min;
		}
};

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

template<typename T, size_t Nm> struct fvec{
	T arr [Nm];
	size_t Ncurr;
	
	fvec():
		Ncurr(0){
	}
	
	void clear(){
		Ncurr = 0;
	}
	
	void push_back(const T& NewEntry){
		arr[Ncurr] = NewEntry;
		Ncurr++;
	}
	
	bool empty() const{
		return (Ncurr == 0);
	}
	
	size_t size() const{
		return Ncurr;
	}
	
	const T& operator[](int ID) const{
		return arr[ID];
	}
	
	T& operator[](int ID){
		return arr[ID];
	}
	
	void erase(int ID){
		Ncurr--;
		arr[ID] = arr[Ncurr];
	}
};

struct Particles {
	private:
	double Positions [DIMENSION*TOTAL_NUMBER_OF_PARTICLES];
	
	ParticleType ParticleTypes [TOTAL_NUMBER_OF_PARTICLES];
	
	vector<int> CellListHead;
	int CellListIndices [TOTAL_NUMBER_OF_PARTICLES];
	int VerletListHead [2*TOTAL_NUMBER_OF_PARTICLES];
	vector<int> VerletIndicesOfNeighbors;

	double ChangeInCoordinates [DIMENSION*TOTAL_NUMBER_OF_PARTICLES];
	double MostTraveledDistances [2];
	int MostTraveledParticleIndex;

	public:
	
	fvec<int, TOTAL_NUMBER_OF_PARTICLES> TypeAParticleIndices;
	fvec<int, TOTAL_NUMBER_OF_PARTICLES> TypeBParticleIndices;
	int NumberOfVerletListBuilds;

	void initialize(double InitialFractionOfAParticles){
		int InitialNumberOfAParticles = static_cast<int>(round(InitialFractionOfAParticles*static_cast<double>(TOTAL_NUMBER_OF_PARTICLES)));
		int InitialNumberOfBParticles = TOTAL_NUMBER_OF_PARTICLES - InitialNumberOfAParticles;

		int NumberOfParticlesInARow(ceil(pow(static_cast<double>(TOTAL_NUMBER_OF_PARTICLES),1.0/static_cast<double>(DIMENSION))));
		double Distance(1.0/static_cast<double>(NumberOfParticlesInARow));
		cerr << "NumberOfParticlesInRow: " << NumberOfParticlesInARow << endl;
		cerr << "Distance: " << Distance*BOX_LENGTH << endl;
		double CurrentPosition [DIMENSION];
		for (int i = 0; i < DIMENSION; i++){
			CurrentPosition[i] = Distance*0.5;
		}
		int NumberOfAParticlesInitialized = 0;
		int NumberOfBParticlesInitialized = 0;
		ParticleType NextTypeToInitialize = ParticleType::A;
		for (int ParticlesInitialized = 0; ParticlesInitialized < TOTAL_NUMBER_OF_PARTICLES; ParticlesInitialized++){
			for (int i = 0; i < DIMENSION; i++){
				Positions[ParticlesInitialized*DIMENSION + i] = CurrentPosition[i];
				ChangeInCoordinates[DIMENSION*ParticlesInitialized+i] = 0.0;
			}
			CurrentPosition[0] += Distance;
			for (int i = 1; i < DIMENSION; i++){
				if (CurrentPosition[i-1] >= 1.0){
					CurrentPosition[i-1] = Distance*0.5;
					CurrentPosition[i] += Distance;
				}
			}
			if (NextTypeToInitialize == ParticleType::A){
				ParticleTypes[ParticlesInitialized] = ParticleType::A;
				NextTypeToInitialize = ParticleType::B;
				TypeAParticleIndices.push_back(ParticlesInitialized);
			}
			else{
				ParticleTypes[ParticlesInitialized] = ParticleType::B;
				NextTypeToInitialize = ParticleType::A;
				TypeBParticleIndices.push_back(ParticlesInitialized);
			}
		}
		for (int i = 0; i < 2; i++){
			MostTraveledDistances[i] = 0.0;
		}
		NumberOfVerletListBuilds = 0;
	}
	
	double getPosition(int ParticleIndex, int Coordinate) const {
		return Positions[DIMENSION*ParticleIndex+Coordinate];
	}

	ParticleType getParticleType(int ParticleIndex) const {
		return ParticleTypes[ParticleIndex];
	}
	
	void switchParticleType(int ParticleID, int IndexInTypeArray, ParticleType TypeBefore) {
		ParticleTypes[ParticleID] = TypeBefore == ParticleType::A ? ParticleType::B : ParticleType::A;
		if (TypeBefore == ParticleType::A){
			TypeAParticleIndices.erase(IndexInTypeArray);
			TypeBParticleIndices.push_back(ParticleID);
		}
		else {
			TypeBParticleIndices.erase(IndexInTypeArray);
			TypeAParticleIndices.push_back(ParticleID);
		}
	}
	
	int getNumberOfAParticles() const {
		return TypeAParticleIndices.size();
	}
	
	int getNumberOfBParticles() const {
		return TypeBParticleIndices.size();
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
	OStream << "#ID\tX       Y       Type | #AParticles:  " << State.getNumberOfAParticles() << "| #BParticles: " << State.getNumberOfBParticles() << "| #Builds: " << State.NumberOfVerletListBuilds << endl;
	for (int i = 0; i < TOTAL_NUMBER_OF_PARTICLES; i++){
		OStream << i << "\t";
		OStream << fixed << setprecision(5) << State.getPosition(i,0) << "\t" << State.getPosition(i,1) << "\t" << State.getParticleType(i) <<  endl;
	}
	return OStream;
}

struct SimulationManager {
	Particles P;
	double Temperature;
	double Beta;
	realRNG RNG;
	double ChemicalPotentialDiff;
	int DistributionNumberOfA [TOTAL_NUMBER_OF_PARTICLES + 1];
	
	void initialize(double InitialFractionOfAParticles, double _Temperature, double _ChemicalPotentialDiff) {
		Temperature = _Temperature;
		Beta = 1.0/Temperature;
		ChemicalPotentialDiff = _ChemicalPotentialDiff;
		P.initialize(InitialFractionOfAParticles);
		P.buildVerletList();
		for (int i = 0; i < TOTAL_NUMBER_OF_PARTICLES + 1; i++){
			DistributionNumberOfA[i] = 0;
		}
	}
	
	void runDisplacementSteps(int StepsPerParticle) {
		for (int i = 0; i < StepsPerParticle; i++){
			for (int j = 0; j < TOTAL_NUMBER_OF_PARTICLES; j++){
				int RandomParticleID = static_cast<int>(RNG.drawRandomNumber()*static_cast<double>(TOTAL_NUMBER_OF_PARTICLES));
				double Deltas [DIMENSION];
				for (int i = 0; i < DIMENSION; i++){
					Deltas[i] = RNG.drawRandomNumber(-MAXIMUM_DISPLACEMENT, MAXIMUM_DISPLACEMENT)*INVERSE_BOX_LENGTH;
				}
				double AcceptanceProbability = exp(-P.computeChangeInPotentialEnergyByMoving(RandomParticleID, Deltas)*Beta);
				if (AcceptanceProbability >= 1.0 || (RNG.drawRandomNumber() < AcceptanceProbability)){
					P.updatePosition(RandomParticleID, Deltas);
				}
			}
		}
	}
	
	void runTypeChanges(int NumberOfTriedChanges) {
		for (int i = 0; i < NumberOfTriedChanges; i++){
			int NumberOfAParticles = P.getNumberOfAParticles();
			int NumberOfBParticles = P.getNumberOfBParticles();
			int RandomParticleIndexInTypeArray;
			int RandomParticleID;
			double ParticleNumbersPrefactor;
			double ChemicalPotentialSign;
			ParticleType TypeOfParticleBefore;
			bool ParticleToSwitchAvailable = false;

			if (RNG.drawRandomNumber() <= 0.5){
				if (NumberOfAParticles > 0){
					ParticleToSwitchAvailable = true;
					RandomParticleIndexInTypeArray = static_cast<int>(RNG.drawRandomNumber()*static_cast<double>(NumberOfAParticles));
					RandomParticleID = P.TypeAParticleIndices[RandomParticleIndexInTypeArray];
					TypeOfParticleBefore = ParticleType::A;
					ParticleNumbersPrefactor = static_cast<double>(NumberOfAParticles)/static_cast<double>(NumberOfBParticles+1);
					ChemicalPotentialSign = 1.0;
				}
			}
			else {
				if (NumberOfBParticles > 0){
					ParticleToSwitchAvailable = true;
					RandomParticleIndexInTypeArray = static_cast<int>(RNG.drawRandomNumber()*static_cast<double>(NumberOfBParticles));
					RandomParticleID = P.TypeBParticleIndices[RandomParticleIndexInTypeArray];
					TypeOfParticleBefore = ParticleType::B;
					ParticleNumbersPrefactor = static_cast<double>(NumberOfBParticles)/static_cast<double>(NumberOfAParticles+1);
					ChemicalPotentialSign = -1.0;
				}
			}
			if (ParticleToSwitchAvailable){
				double AcceptanceProbability = ParticleNumbersPrefactor*exp(-Beta*(P.computeChangeInPotentialEnergyBySwitching(RandomParticleID)  + ChemicalPotentialSign * ChemicalPotentialDiff));
				if (AcceptanceProbability >= 1.0 || (RNG.drawRandomNumber() < AcceptanceProbability)){
					P.switchParticleType(RandomParticleID, RandomParticleIndexInTypeArray, TypeOfParticleBefore);
				}
			}
		}
	}

	void runSimulation(int NumberOfRuns) {
		for (int i = 0; i < NumberOfRuns; i++){
			cerr << i << "|";
			runDisplacementSteps(100);
			runTypeChanges(1);
			DistributionNumberOfA[P.getNumberOfAParticles()]++;
		}
		cerr << endl;
	}

	void writeHistogramToFile() const {
		string FileName("Histogram_T="+to_string(Temperature)+".dat");
		ofstream FileStreamToWrite;
		FileStreamToWrite.open(FileName);
		FileStreamToWrite << "NumberOfAParticles NumberOfOccurences\n";
		for (int i = 0; i < TOTAL_NUMBER_OF_PARTICLES+1; i++){
			FileStreamToWrite << i << " " << DistributionNumberOfA[i] <<"\n";
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

int main(){
	SimulationManager S;
	S.initialize(0.5, 1.0, 1.0);
	cerr << S.P;

	S.runSimulation(100);
	S.writeHistogramToFile();
	S.writeParticleConfigurationToFile("ParticleConfig.dat");

	cerr << S.P;
}

