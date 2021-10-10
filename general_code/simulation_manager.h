#ifndef SIMULATION_MANAGER_INCLUDED
#define SIMULATION_MANAGER_INCLUDED

#include <iostream>
#include <chrono>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <sstream>

#include "global_constants.h"
#include "simulation_config.h"
#include "realRNG.h"
#include "particles.h"

using namespace std;

enum class MCModus{
	SGCMC = 0,
	PressureMC = 1,
	CMC = 2
};

struct SimulationManager {
	Particles P;
	double Temperature;
	double Beta;
	double Pressure;

	double ShearRate;

	realRNG RNG;

	int MinNumberOfA;
	int MaxNumberOfA;

	double NumberOfTriedDisplacements;
	double NumberOfAcceptedDisplacements;

	double NumberOfTriedTypeChanges;
	double NumberOfAcceptedTypeChanges;

	double NumberOfTriedVolumeChanges;
	double NumberOfAcceptedVolumeChanges;

	string FileNameString;
	string DirectoryString;

	SimulationManager(int MinNumberOfA, int MaxNumberOfA, string DirectoryForFreshData):
		MinNumberOfA(MinNumberOfA),
		MaxNumberOfA(MaxNumberOfA),
		NumberOfTriedDisplacements(0.0),
		NumberOfAcceptedDisplacements(0.0),
		NumberOfTriedTypeChanges(0.0),
		NumberOfAcceptedTypeChanges(0.0),
		NumberOfTriedVolumeChanges(0.0),
		NumberOfAcceptedVolumeChanges(0.0),
		DirectoryString(DirectoryForFreshData),
		ShearRate(0.0){
	}

	void initializeParticles(double InitialDensity) {
		P.initialize(MinNumberOfA,InitialDensity);
	}

	void readInParticleState(string FileNameInitialState, int StateCount, double Density) {
		P.readInParticleState(FileNameInitialState, StateCount, Density);
	}

	void setTemperature(double NewTemperature){
		Temperature = NewTemperature;
		Beta = 1.0/Temperature;
	}

	void setPressure(double NewPressure){
		Pressure = NewPressure;
	}

	void setShearRate(double NewShearRate){
		ShearRate = NewShearRate;
	}

	void setFileNameString(int RunNumber, MCModus M, int NumberOfMCSweeps){
		if (M == MCModus::SGCMC || M == MCModus::CMC){
			FileNameString = "N="+to_string(TOTAL_NUMBER_OF_PARTICLES)+"_T="+to_string(Temperature)+"_AvgDens="+to_string(P.getDensity())+"_MCRuns="+to_string(NumberOfMCSweeps)+"_epsAB="+to_string(AB_INTERACTION_STRENGTH)+"_"+to_string(RunNumber);
		}
		else {
			FileNameString = "N="+to_string(TOTAL_NUMBER_OF_PARTICLES)+"_T="+to_string(Temperature)+"_p="+to_string(Pressure)+"_MCRuns="+to_string(NumberOfMCSweeps)+"_epsAB="+to_string(AB_INTERACTION_STRENGTH)+"_"+to_string(RunNumber);
		}
	}

	void resetCountersAndBuffers(){
		NumberOfTriedDisplacements = 0.0;
		NumberOfAcceptedDisplacements = 0.0;
		NumberOfTriedTypeChanges = 0.0;
		NumberOfAcceptedTypeChanges = 0.0;
		NumberOfTriedVolumeChanges = 0.0;
		NumberOfAcceptedVolumeChanges = 0.0;
		P.resetCounters();
	}

	void runDisplacementStep() {
		NumberOfTriedDisplacements++;
		int RandomParticleID = static_cast<int>(RNG.drawRandomNumber()*static_cast<double>(TOTAL_NUMBER_OF_PARTICLES));
		double Deltas [DIMENSION];
		for (int i = 0; i < DIMENSION; i++){
			Deltas[i] = RNG.drawRandomNumber(-MAXIMUM_DISPLACEMENT, MAXIMUM_DISPLACEMENT)*P.getInverseBoxLength();
		}
		double AcceptanceProbability = exp(-P.computeChangeInPotentialEnergyByMoving(RandomParticleID, Deltas)*Beta);
		if (AcceptanceProbability >= 1.0 || (RNG.drawRandomNumber() < AcceptanceProbability)){
			P.updatePosition(RandomParticleID, Deltas);
			NumberOfAcceptedDisplacements++;
		}
	}

	void runDisplacementStepWithShear(double* ChangeInCoordinates, double& AttemptedShearMoves, double& AcceptedShearMoves) {
		NumberOfTriedDisplacements++;
		int RandomParticleID = static_cast<int>(RNG.drawRandomNumber()*static_cast<double>(TOTAL_NUMBER_OF_PARTICLES));
		double Deltas [DIMENSION];
		for (int i = 0; i < DIMENSION; i++){
			Deltas[i] = RNG.drawRandomNumber(-MAXIMUM_DISPLACEMENT, MAXIMUM_DISPLACEMENT)*P.getInverseBoxLength();
		}
		double yCoordinate = P.getPosition(RandomParticleID,1);
		Deltas[0] += ShearRate*2.0*P.getInverseBoxLength()*(yCoordinate-0.5);

		double AcceptanceProbability = exp(-P.computeChangeInPotentialEnergyByMoving(RandomParticleID, Deltas)*Beta);
		if (AcceptanceProbability >= 1.0 || (RNG.drawRandomNumber() < AcceptanceProbability)){
			P.updatePosition(RandomParticleID, Deltas);
			*(ChangeInCoordinates + DIMENSION*RandomParticleID) += Deltas[0];
			*(ChangeInCoordinates + DIMENSION*RandomParticleID+1) += Deltas[1];
			NumberOfAcceptedDisplacements++;
		}
	}

	void runTypeChange() {
		NumberOfTriedTypeChanges++;
		int NumberOfAParticles = P.getNumberOfAParticles();
		int NumberOfBParticles = P.getNumberOfBParticles();
		int RandomParticleIndexInTypeArray;
		double ParticleNumbersPrefactor;

		ParticleType ParticleTypeBefore;
		bool ParticleSwitchAllowed = false;

		if (RNG.drawRandomNumber() <= 0.5){
			if (NumberOfAParticles > MinNumberOfA){
				ParticleSwitchAllowed = true;
				RandomParticleIndexInTypeArray = static_cast<int>(RNG.drawRandomNumber()*static_cast<double>(NumberOfAParticles));
				ParticleTypeBefore = ParticleType::A;
				ParticleNumbersPrefactor = static_cast<double>(NumberOfAParticles)/static_cast<double>(NumberOfBParticles+1);
			}
		}
		else {
			if (NumberOfAParticles < MaxNumberOfA){
				ParticleSwitchAllowed = true;
				RandomParticleIndexInTypeArray = static_cast<int>(RNG.drawRandomNumber()*static_cast<double>(NumberOfBParticles));
				ParticleTypeBefore = ParticleType::B;
				ParticleNumbersPrefactor = static_cast<double>(NumberOfBParticles)/static_cast<double>(NumberOfAParticles+1);
			}
		}
		if (ParticleSwitchAllowed){
			double AcceptanceProbability = ParticleNumbersPrefactor*exp(-Beta*(P.computeChangeInPotentialEnergyBySwitching(RandomParticleIndexInTypeArray, ParticleTypeBefore)));
			if (AcceptanceProbability >= 1.0 || (RNG.drawRandomNumber() < AcceptanceProbability)){
				P.switchParticleType(RandomParticleIndexInTypeArray, ParticleTypeBefore);
				NumberOfAcceptedTypeChanges++;
			}
		}
	}

	void runVolumeChange() {
		NumberOfTriedVolumeChanges++;
		double RandomVolumeChange = RNG.drawRandomNumber(-MAXIMUM_VOLUME_CHANGE,MAXIMUM_VOLUME_CHANGE);
		double AcceptanceProbability = exp(-Beta*(P.computeChangeInPotentialEnergyByChangingVolume(RandomVolumeChange) + Pressure*RandomVolumeChange) + static_cast<double>(TOTAL_NUMBER_OF_PARTICLES) * log(1.0+RandomVolumeChange/P.getVolume())); // when computing the change in pot energy we already rebuild the verlet list and change the box parameters!
		if (AcceptanceProbability >= 1.0 || (RNG.drawRandomNumber() < AcceptanceProbability)){
			NumberOfAcceptedVolumeChanges++;
		}
		else {
			P.changeVolume(-RandomVolumeChange); // if the move got rejected, we have to REVERT the change to the box parameters that were already done when computing the change in pot energy
		}
	}

	void equilibrateFixedSweeps(int NumberOfEquilibrationSteps){
		for (int i = 0; i < NumberOfEquilibrationSteps; i++){
			runSingleSGCMCSweep();
		}
	}

	int equilibrateFixedTime(int MaxRuntimeInMinutes, chrono::time_point<chrono::steady_clock> StartTime){
		int SweepCount = 0;
		for (; chrono::duration_cast<chrono::minutes>(chrono::steady_clock::now()-StartTime).count() < MaxRuntimeInMinutes; SweepCount++){
			runSingleSGCMCSweep();
		}
		return SweepCount;
	}

	void runSingleCMCSweep(){
		for (int j = 0; j < TOTAL_NUMBER_OF_PARTICLES; j++){
			runDisplacementStep();
		}
	}

	void runSingleSGCMCSweep(){
		for (int j = 0; j < TOTAL_NUMBER_OF_PARTICLES; j++){
			if (RNG.drawRandomNumber() <= DISPLACEMENT_PROBABILITY){
				runDisplacementStep();
			}
			else {
				runTypeChange();
			}
		}
	}

	void runSingleSGCMCSweepWithShear(double* ChangeInCoordinates, double& AttemptedShearMoves, double& AcceptedShearMoves){
		for (int j = 0; j < TOTAL_NUMBER_OF_PARTICLES; j++){
			if (RNG.drawRandomNumber() <= DISPLACEMENT_PROBABILITY){
				runDisplacementStepWithShear(ChangeInCoordinates, AttemptedShearMoves, AcceptedShearMoves);
			}
			else {
				runTypeChange();
			}
		}
	}

	void writeSimulationMetaDataToErrorStream(int RunCount, int TimePassedInSeconds, int NumberOfSweeps){
		#pragma omp critical(WRITE_TO_ERROR_STREAM)
		{
			cerr << "Run " << RunCount << ": Simulation of T=" << Temperature << " finished. Simulation-metadata: " << endl;
			cerr << "Tried displacements: " << NumberOfTriedDisplacements << "| Accepted displacements: " << NumberOfAcceptedDisplacements << "| Ratio of accepted displacements: " << static_cast<double>(NumberOfAcceptedDisplacements)/static_cast<double>(NumberOfTriedDisplacements) << endl;
			cerr << "Tried type changes : " << NumberOfTriedTypeChanges   << "| Accepted type changes:  " << NumberOfAcceptedTypeChanges   << "| Ratio of accepted type changes:  " << static_cast<double>(NumberOfAcceptedTypeChanges)/static_cast<double>(NumberOfTriedTypeChanges) << endl;
			cerr << "#VerletListBuilds: " << P.getNumberOfVerletListBuilds() << endl;
			cerr << "Computation time: " << TimePassedInSeconds << " s for " << NumberOfSweeps << " sweeps" << endl << endl;
		}
	}

	void writeSimulationStartInfoToErrorStream(int RunCount){
		#pragma omp critical(WRITE_TO_ERROR_STREAM)
		{
			cerr << "Run " << RunCount <<  ": T=" << Temperature << " | Simulation started." << endl << endl;
		}
	}

	void writeSimulationProgressToErrorStream(int RunCount, int MaxNumberOfSweeps, int SweepCount, double Temperature, int ElapsedTime){
		#pragma omp critical(WRITE_TO_ERROR_STREAM)
		{
			int Progress = MaxNumberOfSweeps >= 100 ? (SweepCount / (MaxNumberOfSweeps/100)) : SweepCount;
			cerr << "Run " << RunCount << ": T = " << Temperature <<  ". Progress: " << Progress << "%| SweepCount: " << SweepCount << "| Time passed: " << ElapsedTime << " s" << endl;
		}
	}

	void runSGCMCSimulationForSingleTemperatureTimeControlled(int RunCount, int MaxRuntimeInMinutes, chrono::time_point<chrono::steady_clock> StartTime, int MaxNumberOfSweeps){

		auto StartOfDataTaking = chrono::steady_clock::now();

		writeSGCMCMetaData();
		writeSimulationStartInfoToErrorStream(RunCount);

		int NextUpdateTime = UPDATE_TIME_INTERVAL;
		int NextPotEnergyComputation = POT_ENERGY_UPDATE_INTERVAL;
		vector<int> NumberOfABuffer;
		vector<double> PotEnergyBuffer;

		int SweepCount = 0;
		for (; 	chrono::duration_cast<chrono::minutes>(chrono::steady_clock::now()-StartTime).count() < MaxRuntimeInMinutes && SweepCount < MaxNumberOfSweeps; SweepCount++){
			runSingleSGCMCSweep();
			NumberOfABuffer.push_back(P.getNumberOfAParticles());
			if (SweepCount == NextPotEnergyComputation){
				PotEnergyBuffer.push_back(P.computePotentialEnergy());
				NextPotEnergyComputation += POT_ENERGY_UPDATE_INTERVAL;
			}
			if (chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() >= NextUpdateTime){
				writeSimulationProgressToErrorStream(RunCount, MaxNumberOfSweeps, SweepCount, Temperature, chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count());
				writeSGCMCResults(NumberOfABuffer, PotEnergyBuffer);
				NumberOfABuffer.clear();
				PotEnergyBuffer.clear();
				NextUpdateTime += UPDATE_TIME_INTERVAL;
			}
		}
		writeSGCMCResults(NumberOfABuffer, PotEnergyBuffer);
		writeParticleStateToFile(DirectoryString+"/State_"+FileNameString+"_"+to_string(0)+".dat");
		writeSimulationMetaDataToErrorStream(RunCount, chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartOfDataTaking).count(), SweepCount);
	}

	void updateAverageTraveledDistances(const Particles& P, const double* ChangeInCoordinates, vector<double>& AverageTraveledDistances, const int NumberOfyValues) const {
		vector<double> NewEntries(NumberOfyValues,0.0);
		vector<int> NumberOfValues(NumberOfyValues,0);
		for (int ParticleID = 0; ParticleID < TOTAL_NUMBER_OF_PARTICLES; ParticleID++){
			int Index = static_cast<int>(P.getPosition(ParticleID,1)*static_cast<double>(NumberOfyValues));
			NewEntries[Index] += ChangeInCoordinates[DIMENSION*ParticleID];
			NumberOfValues[Index]++;
		}
		for (int i = 0; i < NumberOfyValues; i++){
			AverageTraveledDistances.push_back(NewEntries[i] / static_cast<double>(NumberOfValues[i]));
		}
	}

	void writeTraveledDistances(const double* ChangeInCoordinates) const {
		ofstream FileStreamToWrite;
		string FileName = DirectoryString+"/TraveledDistancesvsInity_"+FileNameString+".dat";
		FileStreamToWrite.open(FileName, ios_base::app);
		for (int i = 0; i < TOTAL_NUMBER_OF_PARTICLES; i++){
			FileStreamToWrite << ChangeInCoordinates[DIMENSION*i] << '\t';
		}
		FileStreamToWrite << '\n';
		FileStreamToWrite.close();
	}

	void runSGCMCSimulationForSingleTemperatureWithShear(int RunCount, int MaxRuntimeInMinutes, chrono::time_point<chrono::steady_clock> StartTime, int MaxNumberOfSweeps, const int NumberOfyValues){

		auto StartOfDataTaking = chrono::steady_clock::now();

		writeSGCMCMetaDataWithShear(NumberOfyValues,P);
		writeSimulationStartInfoToErrorStream(RunCount);

		int NextUpdateTime = UPDATE_TIME_INTERVAL;
		int NextPotEnergyComputation = POT_ENERGY_UPDATE_INTERVAL;
		vector<int> NumberOfABuffer;
		vector<double> PotEnergyBuffer;
		vector<double> AverageTraveledDistances;

		double ChangeInCoordinates [DIMENSION*TOTAL_NUMBER_OF_PARTICLES]{};
		double AttemptedShearMoves = 0.0;
		double AcceptedShearMoves = 0.0;

		int SweepCount = 0;
		for (; 	chrono::duration_cast<chrono::minutes>(chrono::steady_clock::now()-StartTime).count() < MaxRuntimeInMinutes && SweepCount < MaxNumberOfSweeps; SweepCount++){
			runSingleSGCMCSweepWithShear(ChangeInCoordinates, AttemptedShearMoves, AcceptedShearMoves);
			NumberOfABuffer.push_back(P.getNumberOfAParticles());

			if (SweepCount == NextPotEnergyComputation){
				PotEnergyBuffer.push_back(P.computePotentialEnergy());
				NextPotEnergyComputation += POT_ENERGY_UPDATE_INTERVAL;

				updateAverageTraveledDistances(P, ChangeInCoordinates, AverageTraveledDistances, NumberOfyValues);
			}
			if (chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() >= NextUpdateTime){
				writeSimulationProgressToErrorStream(RunCount, MaxNumberOfSweeps, SweepCount, Temperature, chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count());
				writeSGCMCResultsWithShear(NumberOfABuffer, PotEnergyBuffer, AverageTraveledDistances, NumberOfyValues);
				NumberOfABuffer.clear();
				PotEnergyBuffer.clear();
				AverageTraveledDistances.clear();
				NextUpdateTime += UPDATE_TIME_INTERVAL;
			}
		}

		PotEnergyBuffer.push_back(P.computePotentialEnergy());
		updateAverageTraveledDistances(P, ChangeInCoordinates, AverageTraveledDistances, NumberOfyValues);

		writeTraveledDistances(ChangeInCoordinates);

		writeSGCMCResultsWithShear(NumberOfABuffer, PotEnergyBuffer, AverageTraveledDistances, NumberOfyValues);
		writeParticleStateToFile(DirectoryString+"/State_"+FileNameString+"_"+to_string(0)+".dat");
		cerr << "Tried shear moves: " << AttemptedShearMoves << " | AcceptedShearMoves: " << AcceptedShearMoves << " | Ratio: " << AcceptedShearMoves/AttemptedShearMoves << endl;
		writeSimulationMetaDataToErrorStream(RunCount, chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartOfDataTaking).count(), SweepCount);
	}

	void runCMCSimulationForSingleTemperature(int RunCount, int MaxRuntimeInMinutes, chrono::time_point<chrono::steady_clock> StartTime, int MaxNumberOfSweeps, int NumberOfSavedStatesPerRun){

		auto StartOfDataTaking = chrono::steady_clock::now();

		writeCMCMetaData();
		writeSimulationStartInfoToErrorStream(RunCount);

		int NextUpdateTime = UPDATE_TIME_INTERVAL;
		int NextPotEnergyComputation = POT_ENERGY_UPDATE_INTERVAL;
		int StateSaveInterval = MaxNumberOfSweeps / NumberOfSavedStatesPerRun;
		int NextStateSave = StateSaveInterval;
		int SavedStateCount = 0;
		vector<double> PotEnergyBuffer;

		int SweepCount = 0;
		for (; chrono::duration_cast<chrono::minutes>(chrono::steady_clock::now()-StartTime).count() < MaxRuntimeInMinutes && SweepCount < MaxNumberOfSweeps; SweepCount++){
			runSingleCMCSweep();
			if (SweepCount == NextPotEnergyComputation){
				PotEnergyBuffer.push_back(P.computePotentialEnergy());
				NextPotEnergyComputation += POT_ENERGY_UPDATE_INTERVAL;
			}
			if (SweepCount == NextStateSave){
				writeParticleStateToFile(DirectoryString+"/State_"+FileNameString+"_"+to_string(SavedStateCount)+".dat");
				NextStateSave += StateSaveInterval;
				SavedStateCount++;
			}
			if (chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() >= NextUpdateTime){
				writeSimulationProgressToErrorStream(RunCount, MaxNumberOfSweeps, SweepCount, Temperature, chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count());
				writeCMCResults(PotEnergyBuffer);
				PotEnergyBuffer.clear();
				NextUpdateTime += UPDATE_TIME_INTERVAL;
			}
		}
		if (SavedStateCount < NumberOfSavedStatesPerRun){
			writeParticleStateToFile(DirectoryString+"/State_"+FileNameString+"_"+to_string(SavedStateCount)+".dat");
		}
		writeCMCResults(PotEnergyBuffer);
		writeSimulationMetaDataToErrorStream(RunCount, chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartOfDataTaking).count(), SweepCount);
	}

	void randomizeInitialPosition(int RunCount, int NumberOfInitialRandomizationSweeps){
		double TemperatureBefore = Temperature;
		setTemperature(AA_INTERACTION_STRENGTH*10.0);
		for (int i = 0; i < NumberOfInitialRandomizationSweeps; i++){
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
			cerr << "#VerletListBuilds: " << P.getNumberOfVerletListBuilds() << endl << endl;
		}
		setTemperature(TemperatureBefore);
	}

	string getMetaDataString(MCModus M) const {
		ostringstream OutputStream;
		OutputStream.precision(numeric_limits<long double>::digits10+1);
		if (M == MCModus::SGCMC){
			OutputStream << fixed << "BOX_LENGTH = " << P.getBoxLength() << " | AA_INTERACTION_STRENGTH = " << AA_INTERACTION_STRENGTH << " | AB_INTERACTION_STRENGTH = " << AB_INTERACTION_STRENGTH << " | MAXIMUM_DISPLACEMENT = " << MAXIMUM_DISPLACEMENT << " | DISPLACEMENT_PROBABILITY = " << DISPLACEMENT_PROBABILITY << " | MinNA = " << MinNumberOfA << " | MaxNA = " << MaxNumberOfA << " | ShearRate = " << ShearRate << endl;
			return OutputStream.str();
		}
		else if (M == MCModus::PressureMC) {
			OutputStream << fixed << "AA_INTERACTION_STRENGTH = " << AA_INTERACTION_STRENGTH << " | AB_INTERACTION_STRENGTH = " << AB_INTERACTION_STRENGTH << " | MAXIMUM_DISPLACEMENT = " << MAXIMUM_DISPLACEMENT << " | VOLUME_CHANGE_PROB = " << VOLUME_CHANGE_PROBABILITY << " | MAX_VOLUME_CHANGE = " << MAXIMUM_VOLUME_CHANGE << " | MinNA = " << MinNumberOfA << " | MaxNA = " << MaxNumberOfA << '\n';
			return OutputStream.str();
		}
		else {
			OutputStream << fixed << "BOX_LENGTH = " << P.getBoxLength() << " | AA_INTERACTION_STRENGTH = " << AA_INTERACTION_STRENGTH << " | AB_INTERACTION_STRENGTH = " << AB_INTERACTION_STRENGTH << " | MAXIMUM_DISPLACEMENT = " << MAXIMUM_DISPLACEMENT << " | NA = " << MinNumberOfA << " | ShearRate = " << ShearRate << '\n';
			return OutputStream.str();
		}
	}

	void writeCMCMetaData() const {
		string FileName(DirectoryString+"/PotEnergySeries_"+FileNameString+".dat");
		ofstream FileStreamToWrite;
		FileStreamToWrite.open(FileName);
		FileStreamToWrite << getMetaDataString(MCModus::CMC);
		FileStreamToWrite.close();
	}

	void writeCMCResults(const vector<double>& PotEnergyBuffer) const {
		string FileName(DirectoryString+"/PotEnergySeries_"+FileNameString+".dat");
		ofstream FileStreamToWrite;
		FileStreamToWrite.open(FileName, ios_base::app);
		for (int i = 0; i < PotEnergyBuffer.size(); i++){
			FileStreamToWrite << PotEnergyBuffer[i] << '\n';
		}
		FileStreamToWrite.close();
	}

	void writeSGCMCMetaData() const {
		string FileName(DirectoryString+"/NA_Series_"+FileNameString+".dat");
		ofstream FileStreamToWrite;
		FileStreamToWrite.open(FileName);
		FileStreamToWrite << getMetaDataString(MCModus::SGCMC);
		FileStreamToWrite.close();
		FileName = DirectoryString+"/PotEnergySeries_"+FileNameString+".dat";
		FileStreamToWrite.open(FileName);
		FileStreamToWrite << getMetaDataString(MCModus::SGCMC);
		FileStreamToWrite.close();
	}

	void writeSGCMCMetaDataWithShear(int NumberOfyValues, const Particles& P) const {
		string FileName(DirectoryString+"/NA_Series_"+FileNameString+".dat");
		ofstream FileStreamToWrite;
		FileStreamToWrite.open(FileName);
		FileStreamToWrite << getMetaDataString(MCModus::SGCMC);
		FileStreamToWrite.close();

		FileName = DirectoryString+"/PotEnergySeries_"+FileNameString+".dat";
		FileStreamToWrite.open(FileName);
		FileStreamToWrite << getMetaDataString(MCModus::SGCMC);
		FileStreamToWrite.close();

		FileName = DirectoryString+"/AvgTraveledDistances_"+FileNameString+".dat";
		FileStreamToWrite.open(FileName);
		FileStreamToWrite << getMetaDataString(MCModus::SGCMC);
		double yDelta = 1.0/static_cast<double>(NumberOfyValues);
		double Currenty = yDelta*0.5;
		for (int i = 0; i < NumberOfyValues-1; i++){
			FileStreamToWrite << Currenty << '\t';
			Currenty += yDelta;
		}
		FileStreamToWrite << Currenty << '\n';
		FileStreamToWrite.close();

		FileName = DirectoryString+"/TraveledDistancesvsInity_"+FileNameString+".dat";
		FileStreamToWrite.open(FileName);
		FileStreamToWrite << getMetaDataString(MCModus::SGCMC);
		for (int i = 0; i < TOTAL_NUMBER_OF_PARTICLES; i++){
			FileStreamToWrite << P.getPosition(i,1) << '\t';
		}
		FileStreamToWrite << '\n';
		FileStreamToWrite.close();
	}

	void writeSGCMCResults(const vector<int>& NumberOfABuffer, const vector<double>& PotEnergyBuffer) const {
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

	void writeSGCMCResultsWithShear(const vector<int>& NumberOfABuffer, const vector<double>& PotEnergyBuffer, const vector<double>& AverageTraveledDistances, const int NumberOfyValues) const {
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

		FileName = DirectoryString+"/AvgTraveledDistances_"+FileNameString+".dat";
		FileStreamToWrite.open(FileName, ios_base::app);

		const int NumberOfRows = AverageTraveledDistances.size()/NumberOfyValues;
		for (int i = 0; i < NumberOfRows; i++){
			for (int j = 0; j < NumberOfyValues-1; j++){
				FileStreamToWrite << AverageTraveledDistances[NumberOfyValues*i+j] << '\t';
			}
			FileStreamToWrite << AverageTraveledDistances[NumberOfyValues*i+NumberOfyValues-1] << '\n';
		}
		FileStreamToWrite.close();
	}

	void writeParticleStateToFile(string FileName) const {
		ofstream FileStreamToWrite;
		FileStreamToWrite.open(FileName);
		FileStreamToWrite << P;
		FileStreamToWrite.close();
	}

	void writePressureMCMetaData() const {
		string FileName(DirectoryString+"/density_series_"+FileNameString+".dat");
		ofstream FileStreamToWrite;
		FileStreamToWrite.open(FileName);
		FileStreamToWrite << getMetaDataString(MCModus::PressureMC);
		FileStreamToWrite.close();
		FileName = DirectoryString+"/PotEnergySeries_"+FileNameString+".dat";
		FileStreamToWrite.open(FileName);
		FileStreamToWrite << getMetaDataString(MCModus::PressureMC);
		FileStreamToWrite.close();
	}

	void writePressureMCResults(const vector<double>& DensityBuffer, const vector<double>& PotEnergyBuffer) const {
		string FileName(DirectoryString+"/density_series_"+FileNameString+".dat");
		ofstream FileStreamToWrite;
		FileStreamToWrite.open(FileName, ios_base::app);
		for (int i = 0; i < DensityBuffer.size(); i++){
			FileStreamToWrite << DensityBuffer[i] << '\n';
		}
		FileStreamToWrite.close();
		FileName = DirectoryString+"/PotEnergySeries_"+FileNameString+".dat";
		FileStreamToWrite.open(FileName, ios_base::app);
		for (int i = 0; i < PotEnergyBuffer.size(); i++){
			FileStreamToWrite << PotEnergyBuffer[i] << '\n';
		}
		FileStreamToWrite.close();
	}

	void runPressureMCForSinglePressure(int RunCount, int NumberOfSavedStatesPerRun, int NumberOfEquilibrationSweeps, int NumberOfMCSweeps) {
		const auto StartTime = chrono::steady_clock::now();
		int NextUpdateTime = UPDATE_TIME_INTERVAL;
		int NextPotEnergyComputation = POT_ENERGY_UPDATE_INTERVAL;
		int NextStateSave = NumberOfMCSweeps;
		int StateSaveInterval = NumberOfMCSweeps / NumberOfSavedStatesPerRun;
		int SavedStateCount = 0;
		writePressureMCMetaData();
		vector<double> DensityBuffer;
		vector<double> PotEnergyBuffer;
		#pragma omp critical(WRITE_TO_ERROR_STREAM)
		{
			cerr << "Run " << RunCount <<  ": p = " << Pressure << " | Simulation started." << endl << endl;
		}
		for (int i = 0; i < NumberOfEquilibrationSweeps; i++){
			for (int j = 0; j < TOTAL_NUMBER_OF_PARTICLES; j++){
				if (RNG.drawRandomNumber() <= VOLUME_CHANGE_PROBABILITY){
					runVolumeChange();
				}
				else {
					runDisplacementStep();
				}
			}
		}
		for (int i = 0; i < NumberOfMCSweeps; i++){
			for (int j = 0; j < TOTAL_NUMBER_OF_PARTICLES; j++){
				if (RNG.drawRandomNumber() <= VOLUME_CHANGE_PROBABILITY){
					runVolumeChange();
				}
				else {
					runDisplacementStep();
				}
			}
			DensityBuffer.push_back(static_cast<double>(TOTAL_NUMBER_OF_PARTICLES)/P.getVolume());
			if (i == NextPotEnergyComputation){
				PotEnergyBuffer.push_back(P.computePotentialEnergy());
				NextPotEnergyComputation += POT_ENERGY_UPDATE_INTERVAL;
			}
			if (i == NextStateSave){
				writeParticleStateToFile(DirectoryString+"/State_"+FileNameString+"_"+to_string(SavedStateCount)+".dat");
				NextStateSave += StateSaveInterval;
				SavedStateCount++;
			}
			if (chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() >= NextUpdateTime){
				#pragma omp critical(WRITE_TO_ERROR_STREAM)
				{
					int Progress = NumberOfMCSweeps >= 100 ? (i / (NumberOfMCSweeps/100)) : i;
					cerr << "Run " << RunCount << ": p = " << Pressure <<  " | Progress: " << Progress << "%| Time passed: " << chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() << " s" << endl;
				}
				writePressureMCResults(DensityBuffer, PotEnergyBuffer);
				DensityBuffer.clear();
				PotEnergyBuffer.clear();
				NextUpdateTime += UPDATE_TIME_INTERVAL;
			}
		}
		writePressureMCResults(DensityBuffer, PotEnergyBuffer);
		#pragma omp critical(WRITE_TO_ERROR_STREAM)
		{
			cerr << "Run " << RunCount << ": Simulation of p =" << Pressure << " finished. Simulation-metadata: " << endl;
			cerr << "Tried displacements: " << NumberOfTriedDisplacements << "| Accepted displacements: " << NumberOfAcceptedDisplacements << "| Ratio of accepted displacements: " << static_cast<double>(NumberOfAcceptedDisplacements)/static_cast<double>(NumberOfTriedDisplacements) << endl;
			cerr << "Tried volume changes: " << NumberOfTriedVolumeChanges << "| Accepted volume changes: " << NumberOfAcceptedVolumeChanges << "| Ratio of accepted volume changes: " << static_cast<double>(NumberOfAcceptedVolumeChanges)/static_cast<double>(NumberOfTriedVolumeChanges) << endl;
			cerr << "#VerletListBuilds: " << P.getNumberOfVerletListBuilds() << endl;
			cerr << "Computation time: " << chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() << " s for " << NumberOfMCSweeps << " MCSweeps." <<  endl << endl;
		}
	}
};

void runPressureMCForMultipleStartStates(double MaxPressure, double MinPressure, double PressureStep, int NumberOfRuns, int RunNumberOffset, double InitialDensity, double Temperature, string DirectoryForFreshData, int NumberOfInitialRandomizationSweeps, int NumberOfEquilibrationSweeps, int NumberOfMCSweeps) {

	const int NumberOfSavedStatesPerRun = NUMBER_OF_SAVED_STATES_PER_TEMPERATURE >= NumberOfRuns ?  NUMBER_OF_SAVED_STATES_PER_TEMPERATURE / NumberOfRuns : 1;

	#pragma omp parallel num_threads(NUMBER_OF_THREADS)
	{
		SimulationManager S(TOTAL_NUMBER_OF_PARTICLES, TOTAL_NUMBER_OF_PARTICLES, DirectoryForFreshData);
		S.setTemperature(Temperature);

		#pragma omp for
		for (int RunCount = RunNumberOffset; RunCount < NumberOfRuns+RunNumberOffset; RunCount++){
			S.initializeParticles(InitialDensity);
			S.randomizeInitialPosition(RunCount, NumberOfInitialRandomizationSweeps);
			for (double CurrentPressure = MaxPressure; CurrentPressure > MinPressure; CurrentPressure -= PressureStep){
				S.setPressure(CurrentPressure);
				S.setFileNameString(RunCount, MCModus::PressureMC, NumberOfMCSweeps);
				S.resetCountersAndBuffers();
				S.runPressureMCForSinglePressure(RunCount, NumberOfSavedStatesPerRun, NumberOfEquilibrationSweeps, NumberOfMCSweeps);
			}
		}
	}
}


#endif //SIMULATION_MANAGER_INCLUDED
