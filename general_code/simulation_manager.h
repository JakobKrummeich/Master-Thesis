#ifndef SIMULATION_MANAGER_INCLUDED
#define SIMULATION_MANAGER_INCLUDED

#include <iostream>
#include <chrono>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "global_variables.h"
#include "realRNG.h"
#include "particles.h"

using namespace std;

constexpr double MAXIMUM_DISPLACEMENT = 0.1;
constexpr double DISPLACEMENT_PROBABILITY = 0.9;

constexpr double MAXIMUM_VOLUME_CHANGE = 0.6;
constexpr double VOLUME_CHANGE_PROBABILITY = 0.01;

constexpr int NUMBER_OF_INITIAL_RANDOMIZATION_SWEEPS = 100;

constexpr int UPDATE_TIME_INTERVAL = 60;
constexpr int POT_ENERGY_UPDATE_INTERVAL = 200;

constexpr int NUMBER_OF_SAVED_STATES_PER_TEMPERATURE = 1;

constexpr int NUMBER_OF_THREADS = 1;

enum class MCModus{
	SGCMC = 0,
	PressureMC = 1
};

struct SimulationManager {
	Particles P;
	double Temperature;
	double Beta;
	double Pressure;

	realRNG RNG;

	int MinNumberOfA;
	int MaxNumberOfA;
	int NumberOfMCSweeps;

	double NumberOfTriedDisplacements;
	double NumberOfAcceptedDisplacements;

	double NumberOfTriedTypeChanges;
	double NumberOfAcceptedTypeChanges;

	double NumberOfTriedVolumeChanges;
	double NumberOfAcceptedVolumeChanges;

	string FileNameString;
	string DirectoryString;

	SimulationManager(int MinNumberOfA, int MaxNumberOfA, int NumberOfMCSweeps):
		MinNumberOfA(MinNumberOfA),
		MaxNumberOfA(MaxNumberOfA),
		NumberOfMCSweeps(NumberOfMCSweeps),
		NumberOfTriedDisplacements(0.0),
		NumberOfAcceptedDisplacements(0.0),
		NumberOfTriedTypeChanges(0.0),
		NumberOfAcceptedTypeChanges(0.0),
		NumberOfTriedVolumeChanges(0.0),
		NumberOfAcceptedVolumeChanges(0.0),
		DirectoryString("fresh_data/N="+to_string(TOTAL_NUMBER_OF_PARTICLES)+"/"){
	}

	void initializeParticles(double InitialDensity) {
		P.initialize(MinNumberOfA,InitialDensity);
	}

	void readInParticleState(string FileNameInitialState) {
		P.readInParticleState(FileNameInitialState);
	}

	void setTemperature(double NewTemperature){
		Temperature = NewTemperature;
		Beta = 1.0/Temperature;
	}
	
	void setPressure(double NewPressure){
		Pressure = NewPressure;
	}

	void setFileNameString(int RunNumber, MCModus M){
		if (M == MCModus::SGCMC){
		FileNameString = "N="+to_string(TOTAL_NUMBER_OF_PARTICLES)+"_T="+to_string(Temperature)+"_AvgDens="+to_string(P.getDensity())+"_MCRuns="+to_string(NumberOfMCSweeps)+"_epsAB="+to_string(AB_INTERACTION_STRENGTH)+"_"+to_string(RunNumber);
		}
		else {
			FileNameString = "N="+to_string(TOTAL_NUMBER_OF_PARTICLES)+"_T="+to_string(Temperature)+"_Pressure="+to_string(Pressure)+"_MCRuns="+to_string(NumberOfMCSweeps)+"_epsAB="+to_string(AB_INTERACTION_STRENGTH)+"_"+to_string(RunNumber);
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

	void runSGCMCSimulationForSingleTemperature(int RunCount, int NumberOfSavedStatesPerRun) {
		const auto StartTime = chrono::steady_clock::now();
		int NextUpdateTime = UPDATE_TIME_INTERVAL;
		int NextPotEnergyComputation = POT_ENERGY_UPDATE_INTERVAL;
		int NextStateSave = NumberOfMCSweeps / 2;
		int StateSaveInterval = (NumberOfMCSweeps / 2) / NumberOfSavedStatesPerRun;
		int SavedStateCount = 0;
		writeSGCMCMetaData();
		vector<int> NumberOfABuffer;
		vector<double> PotEnergyBuffer;
		#pragma omp critical(WRITE_TO_ERROR_STREAM)
		{
			cerr << "Run " << RunCount <<  ": T=" << Temperature << ". Simulation started." << endl << endl;
		}
		for (int i = 0; i < NumberOfMCSweeps; i++){
			for (int j = 0; j < TOTAL_NUMBER_OF_PARTICLES; j++){
				if (RNG.drawRandomNumber() <= DISPLACEMENT_PROBABILITY){
					runDisplacementStep();
				}
				else {
					runTypeChange();
				}
			}
			NumberOfABuffer.push_back(P.getNumberOfAParticles());
			if (i == NextPotEnergyComputation){
				PotEnergyBuffer.push_back(P.computePotentialEnergy());
				NextPotEnergyComputation += POT_ENERGY_UPDATE_INTERVAL;
			}
			if (i == NextStateSave){
				writeParticleStateToFile(DirectoryString+"State_"+FileNameString+"_"+to_string(SavedStateCount)+".dat");
				NextStateSave += StateSaveInterval;
				SavedStateCount++;
			}
			if (chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() >= NextUpdateTime){
				#pragma omp critical(WRITE_TO_ERROR_STREAM)
				{
					int Progress = NumberOfMCSweeps >= 100 ? (i / (NumberOfMCSweeps/100)) : i;
					cerr << "Run " << RunCount << ": T = " << Temperature <<  ". Progress: " << Progress << "%| Time passed: " << chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() << " s" << endl;
				}
				writeSGCMCResults(NumberOfABuffer, PotEnergyBuffer);
				NumberOfABuffer.clear();
				PotEnergyBuffer.clear();
				NextUpdateTime += UPDATE_TIME_INTERVAL;
			}
		}
		writeSGCMCResults(NumberOfABuffer, PotEnergyBuffer);
		#pragma omp critical(WRITE_TO_ERROR_STREAM)
		{
			cerr << "Run " << RunCount << ": Simulation of T=" << Temperature << " finished. Simulation-metadata: " << endl;
			cerr << "Tried displacements: " << NumberOfTriedDisplacements << "| Accepted displacements: " << NumberOfAcceptedDisplacements << "| Ratio of accepted displacements: " << static_cast<double>(NumberOfAcceptedDisplacements)/static_cast<double>(NumberOfTriedDisplacements) << endl;
			cerr << "Tried type changes: " << NumberOfTriedTypeChanges << "| Accepted type changes: " << NumberOfAcceptedTypeChanges << "| Ratio of accepted type changes: " << static_cast<double>(NumberOfAcceptedTypeChanges)/static_cast<double>(NumberOfTriedTypeChanges) << endl;
			cerr << "#VerletListBuilds: " << P.getNumberOfVerletListBuilds() << endl;
			cerr << "Computation time: " << chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() << " s for " << NumberOfMCSweeps << " MCSweeps." <<  endl << endl;
		}
	}

	void randomizeInitialPosition(int RunCount){
		double TemperatureBefore = Temperature;
		setTemperature(AA_INTERACTION_STRENGTH*10.0);
		for (int i = 0; i < NUMBER_OF_INITIAL_RANDOMIZATION_SWEEPS; i++){
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

	void writeSGCMCMetaData() const {
		string FileName(DirectoryString+"/NA_Series_"+FileNameString+".dat");
		string MetaDataString("BOX_LENGTH = "+to_string(P.getBoxLength())+" | AA_INTERACTION_STRENGTH = "+to_string(AA_INTERACTION_STRENGTH)+" | AB_INTERACTION_STRENGTH = "+to_string(AB_INTERACTION_STRENGTH)+" | MAXIMUM_DISPLACEMENT = "+to_string(MAXIMUM_DISPLACEMENT)+" | DISPLACEMENT_PROBABILITY = "+to_string(DISPLACEMENT_PROBABILITY)+" | MinNA = "+to_string(MinNumberOfA)+" | MaxNA = "+to_string(MaxNumberOfA)+'\n');
		ofstream FileStreamToWrite;
		FileStreamToWrite.open(FileName);
		FileStreamToWrite << MetaDataString;
		FileStreamToWrite.close();
		FileName = DirectoryString+"/PotEnergySeries_"+FileNameString+".dat";
		FileStreamToWrite.open(FileName);
		FileStreamToWrite << MetaDataString;
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

	void writeParticleStateToFile(string FileName) const {
		ofstream FileStreamToWrite;
		FileStreamToWrite.open(FileName);
		FileStreamToWrite << P;
		FileStreamToWrite.close();
	}

	void writePressureMCMetaData() const {
		string FileName(DirectoryString+"/density_series_"+FileNameString+".dat");
		string MetaDataString("AA_INTERACTION_STRENGTH = "+to_string(AA_INTERACTION_STRENGTH)+" | AB_INTERACTION_STRENGTH = "+to_string(AB_INTERACTION_STRENGTH)+" | MAXIMUM_DISPLACEMENT = "+to_string(MAXIMUM_DISPLACEMENT)+" | VOLUME_CHANGE_PROB = "+to_string(VOLUME_CHANGE_PROBABILITY)+" | MAX_VOLUME_CHANGE = "+to_string(MAXIMUM_VOLUME_CHANGE)+" | MinNA = "+to_string(MinNumberOfA)+" | MaxNA = "+to_string(MaxNumberOfA)+'\n');
		ofstream FileStreamToWrite;
		FileStreamToWrite.open(FileName);
		FileStreamToWrite << MetaDataString;
		FileStreamToWrite.close();
		FileName = DirectoryString+"/PotEnergySeries_"+FileNameString+".dat";
		FileStreamToWrite.open(FileName);
		FileStreamToWrite << MetaDataString;
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
	
	void runPressureMCForSinglePressure(int RunCount, int NumberOfSavedStatesPerRun) {
		const auto StartTime = chrono::steady_clock::now();
		int NextUpdateTime = UPDATE_TIME_INTERVAL;
		int NextPotEnergyComputation = POT_ENERGY_UPDATE_INTERVAL;
		int NextStateSave = NumberOfMCSweeps / 2;
		int StateSaveInterval = (NumberOfMCSweeps / 2) / NumberOfSavedStatesPerRun;
		int SavedStateCount = 0;
		writePressureMCMetaData();
		vector<double> DensityBuffer;
		vector<double> PotEnergyBuffer;
		#pragma omp critical(WRITE_TO_ERROR_STREAM)
		{
			cerr << "Run " << RunCount <<  ": p=" << Pressure << ". Simulation started." << endl << endl;
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
				writeParticleStateToFile(DirectoryString+"State_"+FileNameString+"_"+to_string(SavedStateCount)+".dat");
				NextStateSave += StateSaveInterval;
				SavedStateCount++;
			}
			if (chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() >= NextUpdateTime){
				#pragma omp critical(WRITE_TO_ERROR_STREAM)
				{
					int Progress = NumberOfMCSweeps >= 100 ? (i / (NumberOfMCSweeps/100)) : i;
					cerr << "Run " << RunCount << ": T = " << Temperature <<  ". Progress: " << Progress << "%| Time passed: " << chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() << " s" << endl;
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
			cerr << "Run " << RunCount << ": Simulation of T=" << Temperature << " finished. Simulation-metadata: " << endl;
			cerr << "Tried displacements: " << NumberOfTriedDisplacements << "| Accepted displacements: " << NumberOfAcceptedDisplacements << "| Ratio of accepted displacements: " << static_cast<double>(NumberOfAcceptedDisplacements)/static_cast<double>(NumberOfTriedDisplacements) << endl;
			cerr << "Tried volume changes: " << NumberOfTriedVolumeChanges << "| Accepted volume changes: " << NumberOfAcceptedVolumeChanges << "| Ratio of accepted volume changes: " << static_cast<double>(NumberOfAcceptedVolumeChanges)/static_cast<double>(NumberOfTriedVolumeChanges) << endl;
			cerr << "#VerletListBuilds: " << P.getNumberOfVerletListBuilds() << endl;
			cerr << "Computation time: " << chrono::duration_cast<chrono::seconds>(chrono::steady_clock::now()-StartTime).count() << " s for " << NumberOfMCSweeps << " MCSweeps." <<  endl << endl;
		}
	}
};

void runSGCMCSimulationForMultipleStartStates(double MaxTemperature, double MinTemperature, double TemperatureStep, int NumberOfRuns, int NumberOfMCSweeps, int RunNumberOffset, double Density) {

	const int NumberOfSavedStatesPerRun = NUMBER_OF_SAVED_STATES_PER_TEMPERATURE >= NumberOfRuns ?  NUMBER_OF_SAVED_STATES_PER_TEMPERATURE / NumberOfRuns : 1;

	#pragma omp parallel num_threads(NUMBER_OF_THREADS)
	{
		SimulationManager S(0, TOTAL_NUMBER_OF_PARTICLES, NumberOfMCSweeps);

		#pragma omp for
		for (int RunCount = RunNumberOffset; RunCount < NumberOfRuns+RunNumberOffset; RunCount++){
			S.initializeParticles(Density);
			S.randomizeInitialPosition(RunCount);
			for (double CurrentTemperature = MaxTemperature; CurrentTemperature > MinTemperature; CurrentTemperature -= TemperatureStep){
				S.setTemperature(CurrentTemperature);
				S.setFileNameString(RunCount, MCModus::SGCMC);
				S.resetCountersAndBuffers();
				S.runSGCMCSimulationForSingleTemperature(RunCount, NumberOfSavedStatesPerRun);
			}
		}
	}
}

void runPressureMCForMultipleStartStates(double MaxPressure, double MinPressure, double PressureStep, int NumberOfRuns, int NumberOfMCSweeps, int RunNumberOffset, double InitialDensity, double Temperature) {

	const int NumberOfSavedStatesPerRun = NUMBER_OF_SAVED_STATES_PER_TEMPERATURE >= NumberOfRuns ?  NUMBER_OF_SAVED_STATES_PER_TEMPERATURE / NumberOfRuns : 1;

	#pragma omp parallel num_threads(NUMBER_OF_THREADS)
	{
		SimulationManager S(TOTAL_NUMBER_OF_PARTICLES, TOTAL_NUMBER_OF_PARTICLES, NumberOfMCSweeps);
		S.setTemperature(Temperature);

		#pragma omp for
		for (int RunCount = RunNumberOffset; RunCount < NumberOfRuns+RunNumberOffset; RunCount++){
			S.initializeParticles(InitialDensity);
			S.randomizeInitialPosition(RunCount);
			for (double CurrentPressure = MaxPressure; CurrentPressure > MinPressure; CurrentPressure -= PressureStep){
				S.setPressure(CurrentPressure);
				S.setFileNameString(RunCount, MCModus::PressureMC);
				S.resetCountersAndBuffers();
				S.runPressureMCForSinglePressure(RunCount, NumberOfSavedStatesPerRun);
			}
		}
	}
}


#endif //SIMULATION_MANAGER_INCLUDED
