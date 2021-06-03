class StructureFactorComputator{
	private:
		const double kMin;
		const double kMax;

		realRNG RNG;

		struct Combination{
			int nx;
			int ny;

			Combination(){
			}

			Combination(int nx, int ny):
				nx(nx),
				ny(ny){
			}

			bool operator==(const Combination& RHS) const {
				return ((nx == RHS.nx && ny == RHS.ny) || (nx == RHS.ny && ny == RHS.nx));
			}
		};

		struct kCombinationMappingEntry{
			double kMagnitude;
			vector<Combination> Combinations;
		};

		struct RowEntry {
			double SAA;
			double SBB;
			double SAB;
			double Scc;
			int NumberOfDataPoints;

			RowEntry():
				SAA(0.0),
				SBB(0.0),
				SAB(0.0),
				Scc(0.0),
				NumberOfDataPoints(0){
			}

			RowEntry& operator+=(const RowEntry& RHS){
				this -> SAA += RHS.SAA;
				this -> SBB += RHS.SBB;
				this -> SAB += RHS.SAB;
				this -> Scc += RHS.Scc;
				this -> NumberOfDataPoints += RHS.NumberOfDataPoints;
				return *this;
			}
		};

		struct SingleTemperatureResults{
			double Temperature;
			vector<RowEntry> StructureFactors;

			SingleTemperatureResults(int NumberOfkEntries, double Temperature):
				StructureFactors(NumberOfkEntries),
				Temperature(Temperature){
			}
		};

		vector<kCombinationMappingEntry> kCombinationMapping;
		vector< SingleTemperatureResults > Results;

		double computeCosSum(const double* const Positions, const fvec<int, TOTAL_NUMBER_OF_PARTICLES>& ParticleIndices, double kx, double ky) const{
			double Sum = 0.0;
			for (int i = 0; i < ParticleIndices.size(); i++){
				Sum += cos( kx* BOX_LENGTH*(*(Positions+DIMENSION*ParticleIndices[i])) + ky * BOX_LENGTH * (*(Positions+DIMENSION*ParticleIndices[i]+1)));
			}
			return Sum;
		}

		double computeSinSum(const double* const Positions, const fvec<int, TOTAL_NUMBER_OF_PARTICLES>& ParticleIndices, double kx, double ky) const{
			double Sum = 0.0;
			for (int i = 0; i < ParticleIndices.size(); i++){
				Sum += sin(kx* BOX_LENGTH * (*(Positions+DIMENSION*ParticleIndices[i])) + ky * BOX_LENGTH * (*(Positions+DIMENSION*ParticleIndices[i]+1)));
			}
			return Sum;
		}

		void findkValuesOnGrid(){
			double kWidth = 0.04;
			double CurrentkMag = kMin;
			int NumberOfAttemptedAveragesPerk = 1000;
			int MaxNumberOfCombinationsPerInterval = 100;
			vector<Combination> CombinationsFound;

			while (CurrentkMag < kMax){
				CombinationsFound.clear();
				kCombinationMappingEntry NewMapping{};
				bool MagnitudeFound = false;
				for (int i = 0; i < NumberOfAttemptedAveragesPerk && CombinationsFound.size() < MaxNumberOfCombinationsPerInterval; i++){
					double RandomAngle = RNG.drawRandomNumber(0.0, 0.5 * M_PI);
					double RandomkMagnitude = RNG.drawRandomNumber(-kWidth,kWidth) + CurrentkMag;
					int nx = round(BOX_LENGTH*RandomkMagnitude*cos(RandomAngle)/(2.0*M_PI));
					int ny = round(BOX_LENGTH*RandomkMagnitude*sin(RandomAngle)/(2.0*M_PI));
					double GridkMagnitude = 2.0*M_PI/BOX_LENGTH*sqrt(static_cast<double>(nx*nx+ny*ny));
					if (GridkMagnitude <= (CurrentkMag + kWidth) && GridkMagnitude >= (CurrentkMag - kWidth) && GridkMagnitude >= kMin){
						bool sameCombinationAlreadyFoundBefore = false;
						Combination NewCombination(nx,ny);
						for (int j = 0; j < CombinationsFound.size() && !sameCombinationAlreadyFoundBefore; j++){
							if (CombinationsFound[j] == NewCombination){
								sameCombinationAlreadyFoundBefore = true;
							}
						}
						if (!sameCombinationAlreadyFoundBefore){
							CombinationsFound.push_back(NewCombination);
							if (!MagnitudeFound){
								NewMapping.kMagnitude = GridkMagnitude;
								NewMapping.Combinations.push_back(NewCombination);
								MagnitudeFound = true;
							}
							else {
								NewMapping.kMagnitude = CurrentkMag;
								NewMapping.Combinations.push_back(NewCombination);
							}
						}
					}
				}
				if (MagnitudeFound){
					kCombinationMapping.push_back(NewMapping);
				}
				CurrentkMag += 2.0*kWidth;
			}
		}

	public:
		StructureFactorComputator(double MaxTemperature, double MinTemperature, double TemperatureStep, double kMax):
			kMin(2.0*M_PI/BOX_LENGTH),
			kMax(kMax){
			findkValuesOnGrid();
			for (double CurrentTemperature = MaxTemperature; CurrentTemperature > MinTemperature; CurrentTemperature -= TemperatureStep){
				Results.push_back(SingleTemperatureResults(kCombinationMapping.size(), CurrentTemperature));
			}
		}

		void computeNewStructureFactorValues(const double* const Positions, const fvec<int, TOTAL_NUMBER_OF_PARTICLES>& TypeAParticleIndices, const fvec<int, TOTAL_NUMBER_OF_PARTICLES>& TypeBParticleIndices, int TemperatureIndex){
			if (TemperatureIndex >= Results.size()){
				cerr << "Invalid temperature index in computeNewStructureFactorValues! Size of Results:  " << Results.size() << " , TemperatureIndex given: " << TemperatureIndex << endl;
				return;
			}
			double xA = static_cast<double>(TypeAParticleIndices.size())/static_cast<double>(TOTAL_NUMBER_OF_PARTICLES);
			double xB = static_cast<double>(TypeBParticleIndices.size())/static_cast<double>(TOTAL_NUMBER_OF_PARTICLES);
			for (int i = 0; i < kCombinationMapping.size(); i++){
				for (int j = 0; j < kCombinationMapping[i].Combinations.size(); j++){
					int nx = kCombinationMapping[i].Combinations[j].nx;
					int ny = kCombinationMapping[i].Combinations[j].ny;
					for (int k = 0; k < (abs(nx) == abs(ny) ? 1 : 2); k++){
						for (int Sign = 1; Sign >= ((nx == 0 || ny == 0)? 1: -1); Sign -= 2){
							double kx = Sign * nx * 2.0*M_PI/BOX_LENGTH;
							double ky = ny * 2.0*M_PI/BOX_LENGTH;
							double cosSumA = computeCosSum(Positions, TypeAParticleIndices, kx, ky);
							double sinSumA = computeSinSum(Positions, TypeAParticleIndices, kx, ky);
							double cosSumB = computeCosSum(Positions, TypeBParticleIndices, kx, ky);
							double sinSumB = computeSinSum(Positions, TypeBParticleIndices, kx, ky);

							double SAA = cosSumA*cosSumA + sinSumA * sinSumA;
							double SBB = cosSumB*cosSumB + sinSumB * sinSumB;
							double SAB = cosSumA*cosSumB + sinSumA * sinSumB;

							Results[TemperatureIndex].StructureFactors[i].SAA += SAA;
							Results[TemperatureIndex].StructureFactors[i].SBB += SBB;
							Results[TemperatureIndex].StructureFactors[i].SAB += SAB;
							Results[TemperatureIndex].StructureFactors[i].Scc += xB*xB*SAA + xA*xA*SBB - 2.0 * xA * xB * SAB;
							Results[TemperatureIndex].StructureFactors[i].NumberOfDataPoints++;
						}
						swap(nx,ny);
					}
				}
			}
		}

		void writeResultsToFile(string FileName) {
			for (int TemperatureIndex = 0; TemperatureIndex < Results.size(); TemperatureIndex++){
				for (int kIndex = 0; kIndex < kCombinationMapping.size(); kIndex++){
					Results[TemperatureIndex].StructureFactors[kIndex].SAA /= (static_cast<double>(Results[TemperatureIndex].StructureFactors[kIndex].NumberOfDataPoints) * static_cast<double>(TOTAL_NUMBER_OF_PARTICLES));
					Results[TemperatureIndex].StructureFactors[kIndex].SBB /= (static_cast<double>(Results[TemperatureIndex].StructureFactors[kIndex].NumberOfDataPoints) * static_cast<double>(TOTAL_NUMBER_OF_PARTICLES));
					Results[TemperatureIndex].StructureFactors[kIndex].SAB /= (static_cast<double>(Results[TemperatureIndex].StructureFactors[kIndex].NumberOfDataPoints) * static_cast<double>(TOTAL_NUMBER_OF_PARTICLES));
					Results[TemperatureIndex].StructureFactors[kIndex].Scc /= (static_cast<double>(Results[TemperatureIndex].StructureFactors[kIndex].NumberOfDataPoints) * static_cast<double>(TOTAL_NUMBER_OF_PARTICLES));
				}
				ofstream FileStreamToWriteTo;
				FileStreamToWriteTo.open(FileName+"_T="+to_string(Results[TemperatureIndex].Temperature)+".dat");
				FileStreamToWriteTo << "k\tAAStructureFactor\tBBStructureFactor\tABStructureFactor\tConcentrationStructureFactor\n";
				for (int i = 0; i < kCombinationMapping.size(); i++){
					FileStreamToWriteTo << setprecision(10) << kCombinationMapping[i].kMagnitude << '\t' << Results[TemperatureIndex].StructureFactors[i].SAA << '\t' << Results[TemperatureIndex].StructureFactors[i].SBB << '\t' << Results[TemperatureIndex].StructureFactors[i].SAB << '\t' << Results[TemperatureIndex].StructureFactors[i].Scc << '\n';
				}
				FileStreamToWriteTo.close();
			}
		}

		void addComputationResultsFromOtherComputator(const StructureFactorComputator& OtherComputator){
			if (Results.size() != OtherComputator.Results.size()){
				cerr << "WARNING: The Results of the computators have different sizes! Can't add them. Returning without doing anything." << endl;
				return;
			}
			for (int i = 0; i < Results.size(); i++){
				for (int j = 0; j < kCombinationMapping.size(); j++){
					Results[i].StructureFactors[j] += OtherComputator.Results[i].StructureFactors[j];
				}
			}
		}
};
