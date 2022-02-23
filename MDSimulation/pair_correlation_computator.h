#ifndef PAIR_COMPUTATOR_INCLUDED
#define PAIR_COMPUTATOR_INCLUDED

#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <thread>

#include "../particle_type.h"
#include "particles.h"

using namespace std;

class PairCorrelationComputator {

	private:

		const double deltar; //dimensionless
		const double boxLength;
		const double maxLength;

		vector<double> img22AA;
		vector<double> img22BB;
		vector<double> img22AB;

		int numberOfValues;
		int numberOfAveragedPositions;

	public:

		PairCorrelationComputator(int numberOfValues, double maxLength, double boxLength): //maxLength and boxLength in units of sigma
			numberOfValues(numberOfValues),
			boxLength(boxLength),
			maxLength(maxLength),
			deltar(maxLength/(boxLength*static_cast<double>(numberOfValues))),
			numberOfAveragedPositions(0),
			img22AA(numberOfValues, 0.0),
			img22BB(numberOfValues, 0.0),
			img22AB(numberOfValues, 0.0)
		{
		}

		void computeImg22 (const Particles& P) {
			numberOfAveragedPositions++;
			const int numberOfAParticles = P.getNumberOfAParticles();
			const int numberOfBParticles = P.getNumberOfBParticles();
			const int totalNumber = numberOfAParticles + numberOfBParticles;

			for (int particleIndex = 0; particleIndex < totalNumber; particleIndex++){
				for (int otherParticleIndex = particleIndex+1; otherParticleIndex < totalNumber; otherParticleIndex++){
					double x0 = P.getPosition(particleIndex,0);
					double y0 = P.getPosition(particleIndex,1);

					double x1 = P.getPosition(otherParticleIndex,0);
					double y1 = P.getPosition(otherParticleIndex,1);

					double deltax = x1 - x0;
					double deltay = y1 - y0;

					if (deltay > 0.5){
						deltay -= 1.0;
						deltax -= P.getxDisplacement();
					}
					else if (deltay <= -0.5){
						deltay += 1.0;
						deltax += P.getxDisplacement();
					}
					while (deltax > 0.5){
						deltax -= 1.0;
					}
					while (deltax <= -0.5){
						deltax += 1.0;
					}
					double r = sqrt(deltax*deltax + deltay*deltay);
					if (r*boxLength < maxLength){

						double inverserCubed = 1.0/(r*r*r);
						int binIndex = static_cast<int>(r*boxLength*static_cast<double>(numberOfValues)/maxLength);

						ParticleType type0 = P.getParticleType(particleIndex);
						ParticleType type1 = P.getParticleType(otherParticleIndex);

						const double newValue = sqrt(7.5/M_PI) * (deltax * deltay) * inverserCubed / (deltar);

						if (type0 != type1){
							img22AB[binIndex] += newValue / static_cast<double>(numberOfAParticles*numberOfBParticles);
						}
						else if (type0 == ParticleType::A){
							img22AA[binIndex] += newValue / static_cast<double>(numberOfAParticles*numberOfAParticles);
						}
						else {
							img22BB[binIndex] += newValue / static_cast<double>(numberOfBParticles*numberOfBParticles);
						}
					}
				}
			}
		}

		void writeResults(string filePath) const {
			ofstream ofs(filePath+"Img22.dat");

			while (!ofs.is_open()){
				ofs.open(filePath+"Img22.dat");
				this_thread::sleep_for(std::chrono::milliseconds(100));
			}

			ofs << "r [sigma]\tImg22AA\tImg22BB\tImg22AB" << endl;
			double r = deltar*0.5;
			for (int i = 0; i < numberOfValues; i++){
				double inverseNumberOfAverages = 1.0/(static_cast<double>(numberOfAveragedPositions));
				ofs << r*boxLength << '\t' << img22AA[i]*inverseNumberOfAverages << '\t' << img22BB[i]*inverseNumberOfAverages << '\t' << img22AB[i]*inverseNumberOfAverages << '\n';
				r += deltar;
			}
		}
};

#endif //PAIR_COMPUTATOR_INCLUDED
