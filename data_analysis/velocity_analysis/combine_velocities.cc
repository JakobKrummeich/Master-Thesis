#include <iostream>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <vector>

using namespace std;

void readNewVelocityFile(string filePath, vector<double>& yValues, vector<double>& xVelocities, vector<double>& yVelocities) {

	ifstream fS;
	fS.open(filePath);

	string currentString;
	getline(fS, currentString);

	istringstream iss(currentString);
	double entry;
	while (iss >> entry){
		yValues.push_back(entry);
	}
	iss.clear();

	getline(fS, currentString);
	iss.str(currentString);
	while (iss >> entry){
		xVelocities.push_back(entry);
	}
	iss.clear();
	getline(fS, currentString);

	iss.str(currentString);
	while (iss >> entry){
		yVelocities.push_back(entry);
	}
	fS.close();
}

using namespace std;

int main(int argc, char* argv[]){
	if (argc < 2){
		cerr << "velocity-combiner failed as we need an output-directory!" << endl;
		return 0;
	}

	string outputFilename = argv[1];
	outputFilename += "avg_combined_velocities.dat";

	vector<double> globalyValues;
	vector<double> avgxVelocities;
	vector<double> avgyVelocities;

	int numberOfFilesRead = 0;
	int numberOfColumns = 0;

	string inputFilenameVelocity;
	while (cin >> inputFilenameVelocity){
		cerr << "new velocity file: " << inputFilenameVelocity << endl;
		vector<double> yValues;
		vector<double> xVelocities;
		vector<double> yVelocities;
		readNewVelocityFile(inputFilenameVelocity, yValues, xVelocities, yVelocities);

		if (numberOfFilesRead == 0){
			globalyValues = yValues;
			avgxVelocities = xVelocities;
			avgyVelocities = yVelocities;
			numberOfColumns = xVelocities.size();
		}
		else {
			for (int i = 0; i < numberOfColumns; i++){
				avgxVelocities[i] += xVelocities[i];
				avgyVelocities[i] += yVelocities[i];
			}
		}
		numberOfFilesRead ++;
	}
	for (int i = 0; i < numberOfColumns; i++){
		avgxVelocities[i] /= static_cast<double>(numberOfFilesRead);
		avgyVelocities[i] /= static_cast<double>(numberOfFilesRead);
	}

	ofstream ofs;
	ofs.open(outputFilename);
	ofs << fixed << setprecision(numeric_limits<long double>::digits10+1);
	for (int i = 0; i < numberOfColumns; i++){
		ofs << globalyValues[i] << '\t';
	}
	ofs << '\n';
	for (int i = 0; i < numberOfColumns; i++){
		ofs << avgxVelocities[i] << '\t';
	}
	ofs << '\n';
	for (int i = 0; i < numberOfColumns; i++){
		ofs << avgyVelocities[i] << '\t';
	}
	ofs << '\n';
	ofs.close();
}

