#include <iostream>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <vector>

//assumes following format: n lines comments/description, then data of the form kxm, i.e. k rows and m columns

using namespace std;

void readNewFile(string filePath, string& header, vector<double>& values, int numberOfHeaderLines) {

	header = "";
	values.clear();

	ifstream fS;
	fS.open(filePath);

	string currentString;
	for (int i = 0; i < numberOfHeaderLines; i++){
		getline(fS, currentString);
		header += currentString;
		header += "\n";
	}

	double entry;
	while (fS >> entry){
		values.push_back(entry);
	}
	fS.close();
}

using namespace std;

int main(int argc, char* argv[]){
	if (argc < 5){
		cerr << "combiner failed as we need an output-directory, input-filename, numberOfHeaderLines and numberOfColumns!" << endl;
		return 0;
	}

	string outputFilename = argv[1];
	outputFilename += "avg_combined_";
	outputFilename += argv[2];
	
	const int numberOfHeaderLines = stoi(argv[3]);
	const int numberOfColumns = stoi(argv[4]);

	vector<double> avgValues;
	string header;

	int numberOfFilesRead = 0;

	string inputFilename;
	while (cin >> inputFilename){
		cerr << "new file: " << inputFilename << endl;

		vector<double> newValues;

		readNewFile(inputFilename, header, newValues, numberOfHeaderLines);

		if (numberOfFilesRead == 0){
			avgValues = newValues;
		}
		else {
			for (int i = 0; i < newValues.size(); i++){
				avgValues[i] += newValues[i];
			}
		}
		numberOfFilesRead++;
	}
	for (int i = 0; i < avgValues.size(); i++){
		avgValues[i] /= static_cast<double>(numberOfFilesRead);
	}

	ofstream ofs;
	ofs.open(outputFilename);
	ofs << fixed << setprecision(numeric_limits<long double>::digits10+1);
	ofs << header;
	for (int i = 0, columnCounter = 0; i < avgValues.size(); i++){
		ofs << avgValues[i];
		char separator;
		if (columnCounter < numberOfColumns-1){
			separator = '\t';
			columnCounter++;
		}
		else {
			separator = '\n';
			columnCounter = 0;
		}
		ofs << separator;
	}
	ofs.close();
}

