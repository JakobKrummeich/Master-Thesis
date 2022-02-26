#ifndef UTILITY_FUNCTIONS_INCLUDED
#define UTILITY_FUNCTIONS_INCLUDED

#include <string>
#include <fstream>

using namespace std;

void skipLines(ifstream& FileStreamToReadIn, int NumberOfLinesToSkip){
	string s;
	for (int i = 0; i < NumberOfLinesToSkip; i++){
		getline(FileStreamToReadIn, s, '\n');
	}
}

#endif //UTILITY_FUNCTIONS_INCLUDED
