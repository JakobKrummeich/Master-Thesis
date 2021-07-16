#ifndef UTILITY_FUNCTIONS_INCLUDED
#define UTILITY_FUNCTIONS_INCLUDED

#include <string>
#include <fstream>

void skipLines(ifstream& FileStreamToReadIn, int NumberOfLinesToSkip){
	string s;
	for (int i = 0; i < NumberOfLinesToSkip; i++){
		getline(FileStreamToReadIn, s, '\n');
	}
}

#endif //UTILITY_FUNCTIONS_INCLUDED
