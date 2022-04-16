#ifndef INTEGRATION_INCLUDED
#define INTEGRATION_INCLUDED

#include <vector>
#include <fstream>
#include <limits>
#include <iomanip>
#include <cmath>
#include "value_pair.h"

using namespace std;

double computeIntegral(const vector<ValuePair>& Distribution) {
	double h = Distribution[1].xValue - Distribution[0].xValue;
	double Integral = 0.5 * (Distribution[0].yValue + Distribution.back().yValue);
	for (unsigned long i = 1; i < Distribution.size() - 1; i++){
		Integral += Distribution[i].yValue;
	}
	return h*Integral;
}

double computeFirstMoment(const vector<ValuePair>& Distribution) {
	vector<ValuePair> FirstMomentDistribution(Distribution);
	for (unsigned long i = 0; i < Distribution.size(); i++){
		FirstMomentDistribution[i].yValue = Distribution[i].yValue * Distribution[i].xValue;
	}
	return computeIntegral(FirstMomentDistribution);
}

double computeOrderParameter(const vector<ValuePair>& Distribution) {
	vector<ValuePair> orderParameterDistribution(Distribution);
	for (int i = 0; i < Distribution.size(); i++){
		orderParameterDistribution[i].yValue = Distribution[i].yValue * abs(Distribution[i].xValue - 0.5);
	}
	return computeIntegral(orderParameterDistribution);
}

double computeSecondCentralMoment(const vector<ValuePair>& Distribution, double Mean) {
	vector<ValuePair> SecondCentralMomentDistribution(Distribution);
	for (unsigned long i = 0; i < Distribution.size(); i++){
		double Difference = (Distribution[i].xValue - Mean);
		SecondCentralMomentDistribution[i].yValue = Distribution[i].yValue * Difference * Difference;
	}
	return computeIntegral(SecondCentralMomentDistribution);
}

double computeFourthCentralMoment(const vector<ValuePair>& Distribution, double Mean) {
	vector<ValuePair> FourthCentralMomentDistribution(Distribution);
	for (unsigned long i = 0; i < Distribution.size(); i++){
		double Difference = (Distribution[i].xValue - Mean);
		FourthCentralMomentDistribution[i].yValue = Distribution[i].yValue * Difference * Difference * Difference * Difference;
	}
	return computeIntegral(FourthCentralMomentDistribution);
}

void normalizeDistribution(vector<ValuePair>& Distribution) {
	double TotalIntegral = computeIntegral(Distribution);
	for (unsigned long i = 0; i < Distribution.size(); i++){
		Distribution[i].yValue /= TotalIntegral;
	}
}
		
#endif /* INTEGRATION_INCLUDED */
