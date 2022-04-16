#ifndef VALUE_PAIR_INCLUDED
#define VALUE_PAIR_INCLUDED

struct ValuePair {
	double xValue;
	double yValue;

	ValuePair()
	{
	}

	ValuePair(double xValue, double yValue):
		xValue(xValue),
		yValue(yValue){
	}
};

#endif /* VALUE_PAIR_INCLUDED */
