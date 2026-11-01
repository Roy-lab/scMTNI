#ifndef _SLIM_FACTOR
#define _SLIM_FACTOR

#include <iostream>
#include<unordered_set>
#include<unordered_map>
#include "CommonTypes.H"

using namespace std;

class Variable;

class SlimFactor
{
	public:
		SlimFactor();
		~SlimFactor();

		int fId;
		double mbScore;
        unordered_set<int> mergedMB;
};
#endif
