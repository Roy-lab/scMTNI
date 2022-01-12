#ifndef _SLIM_FACTOR
#define _SLIM_FACTOR
#include <iostream>
#include<unordered_set>
#include<unordered_map>
using namespace std;

#include "CommonTypes.H"
class Variable;
class Potential;

class SlimFactor
{
	public:
		SlimFactor();
		SlimFactor(int);
		~SlimFactor();
		//This returns the ids of the maximal subsets of this factor. That is
		//all subsets with vCnt-1 factors
		//int generateMaximalSubsets(int**);
		//int getSetDiff(SlimFactor* aFactor,int* diffVIds,int& diffSize);
		//bool isMemberVariable(int v);
		//int genMBSubsets(int);
		//int thresholdToOne(double);
		//bool allEntriesInsignificant();
		//int showFactor(ostream&, VSET& variableSet, bool newLine=true);
		int setMBWts(unordered_map<int,double>&);

		int fId;
		//int* vIds;
		//int vCnt;
		//int secondPId;
		//double mutualInfo;
		//double jointEntropy;
		//double marginalEntropy;
		double mbScore;
		//double moveScore;
		//double marginalLL;
	//	double confidence;
		//Number of subsets that might be referring to this factor
		//int refCnt;
		//The collection of all Markov blankets that are good for this factor
		//INTDBLMAP goodMBIDs;
		//The merged Markov blanket 
		//INTINTMAP mergedMB;
        unordered_set<int> mergedMB; //unordered_map<int,int> mergedMB;
		//INTINTMAP tabulist;
		//INTDBLMAP candidateNeighbours;
		//vector<int> candidateNeighbours_vect;
		//vector<INTINTMAP*> mbSubsets;
		//Store the start indices for Markov blanket subsets of different sizes
		//map<int,int> mbSubsetStartInd;
		//Potential* canonicalParams;
		unordered_map<int,double> mbWts; //INTDBLMAP mbWts;
};
#endif
