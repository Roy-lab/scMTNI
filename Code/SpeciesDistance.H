#ifndef _SPECIES_DATA_
#define _SPECIES_DATA_

#include <string>
#include <vector>
#include <unordered_map>
using namespace std;
class SpeciesDistance
{
	public:
	SpeciesDistance();
	~SpeciesDistance();
	int readSpeciesTree(const char*);
	double getEdgeStatusProb(vector<int>& edgeStatus);
	double getProbGainGain(string&);
	double getProbGainLoss(string&);
	double getProbLossGain(string&);
	double getProbLossLoss(string&);
    int setProot(double aVal);
    int setspeciesNameIDMap(unordered_map<string,int>& spNameIDMap);
	struct Species;

	struct Species 
	{
		string name;
		Species* parent;
        //update by Shilu
        //Species* leftchild;
        //Species* rightchild;
        vector<Species*> children;
        double p_maintain_edge;
		double p_gain;	
	};
	
	Species* getRoot();
	double getSubTree(bool,Species*,vector<int>&);
	private:
	unordered_map<string,Species*> speciesSet;
    unordered_map<string,int> speciesNameIDMap;
	unordered_map<int,double> binEdgeKeyProbMap;
	Species* root;
	double proot;
};
#endif
