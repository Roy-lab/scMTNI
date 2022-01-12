#ifndef _SPECIES_DATA_MANAGER_
#define _SPECIES_DATA_MANAGER_

#include <map>
#include <string>
using namespace std;

class VariableManager;
class FactorGraph;
class EvidenceManager;
class PotentialManager;
class SpeciesDataManager
{
	public:
	SpeciesDataManager();
	~SpeciesDataManager();
	int setVariableManager(VariableManager*);
	int createFactorGraph();
	int setEvidenceManager(EvidenceManager*);
	int setPotentialManager(PotentialManager*);
	int setOutputLoc(const char*);
	//int readRegulators(const char*);
	//int readTargets(const char*);
	int setMotifNetwork(const char*);
	VariableManager* getVariableManager();
	FactorGraph* getFactorGraph();
	EvidenceManager* getEvidenceManager();
	PotentialManager* getPotentialManager();
	const char* getOutputLoc();
	//map<string,int>& getRegulators();
	//map<string,int>& getTargets();
    unordered_map<int,unordered_map<int,double>*>& getMotifNetwork();
	private:
	VariableManager* varMgr;
	FactorGraph* fgraph;
	EvidenceManager* evMgr;
	PotentialManager* potMgr;
	char outputLoc[1024];
	//map<string,int> regulatorSet;
	//map<string,int> targetSet;
    unordered_map<int,unordered_map<int,double>*> motifNetwork;
};
#endif
