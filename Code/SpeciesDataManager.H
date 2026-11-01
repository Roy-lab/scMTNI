#ifndef _SPECIES_DATA_MANAGER_
#define _SPECIES_DATA_MANAGER_

#include <map>
#include <string>
using namespace std;

class VariableManager;
class PotentialManager;

class SpeciesDataManager
{
	public:
		~SpeciesDataManager();

		void setPotentialManager(PotentialManager*);
		void setOutputLoc(const char*);
		void setMotifNetwork(const char*);

		VariableManager* getVariableManager();
		PotentialManager* getPotentialManager();
		const char* getOutputLoc();
		unordered_map<int,unordered_map<int,double>*>& getMotifNetwork();
		SlimFactor* getFactor(int varID);

	private:
		PotentialManager* potMgr;
		char outputLoc[1024];
		unordered_map<int,unordered_map<int,double>*> motifNetwork;
		unordered_map<int, SlimFactor*> factors;
};
#endif
