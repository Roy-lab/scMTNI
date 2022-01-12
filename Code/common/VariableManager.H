#ifndef _VARIABLE_MANAGER
#define _VARIABLE_MANAGER

#include <unordered_map>
#include <string>
#include <vector>

using namespace std;

class Error;
class Variable;
class VariableManager
{
	public:
		VariableManager();
		~VariableManager();
		Error::ErrorCode readVariables(const char*);
		int getVarID(const char*);
		bool isValid(int,int);
		vector<Variable*>& getVariableSet(); //map<int,Variable*>& getVariableSet();
		Variable* getVariableAt(int);

		Error::ErrorCode readVariablesFromTable(vector<string>&);
        Error::ErrorCode readVariablesFromData(string,int);
	private:
		unordered_map<string,int> varNameIDMap;
		vector<Variable*> variableSet;//map<int,Variable*> variableSet;
		//int globalVarId;
};
#endif
