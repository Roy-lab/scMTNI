#ifndef _VARIABLE_MANAGER
#define _VARIABLE_MANAGER

#include <string>
#include <vector>
#include <unordered_map>

using namespace std;

class Variable;

class VariableManager
{
	public:
		Variable* getVariable(string varName);
		Variable* getVariable(int varID);
		vector<Variable*>& getVariableSet();
        void addVariable(string varName, int varID);

	private:
		vector<Variable*> variables;
		unordered_map<int, Variable*> varMap;
};
#endif
