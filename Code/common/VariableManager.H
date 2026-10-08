#ifndef _VARIABLE_MANAGER
#define _VARIABLE_MANAGER

#include <string>
#include <vector>

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
};
#endif
