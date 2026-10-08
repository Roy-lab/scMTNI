#ifndef _VARIABLE
#define _VARIABLE

#include <map>
#include <string>
#include <vector>
using namespace std;
#include "CommonTypes.H"

class Variable;

typedef vector<Variable*> VSET;
typedef vector<Variable*>::iterator VSET_ITER;

class Variable
{
	public:
		// This is how the variable will be used throughout the system
		void setID(int);
		int getID();

		// This is the name assigned to the variable by the user
		void setName(string name);
		const string& getName();

	private:
		string name;
		int vId;
};
#endif
