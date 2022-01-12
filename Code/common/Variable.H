#ifndef _VARIABLE
#define _VARIABLE

#include <map>
#include <string>
#include <vector>
using namespace std;
#include "CommonTypes.H"

class Variable;

//typedef map<int,Variable*> VSET;
//typedef map<int,Variable*>::iterator VSET_ITER;
typedef vector<Variable*> VSET;
typedef vector<Variable*>::iterator VSET_ITER;

class Variable
{
	public:
		Variable();
		~Variable();
		//This is how the variable will be used throughout
		//the system
		int setID(int);
		int getID();
		//This is the name assigned to the variable by the user
		int setName(const char*);
		const string& getName();
		
		int setNextValue(int);
		int setValues(INTVECT&);
		INTVECT& getValues();
		int getValueCnt();
		bool isValidValue(int);
		int initEvidence(INTDBLMAP&);
	private:
		//The set of values this variable can take
		INTVECT valueSet;
		string name;
		int vId;
};
#endif
