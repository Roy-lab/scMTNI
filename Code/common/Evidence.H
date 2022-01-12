#ifndef _EVIDENCE
#define _EVIDENCE
//This is the evidence for a particular variable
//This can take the place of the likelihood

#include <fstream>
using namespace std;
#include "CommonTypes.H"

class Evidence
{
	public:
		Evidence();
		~Evidence();
		typedef enum
		{
			HARD,
			SOFT
		}EvidenceType;

		int assocVariable(int);
		int getAssocVariable();
		int setData(INTDBLMAP&);
		int addToData(int,double);
		INTDBLMAP& getData();
		int setType(EvidenceType);
		EvidenceType getType();
		int dumpEvidence(ostream&);
		int makeHidden();
		bool isHidden();
		//Get the value which is the most likely, i.e., the one with the max
		//value in the evidData
		int getMLVal();
		//This will prevent the loop over the values of the associated variable
		//in case of hard evidence
		int getHardEvidVal();
		int setEvidVal(double);
		double getEvidVal();
	private:
		INTDBLMAP evidData;
		int varID;
		EvidenceType eType;
		bool hiddenFlag;
		//Store the value of the hard evidence
		int hardValue;
		double contValue;
};

#endif
