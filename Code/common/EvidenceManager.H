#ifndef _EVIDENCE_MANAGER
#define _EVIDENCE_MANAGER

#include <fstream>
#include <map>
#include <string>
#include <vector>
#include "Matrix.H"
#include "gsl/gsl_randist.h"
using namespace std;

class Error;
class Evidence;

//EMAP stores the map of variable ids and their evidence. This corresponds
//to one line in the datafile
typedef map<int,Evidence*> EMAP;
typedef map<int,Evidence*>::iterator EMAP_ITER;

typedef vector<EMAP*> EVIDENCE_SET;

class VariableManager;
class EvidenceManager
{
	public:
		EvidenceManager();
		~EvidenceManager();
		//This associates the variable information with the variable id
		int setVariableManager(VariableManager*);
		//Read the evidence file
		Error::ErrorCode loadEvidenceFromFile(const char*);
		Error::ErrorCode loadEvidenceFromFile_Continuous(const char*);
		Error::ErrorCode loadEvidenceFromTable(vector<string>& inputLines);
        int deleteData();
		int normalize();
		int randomizeEvidence(gsl_rng*);
		int getNumberOfEvidences();
		//Get the evidence set at the specified id
		EMAP* getEvidenceAt(int);
		EMAP* getRandomEvidenceAt(int);
		int addEvidence(EMAP*);
		//This is used during EM, updating missing/hidden values using 
		//expected counts. eSetId corresponds to the data entry, vId is the
		//variable whose evidence we are updating and evidData is the set
		//of values and their evidence values (soft evidence)
		int addToEvidence(int eSetId,int vId, INTDBLMAP& evidData);	
		int dumpEvidenceSet(ostream&);
		int dumpSummaryStat(ostream&);
		int getMLSettings(ostream&);

		int setFoldCnt(int);
		int generateValidationSet(const char*,int vSetSize,gsl_rng* r);
		int setPreRandomizeSplit();
		int setPreRandomizeSplitSeed(int);
		int splitData(int);
		vector<int>& getTrainingSet();
        int getTrainSetSize();
        int getTestSetSize();
		//INTINTMAP& getTestSet();
        vector<int>& getTestSet();
		INTINTMAP& getValidationSet();
		int populateRandIntegers(gsl_rng*, INTINTMAP&,int,int);
        //added by shilu:
        Matrix* getData();
        VariableManager* getVariableManager();
    
	private:
		int populateEvidence(Evidence**,const char*);
		int populateEvidence_Continuous(Evidence**,const char*);
		int populateRandIntegers(gsl_rng*, int*,int);
		EVIDENCE_SET evidenceSet;
		EVIDENCE_SET randEvidenceSet;
		VariableManager* vMgr;
		//INTINTMAP trainIndex;
		//INTINTMAP testIndex;
        vector<int> testIndex;
        vector<int> trainIndex;
        int traindataSize;
        int testdataSize;
		INTINTMAP validationIndex;
		int foldCnt;
		INTDBLMAP maxAbsValues;
		bool preRandomizeSplit;
		int randseed;
        //added by shilu:
        Matrix* dataMat;
};
#endif
