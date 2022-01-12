#ifndef _META_LEARNER_
#define _META_LEARNER_

#include <unordered_map>
#include <unordered_set>
#include <map>
#include <string>
#include <vector>
#include "CommonTypes.H"
#include <utility>
using namespace std;


class VariableManager;
class EvidenceManager;
class PotentialManager;
class MetaMove;
class SpeciesDistance;
class MappedOrthogroupReader;
class MetaLearner
{
	public:
	MetaLearner();
	~MetaLearner();
	int setInputFName(const char*);
	int setMaxFactorSize(int);
	int setMaxFactorSize_Approx(int);
	int setMBCntForApproxGreedy(int);
	int setPenalty(double);
	int setConvergenceThreshold(double);
	int setRestrictedList(const char*);
	int setLambdaLearnSize(int);
	int setLambdaLearn();
	int setTrueGraph(const char* aGraphName);
	int setScoreOpt_Entropy();
	int setScoreOpt_PLL();
	int setInputOGList(const char*);
        
	int setOrthogroupReader(MappedOrthogroupReader*);
	int setSpeciesDistances(SpeciesDistance*);
	int setBeta1(double);
	int setBeta2(double);
    int setEdgePrior_PerSpecies();
	int init();
	int initOld();
	int doCrossValidation(int);
    int doOneFold();
	int start(int);
	int showModelParameters(int);
	int setPreRandomizeSplit();
	int setPreRandomizeSplitSeed(int seed);
    //added by Shilu
    int init_onefold();
    int setINDEP();
    int setsplitGenes();
    void process_mem_usage(double& vm_usage, double& resident_set);
	private:
		
	double getInitScore();
	double getScore();
	double getPriorDelta();
	int clearFoldSpecData();
	int initEdgeSet();
    int initEdgeSet_onefold();
	double precomputeEmptyGraphPrior();
	double precomputePerSpeciesPrior(int specID, int targetID, Variable* target, SpeciesDataManager* sdm, vector<vector<double>> &edgePresenceProb,map<int,MappedOrthogroup*>& orthogroupSet);
    //int precomputePerSpeciesPrior();
	int initCondsetMap();
	int initCondsetMap_Nopool();
	INTINTMAP* getConditionSet(int cind);
    int setinputVariableNames();
	int collectMoves(int);
	int collectMoves_Orthogroups(int);
    int initCondsetMap_TreeAll(SpeciesDistance::Species* node);
    int collectMoves_Orthogroups_INDEP(int);
	int collectMoves_Deletions();
	int sortMoves();
	double makeMoves(int& successMove);

	int makeDelMoves();
    int attemptMove(MetaMove*); //int attemptMove(MetaMove*,map<int,INTINTMAP*>&);
	int attemptDelMove(MetaMove* move,map<int,INTINTMAP*>& affectedVars);
	//int getNewPLLScore(int,INTINTMAP& conditionSet,Variable* u, Variable* v,double&,double&,int,INTDBLMAP&);
    int getNewPLLScore(int,Variable* u, Variable* v,double&,double&,int);//unordered_map<int,double>&);
	double getPLLScore_Condition(int,SlimFactor* sFactor);
    //double getPLLScore_Condition_onefold(string& condKey,SlimFactor* sFactor);
	double getPLLScore_Condition(string& condKey,SlimFactor* sFactor,int&,INTDBLMAP&);
	//double getPLLScore_Condition_Tracetrick(string& condKey,SlimFactor* sFactor,int&,INTDBLMAP&);
    double getPLLScore_Condition_Tracetrick(int,SlimFactor* sFactor,int&);//,unordered_map<int,double>&);
	double getValidationPLLScore_Condition(int,int);
	double getEdgePrior(int specID,Variable* u, Variable* v,int orthoGrpNo);
	double getEdgePrior(int specID,Variable* u, Variable* v,int orthoGrpNo,int);
	double getEdgePrior_PerSpecies(int tfID, int targetID,SpeciesDataManager*);
	int genCondSetKey(INTINTMAP& condSet, string& aKey);
	int cleanUp();
	int dumpAllGraphs(int,int);
	bool checkMBSize(int,int,int);

	int generateData(int sampleCnt, int burnin,int foldid);
	int generateInitSample(INTDBLMAP& initialSample,gsl_rng* rndgen,map<int,Potential*>*);
	int writeEvidenceTo(ofstream* oFile,INTDBLMAP& sample);
	int initCondsetMap_Tree(SpeciesDistance::Species* node);
    vector<unordered_map<int,int> > initCondsetMap_Tree_backtrack(SpeciesDistance::Species* node, vector<unordered_map<int,int> > mycondition, int parentstatus, int ntransition);
	int getLeaves(SpeciesDistance::Species* node, map<string,int>& leaves);

	void readEvidenceTable(string fileName);

    vector<string> speciesIDNameMap; //map<int,string> speciesIDNameMap;
	unordered_map<string,int> speciesNameIDMap;
	vector<SpeciesDataManager*> speciesDataSet;
    unordered_set<int> intermediatecells;
    //map<string,SpeciesDataManager*> speciesDataSet;
	map<string,PotentialManager*> pooledPotentials;
	map<string,EvidenceManager*> pooledData;
	//map<string,int> pooledKeysForDel;
	int maxFactorSize;
	int maxFactorSizeApprox;
	int maxMBCntGreedy;
	//char inputFName[1024];
	string inputFName;
	double penalty;
	double convThreshold;
	char restrictedFName[1024];
	int lambdaLearnSize;
	bool learnLambda;

	//map<string,INTINTMAP*> edgeConditionMap;
    unordered_map<string,vector<int>*> affectedOGPairs; //map<pair<int,int>,STRINTMAP*> affectedOGPairs;
	vector<vector<double> > ogpairPrior; //map<string,double> ogpairPrior;
	FactorGraph* fg;
	vector<MetaMove*> moveSet;
	char sharedEdgesFName[1024];
	map<int,INTINTMAP*> condsetMap;
    vector<vector<int>> condsetMap0;
    vector<vector<int>> condsetMap_Tree;
    //map<int,INTINTMAP*> condsetMap_Tree;
	//map<string,int> condsetKeyIDMap;
	//map<int,string> condsetIDKeyMap;
	//map<int,int> restrictedOrthoList;
	bool dataPool;
	bool usePLL;
	bool useEntropy;
	map<int,INTDBLMAP*> validationPLLs;
	char trueGraphFName[1024];
	map<int,double> finalScores;
	map<string,int> subgraphVarSet;
	map<string,int> bottomClip;
	map<string,int> topClip;
	int totalFoldCnt;
	bool preRandomizeSplit;
	int preRandSeed;
	double globalVar;
	MappedOrthogroupReader* ogr;
	SpeciesDistance* speciesData;
	char orthoMapFName[1024];
	char speciesOrder[1024];
	vector<int> inputOGList;  //map<int,int> inputOGList;
    vector<int>inputRegulatorOGs; //map<int,int> inputRegulatorOGs;
    unordered_set<int> inputVariables;
    unordered_map<int,unordered_set<string>> inputVariableNames;
    vector<vector<double>> varNeighborhoodPrior_PerSpecies;
    vector<vector<vector<double>>> edgePresenceProb_PerSpecies;
    //map<string,map<int,double>*> varNeighborhoodPrior_PerSpecies0;//This stores the empty graph
	//map<string,map<string,double>*> edgePresenceProb_PerSpecies0;//This stores the edge probability
	double beta1;
	double beta2;
	vector<string> inputTable;
	//unordered_map<string,unordered_map<string,double>*> likelihood;
	bool INDEP;
    bool splitGenes;
    double initGlobalScore;
    double edgePresenceProb;
};
#endif
