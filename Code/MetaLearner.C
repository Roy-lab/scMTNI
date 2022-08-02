#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include <utility>
#include "Error.H"
#include "Variable.H"
#include "VariableManager.H"

#include "Evidence.H"
#include "EvidenceManager.H"

#include "Potential.H"
#include "SlimFactor.H"
#include "PotentialManager.H"

#include "FactorGraph.H"
#include "MetaMove.H"

#include "Utils.H"

#include "SpeciesDistance.H"
#include "SpeciesDataManager.H"
#include "GeneMap.H"
#include "MappedOrthogroup.H"
#include "MappedOrthogroupReader.H"
#include "MetaLearner.H"
#include <chrono>
#include <unistd.h>
using namespace std::chrono;
using namespace std;

MetaLearner::MetaLearner()
{
    restrictedFName[0]='\0';
    lambdaLearnSize=0;
    learnLambda=false;
    dataPool=false;
    usePLL=true;
    useEntropy=true;
    trueGraphFName[0]='\0';
    preRandomizeSplit=false;
    preRandSeed=-1;
    convThreshold=1e-4;
    beta1=-0.9;
    beta2=4.0;
    INDEP=false;
    initGlobalScore=0;
    splitGenes=false;
}

MetaLearner::~MetaLearner()
{
    inputVariables.clear();
    inputOGList.clear();
    inputRegulatorOGs.clear();
}

int
MetaLearner::setInputFName(const char* aFName)
{
    //strcpy(inputFName,aFName);
    inputFName = aFName;
    return 0;
}

int
MetaLearner::setMaxFactorSize(int aVal)
{
    maxFactorSize=aVal;
    return 0;
}

int
MetaLearner::setMaxFactorSize_Approx(int aVal)
{
    maxFactorSizeApprox=aVal;
    return 0;
}

int
MetaLearner::setMBCntForApproxGreedy(int aCnt)
{
    maxMBCntGreedy=aCnt;
    return 0;
}

int
MetaLearner::setPenalty(double aVal)
{
    penalty=aVal;
    return 0;
}

int
MetaLearner::setINDEP()
{
    INDEP=true;
    return 0;
}

int
MetaLearner::setsplitGenes()
{
    splitGenes=true;
    return 0;
}

int
MetaLearner::setConvergenceThreshold(double aVal)
{
    convThreshold=aVal;
    return 0;
}

int
MetaLearner::setRestrictedList(const char* aFName)
{
    strcpy(restrictedFName,aFName);
    ifstream inFile(restrictedFName);
    string buffer;
    while(inFile.good())
    {
        getline(inFile,buffer);
        if(buffer.length()<=0)
        {
            continue;
        }
        int ogid=atoi(buffer.c_str());
        inputRegulatorOGs.push_back(ogid); //inputRegulatorOGs[ogid]=0;
        inputVariables.insert(ogid);
    }
    inFile.close();
    /*for(int i=0;i<inputRegulatorOGs.size();i++){
        cout <<"inputRegulatorOGs[" << i << "]=" << inputRegulatorOGs[i] << endl;
    }*/
    return 0;
}

int
MetaLearner::setLambdaLearnSize(int vSize)
{
    lambdaLearnSize=vSize;
    return 0;
}

int
MetaLearner::setLambdaLearn()
{
    learnLambda=true;
    return 0;
}

int
MetaLearner::setTrueGraph(const char* aGraphName)
{
    strcpy(trueGraphFName,aGraphName);
    return 0;
}


int
MetaLearner::setPreRandomizeSplit()
{
    preRandomizeSplit=true;
    return 0;
}

int
MetaLearner::setPreRandomizeSplitSeed(int seed)
{
    preRandSeed=seed;
    return 0;
}


int
MetaLearner::setScoreOpt_Entropy()
{
    usePLL=false;
    useEntropy=true;
    return 0;
}

int
MetaLearner::setScoreOpt_PLL()
{
    usePLL=true;
    useEntropy=false;
    return 0;
}


int
MetaLearner::setInputOGList(const char* aFName)
{
    ifstream inFile(aFName);
    char buffer[1024];
    while(inFile.good())
    {
        inFile.getline(buffer,1023);
        if(strlen(buffer)<=0)
        {
            continue;
        }
        int og=atoi(buffer);
        inputOGList.push_back(og); //inputOGList[og]=0;
        inputVariables.insert(og);
    }
    inFile.close();
    /*for(int i=0;i<inputOGList.size();i++){
        cout <<"inputOGList[" << i << "]=" << inputOGList[i] << endl;
    }*/
    return 0;
}

int
MetaLearner::setOrthogroupReader(MappedOrthogroupReader* aPtr)
{
    ogr=aPtr;
    return 0;
}

int
MetaLearner::setSpeciesDistances(SpeciesDistance* aPtr)
{
    speciesData=aPtr;
    speciesData->setspeciesNameIDMap(speciesNameIDMap);
    return 0;
}


int
MetaLearner::setBeta1(double b1)
{
    beta1 = b1;
    return 0;
}

/*int
MetaLearner::setEdgePrior_PerSpecies()
{
    edgePresenceProb=1/(1+exp(-1*beta1));
    if(edgePresenceProb<1e-50){
        edgePresenceProb=1e-50;
    }else if(edgePresenceProb>1-1e-50){
        edgePresenceProb=1-1e-50;
    }
    cout <<"Edge Specific Prior=" << edgePresenceProb << endl;
    return 0;
}*/


int
MetaLearner::setBeta2(double b2)
{
    beta2 = b2;
    return 0;
}
 
void
MetaLearner::process_mem_usage(double& vm_usage, double& resident_set)
{
    vm_usage     = 0.0;
    resident_set = 0.0;

    // the two fields we want
    unsigned long vsize;
    long rss;
    {
        std::string ignore;
        std::ifstream ifs("/proc/self/stat", std::ios_base::in);
        ifs >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
                >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore >> ignore
                >> ignore >> ignore >> vsize >> rss;
    }

    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    vm_usage = vsize/1024.0/1024.0;  //MB
    resident_set = rss * page_size_kb;
}



int
MetaLearner::init_onefold()
{
    ifstream inFile(inputFName);
    cout << "MetaLearner::init_onefold() read:" << inputFName << endl;
    string buffer;
    //char buffer[1024];
    int datasetId = 0;
    if(splitGenes){
        setinputVariableNames();
    }
    // set sparsity prior: edgePresenceProb
    //setEdgePrior_PerSpecies();
    map<int,MappedOrthogroup*>& orthogroupSet=ogr->getMappedOrthogroups();
    while (inFile.good())
    {
        getline(inFile, buffer);
        if (buffer.empty() || buffer.find("#") == 0)
        {
            continue;
        }
        vector<string> strs;
        
        // strip trailing newlines
        buffer.erase(buffer.find_last_not_of(" \n\r\t") + 1);
        
        strs = Utils::split(buffer, '\t');
        assert(strs.size() >= 3);
        string specName = strs[0];
        string datasetSuff = strs[1];
        string outputLoc = strs[2];
        //string regulatorFName = strs[3];
        //string targetFName = strs[4];
        string motifNetwork = strs[5];

        PotentialManager* potMgr = new PotentialManager;
        /*EvidenceManager* evMgr = new EvidenceManager;
        if (preRandomizeSplit)
        {
            evMgr->setPreRandomizeSplit();
            evMgr->setPreRandomizeSplitSeed(preRandSeed);
        }
        potMgr->setEvidenceManager(evMgr);*/
        string tableFileName = datasetSuff+ ".table";
        readEvidenceTable(tableFileName);
        Error::ErrorCode eCode;
        if(splitGenes){
            eCode=potMgr->loadEvidenceFromTable(inputTable,inputVariableNames[datasetId]);
        }else{
            eCode=potMgr->loadEvidenceFromTable(inputTable);
        }
        if(eCode!=Error::SUCCESS)
        {
            cout <<Error::getErrorString(eCode) << endl;
            return -1;
        }
        VariableManager* varMgr=potMgr->getVariableManager();//evMgr->getVariableManager();
        potMgr->setOutputDir(outputLoc.c_str());
        SpeciesDataManager* spMgr=new SpeciesDataManager;
        speciesDataSet.push_back(spMgr); //speciesDataSet[specName]=spMgr;
        spMgr->setVariableManager(varMgr);
        //spMgr->setEvidenceManager(evMgr);
        spMgr->setPotentialManager(potMgr);
        spMgr->createFactorGraph();
        spMgr->setOutputLoc(outputLoc.c_str());
        spMgr->setMotifNetwork(motifNetwork.c_str()); //add motif prior network in CVN for single cell
        speciesIDNameMap.push_back(specName); //speciesIDNameMap[datasetId]=specName;
        speciesNameIDMap[specName]=datasetId;
        cout << datasetId << "=" << specName << " motifNetwork=" << motifNetwork<<endl;
        
        //move one fold cross-validation here:
        potMgr->reset();
        potMgr->init();
        char foldOutputDirCmd[1024];
        sprintf(foldOutputDirCmd,"mkdir -p %s/fold0",outputLoc.c_str());
        system(foldOutputDirCmd);
        FactorGraph* condspecGraph=spMgr->getFactorGraph();
        //initialize factor graph with pll score:
        //sparsity prior for each target:
        vector<double>varNeighborhoodPrior;
        vector<Variable*>& varSet=varMgr->getVariableSet();
        vector<vector<double>> edgePresenceProb(condspecGraph->getFactorCnt(),vector<double> (condspecGraph->getFactorCnt(),0));
        for(int f=0;f<condspecGraph->getFactorCnt();f++)
        {
            SlimFactor* sFactor=condspecGraph->getFactorAt(f);
            //double pll=getPLLScore_Condition_onefold((string&)eIter->first,sFactor);
            double pll=potMgr->computeMeanVarPseudoLikelihood_onefold(sFactor->fId);
            Variable* target=varSet[sFactor->fId];
            double priorScore=precomputePerSpeciesPrior(datasetId,sFactor->fId,target,spMgr,edgePresenceProb,orthogroupSet);
            varNeighborhoodPrior.push_back(priorScore);
            //cout << "cell=" << datasetId <<" f=" << f << " sFactor->fId=" << sFactor->fId<<" pll=" << pll << " priorScore=" << priorScore << endl;
            sFactor->mbScore=pll+priorScore;
            initGlobalScore=initGlobalScore+pll+priorScore; //old getInitScore();
            //sFactor->marginalLL=pll;
        }
        varNeighborhoodPrior_PerSpecies.push_back(varNeighborhoodPrior);
        edgePresenceProb_PerSpecies.push_back(edgePresenceProb);
        potMgr->deleteData();
        datasetId+=1; //=datasetId*2;
        strs.clear();
        cout <<"variable size=" << condspecGraph->getFactorCnt() << endl;
    }
    inFile.close();
    inputTable.clear();
    cout << "-------------------------------------------------------------" << endl;
    return 0;
}

int
MetaLearner::init()
{
    ifstream inFile(inputFName);
    cout << "MetaLearner::init() read:" << inputFName << endl;
    string buffer;
    //char buffer[1024];
    int datasetId = 0;
    while (inFile.good())
    {
        getline(inFile, buffer);
        if (buffer.empty() || buffer.find("#") == 0)
        {
            continue;
        }
        /*
         inFile.getline(buffer,1023);
         if(strlen(buffer)<=0)
         {
         continue;
         }
         if(strstr(buffer,"#")!=NULL)
         {
         continue;
         }
         */
        vector<string> strs;
        
        // strip trailing newlines
        buffer.erase(buffer.find_last_not_of(" \n\r\t") + 1);
        
        strs = Utils::split(buffer, '\t');
        assert(strs.size() == 6);
        string specName = strs[0];
        string datasetSuff = strs[1];
        string outputLoc = strs[2];
        //string regulatorFName = strs[3];
        //string targetFName = strs[4];
        string motifNetwork = strs[5];
        
        PotentialManager* potMgr = new PotentialManager;
        EvidenceManager* evMgr = new EvidenceManager;
        if (preRandomizeSplit)
        {
            evMgr->setPreRandomizeSplit();
            evMgr->setPreRandomizeSplitSeed(preRandSeed);
        }
        potMgr->setEvidenceManager(evMgr);
        string tableFileName = datasetSuff; //+ ".table";
        readEvidenceTable(tableFileName);
        
        //merge varMgr->readVariablesFromTable with evMgr->loadEvidenceFromTable(inputTable);
        /*VariableManager* varMgr = new VariableManager;
        Error::ErrorCode eCode = varMgr->readVariablesFromTable(inputTable);
        if(eCode != Error::SUCCESS)
        {
            cout << Error::getErrorString(eCode) << endl;
            return -1;
        }
        evMgr->setVariableManager(varMgr);*/
        Error::ErrorCode eCode=evMgr->loadEvidenceFromTable(inputTable);
        if(eCode!=Error::SUCCESS)
        {
            cout <<Error::getErrorString(eCode) << endl;
            return -1;
        }
        VariableManager* varMgr=evMgr->getVariableManager();
        potMgr->setOutputDir(outputLoc.c_str());
        SpeciesDataManager* spMgr=new SpeciesDataManager;
        speciesDataSet.push_back(spMgr); //speciesDataSet[specName]=spMgr;
        spMgr->setVariableManager(varMgr);
        spMgr->setEvidenceManager(evMgr);
        spMgr->setPotentialManager(potMgr);
        spMgr->createFactorGraph();
        spMgr->setOutputLoc(outputLoc.c_str());
        //spMgr->readRegulators(regulatorFName.c_str());  //comment by Shilu July 5
        //spMgr->readTargets(targetFName.c_str());
        spMgr->setMotifNetwork(motifNetwork.c_str());  //add motif network in CVN for single cell
        speciesIDNameMap.push_back(specName); //speciesIDNameMap[datasetId]=specName;
        speciesNameIDMap[specName]=datasetId;
        //unordered_map<string,double>* likelihoodPerSpecies=new unordered_map<string,double>;
        //likelihood[specName]=likelihoodPerSpecies;
        cout << datasetId << "=" << specName << endl;
        datasetId+=1; //=datasetId*2;
        strs.clear();
    }
    inFile.close();
    inputTable.clear();
    cout << "-------------------------------------------------------------" << endl;
    return 0;
}

void
MetaLearner::readEvidenceTable(string fileName)
{
    /* Reads gene names and expression levels from a tab-separated file, where the first row is assumed (for now)
     to be headers. In subsequent rows, the first column is a gene name and remaining columns are expression levels.
     This method just loads the meaningful lines into a vector that is then passed to VariableManager::readVariablesFromTable()
     and EvidenceManager::loadEvidenceFromTable(). It's done this way so as not to break encapsulation of private members
     of VariableManager and EvidenceManager.
     */
    cout << "readEvidenceTable:" << fileName << endl;
    ifstream inFile(fileName);
    assert(inFile);
    
    bool headerLine = true;
    string inputLine;
    inputTable.clear();
    
    // Read all the lines with meaningful data
    while (getline(inFile, inputLine))
    {
        if (inputLine.empty() || inputLine.find('#') == 0)
        {
            continue;
        }
        if (headerLine)
        {
            headerLine = false;
            continue;
        }
        inputTable.push_back(inputLine);
    }
}

int
MetaLearner::doCrossValidation(int foldCnt)
{
    totalFoldCnt=foldCnt;
    /*gsl_rng* r=gsl_rng_alloc(gsl_rng_default);
    for(map<string,SpeciesDataManager*>::iterator eIter=speciesDataSet.begin();eIter!=speciesDataSet.end();eIter++)
    {
        EvidenceManager* evMgr=eIter->second->getEvidenceManager();
        evMgr->setFoldCnt(foldCnt);
        evMgr->splitData(0);
    }
    gsl_rng_free(r);*/
    //The first key is for the fold number
    //For each fold we have a trained model. For each trained model we have the likelihood on
    //all the test sets, including the self test.

    for(int f=0;f<foldCnt;f++)
    {
        //for(map<string,SpeciesDataManager*>::iterator eIter=speciesDataSet.begin();eIter!=speciesDataSet.end();eIter++)
        for(int i=0;i<speciesDataSet.size();i++){
            EvidenceManager* evMgr=speciesDataSet[i]->getEvidenceManager();
            evMgr->setFoldCnt(foldCnt);
            evMgr->splitData(f);
            PotentialManager* potMgr=speciesDataSet[i]->getPotentialManager();
            potMgr->reset();
            potMgr->init(f);
            char outputDir[1024];
            sprintf(outputDir,"%s/fold%d",speciesDataSet[i]->getOutputLoc(),f);
            char foldOutputDirCmd[1024];
            sprintf(foldOutputDirCmd,"mkdir -p %s",outputDir);
            system(foldOutputDirCmd);
        }
        //clearFoldSpecData();
        initEdgeSet(); //need to double check for cross-validation with more folds (moved from start(f));
        initGlobalScore=getInitScore();
        start(f);
        //Now get all the PLL scores on the validation sets
        //for(map<string,SpeciesDataManager*>::iterator eIter=speciesDataSet.begin();eIter!=speciesDataSet.end();eIter++)
        for(int i=0;i<speciesDataSet.size();i++){
            double pll=0;
            VariableManager* varMgr=speciesDataSet[i]->getVariableManager();
            vector<Variable*>& varSet=varMgr->getVariableSet();
            int specID=i;//speciesNameIDMap[eIter->first];
            EvidenceManager* evMgr=speciesDataSet[i]->getEvidenceManager();
            INTINTMAP& vSet=evMgr->getValidationSet();
            if(vSet.size()==0){
                cout << "No validation data!" << endl;
                break;
            }
            for(int j=0;j<varSet.size();j++)
            {
                pll=pll+getValidationPLLScore_Condition(specID,j);
            }
            INTDBLMAP* pllSet=NULL;
            if(validationPLLs.find(specID)==validationPLLs.end())
            {
                pllSet=new INTDBLMAP;
                validationPLLs[specID]=pllSet;
            }
            else
            {
                pllSet=validationPLLs[specID];
            }
            (*pllSet)[f]=pll;
        }
    }

    char scoreFName[1024];
    sprintf(scoreFName,"%s/scoreFile.txt",speciesDataSet[0]->getOutputLoc()); //speciesDataSet.begin()->second->getOutputLoc()
    ofstream sFile(scoreFName);
    for(map<int,double>::iterator dIter=finalScores.begin();dIter!=finalScores.end();dIter++)
    {
        sFile<< dIter->second << endl;
    }
    sFile.close();
    return 0;
}

int
MetaLearner::doOneFold()
{
    totalFoldCnt=1;
    start(0);
    char scoreFName[1024];
    sprintf(scoreFName,"%s/scoreFile.txt",speciesDataSet[0]->getOutputLoc()); //speciesDataSet.begin()->second->getOutputLoc()
    ofstream sFile(scoreFName);
    for(map<int,double>::iterator dIter=finalScores.begin();dIter!=finalScores.end();dIter++)
    {
        sFile<< dIter->second << endl;
    }
    sFile.close();
    return 0;
}

int
MetaLearner::start(int f)
{
    auto start = high_resolution_clock::now();
    cout << "MetaLearner::start" << endl;
    //Repeat until convergence
    double currGlobalScore=initGlobalScore;//currGlobalScore=-1; currGlobalScore=getInitScore();
    int maxMBSizeApprox=maxFactorSizeApprox;//maxFactorSizeApprox-1; by shilu
    //int currK=1;
    int currK=maxMBSizeApprox;
    double initPrior=0;
    if(!INDEP){
        cout << "precomputeEmptyGraphPrior" << endl;
        initPrior=precomputeEmptyGraphPrior();
        initCondsetMap_Tree(speciesData->getRoot());
        cout <<"condsetMap_Tree" << endl;
        for(int setIter=0;setIter<condsetMap_Tree.size();setIter++)
        {
            //INTINTMAP* cset=setIter->second;
            vector<int> & cset=condsetMap_Tree[setIter];
            cout << setIter;
            for(auto cIter=0;cIter<cset.size();cIter++)
            {
                cout << " "<<speciesIDNameMap[cIter]<<"=" << cset[cIter];
            }
            cout << endl;
        }
    }
    //precomputePerSpeciesPrior();
    //initEdgeSet();
    //initEdgeSet_onefold();
    cout << "inputRegulatorOGs.size() = " << inputRegulatorOGs.size() << " inputOGList.size() = " << inputOGList.size() << endl;

    if(strlen(trueGraphFName)==0)
    {
        //while(currK<=maxMBSizeApprox)
        //{
        int iter = 0; // vperiyasamy added
        bool notConverged=true;
        while(notConverged && iter < 100) // hardcode 100 iterations for now
        {
            //cout << "ITERATION " << iter << endl;
            //collect the candidate edges
            if(!INDEP){
                collectMoves_Orthogroups(currK);
            }else{
                collectMoves_Orthogroups_INDEP(currK);
            }
            //sortMoves();  //not needed
            int successMove=0;
            double diff=makeMoves(successMove);  //diff=sum of move->getScoreImprovement(), which is equal to getScore()-currGlobalScore;
            /*if(currGlobalScore==-1)
            {
                currGlobalScore=getInitScore();
                //	currGlobalScore=currGlobalScore+initPrior;
            }*/
            double priorChange=0;
            if(!INDEP){
                priorChange=getPriorDelta();
            }
            //collectMoves_Deletions();
            //sortMoves();
            //makeDelMoves();
            //priorChange=getPriorDelta();
            double newScore=currGlobalScore+diff;//getScore();
            //newScore=newScore+initPrior+priorChange;
            //double diff=newScore-currGlobalScore;
            if(diff<=convThreshold)
            {
                notConverged=false;
            }
            //dumpAllGraphs(currK,f);
            currGlobalScore=newScore;
            cout << "ITERATION " << iter <<" newScore=" << newScore << " diffscore=" <<diff  <<" priorChange=" << priorChange << " successMove=" <<successMove <<endl;
            iter++;
        }
        //currK++;
        //}
        dumpAllGraphs(currK,f);
        cout <<"Final Score " << currGlobalScore << endl;
        finalScores[f]=currGlobalScore;
        showModelParameters(f);
    }
    else
    {
        //populateGraphsFromFile();
        double globalScore=getScore();
        generateData(10000,30000,f);
        showModelParameters(f);
        //cout <<"Final Score " << globalScore<< endl;
        cout <<"Final Score\tTOTALFOLD" << totalFoldCnt <<"_fold"<< f<<"\t" << globalScore<< endl;
    }
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);    
    cout << "MetaLearner::start runtime: " <<duration.count() << " ms" <<endl;
    return 0;
}


double
MetaLearner::getInitScore()
{
    cout <<"getInitScore" << endl;
    double initScore=0;
    //for(map<string,SpeciesDataManager*>::iterator gIter=speciesDataSet.begin();gIter!=speciesDataSet.end();gIter++)
    for(int i=0;i<speciesDataSet.size();i++){
        FactorGraph* fg=speciesDataSet[i]->getFactorGraph();//gIter->second->getFactorGraph();
        vector<SlimFactor*>& factorSet=fg->getAllFactors();
        for(int fIter=0;fIter<factorSet.size();fIter++)
        {
            SlimFactor* sFactor=factorSet[fIter];
            initScore=initScore+sFactor->mbScore;//same as sFactor->marginalLL;
        }
    }
    return initScore;
}

double
MetaLearner::getScore()
{
    double gScore=0;
    
    //for(map<string,SpeciesDataManager*>::iterator gIter=speciesDataSet.begin();gIter!=speciesDataSet.end();gIter++)
    for(int i=0;i<speciesDataSet.size();i++){
        FactorGraph* fg=speciesDataSet[i]->getFactorGraph();
        vector<SlimFactor*>& factorSet=fg->getAllFactors();
        for(int fIter=0;fIter<factorSet.size();fIter++)
        {
            SlimFactor* sFactor=factorSet[fIter];
            gScore=gScore+sFactor->mbScore;
        }
    }
    return gScore;
}
//Shilu: update ogpairPrior, edgeConditionMap not needed
double
MetaLearner::getPriorDelta()
{
    double oldStructPrior=0;
    double newStructPrior=0;
    //cout << "MetaLearner::getPriorDelta() Update the OldpriorScore" << endl;
    //Need to consider the old contribution of the edges, delete that from the overall prior and add the new contribution
    for(auto edgeIter=affectedOGPairs.begin();edgeIter!=affectedOGPairs.end();edgeIter++)
    {
        vector<string> keyid = Utils::split(edgeIter->first,'-');
        int regi=stoi(keyid[0]);//edgeIter->first.first;
        int targeti=stoi(keyid[1]); //edgeIter->first.second;
        double aval=speciesData->getEdgeStatusProb(*(edgeIter->second));
        double edgePrior=log(aval);
        double oldEdgePrior=ogpairPrior[regi][targeti];  //oldEdgePrior=ogpairPrior[edgeIter->first];
        oldStructPrior+=oldEdgePrior;
        newStructPrior+=edgePrior;
        /*INTINTMAP* currEdgeStatus=edgeConditionMap[edgeIter->first];
        STRINTMAP* newEdgeStatus=edgeIter->second;
        for(STRINTMAP_ITER sIter=newEdgeStatus->begin();sIter!=newEdgeStatus->end();sIter++)
        {
            int specID=speciesNameIDMap[sIter->first];
            (*currEdgeStatus)[specID]=sIter->second;  //edgeConditionMap[edgeIter->first][specID]=sIter->second
        }*/
        ogpairPrior[regi][targeti]=edgePrior; //ogpairPrior[edgeIter->first]=edgePrior;
        //cout <<" ogpairPrior[" << regi<<"][" << targeti << "]=" <<ogpairPrior[regi][targeti] << endl;
        keyid.clear();
    }
    for(auto edgeIter=affectedOGPairs.begin();edgeIter!=affectedOGPairs.end();edgeIter++)
    {
        edgeIter->second->clear();
        delete edgeIter->second;
    }
    affectedOGPairs.clear();
    double priorDelta=newStructPrior-oldStructPrior;
    return priorDelta;
}

int
MetaLearner::clearFoldSpecData()
{
    //Clear existing graphs
    //for(map<string,SpeciesDataManager*>::iterator fIter=speciesDataSet.begin();fIter!=speciesDataSet.end();fIter++)
    for(int i=0;i<speciesDataSet.size();i++){
        delete speciesDataSet[i]; //fIter->second;
    }
    speciesDataSet.clear();
    return 0;
}

int
MetaLearner::initEdgeSet()
{
    
    if(condsetMap.size()==0)
    {
        /*if(dataPool)
        {
            initCondsetMap();
        }
        else
        {*/
        initCondsetMap_Nopool();
        //}
    }
    //initCondsetMap_Tree(speciesData->getRoot());
    
    //for(map<string,SpeciesDataManager*>::iterator eIter=speciesDataSet.begin();eIter!=speciesDataSet.end();eIter++)
    for(int i=0;i<speciesDataSet.size();i++){ //i is specID
        FactorGraph* condspecGraph=speciesDataSet[i]->getFactorGraph();
        //string& species=speciesIDNameMap[i];
        //int specID=speciesNameIDMap[eIter->first];
        //INTINTMAP* cset=getConditionSet(specID);
        //string condKey;
        //genCondSetKey(*cset,condKey);
        //map<int,double>* varNeighborhoodPrior=varNeighborhoodPrior_PerSpecies[eIter->first];
        for(int f=0;f<condspecGraph->getFactorCnt();f++)
        {
            SlimFactor* sFactor=condspecGraph->getFactorAt(f);
            /*if(sFactor->fId==111)
            {
                cout <<"found the target of interest" << endl;
            }*/
            double pll=getPLLScore_Condition(i,sFactor);
            //double priorScore=(*varNeighborhoodPrior)[sFactor->fId];
            sFactor->mbScore=pll;//pll+priorScore;
            //sFactor->marginalLL=pll;//pll+priorScore;
        }
    }
    return 0;
}

int
MetaLearner::initEdgeSet_onefold()
{
    cout <<"initEdgeSet_onefold" << endl;
    initCondsetMap_Nopool();
    initCondsetMap_Tree(speciesData->getRoot());
    cout <<"condsetMap_Tree" << endl;
    for(int setIter=0;setIter<condsetMap_Tree.size();setIter++)
    {
        //INTINTMAP* cset=setIter->second;
        vector<int> & cset=condsetMap_Tree[setIter];
        cout << setIter;
        for(auto cIter=0;cIter<cset.size();cIter++)
        {
            cout << " "<<cIter<<"=" << cset[cIter];
        }
        cout << endl;
    }
    //for(map<string,SpeciesDataManager*>::iterator eIter=speciesDataSet.begin();eIter!=speciesDataSet.end();eIter++)
    for(int i=0;i<speciesDataSet.size();i++){
        SpeciesDataManager* spd=speciesDataSet[i];
        FactorGraph* condspecGraph=spd->getFactorGraph();
        PotentialManager* potMgr=spd->getPotentialManager();
        //int specID=speciesNameIDMap[eIter->first];
        for(int f=0;f<condspecGraph->getFactorCnt();f++)
        {
            SlimFactor* sFactor=condspecGraph->getFactorAt(f);
            //double pll=getPLLScore_Condition_onefold((string&)eIter->first,sFactor);
            //VariableManager* varMgr=spd->getVariableManager();
            //VSET& varSet=varMgr->getVariableSet();
            double pll=potMgr->computeMeanVarPseudoLikelihood_onefold(sFactor->fId);
            sFactor->mbScore=pll;//pll+priorScore;
            //sFactor->marginalLL=pll;//pll+priorScore;
        }
        //EvidenceManager* evMgr=spd->getEvidenceManager();
        //evMgr->deleteData();
        //potMgr->deleteData();
    }
    return 0;
}


// Update by shilu: compute ogpairPrior, edgeConditionMap not needed
double
MetaLearner::precomputeEmptyGraphPrior()
{
    //map<string,int> ogpairEdgeCntMap;
    //map<int,int> ogMaxSpecCnt_TF;
    //map<int,int> ogMaxSpecCnt_Tgt;
    vector<int> edgeStatus(speciesIDNameMap.size(),0);
    /*for(int specIter=0;specIter<speciesIDNameMap.size();specIter++)
    {
        edgeStatus[speciesIDNameMap[specIter]]=0;
    }*/
    double prior=speciesData->getEdgeStatusProb(edgeStatus);
    double logPrior=log(prior);
    int edgeCnt=0;
    for(int regi=0;regi<inputRegulatorOGs.size();regi++)
    {
        int ogno_tf=inputRegulatorOGs[regi];
        vector<double> targetPrior;
        for(int targeti=0;targeti<inputOGList.size();targeti++)
        {
            int ogno_tgt=inputOGList[targeti];
            if(ogno_tf==ogno_tgt)
            {
                targetPrior.push_back(0);
                continue;
            }
            /*char key[1024];
            sprintf(key,"%d-%d",ogno_tf,ogno_tgt);
            string keystr(key);
            INTINTMAP* edgeSpecAssign=new INTINTMAP;
            for(map<string,int>::iterator tIter=speciesNameIDMap.begin();tIter!=speciesNameIDMap.end();tIter++)
            {
                (*edgeSpecAssign)[tIter->second]=0;
            }
            edgeConditionMap[keystr]=edgeSpecAssign;
            ogpairPrior[keystr]=logPrior;*/
            targetPrior.push_back(logPrior); //ogpairPrior[regi][targeti]=logPrior;
            edgeCnt++;
            //cout << "regi=" << regi << " targeti=" << targeti << " ogpairPrior[][]=" << logPrior << endl;
        }
        ogpairPrior.push_back(targetPrior);
    }
    double emptyGraphPrior=edgeCnt*logPrior;
    return emptyGraphPrior;
}

/*
// compute ogpairPrior
double
MetaLearner::precomputeEmptyGraphPrior()
{
    map<string,int> ogpairEdgeCntMap;
    map<int,int> ogMaxSpecCnt_TF;
    map<int,int> ogMaxSpecCnt_Tgt;
    STRINTMAP edgeStatus;
    for(map<string,SpeciesDataManager*>::iterator sIter=speciesDataSet.begin();sIter!=speciesDataSet.end();sIter++)
    {
        SpeciesDataManager* sdm=sIter->second;
        VariableManager* vMgr=sdm->getVariableManager();
        VSET& varSet=vMgr->getVariableSet();
        map<string,int>& regulators=sdm->getRegulators();
        map<string,int>& targets=sdm->getTargets();
        for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
        {
            Variable* v=vIter->second;
            if(targets.size()>0 && (targets.find(v->getName())==targets.end()))
            {
                continue;
            }
            int ogno_tgt=ogr->getMappedOrthogroupID(v->getName().c_str(),sIter->first.c_str());
            if(inputOGList.find(ogno_tgt) == inputOGList.end())
            {
                continue;
            }
            if(ogno_tgt==-1)
            {
            //cout <<"No Orthogroup for " << v->getName() << " species " << sIter->first << endl;
                continue;
            }
            MappedOrthogroup* ogrp_Tgt=ogr->getMappedOrthogroup(v->getName().c_str(),sIter->first.c_str());
            if(ogrp_Tgt==NULL)
            {
                cout <<"No target by name " << v->getName() << " in " << sIter->first << endl;
            }
            GeneMap* gmap=ogrp_Tgt->getSpeciesHits(sIter->first.c_str());
            map<string,map<string,STRINTMAP*>*>& targetsInSpec=gmap->getGeneSet();
            int targetsWithVals=0;
            for(map<string,map<string,STRINTMAP*>*>::iterator tIter=targetsInSpec.begin();tIter!=targetsInSpec.end();tIter++)
            {
                if(vMgr->getVarID(tIter->first.c_str())==-1)
                {
                    continue;
                }
                targetsWithVals++;
            }
            if(ogMaxSpecCnt_Tgt.find(ogno_tgt)==ogMaxSpecCnt_Tgt.end())
            {
                ogMaxSpecCnt_Tgt[ogno_tgt]=targetsWithVals;
            }
            else
            {
                int currval=ogMaxSpecCnt_Tgt[ogno_tgt];
                if(currval<targetsWithVals)
                {
                    ogMaxSpecCnt_Tgt[ogno_tgt]=targetsWithVals;
                }
            }
            //Now we need to figure out how many edges this ogno_tgt-ogno_tf pair is going to contribute to
            //This is essentially the max number of active targets
            for(map<string,int>::iterator rIter=regulators.begin();rIter!=regulators.end();rIter++)
            {
                int regID=vMgr->getVarID(rIter->first.c_str());
                if(regID==-1)
                {
                    continue;
                }
                int ogno_tf=ogr->getMappedOrthogroupID(rIter->first.c_str(),sIter->first.c_str());
                if(ogno_tf==-1)
                {
                    cout <<"No Orthogroup for " <<rIter->first << " species " << sIter->first << endl;
                    continue;
                }
                MappedOrthogroup* ogrp_TF=ogr->getMappedOrthogroup(rIter->first.c_str(),sIter->first.c_str());
                if(ogrp_TF==NULL)
                {
                    cout <<"No gene " << rIter->first << " in " << sIter->first << endl;
                    continue;
                }
                char key[1024];
                sprintf(key,"%d-%d",ogno_tf,ogno_tgt);
                string keystr(key);
                if(edgeConditionMap.find(keystr)==edgeConditionMap.end())
                {
                    INTINTMAP* edgeSpecAssign=new INTINTMAP;
                    for(map<string,int>::iterator tIter=speciesNameIDMap.begin();tIter!=speciesNameIDMap.end();tIter++)
                    {
                        (*edgeSpecAssign)[tIter->second]=0;
                    }
                    edgeConditionMap[keystr]=edgeSpecAssign;
                }
                GeneMap* tfgenemap=ogrp_TF->getSpeciesHits(sIter->first.c_str());
                map<string,map<string,STRINTMAP*>*>& tfsInSpec=tfgenemap->getGeneSet();
                int tfsWithVals=0;
                for(map<string,map<string,STRINTMAP*>*>::iterator tIter=tfsInSpec.begin();tIter!=tfsInSpec.end();tIter++)
                {
                    if(vMgr->getVarID(tIter->first.c_str())==-1)
                    {
                        continue;
                    }
                    tfsWithVals++;
                }
                if(ogMaxSpecCnt_TF.find(ogno_tf)==ogMaxSpecCnt_TF.end())
                {
                    ogMaxSpecCnt_TF[ogno_tf]=tfsWithVals;
                }
                else
                {
                    int currval=ogMaxSpecCnt_TF[ogno_tf];
                    if(currval<tfsWithVals)
                    {
                        ogMaxSpecCnt_TF[ogno_tf]=tfsWithVals;
                    }
                }
            }

        }
        edgeStatus[sIter->first]=0;
    }
    double prior=speciesData->getEdgeStatusProb(edgeStatus);
    double logPrior=log(prior);
    int edgeCnt=0;
    cout << "Compute ogpairPrior! ogMaxSpecCnt_TF.size=" << ogMaxSpecCnt_TF.size() << " ogMaxSpecCnt_Tgt.size="<< ogMaxSpecCnt_Tgt.size() << endl;
    for(map<int,int>::iterator tfIter=ogMaxSpecCnt_TF.begin();tfIter!=ogMaxSpecCnt_TF.end();tfIter++)
    {
        for(map<int,int>::iterator tgtIter=ogMaxSpecCnt_Tgt.begin();tgtIter!=ogMaxSpecCnt_Tgt.end();tgtIter++)
        {
            edgeCnt=edgeCnt+(tfIter->second*tgtIter->second);
            char key[1024];
            sprintf(key,"%d-%d",tfIter->first,tgtIter->first);
            string keystr(key);
            ogpairPrior[keystr]=logPrior;
            cout << "ogpairPrior["<< keystr <<"]=" << logPrior <<endl;
        //cout <<"ogMaxSpecCnt_TF[" <<tfIter->first << "]=" << tfIter->second << " ogMaxSpecCnt_Tgt[" <<tgtIter->first << "]=" <<tgtIter->second << endl;
        }
    }
    double emptyGraphPrior=edgeCnt*logPrior;
    return emptyGraphPrior;
}*/


// Shilu: compute varNeighborhoodPrior_PerSpecies[specID][targeti]: set to 0
// edgePresenceProb is the same for CVN for each edge per cell, so no need to store into a map
//targetID is variable ID
double
MetaLearner::precomputePerSpeciesPrior(int specID, int targetID, Variable* target, SpeciesDataManager* sdm, vector<vector<double>> &edgePresenceProb, map<int,MappedOrthogroup*>& orthogroupSet)
{
    //double initPrior=edgePresenceProb;
    string spec=speciesIDNameMap[specID];
    int ogno_tgt=ogr->getMappedOrthogroupID(target->getName().c_str(),spec.c_str());
    double NeighborhoodPrior=0;
    for(int regi=0;regi<inputRegulatorOGs.size();regi++)
    {
        int ogno_tf=inputRegulatorOGs[regi];
        if(ogno_tf==ogno_tgt)
        {
            continue;
        }
        MappedOrthogroup* tfogrp=orthogroupSet[ogno_tf]; //oIter->second; four species
        vector<string>& tfgrpMembers=tfogrp->getOrthoMembers();
        VariableManager* vMgr=sdm->getVariableManager();
        int regID=vMgr->getVarID(tfgrpMembers[specID].c_str());
        double initPrior=getEdgePrior_PerSpecies(regID,targetID,sdm);
        edgePresenceProb[regID][targetID]=initPrior;
        NeighborhoodPrior+=log(1-initPrior);
        //cout << spec <<" regi=" << regi << " TF_OGID=" <<ogno_tf << " regVarID=" <<regID << " targetID=" <<targetID << " initPrior=" << initPrior << endl;
    }
    //cout << spec << "target=" << target->getName() <<" NeighborhoodPrior=" << NeighborhoodPrior << endl;
    return NeighborhoodPrior;
}

/*int
MetaLearner::precomputePerSpeciesPrior()
{
    for(map<string,SpeciesDataManager*>::iterator sIter=speciesDataSet.begin();sIter!=speciesDataSet.end();sIter++)
    {
        SpeciesDataManager* sdm=sIter->second;
        VariableManager* vMgr=sdm->getVariableManager();
        VSET& varSet=vMgr->getVariableSet();
        map<string,int>& regulators=sdm->getRegulators();
        map<string,int>& targets=sdm->getTargets();
        map<string,double>* edgePresenceProb=new map<string,double>;
        edgePresenceProb_PerSpecies[sIter->first]=edgePresenceProb;
        map<int,double>* varNeighborhoodPrior=new map<int,double>;
        varNeighborhoodPrior_PerSpecies[sIter->first]=varNeighborhoodPrior;

        for(VSET_ITER vIter=varSet.begin();vIter!=varSet.end();vIter++)
        {
            Variable* v=vIter->second;
            if(targets.size()>0 && (targets.find(v->getName())==targets.end()))
            {
                continue;
            }
            //Now we need to figure out how many edges this ogno_tgt-ogno_tf pair is going to contribute to
            //This is essentially the max number of active targets
            for(map<string,int>::iterator rIter=regulators.begin();rIter!=regulators.end();rIter++)
            {
                int regID=vMgr->getVarID(rIter->first.c_str());
                if(regID==-1)
                {
                    continue;
                }
                string edgeKey;
                edgeKey.append(rIter->first.c_str());
                edgeKey.append("\t");
                edgeKey.append(v->getName().c_str());
                double initPrior=getEdgePrior_PerSpecies(regID,vIter->first,sdm);
                (*edgePresenceProb)[edgeKey]=initPrior;
                if(varNeighborhoodPrior->find(vIter->first)==varNeighborhoodPrior->end())
                {
                    (*varNeighborhoodPrior)[vIter->first]=log(1-initPrior);
                }
                else
                {
                    (*varNeighborhoodPrior)[vIter->first]=(*varNeighborhoodPrior)[vIter->first]+log(1-initPrior);
                }
            }
        }
    }
    
    return 0;
}*/

/*
int
MetaLearner::initCondsetMap()
{
    int lastcondind=(int)pow(2,speciesDataSet.size());
    for(int i=1;i<lastcondind;i++)
    {
        INTINTMAP* cset=new INTINTMAP;
        int currind=i;
        for(map<int,string>::iterator eIter=speciesIDNameMap.begin();eIter!=speciesIDNameMap.end();eIter++)
        {
            if(currind==0)
            {
                (*cset)[eIter->first]=0;
                continue;
            }
            if((currind%2)==1)
            {
                (*cset)[eIter->first]=1;
            }
            else
            {
                (*cset)[eIter->first]=0;
            }
            currind=currind/2;
        }
        condsetMap[i]=cset;
        string condKey;
        genCondSetKey(*cset,condKey);
        condsetKeyIDMap[condKey]=i;
        condsetIDKeyMap[i]=condKey;
    }
    return 0;
}*/

/*
// update by Shilu Oct4 2019: consider all combinations, no restrictions.
int
MetaLearner::initCondsetMap_TreeAll(SpeciesDistance::Species* node)
{
    int npos=pow(2,speciesIDNameMap.size());
    vector<vector<int> > mycondition;
    for(int j=0;j<speciesIDNameMap.size();j++)  //(map<string,int> speciesNameIDMap
    {
        mycondition.resize(pow(2,j+1));
        if(j==0)
        {
            mycondition[0].push_back(0);
            mycondition[1].push_back(1);
        }else{
            vector<vector<int> > temp;
            for (int r = 0; r < mycondition.size(); r++) {
                vector<int> temp1(mycondition[r]);
                temp1.push_back(0);
                temp.push_back(temp1);
                vector<int> temp2(mycondition[r]);
                temp2.push_back(1);
                temp.push_back(temp2);
            }
            for (int r = 0; r < mycondition.size(); r++) {
                mycondition[r]=temp[r];
            }
        }
    }
c
    return 0;
}*/

/*
int
MetaLearner::initCondsetMap_Tree(SpeciesDistance::Species* node)
{
    //reach leaves
    if(node->children.empty())
    {
        vector<int> cset;
        //map<int,int>* cset=new map<int,int>;
        int specID=speciesNameIDMap[node->name];
        //for(map<string,int>::iterator sIter=speciesNameIDMap.begin();sIter!=speciesNameIDMap.end();sIter++)
        for(int j=0;j<speciesIDNameMap.size();j++)
        {
            if(j==specID)
            {
                cset.push_back(1); //(*cset)[j]=1;
            }
            else
            {
                cset.push_back(0); //(*cset)[j]=0;
            }
        }
        //int id=condsetMap_Tree.size();
        condsetMap_Tree.push_back(cset);//condsetMap_Tree[id]=cset;
        
        vector<int> cset1;
        //cset=new map<int,int>;
        //for(map<string,int>::iterator sIter=speciesNameIDMap.begin();sIter!=speciesNameIDMap.end();sIter++)
        for(int j=0;j<speciesIDNameMap.size();j++)
        {
            if(j==specID)
            {
                cset1.push_back(0); //(*cset)[j]=0;
            }
            else
            {
                cset1.push_back(1); //(*cset)[j]=1;
            }
        }
        //id=condsetMap_Tree.size();
        condsetMap_Tree.push_back(cset1); //condsetMap_Tree[id]=cset;
        
        return 0;
    }
    
    //Otherwise get me all children under me
    map<string,int> leaves;
    getLeaves(node,leaves);
    cout << "Descendants of " << node->name << " : ";
    for(map<string,int>::iterator lIter=leaves.begin();lIter!=leaves.end();lIter++)
    {
        cout << lIter->first  << " ";
    }
    cout << endl;
    
    // if it is root, first set everything to 1
    if(node->parent==NULL)
    {
        //map<int,int>* cset=new map<int,int>;
        vector<int> cset;
        //for(map<string,int>::iterator sIter=speciesNameIDMap.begin();sIter!=speciesNameIDMap.end();sIter++)
        for(int j=0;j<speciesIDNameMap.size();j++)
        {
            
            cset.push_back(1); //(*cset)[j]=1;
        }
        //int id=condsetMap_Tree.size();
        //condsetMap_Tree[id]=cset;
        condsetMap_Tree.push_back(cset);
    }
    if(leaves.size()>1)
    {
        //Get me all children under me and set their coordinates to 1
        //map<int,int>* cset=new map<int,int>;
        cout<<"leaves.size()>1" << endl;
        vector<int> cset;
        //for(map<string,int>::iterator sIter=speciesNameIDMap.begin();sIter!=speciesNameIDMap.end();sIter++)
        for(int j=0;j<speciesIDNameMap.size();j++)
        {
            string specname=speciesIDNameMap[j];
            if(leaves.find(specname)==leaves.end())
            {
                cset.push_back(0); //(*cset)[j]=0;
            }
            else
            {
                cset.push_back(1); //(*cset)[j]=1;
            }
            cout <<specname<<"="<<cset[j]<<endl;
        }
        //int id=condsetMap_Tree.size();
        //condsetMap_Tree[id]=cset;
        condsetMap_Tree.push_back(cset);
        
        //Get me all children under me and set their coordinates to 0
        //cset=new map<int,int>;
        vector<int> cset1;
        //for(map<string,int>::iterator sIter=speciesNameIDMap.begin();sIter!=speciesNameIDMap.end();sIter++)
        for(int j=0;j<speciesIDNameMap.size();j++)
        {
            string specname=speciesIDNameMap[j];
            if(leaves.find(specname)==leaves.end())
            {
                cset1.push_back(1); //(*cset)[j]=1;
            }
            else
            {
                cset1.push_back(0); //(*cset)[j]=0;
            }
        }
        //id=condsetMap_Tree.size();
        //condsetMap_Tree[id]=cset;
        condsetMap_Tree.push_back(cset1);
    }
    
    leaves.clear();
    if((!node->children.empty()) && node->parent!=NULL)
    {
        cout << "intermediate state: " <<node->name << endl;
        vector<int> cset;
        for(int j=0;j<speciesIDNameMap.size();j++)
        {
            string specname=speciesIDNameMap[j];
            if(node->name==specname){
                cset.push_back(1);
                cout << specname << "=1\t" ;
            }else{
                cset.push_back(0);
                cout << specname << "=0\t" ;
            }
    }
        cout << endl;
        condsetMap_Tree.push_back(cset);
    }
    if(!node->children.empty())
    {
        for(int i=0;i<node->children.size();i++)
        {
            initCondsetMap_Tree(node->children[i]);
        }
    }
    return 0;
}*/

int
MetaLearner::initCondsetMap_Tree(SpeciesDistance::Species* node){
    //cout <<"rootnode=" << node->name << endl;
    // this is root, it has 0/1 two options:
    int npos=pow(2,speciesIDNameMap.size());
    vector<unordered_map<int,int> > mycondition,mycondition1;
    int id=speciesNameIDMap[node->name];
    unordered_map<int,int> map0,map1;
    map0[id]=0;
    mycondition.push_back(map0);
    map1[id]=1;
    mycondition1.push_back(map1);
    map0.clear();
    map1.clear();
    if(!node->children.empty())
    {
        for(int i=0;i<node->children.size();i++)
        {
            mycondition=initCondsetMap_Tree_backtrack(node->children[i],mycondition,0,0);
            mycondition1=initCondsetMap_Tree_backtrack(node->children[i],mycondition1,1,0);
        }
    }
    mycondition.insert(mycondition.end(), mycondition1.begin(), mycondition1.end());

    //copy mycondition to condsetMap_Tree
    //remove all 0s combination from mycondition:
    for(int i=0;i<mycondition.size();i++)
    {
        //cout<<"condition " << i << ": ";
        vector<int> cset(speciesIDNameMap.size(),0);
        int zerocount=0;
        for(auto sIter=mycondition[i].begin();sIter!=mycondition[i].end();sIter++)
        {
            cset[sIter->first]=sIter->second; //cell id: sIter->first, edge status: sIter->second
            //cout << speciesIDNameMap[sIter->first] << "=" <<sIter->second <<" ";
            zerocount+=1-sIter->second; //count number of 0s in the condition
        }
        if(zerocount<speciesIDNameMap.size()){
            condsetMap_Tree.push_back(cset);
        }
        mycondition[i].clear();
        cset.clear();
        //cout << endl;
    }
    // add cell-specific condition for intermediate cells (already included)
    for(auto it=intermediatecells.begin();it!=intermediatecells.end();it++)
    {
        string species=speciesIDNameMap[*it];
        //cout << "add cell-specific condition for intermediate cell: " <<species << endl;
        // only in this cell:
        vector<int> cset1;//cset0;
        for(int j=0;j<speciesIDNameMap.size();j++)
        {
            string specname=speciesIDNameMap[j];
            if(j==*it){
                //cset0.push_back(0);
                cset1.push_back(1);
            }else{
                //cset0.push_back(1);
                cset1.push_back(0);
            }
        }
        //condsetMap_Tree.push_back(cset0);
        condsetMap_Tree.push_back(cset1);
        //cset0.clear();
        cset1.clear();
    }
    
    mycondition.clear();
    for(int i=0;i<mycondition1.size();i++)
    {
        mycondition1[i].clear();
    }
    mycondition1.clear();
    return 0;
}
//shilu: use backtracking to set up constrained conditions:
vector<unordered_map<int,int> >
MetaLearner::initCondsetMap_Tree_backtrack(SpeciesDistance::Species* node, vector<unordered_map<int,int> > mycondition, int parentstatus, int ntransition)
{
    
    int id=speciesNameIDMap[node->name];
    //cout <<"node=" << node->name << " id=" << id <<" parentstatus=" <<parentstatus<<" ntransition=" << ntransition << endl;
    // intermediate cell:
    if((!node->children.empty()) && node->parent!=NULL)
    {
        //cout << "intermediate cell: " <<node->name << endl;
        intermediatecells.insert(id);
        if(ntransition==0){
            //not reach max transition yet in this branch, it has 0/1 two options:
            vector<unordered_map<int,int> > outcondition, cond1;
            for (int r = 0; r < mycondition.size(); r++) {
                unordered_map<int,int> temp1(mycondition[r]);
                temp1[id]=parentstatus;
                outcondition.push_back(temp1);
                unordered_map<int,int> temp2(mycondition[r]);
                temp2[id]=1-parentstatus;
                cond1.push_back(temp2);
                temp1.clear();
                temp2.clear();
            }
            if(!node->children.empty())
            {
                for(int i=0;i<node->children.size();i++)
                {
                    outcondition=initCondsetMap_Tree_backtrack(node->children[i],outcondition,parentstatus,ntransition);
                    cond1=initCondsetMap_Tree_backtrack(node->children[i],cond1,1-parentstatus,ntransition+1);
                }
            }
            outcondition.insert(outcondition.end(), cond1.begin(), cond1.end());
            
            for (int r = 0; r < mycondition.size(); r++)
            {
                mycondition[r].clear();
                cond1[r].clear();
            }
            cond1.clear();
            mycondition.clear();
            return outcondition;
        }else if(ntransition==1){
            //transition has occurred, it needs to be the same as parent:
            for (int r = 0; r < mycondition.size(); r++) {
                mycondition[r][id]=parentstatus;
            }
            if(!node->children.empty())
            {
                for(int i=0;i<node->children.size();i++)
                {
                    mycondition=initCondsetMap_Tree_backtrack(node->children[i],mycondition,parentstatus,ntransition);
                }
            }
            return mycondition;
        }
            
    }
    
    //reach leaves
    if(node->children.empty())
    {
        if(ntransition==0){
            //not reach max transition yet in this branch, it has 0/1 two options:
            vector<unordered_map<int,int> > outcondition;
            for (int r = 0; r < mycondition.size(); r++) {
                unordered_map<int,int> temp1(mycondition[r]);
                temp1[id]=parentstatus;
                outcondition.push_back(temp1);
                unordered_map<int,int> temp2(mycondition[r]);
                temp2[id]=1-parentstatus;
                outcondition.push_back(temp2);
                temp1.clear();
                temp2.clear();
                mycondition[r].clear();
            }
            mycondition.clear();
            return outcondition;
        }else if(ntransition==1){
            //transit twice, it needs to be the same as parent:
            for (int r = 0; r < mycondition.size(); r++) {
                mycondition[r][id]=parentstatus;
            }
            return mycondition;
        }
    }
}



// update by Shilu get all its children, not just leaves!
int
MetaLearner::getLeaves(SpeciesDistance::Species* node, map<string,int>& leaves)
{
    if(node->children.empty())
    {
        leaves[node->name]=0;
        return 0;
    }
    if(!node->children.empty())
    {
        for(int i=0;i<node->children.size();i++)
        {
            leaves[node->children[i]->name]=0;
            getLeaves(node->children[i],leaves);
        }
    }
    return 0;
}

// fill in condsetMap
int
MetaLearner::initCondsetMap_Nopool()
{
    //int ind=0;
    cout << "MetaLearner::initCondsetMap_Nopool" << endl;
    for(int ind=0;ind<speciesIDNameMap.size();ind++)
    {
        INTINTMAP* cset=new INTINTMAP;
        //cout << "ind=" << ind << " species=" << speciesIDNameMap[ind];
        //vector<int> cset0(speciesIDNameMap.size(),0);
        for(int dIter=0;dIter<speciesIDNameMap.size();dIter++)
        {
            if(ind==dIter)
            {
                (*cset)[dIter]=1;
                //cset0[dIter]=1;
            }
            else
            {
                (*cset)[dIter]=0;
            }
            cout << " cell=" << dIter << " status=" <<(*cset)[dIter];
        }
        //cout << " condsetMap[ind]=cset" << endl;
        //int currind=(int)pow(2.0,ind);
        condsetMap[ind]=cset;
        //condsetMap0.push_back(cset0);
        //string condKey;
        //genCondSetKey(*cset,condKey);
        //condsetKeyIDMap[condKey]=currind;
        //condsetIDKeyMap[currind]=condKey;
        //ind=ind+1;
    }
    return 0;
}

INTINTMAP*
MetaLearner::getConditionSet(int cind)
{
    if(condsetMap.find(cind)==condsetMap.end())
    {
        cout <<"Did not find any condition sets associated with " << cind << endl;
        exit(0);
    }
    return condsetMap[cind];
}

int
MetaLearner::setinputVariableNames()
{
    map<int,MappedOrthogroup*>& orthogroupSet=ogr->getMappedOrthogroups();
    for(auto og=inputVariables.begin(); og!=inputVariables.end();og++)
    {
        MappedOrthogroup* targetogrp=orthogroupSet[*og];
        vector<string>& targetgrpMembers=targetogrp->getOrthoMembers();
        for(int specID=0;specID<targetgrpMembers.size();specID++)
        {
            inputVariableNames[specID].insert(targetgrpMembers[specID]);
            //cout <<"og=" << *og <<" specID=" << specID << " varName=" << targetgrpMembers[specID] << endl;
        }
    }
    cout << "MetaLearner::setinputVariableNames number of variables: " <<inputVariables.size() << endl;
    return 0;
}
    
int
MetaLearner::collectMoves_Orthogroups(int currK)
{
    //auto start = high_resolution_clock::now();
    for(int i=0;i<moveSet.size();i++)
    {
        delete moveSet[i];
    }
    moveSet.clear();
    map<int,MappedOrthogroup*>& orthogroupSet=ogr->getMappedOrthogroups();
    
    //cout << "inputOGList.size() = " << inputOGList.size() << endl;
    //cout << "inputRegulatorOGs.size() = " << inputRegulatorOGs.size() << endl;
    //Now we will have a move for one orthogroup at a time
    for(int targeti=0; targeti<inputOGList.size();targeti++)
    {
        int targetOGIter=inputOGList[targeti];
        MappedOrthogroup* targetogrp=orthogroupSet[targetOGIter]; //oIter->second; four species
        //map<string,GeneMap*>& targetgrpMembers=targetogrp->getOrthoMembers();
        vector<string>& targetgrpMembers=targetogrp->getOrthoMembers();
        int n=speciesIDNameMap.size();
        
        //score the score of best regulator:
        vector<double> bestscore_PerSpecies(n,0);
        vector<double> bestscoreImprovement_PerSpecies(n,0);
        vector<int> besttarget_PerSpecies(n,-1);
        vector<int> besttf_PerSpecies(n,-1);
        //unordered_map<int,unordered_map<int,double>*> bestregWt_PerSpecies; only need for each species!
        int bestcsetid=-1;
        int bestregi=-1; //inputRegulatorOGs vector index
        //vector<double> bestScoreImprovementTF(n,0);
        double bestScoreImprovement_TF=0;
        //The logic of this is we will basically search for the utility of each regulator across every species. The datalikelihood
        //term is computed separately from the prior. Then we will consider what will happen if were to make moves for all species.
        
        //for each regulator
        for(int regi=0;regi<inputRegulatorOGs.size();regi++)
        {
            int regOGIter=inputRegulatorOGs[regi];
            if(regOGIter==targetOGIter)  //if(regOGIter->first==oIter->first)
            {
                continue;
            }
            MappedOrthogroup* tfogrp=orthogroupSet[regOGIter]; //orthogroupSet[regOGIter->first];
            //map<string,GeneMap*>& tfgrpMembers=tfogrp->getOrthoMembers();
            vector<string>& tfgrpMembers=tfogrp->getOrthoMembers();
            double oldpriorScore=ogpairPrior[regi][targeti]; //double oldpriorScore=ogpairPrior[ogPair];
            
            //only store current tf for each species!
            vector<double> score_PerSpecies(n,0);
            vector<double> scoreImprovement_PerSpecies(n,0);
            vector<int> target_PerSpecies(n,0);  //store variable index
            vector<int> tf_PerSpecies(n,0);
            //map<int,INTDBLMAP*> regWt_PerSpecies;
            int nscoreImp=0;            
            
            //for each species:
            for(int specID=0;specID<speciesIDNameMap.size();specID++)
            {
                //unordered_map<int,double>* bestregWt_PerSpecies=new unordered_map<int,double>; no need to store weights
                //int specID=specIter->first;
                string spec=speciesIDNameMap[specID];//specIter->second;
                SpeciesDataManager* sdm=speciesDataSet[specID];
                //map<string,int>& candidateregulators_Species=sdm->getRegulators();  //comment out
                FactorGraph* speciesGraph=sdm->getFactorGraph();
                VariableManager* vMgr=sdm->getVariableManager();
                vector<Variable*>& varSet=vMgr->getVariableSet();
                int targetID=vMgr->getVarID(targetgrpMembers[specID].c_str());
                int regID=vMgr->getVarID(tfgrpMembers[specID].c_str());
                //cout << "regi=" << regi <<  " regOGID=" << regOGIter  <<" targeti=" << targeti << " targetOGID=" <<targetOGIter << " TFname=" << tfgrpMembers[specID] <<" targetgene=" << targetgrpMembers[specID] << " regvaribleID=" << regID << " targetvaribleID=" <<targetID << endl;

                /*
                //GeneMap* geneMap_Tgts=specIter->second;
                GeneMap* geneMap_Tgts=targetgrpMembers[spec]; //GeneMap* geneMap_Tgts=targetgrpMembers[specIter->first];
                map<string,map<string,STRINTMAP*>*>& speciesTargetSet=geneMap_Tgts->getGeneSet();
                GeneMap* geneMap_TFs=tfgrpMembers[spec]; //GeneMap* geneMap_TFs=tfgrpMembers[specIter->first];
                map<string,map<string,STRINTMAP*>*>& speciesTFSet=geneMap_TFs->getGeneSet();
                 */
                //If the species has multiple genes in it's list, we consider the member which will give the max improvement, that is
                //highest data likelihood. Similarly, we will extract the highest likelihood if there are multiple regulators
                double maxScore=-999999;
                double maxScoreImprovement=-999999;
                //Best here makes sense only if there are multiple TFs and targets
                int bestTarget=-1;
                int bestTF=-1;
                //speciesTargetSet.size()=1 and speciesTFSet.size()=1
                //for(map<string,map<string,STRINTMAP*>*>::iterator vIter=speciesTargetSet.begin();vIter!=speciesTargetSet.end();vIter++)
                //{
                //int targetID=vMgr->getVarID(vIter->first.c_str());
                Variable* target=varSet[targetID];
                SlimFactor* sFactor=speciesGraph->getFactorAt(targetID);
                //If the edge already exists in the MB of sFactor continue
                if(sFactor->mergedMB.find(regID)!=sFactor->mergedMB.end())
                {
                    //cout <<"Edge exists in the MB of sFactor" << endl;
                    continue;
                }
                if(sFactor->mergedMB.size()>=currK)
                {
                    continue;
                }
                Variable* reg=varSet[regID];
                
                //Otherwise get the score delta of adding this reguluator in sFactor's MB.
                double scoreImprovement=0;
                double newScore=0;
                double newTargetScore=0;
                //unordered_map<int,double> regwt;
                //INTINTMAP* cset=getConditionSet(specID);
                //getNewPLLScore(specID,*cset,reg,target,newScore,scoreImprovement,targetOGIter,regwt); //*cset not needed in CVN
                getNewPLLScore(specID,reg,target,newScore,scoreImprovement,targetOGIter);
                cout << "specID=" << specID <<" regi=" << regi <<  " regOGID=" << regOGIter  <<" targeti=" << targeti << " targetOGID=" <<targetOGIter << " TFname=" << tfgrpMembers[specID] <<" targetgene=" << targetgrpMembers[specID] << " regvaribleID=" << regID << " targetvaribleID=" <<targetID << " scoreImprovement=" << scoreImprovement << " newScore=" << newScore << endl;
                if(scoreImprovement>0) //&& (scoreImprovement>maxScoreImprovement))
                {
                    bestTarget=targetID;  //variable index
                    bestTF=regID;   //variable index
                    maxScoreImprovement=scoreImprovement;
                    maxScore=newScore;
                }else{
                    //regwt.clear();
                    continue;
                }
                //At this stage we are done with this species, and if maxDLL >0 we proceed with updating the information for this species
                scoreImprovement_PerSpecies[specID]=maxScoreImprovement;
                score_PerSpecies[specID]=maxScore;
                target_PerSpecies[specID]=bestTarget;
                tf_PerSpecies[specID]=bestTF;
                nscoreImp++;
                /*unordered_map<int,double>* regwtforspecies=new unordered_map<int,double>;  //copy of wts/regwt
                for(auto wIter=regwt.begin();wIter!=regwt.end();wIter++)
                {
                    (*regwtforspecies)[wIter->first]=wIter->second;
                }
                //regWt_PerSpecies[specID]=regwtforspecies;
                regwt.clear();*/
                
            } //CVNvariant: species dataset end
            if(nscoreImp==0){
                continue;
                //cout <<"Skipping no hit regulator " <<  regOGIter << endl;
            } 
            //Now we wish to see how good it would be add these edges in different conditions for all cell types:
            double bestImprovement=0;  //adding edge or not, which improvement is higher
            int csetid=-1;
            //cout <<"Compute for each condsetMap_Tree" << endl;
            for(int setIter=0;setIter<condsetMap_Tree.size();setIter++)
            {
                vector<int> & cset=condsetMap_Tree[setIter];
                vector<int> speciesEdgeStat(cset.size(),0);
                int valid=1;
                //compute score improvement+prior for each condition:
                double netImprovement=0;
                for(int i=0;i<cset.size();i++)
                {
                    if(cset[i]==1){
                        speciesEdgeStat[i]=1;
                        if(scoreImprovement_PerSpecies[i]>0){
                            netImprovement += scoreImprovement_PerSpecies[i];
                        }else{
                            valid=0;
                            break;
                        }
                    }
                    else{
                        speciesEdgeStat[i]=0;
                    }
                }
		if(valid==0){
		    // This configuration is not valid
                    continue;
                }
                //compute the prior
                double ePrior=log(speciesData->getEdgeStatusProb(speciesEdgeStat));
                netImprovement=netImprovement+ePrior-oldpriorScore;
                
                if(netImprovement>bestImprovement)
                {
                    bestImprovement=netImprovement;
                    csetid=setIter;
                }
                cout << "condition"<<setIter <<": ePrior=" << ePrior << " oldpriorScore=" << oldpriorScore << " netImprovementwPrior=" << netImprovement << endl;
                speciesEdgeStat.clear();
            }
            cout << "Best condition: "<<csetid << " bestImprovement=" << bestImprovement << endl;                
            if(csetid==-1) //add nothing
            {
                continue;
            }
            if(bestImprovement>bestScoreImprovement_TF) // ePrior is added
            {
                //best condition:
                vector<int> & cset=condsetMap_Tree[csetid];
                bestscore_PerSpecies.clear();
                bestscoreImprovement_PerSpecies.clear();
                besttarget_PerSpecies.clear();
                besttf_PerSpecies.clear();
                for(int i=0;i<cset.size();i++)
                {
                    if(cset[i]==0){
                        continue;
                    }
                    bestscore_PerSpecies[i]=score_PerSpecies[i];
                    bestscoreImprovement_PerSpecies[i]=scoreImprovement_PerSpecies[i];
                    besttarget_PerSpecies[i]=target_PerSpecies[i];
                    besttf_PerSpecies[i]=tf_PerSpecies[i];
                    //bestregWt_PerSpecies->clear(); //pointer to old unordered_map<int,double>* regwtforspecies
                    //bestregWt_PerSpecies=regwtforspecies;
                }
                bestcsetid=csetid;
                bestScoreImprovement_TF=bestImprovement;
                bestregi=regi; // inputRegulatorOGs vector index
            }
        }//all inputRegulatorOGs(regulators) end

        // add the best regulator for this target:
        if(bestcsetid==-1)  //edge status cset[specID]
        {
            continue;
        }
        vector<int> & cset=condsetMap_Tree[bestcsetid];
        for(int i=0;i<cset.size();i++)
        {
            if(cset[i]==0){
                continue;
            }
            MetaMove* move=new MetaMove;
            //int tfid=besttf_PerSpecies[i];
            //int tgtid=besttarget_PerSpecies[i];
            move->setSrcVertex(besttf_PerSpecies[i]); //TF variable id
            move->setTFID(bestregi); //inputRegulatorOGs vector index
            move->setConditionSetInd(i);
            move->setTargetVertex(besttarget_PerSpecies[i]); //target variable id
            move->setTargetID(targeti); //inputOGList vector index
            move->setTargetMBScore(bestscore_PerSpecies[i]);
            move->setScoreImprovement(bestscoreImprovement_PerSpecies[i]);
            //move->setSrcWeight(*bestregWt_PerSpecies); //old *bestregWt_PerSpecies[specID]
            moveSet.push_back(move);
            //cout <<"Found edge for "<<speciesIDNameMap[i]<<" TFvar=" << besttf_PerSpecies[i] << " Targetvar=" << besttarget_PerSpecies[i] <<  " regi=" << bestregi << " targeti="<<targeti << " score delta=" << bestscoreImprovement_PerSpecies[i]<<endl;
            //bestregWt_PerSpecies->clear(); //pointer to old unordered_map<int,double>* regwtforspecies
            //delete bestregWt_PerSpecies;
        } //species dataset end
        /*double vm, rss;
        process_mem_usage(vm, rss);
        cout << "MetaLearner::collectMoves_Orthogroups VM: " << vm << " MB; RSS: " << rss << endl;*/
        
        /*for(auto wtIter=bestregWt_PerSpecies.begin();wtIter!=bestregWt_PerSpecies.end();wtIter++)
        {
            wtIter->second->clear();
            delete wtIter->second;
        }
        bestregWt_PerSpecies.clear();*/
        besttarget_PerSpecies.clear();
        bestscore_PerSpecies.clear();
        besttf_PerSpecies.clear();
        bestscoreImprovement_PerSpecies.clear();
        //bestScoreImprovementTF.clear();
    }
    //auto stop = high_resolution_clock::now();
    //auto duration = duration_cast<microseconds>(stop - start);
    //cout << "collectMoves_Orthogroups time: " << duration.count() << " microseconds" <<endl;
    return 0;
} //collectMoves_Orthogroups


int
MetaLearner::collectMoves_Orthogroups_INDEP(int currK)
{
    cout << "MetaLearner::collectMoves_Orthogroups_INDEP" << endl;
    //auto start = high_resolution_clock::now();
    for(int i=0;i<moveSet.size();i++)
    {
        delete moveSet[i];
    }
    moveSet.clear();
    map<int,MappedOrthogroup*>& orthogroupSet=ogr->getMappedOrthogroups();
    
    //cout << "inputOGList.size() = " << inputOGList.size() << endl;
    //cout << "inputRegulatorOGs.size() = " << inputRegulatorOGs.size() << endl;
    //Now we will have a move for one orthogroup at a time
    for(int targeti=0; targeti<inputOGList.size();targeti++)
    {
        int targetOGIter=inputOGList[targeti];
        MappedOrthogroup* targetogrp=orthogroupSet[targetOGIter]; //oIter->second; four species
        //map<string,GeneMap*>& targetgrpMembers=targetogrp->getOrthoMembers();
        vector<string>& targetgrpMembers=targetogrp->getOrthoMembers();
        int n=speciesIDNameMap.size();
        vector<double> bestscore_PerSpecies(n,0);
        vector<double> bestscoreImprovement_PerSpecies(n,0);
        vector<int> besttarget_PerSpecies(n,-1);
        vector<int> besttf_PerSpecies(n,-1);
        //unordered_map<int,unordered_map<int,double>*> bestregWt_PerSpecies;
        int bestcsetid=-1;
        vector<double> bestScoreImprovementTF(n,0);
        //The logic of this is we will basically search for the utility of each regulator across every species. The datalikelihood
        //term is computed separately from the prior. Then we will consider what will happen if were to make moves for all species.
        
        /*//only store current tf for each species!
        map<int,double> score_PerSpecies;
        map<int,double> scoreImprovement_PerSpecies;
        map<int,int> target_PerSpecies;
        map<int,int> tf_PerSpecies;
        map<int,INTDBLMAP*> regWt_PerSpecies;*/
        
        // we need to start from children to parents:loop from end of speciesIDNameMap to start of it!
        //for(map<int,string>::iterator specIter=speciesIDNameMap.begin();specIter!=speciesIDNameMap.end();specIter++)
        for(int specID=0;specID<speciesIDNameMap.size();specID++)
            //for(map<string,SpeciesDataManager*>::iterator specIter=speciesDataSet.begin();specIter!=speciesDataSet.end();specIter++)
        {
            //unordered_map<int,double>* bestregWt_PerSpecies=new unordered_map<int,double>; no need to store weights
            //int specID=specIter->first;
            string spec=speciesIDNameMap[specID];//specIter->second;
            SpeciesDataManager* sdm=speciesDataSet[specID];
            //map<string,int>& candidateregulators_Species=sdm->getRegulators();  //comment out
            FactorGraph* speciesGraph=sdm->getFactorGraph();
            VariableManager* vMgr=sdm->getVariableManager();
            vector<Variable*>& varSet=vMgr->getVariableSet();
            int targetID=vMgr->getVarID(targetgrpMembers[specID].c_str());
            int bestregi=-1;
            //for each regulator
            //for(map<int,int>::iterator regOGIter=inputRegulatorOGs.begin();regOGIter!=inputRegulatorOGs.end();regOGIter++)
            for(int regi=0;regi<inputRegulatorOGs.size();regi++)
            {
                int regOGIter=inputRegulatorOGs[regi];
                if(regOGIter==targetOGIter)  //if(regOGIter->first==oIter->first)
                {
                    continue;
                }
                MappedOrthogroup* tfogrp=orthogroupSet[regOGIter]; //orthogroupSet[regOGIter->first];
                //map<string,GeneMap*>& tfgrpMembers=tfogrp->getOrthoMembers();
                vector<string>& tfgrpMembers=tfogrp->getOrthoMembers();
                int regID=vMgr->getVarID(tfgrpMembers[specID].c_str());
                //If the species has multiple genes in it's list, we consider the member which will give the max improvement, that is
                //highest data likelihood. Similarly, we will extract the highest likelihood if there are multiple regulators
                double maxScore=-999999;
                double maxScoreImprovement=-999999;
                //Best here makes sense only if there are multiple TFs and targets
                int bestTarget=-1;
                int bestTF=-1;
                //speciesTargetSet.size()=1 and speciesTFSet.size()=1
                //for(map<string,map<string,STRINTMAP*>*>::iterator vIter=speciesTargetSet.begin();vIter!=speciesTargetSet.end();vIter++)
                //{
                //int targetID=vMgr->getVarID(vIter->first.c_str());
                Variable* target=varSet[targetID];
                SlimFactor* sFactor=speciesGraph->getFactorAt(targetID);
                //If the edge already exists in the MB of sFactor continue
                if(sFactor->mergedMB.find(regID)!=sFactor->mergedMB.end())
                {
                    //cout <<"Edge exists in the MB of sFactor" << endl;
                    continue;
                }
                if(sFactor->mergedMB.size()>=currK)
                {
                    continue;
                }
                Variable* reg=varSet[regID];
                
                //Otherwise get the score delta of adding this reguluator in sFactor's MB.
                double scoreImprovement=0;
                double newScore=0;
                double newTargetScore=0;
                //unordered_map<int,double> regwt;
                //INTINTMAP* cset=getConditionSet(specID);
                //getNewPLLScore(specID,*cset,reg,target,newScore,scoreImprovement,targetOGIter,regwt); //*cset not needed in CVN
                getNewPLLScore(specID,reg,target,newScore,scoreImprovement,targetOGIter);
                //cout << "specID=" << specID <<" regi=" << regi <<  " regOGID=" << regOGIter  <<" targeti=" << targeti << " targetOGID=" <<targetOGIter << " TFname=" << tfgrpMembers[specID] <<" targetgene=" << targetgrpMembers[specID] << " regvaribleID=" << regID << " targetvaribleID=" <<targetID << " scoreImprovement=" << scoreImprovement << " newScore=" << newScore << endl;
                //Now we wish to see how good it would be add these edges in different groups
                // create my own edge status: add this one or not
                int regAdd=-1;
                if(scoreImprovement>0) //&& (scoreImprovement>maxScoreImprovement))
                {
                    bestTarget=targetID;  //variable index
                    bestTF=regID;   //variable index
                    maxScoreImprovement=scoreImprovement;
                    maxScore=newScore;
                    regAdd=1;
                }else{
                    //regwt.clear();
                    continue;
                }
                /*unordered_map<int,double>* regwtforspecies=new unordered_map<int,double>;  //copy of wts/regwt
                for(auto wIter=regwt.begin();wIter!=regwt.end();wIter++)
                {
                    (*regwtforspecies)[wIter->first]=wIter->second;
                }
                //regWt_PerSpecies[specID]=regwtforspecies;
                regwt.clear();*/
                
                // add this one or not
                if(maxScoreImprovement>bestScoreImprovementTF[specID]) // ePrior is added
                {
                    bestscore_PerSpecies[specID]=maxScore;
                    bestscoreImprovement_PerSpecies[specID]=maxScoreImprovement;
                    besttarget_PerSpecies[specID]=bestTarget;
                    besttf_PerSpecies[specID]=bestTF;
                    //bestregWt_PerSpecies->clear(); //pointer to old unordered_map<int,double>* regwtforspecies
                    //bestregWt_PerSpecies=regwtforspecies;
                    bestcsetid=regAdd;
                    bestScoreImprovementTF[specID]=maxScoreImprovement;
                    bestregi=regi;
                }/*else
                {
                    regwtforspecies->clear();
                    delete regwtforspecies;
                }*/
            }//all inputRegulatorOGs(regulators) end
            
            // add the regulator: bestcsetid=1 add bestcsetid=0 continue
            if(bestcsetid!=1 || bestscoreImprovement_PerSpecies[specID]<=0)  //edge status cset[specID]
            {
                continue;
            }
            MetaMove* move=new MetaMove;
            //int tfid=besttf_PerSpecies[specID];
            //int tgtid=besttarget_PerSpecies[specID];
            move->setSrcVertex(besttf_PerSpecies[specID]);
            move->setTFID(bestregi);
            move->setConditionSetInd(specID);
            move->setTargetVertex(besttarget_PerSpecies[specID]);
            move->setTargetID(targeti);
            move->setTargetMBScore(bestscore_PerSpecies[specID]);
            move->setScoreImprovement(bestscoreImprovement_PerSpecies[specID]);
            //move->setSrcWeight(*bestregWt_PerSpecies);
            moveSet.push_back(move);
            //bestregWt_PerSpecies->clear(); //pointer to old unordered_map<int,double>* regwtforspecies
            //delete bestregWt_PerSpecies;
            
        } //species dataset end
        besttarget_PerSpecies.clear();
        bestscore_PerSpecies.clear();
        besttf_PerSpecies.clear();
        bestscoreImprovement_PerSpecies.clear();
        bestScoreImprovementTF.clear();
        /*for(auto wtIter=bestregWt_PerSpecies.begin();wtIter!=bestregWt_PerSpecies.end();wtIter++)
        {
            wtIter->second->clear();
            delete wtIter->second;
        }
        bestregWt_PerSpecies.clear();*/
    }
    //auto stop = high_resolution_clock::now();
    //auto duration = duration_cast<microseconds>(stop - start);
    //cout << "collectMoves_Orthogroups time: " << duration.count() << " microseconds" <<endl;
    return 0;
}




// currPrior is the species-specific prior, which we do not need
// priorScore is computed by not counted
// getNewPLLScore(specID,*cset,reg,target,newScore,scoreImprovement,oIter->first,regwt);
// shilu: removed INTINTMAP& conditionSet
// u is reg and v is target
// species-specific prior: sum_reg[log(p)]+sum_nonreg[log(1-p)]
// score=likelihood+species-specific prior:
int
MetaLearner::getNewPLLScore(int cid, Variable* u, Variable* v, double& targetmbScore, double& scoreImprovement,int orthoGrpNo) //unordered_map<int,double>& regwt
{
    //Each condition set has an evidence manager which pools the data from the specific conditions.
    //string condKey;
    //genCondSetKey(conditionSet,condKey);
    SpeciesDataManager* sdm=speciesDataSet[cid];
    //PotentialManager* potMgr=sdm->getPotentialManager();
    //VSET& varSet=sdm->getVariableManager()->getVariableSet();
    FactorGraph* fg=sdm->getFactorGraph();
    //SlimFactor* sFactor=fg->getFactorAt(u->getID());  //regulator
    SlimFactor* dFactor=fg->getFactorAt(v->getID());  //target
    //map<int,double>* varNeighborhoodPrior=varNeighborhoodPrior_PerSpecies[speciesIDNameMap[cid]];
    vector<double>& varNeighborhoodPrior=varNeighborhoodPrior_PerSpecies[cid];
    vector<vector<double>>& edgePresenceProb=edgePresenceProb_PerSpecies[cid];
    //map<string,double>* edgePresenceProb=edgePresenceProb_PerSpecies[speciesIDNameMap[cid]];
    double currPrior=varNeighborhoodPrior[v->getID()]; //target
    bool toDel_d=true;
    /*if(dFactor->mergedMB.find(u->getID())!=dFactor->mergedMB.end())
    {
        toDel_d=false;
    }*/  //already checked mergedMB before computing getNewPLLScore
    double plus=0;
    double minus=0;
    dFactor->mergedMB.insert(u->getID()); //Aug 23: dFactor->mergedMB[u->getID()]=0;
    int status=0;
    for(auto mIter=dFactor->mergedMB.begin();mIter!=dFactor->mergedMB.end();mIter++)
    {
        /*Variable* aVar=varSet[mIter->first];
        string regulatorKey(aVar->getName().c_str());
        regulatorKey.append("\t");
        regulatorKey.append(v->getName().c_str());*/
        double p=edgePresenceProb[*mIter][v->getID()]; //(*edgePresenceProb)[regulatorKey];
        // update by Shilu to avoid inf
        if(p==0 or p==1)
        {
            continue;
        }else{
            minus=minus+log(1-p);
            plus=plus+log(p);
        }
    }
    double pll_d=getPLLScore_Condition_Tracetrick(cid,dFactor,status);
    //double pll_d=getPLLScore_Condition(speciesIDNameMap[cid],dFactor,status,regwt);
    if(status==-1)
    {
        scoreImprovement=-1;
        if(toDel_d)
        {
            auto dIter=dFactor->mergedMB.find(u->getID());
            dFactor->mergedMB.erase(dIter);
        }
        return 0;
    }
    currPrior=currPrior+plus-minus;
    cout << "cell=" << cid << " regID=" <<u->getID() <<" varID=" << v->getID() << " oldsparsityPrior=" << varNeighborhoodPrior[v->getID()] <<" sparsityPriorchanged=" << plus-minus << " currsparsityPrior=" << currPrior << " pll_d=" << pll_d << " newscore=" << pll_d+currPrior << " oldscore=" << dFactor->mbScore << endl;
    pll_d=pll_d+currPrior;
    /* comment by shilu
     double priorScore=0;
     double oldpriorScore=0;
     if(speciesData->getRoot()!=NULL)
     {
     priorScore=getEdgePrior(cid,u,v,orthoGrpNo);
     int uogno=ogr->getMappedOrthogroupID(u->getName().c_str(),speciesIDNameMap[cid].c_str());
     int vogno=ogr->getMappedOrthogroupID(v->getName().c_str(),speciesIDNameMap[cid].c_str());
     if(vogno!=-1)
     {
     char ogPairKey[256];
     sprintf(ogPairKey,"%d-%d",uogno,vogno);
     string ogPair(ogPairKey);
     oldpriorScore=ogpairPrior[ogPair];
     }
     }*/
    targetmbScore=pll_d;
    //scoreImprovement=0;
    //double dImpr=(targetmbScore+priorScore)-(dFactor->mbScore+oldpriorScore);
    //Don't include the prior. Just use the data likelihood improvement
    double dImpr=targetmbScore-dFactor->mbScore;
    //

    if(dImpr<=0)
    {
        scoreImprovement=-1;
    }
    else
    {
        scoreImprovement=dImpr;
    }
    if(toDel_d)
    {
        auto dIter=dFactor->mergedMB.find(u->getID());
        dFactor->mergedMB.erase(dIter);
    }
    cout << "getPLLScore_Condition_Tracetrick likelihood=" <<pll_d << " dFactor->mbScore="<<dFactor->mbScore<<" in " <<speciesIDNameMap[cid] << " for target " << dFactor << " status=" <<status << " toDel_d=" << toDel_d<< endl;

    return 0;
}
/*
double
MetaLearner::getPLLScore_Condition_onefold(string& speciesName,SlimFactor* sFactor)
{
    SpeciesDataManager* spd=speciesDataSet[speciesName];
    PotentialManager* potMgr=spd->getPotentialManager();
    //EvidenceManager* evMgr=spd->getEvidenceManager();
    VariableManager* varMgr=spd->getVariableManager();
    VSET& varSet=varMgr->getVariableSet();
    double unreg_pll=potMgr->getPseudoLikelihood_onefold(sFactor,varSet);
    return unreg_pll;
}*/

double
MetaLearner::getPLLScore_Condition(int specID,SlimFactor* sFactor)
{
    SpeciesDataManager* spd=speciesDataSet[specID];
    PotentialManager* potMgr=spd->getPotentialManager();
    //EvidenceManager* evMgr=spd->getEvidenceManager();
    VariableManager* varMgr=spd->getVariableManager();
    vector<Variable*>& varSet=varMgr->getVariableSet();
    /*string mbkey;
    char mbchar[256];
    for(INTINTMAP_ITER mIter=sFactor->mergedMB.begin();mIter!=sFactor->mergedMB.end();mIter++)
    {
        sprintf(mbchar,"-%d",mIter->first);
        mbkey.append(mbchar);
    }
    if(strcmp(mbkey.c_str(),"-1-35-39-45")==0)
    {
        cout <<"Found " << mbkey.c_str() << endl;
    }*/
    double unreg_pll=potMgr->getPseudoLikelihood(sFactor,varSet,false);
    if((isnan(unreg_pll)) || (isinf(unreg_pll)))
    {
        cout <<"Found nan for " << sFactor->fId<< ": MB: ";
        for(auto mIter=sFactor->mergedMB.begin();mIter!=sFactor->mergedMB.end();mIter++)
        {
            cout <<"-"<<*mIter;
        }
        cout<< endl;
    }
    /* //comment by shilu
    //double varCnt=sFactor->mergedMB.size()+1;
    double varCnt=sFactor->mergedMB.size();
    double paramCnt=2*varCnt;
    //double paramCnt=varCnt;
    paramCnt=paramCnt+((varCnt*(varCnt-1))/2);
    double pll=unreg_pll-((paramCnt/2)*log(evMgr->getTestSet().size()));
    //pll=unreg_pll-(0.5*log(evMgr->getTestSet().size()));
    pll=unreg_pll;*/
    return unreg_pll;
}

/*
// commented by vperiyasamy
double
MetaLearner::getPLLScore_Condition(string& speciesName,SlimFactor* sFactor,int& status,INTDBLMAP& regWts)
{
    PotentialManager* potMgr=speciesDataSet[speciesName]->getPotentialManager(); // grab managers
    EvidenceManager* evMgr=speciesDataSet[speciesName]->getEvidenceManager(); // for this
    VariableManager* varMgr=speciesDataSet[speciesName]->getVariableManager(); // species
    VSET& varSet=varMgr->getVariableSet();
    
    string mbkey;
    char mbchar[256];
    if(sFactor->mergedMB.size()>0) // sFactor here is for vertex v, if no parents then what?
    {
        //	cout <<"Stop here" << endl;
    }
    int statuspot=0;
    Potential* apot=new Potential;
    
    // get unregularized PLL score?
    double unreg_pll=potMgr->getPseudoLikelihood(sFactor,varSet,false,statuspot,apot);
    if(statuspot==-1)
    {
        status=-1;
        delete apot;
        return 0;
    }
    
    // get condition weights from filled in Potential
    INTDBLMAP& potwts=apot->getCondWeight();
    for(INTDBLMAP_ITER wIter=potwts.begin();wIter!=potwts.end();wIter++)
    {
        regWts[wIter->first]=wIter->second;
    }
    
    if((isnan(unreg_pll)) || (isinf(unreg_pll)))
    {
        cout <<"Found nan for " << sFactor->fId<< ": MB: ";
        for(auto mIter=sFactor->mergedMB.begin();mIter!=sFactor->mergedMB.end();mIter++)
        {
            cout <<"-"<<*mIter;
        }
        cout<< endl;
    }
    
    // increase parent count because we have a new edge
    double varCnt=sFactor->mergedMB.size()+1.0;
    //double varCnt=sFactor->mergedMB.size();
    
    double paramCnt=2*varCnt;
    paramCnt=paramCnt+((varCnt*(varCnt-1))/2);
    //paramCnt=paramCnt+((varCnt*(varCnt-1))/2);
    double pll=unreg_pll-(0.5*paramCnt*log(evMgr->getTestSetSize()));//getTestSet().size())); // MERLIN uses training set size? does this make a difference
    //double pll=unreg_pll;
    pll=unreg_pll; // ?????? why are we just putting it back to unreg_pll
    delete apot;
    return pll;
}*/

// shilu: more efficient version
double
//MetaLearner::getNewPLLScore_Condition_Tracetrick(int csetId, int vId, int uId, Potential* newPot)
MetaLearner::getPLLScore_Condition_Tracetrick(int specID,SlimFactor* sFactor,int& status)  //sFactor is target, unordered_map<int,double>& regWts
{
    
    //cout << "MetaLearner::getPLLScore_Condition_Tracetrick() for " << speciesName ;
    PotentialManager* potMgr=speciesDataSet[specID]->getPotentialManager(); // grab managers
    EvidenceManager* evMgr=speciesDataSet[specID]->getEvidenceManager(); // for this
    //VariableManager* varMgr=speciesDataSet[speciesName]->getVariableManager(); // species
    //VSET& varSet=varMgr->getVariableSet();

    // NOTE: use test set in order to get same results as non tracetrick
    //int datasize = evMgr->getTestSetSize(); //tSet->size();
    double pll=potMgr->computePotentialMBCovMean(sFactor,status);
    /*if(status==-1) {
        return 0;
    }*/
    //cout << " data likelihood is " << pll << " getPLLScore_Condition_Tracetrick runtime: " <<duration.count() << " ms" <<endl;
    return pll;
}

/*
// added by vperiyasamy
double
MetaLearner::getPLLScore_Condition_Tracetrick(string& speciesName,SlimFactor* sFactor,int& status,unordered_map<int,double>& regWts)  //sFactor is target
{
    
    //cout << "MetaLearner::getPLLScore_Condition_Tracetrick() for " << speciesName ;
    PotentialManager* potMgr=speciesDataSet[speciesName]->getPotentialManager(); // grab managers
    EvidenceManager* evMgr=speciesDataSet[speciesName]->getEvidenceManager(); // for this
    VariableManager* varMgr=speciesDataSet[speciesName]->getVariableManager(); // species
    VSET& varSet=varMgr->getVariableSet();
    
    // first define newPot and set it up
    Potential* newPot = new Potential;
    newPot->setAssocVariable(varSet[sFactor->fId],Potential::FACTOR); // set v as the target for this potential
    //vector <int> children;
    //cout << " Target is " << varSet[sFactor->fId]->getName() << " id=" << sFactor->fId<< " market blank is: ";
    //children.push_back(sFactor->fId);
    // second define parentPot and set it up (store regulator of sFactor[target])
    Potential* parentPot=new Potential;
    for(auto aIter=sFactor->mergedMB.begin();aIter!=sFactor->mergedMB.end();aIter++) {
        Variable* aVar=varSet[*aIter];
        newPot->setAssocVariable(aVar,Potential::MARKOV_BNKT);  //set sFactor->mergedMB as the MARKOV_BNKT
        //cout << aVar->getName() <<" id=" <<aVar->getID() << " " ;
        parentPot->setAssocVariable(aVar,Potential::MARKOV_BNKT); //put sFactor->mergedMB into parent
        //children.push_back(aVar->getID());
    }
    //cout << endl;
    
    newPot->potZeroInit(); // initialize meanWrap[INTDBLMAP], covariance[matrix], mean[1d matrix]
    //update mean and covariance
    potMgr->populatePotential_Eff(newPot,false);
    //This function creates a submatrix of the covariance matrix and inverts it
    newPot->initMBCovMean(); // compute mbcondVar[CondVariance], mbcondMean_Part[CondBias], mbcondMean_Vect[CondWeight]
    if(newPot->getCondVariance() < 0) {
        status=-1;
        delete newPot;
        return 0;
    }
    // get condition weights from filled in Potential
    unordered_map<int,double>& potwts = newPot->getCondWeight();
    for(auto wIter = potwts.begin(); wIter != potwts.end(); wIter++) {
        regWts[wIter->first] = wIter->second;
    }
    
    //double pll=0;
    //Need to fix this to be set automatically
    //get parameter prior
    //double paramPrior=0;
    // Initialize parentPot:
    parentPot->potZeroInit();
    //update mean and covariance
    potMgr->populatePotential_Eff(parentPot,false);
    
    // NOTE: use test set in order to get same results as non tracetrick
    //vector<int>* tSet = &evMgr->getTestSet(); //INTINTMAP* tSet = &evMgr->getTestSet();
    int datasize = evMgr->getTestSetSize(); //tSet->size();
    double jointll1 = newPot->computeLL_Tracetrick(datasize);
    double jointll2 = parentPot->computeLL_Tracetrick(datasize);
    //double vCnt=(double)newPot->getAssocVariables().size();
    //double paramCnt=paramCnt+(2*vCnt)+((vCnt*(vCnt-1))/2);
    double pll = jointll1 - jointll2;
    //pll=pll-(lambda*paramCnt*log(datasize));
    //pll = pll - (0.5 * paramCnt * log(datasize));
    delete newPot;
    delete parentPot;
    //cout << " data likelihood is " << pll << " getPLLScore_Condition_Tracetrick runtime: " <<duration.count() << " ms" <<endl;
    return pll;
}*/

double
MetaLearner::getValidationPLLScore_Condition(int cId, int vId)
{
    SpeciesDataManager* sdm=speciesDataSet[cId];
    EvidenceManager* evMgr=sdm->getEvidenceManager();
    VariableManager* varMgr=sdm->getVariableManager();
    INTINTMAP& vSet=evMgr->getValidationSet();
    if(vSet.size()==0){
        return 0;
    }
    double pll=0;
    //Need to fix this to be set automatically
    double wt=0.5;
    if(!dataPool)
    {
        wt=1;
    }
    double paramCnt=0;
    vector<Variable*>& varSet=varMgr->getVariableSet();
    map<int,Potential*> potSet;
    int thresholded=0;
    for(INTINTMAP_ITER eIter=vSet.begin();eIter!=vSet.end();eIter++)
    {
        EMAP* evidMap=evMgr->getEvidenceAt(eIter->first);
        //Go over all condition sets that include cInd
        double cll=0;
        for(map<int,INTINTMAP*>::iterator csIter=condsetMap.begin();csIter!=condsetMap.end();csIter++)
        {
            INTINTMAP* cset=csIter->second;
            if((*cset)[cId]==0)
            {
                continue;
            }
            SpeciesDataManager* localSdm=speciesDataSet[csIter->first];
            FactorGraph* fg=localSdm->getFactorGraph();
            SlimFactor* sFactor=fg->getFactorAt(vId);
            Potential* sPot=NULL;
            if(potSet.find(csIter->first)==potSet.end())
            {
                sPot=new Potential;
                potSet[csIter->first]=sPot;
                sPot->setAssocVariable(varSet[sFactor->fId],Potential::FACTOR);
                for(auto mIter=sFactor->mergedMB.begin();mIter!=sFactor->mergedMB.end();mIter++)
                {
                    Variable* aVar=varSet[*mIter];
                    sPot->setAssocVariable(aVar,Potential::MARKOV_BNKT);
                }
                sPot->potZeroInit();
                string condKey;
                genCondSetKey(*cset,condKey);
                PotentialManager* potMgr=pooledPotentials[condKey];
                //potMgr->populatePotential(sPot,false);
                potMgr->populatePotential_Eff(sPot,false);
                sPot->initMBCovMean();
            }
            else
            {
                sPot=potSet[csIter->first];
            }
            double pval=sPot->getCondPotValueFor(evidMap);
            
            if(pval<1e-50)
            {
                pval=1e-50;
                thresholded++;
            }
            if((pval==0) || (isnan(pval)) || (isinf(pval)))
            {
                //	cout <<"Stop here" << endl;
            }
            cll=cll+(pval*wt);
            if(eIter==vSet.begin())
            {
                double vCnt=(double)sPot->getAssocVariables().size();
                paramCnt=paramCnt+(2*vCnt)+((vCnt*(vCnt-1))/2);
            }
        }
        //Check here for really small loglikelihoods
        pll=pll+log(cll);
        
    }
    pll=pll-(0.5*paramCnt*log(vSet.size()));
    for(map<int,Potential*>::iterator pIter=potSet.begin();pIter!=potSet.end();pIter++)
    {
        delete pIter->second;
    }
    potSet.clear();
    return pll;
}

/*
//This function computes the prior score of the edge between u and v in the species specified by conditionset
double
MetaLearner::getEdgePrior(int specID,Variable* u, Variable* v,int orthoGrpNo)
{
    double prior=0;
    if(speciesDataSet.size()==1)
    {
        return 0;
    }
    STRINTMAP edgeStatus;
    string& species=speciesIDNameMap[specID];
    for(map<string,int>::iterator sIter=speciesNameIDMap.begin();sIter!=speciesNameIDMap.end();sIter++)
    {
        if(sIter->second==specID)
        {
            continue;
        }
        int hit=0;
        STRINTMAP* specTFs=ogr->getOrtholog(species.c_str(),u->getName().c_str(),sIter->first.c_str());
        STRINTMAP* specTargets=ogr->getOrtholog(species.c_str(),v->getName().c_str(),sIter->first.c_str());
        if((specTFs==NULL) || (specTargets==NULL))
        {
            edgeStatus[sIter->first]=0;
        }
        else
        {
            SpeciesDataManager* sdm=speciesDataSet[sIter->first];
            FactorGraph* speciesGraph=sdm->getFactorGraph();
            VariableManager* vMgr=sdm->getVariableManager();
            VSET& varSet=vMgr->getVariableSet();
            //Do any of the genes in this species have data and if is anybody connected to u already
            for(STRINTMAP_ITER gIter=specTargets->begin();gIter!=specTargets->end();gIter++)
            {
                int tgtID=vMgr->getVarID(gIter->first.c_str());
                if(tgtID==-1)
                {
                    continue;
                }
                SlimFactor* sFactor=speciesGraph->getFactorAt(tgtID);
                for(STRINTMAP_ITER fIter=specTFs->begin();fIter!=specTFs->end();fIter++)
                {
                    int tfID=vMgr->getVarID(fIter->first.c_str());
                    if(tfID==-1)
                    {
                        continue;
                    }
                    if(sFactor->mergedMB.find(tfID)!=sFactor->mergedMB.end())
                    {
                        hit++;
                    }
                }
            }
        }
        if(hit>0)
        {
            edgeStatus[sIter->first]=1;
        }
        else
        {
            edgeStatus[sIter->first]=0;;
        }
    }
    edgeStatus[species]=1;
    prior=speciesData->getEdgeStatusProb(edgeStatus);
    double logPrior=log(prior);
    edgeStatus.clear();
    return logPrior;
}
*/

// prior probability: between regulator j and target k as a logistic function, remove motifweight
// tfID is variable ID for regulator, targetID is variable ID for target.
double
MetaLearner::getEdgePrior_PerSpecies(int tfID, int targetID,SpeciesDataManager* sdm)
{
    INTDBLMAP* regPriors=NULL;
    //double prior=1/(1+exp(-1*beta1));
    double motifweight=0;
    //double chipweight=0;
    unordered_map<int,unordered_map<int,double>*>& motifNetwork=sdm->getMotifNetwork();
    if(motifNetwork.find(tfID)!=motifNetwork.end())
    {
        unordered_map<int,double>* values=motifNetwork[tfID];
        if(values->find(targetID)!=values->end())
        {
            motifweight=(*values)[targetID];
            //cout <<"tfID=" << tfID <<" targetID=" <<targetID<<" motifweight=" << motifweight <<  " prior=" << 1/(1+exp(-1*(beta1+motifweight*beta2)))<<endl;
        }
    }
    double fwt=motifweight*beta2;
    double prior=1/(1+exp(-1*(beta1+fwt)));
    if(prior<1e-6)
    {
        prior=1e-6;
    }
    if(prior==1)
    {
        prior=1-1e-6;
    }
    return prior;
}

/*
// prior probability: between regulator j and target k as a logistic function:
double
MetaLearner::getEdgePrior_PerSpecies(int tfID, int targetID,SpeciesDataManager* sdm)
{
    INTDBLMAP* regPriors=NULL;
    double prior=1/(1+exp(-1*beta1));
    double motifweight=0;
    double chipweight=0;
    map<int,map<int,double>*>& motifNetwork=sdm->getMotifNetwork();
    if(motifNetwork.find(tfID)!=motifNetwork.end())
    {
        map<int,double>* values=motifNetwork[tfID];
        if(values->find(targetID)!=values->end())
        {
            motifweight=(*values)[targetID];
        }
    }
    double fwt=motifweight*beta2;
    prior=1/(1+exp(-1*(beta1+fwt)));
    if(prior<1e-6)
    {
        prior=1e-6;
    }
    if(prior==1)
    {
        prior=1-1e-6;
    }
    return prior;
}

//This function computes the prior score of the edge between u and v in the species specified by conditionset
double
MetaLearner::getEdgePrior(int specID,Variable* u, Variable* v,int orthoGrpNo, int addRemStat)
{
    double prior=0;
    STRINTMAP edgeStatus;
    string& species=speciesIDNameMap[specID];
    for(map<string,int>::iterator sIter=speciesNameIDMap.begin();sIter!=speciesNameIDMap.end();sIter++)
    {
        int hit=0;
        STRINTMAP* specTargets=ogr->getOrtholog(species.c_str(),v->getName().c_str(),sIter->first.c_str());
        STRINTMAP* specTFs=ogr->getOrtholog(species.c_str(),u->getName().c_str(),sIter->first.c_str());
        if((specTargets==NULL)||(specTFs==NULL))
        {
            edgeStatus[sIter->first]=0;
        }
        else
        {
            SpeciesDataManager* sdm=speciesDataSet[sIter->first];
            FactorGraph* speciesGraph=sdm->getFactorGraph();
            VariableManager* vMgr=sdm->getVariableManager();
            VSET& varSet=vMgr->getVariableSet();
            //Do any of the genes in this species have data and if is anybody connected to u already
            for(STRINTMAP_ITER gIter=specTargets->begin();gIter!=specTargets->end();gIter++)
            {
                int tgtID=vMgr->getVarID(gIter->first.c_str());
                if(tgtID==-1)
                {
                    continue;
                }
                SlimFactor* sFactor=speciesGraph->getFactorAt(tgtID);
                for(STRINTMAP_ITER fIter=specTFs->begin();fIter!=specTFs->end();fIter++)
                {
                    int tfID=vMgr->getVarID(fIter->first.c_str());
                    if(tfID==-1)
                    {
                        continue;
                    }
                    if(sFactor->mergedMB.find(tfID)!=sFactor->mergedMB.end())
                    {
                        hit++;
                    }
                }
            }
        }
        if(hit>0)
        {
            edgeStatus[sIter->first]=1;
        }
    }
    //Now we need to figure out what the status of the edge is if v has duplicates and it is present in the dataset. Also if u has
    //duplicates and is present in the dataset.
    SpeciesDataManager* sdm=speciesDataSet[species];
    FactorGraph* speciesGraph=sdm->getFactorGraph();
    VariableManager* vMgr=sdm->getVariableManager();
    MappedOrthogroup* ogrp_TF=ogr->getMappedOrthogroup(u->getName().c_str(),species.c_str());
    MappedOrthogroup* ogrp_Tgt=ogr->getMappedOrthogroup(v->getName().c_str(),species.c_str());
    map<string,map<string,STRINTMAP*>*>& targets=ogrp_Tgt->getSpeciesHits(species.c_str())->getGeneSet();
    map<string,map<string,STRINTMAP*>*>& tfs=ogrp_TF->getSpeciesHits(species.c_str())->getGeneSet();
    int hitinspecies=0;
    for(map<string,map<string,STRINTMAP*>*>::iterator tgtIter=targets.begin();tgtIter!=targets.end();tgtIter++)
    {
        int vId=vMgr->getVarID(tgtIter->first.c_str());
        if(vId==-1)
        {
            continue;
        }
        SlimFactor* sFactor=speciesGraph->getFactorAt(vId);
        for(map<string,map<string,STRINTMAP*>*>::iterator tfIter=tfs.begin();tfIter!=tfs.end();tfIter++)
        {
            int uId=vMgr->getVarID(tfIter->first.c_str());
            if(uId==-1)
            {
                continue;
            }
            if(sFactor->mergedMB.find(uId)!=sFactor->mergedMB.end())
            {
                hitinspecies++;
            }
        }
    }
    if(hitinspecies>0)
    {
        edgeStatus[species]=1;
    }
    else
    {
        edgeStatus[species]=0;
    }
    prior=speciesData->getEdgeStatusProb(edgeStatus);
    double logPrior=log(prior);
    return logPrior;
}*/

int
MetaLearner::genCondSetKey(INTINTMAP& condSet, string& aKey)
{
    char keypair[256];
    for(INTINTMAP_ITER cIter=condSet.begin();cIter!=condSet.end();cIter++)
    {
        sprintf(keypair,"-%d=%d",cIter->first,cIter->second);
        aKey.append(keypair);
    }
    return 0;
}

int
MetaLearner::sortMoves()
{
    for(int m=0;m<moveSet.size();m++)
    {
        for(int n=m+1;n<moveSet.size();n++)
        {
            MetaMove* m1=moveSet[m];
            MetaMove* m2=moveSet[n];
            if(m1->getScoreImprovement()<m2->getScoreImprovement())
            {
                moveSet[m]=m2;
                moveSet[n]=m1;
            }
        }
    }
    return 0;
}

double
MetaLearner::makeMoves(int& successMove)
{
    /*map<int,INTINTMAP*> affectedVariables;
    for(map<string,SpeciesDataManager*>::iterator gIter=speciesDataSet.begin();gIter!=speciesDataSet.end();gIter++)
    {
        int cind=speciesNameIDMap[gIter->first];
        INTINTMAP* csVars=new INTINTMAP;
        affectedVariables[cind]=csVars;
    }*/
    //int successMove=0;
    double netScoreDelta=0;
    //int net1Move=0;
    //int net2Move=0;
    for(int m=0;m<moveSet.size();m++)
    {
        MetaMove* move=moveSet[m];
        if(attemptMove(move)==0) //if(attemptMove(move,affectedVariables)==0)
        {
            /*if(move->getConditionSetInd()==1)
            {
                net1Move++;
            }
            else if(move->getConditionSetInd()==2)
            {
                net2Move++;
            }*/
            successMove++;
            netScoreDelta=netScoreDelta+move->getScoreImprovement();
        }
        else
        {
            cout <<"Move not valid " << endl;
        }
    }
    /*double vm, rss;
    process_mem_usage(vm, rss);
    cout << "MetaLearner::makeMoves VM: " << vm << " MB; RSS: " << rss << endl;*/
    /*for(map<int,INTINTMAP*>::iterator cIter=affectedVariables.begin();cIter!=affectedVariables.end();cIter++)
    {
        cIter->second->clear();
        delete cIter->second;
    }
    affectedVariables.clear();*/
    //cout <<"Total successful moves " << successMove << " out of total " << moveSet.size() << " with net score improvement " << netScoreDelta << endl;
    return netScoreDelta;
}

/*
int
MetaLearner::makeDelMoves()
{
    map<int,INTINTMAP*> affectedVariables;
    for(map<string,SpeciesDataManager*>::iterator gIter=speciesDataSet.begin();gIter!=speciesDataSet.end();gIter++)
    {
        int cind=speciesNameIDMap[gIter->first];
        INTINTMAP* csVars=new INTINTMAP;
        affectedVariables[cind]=csVars;
    }
    int successMove=0;
    double netScoreDelta=0;
    int net1Move=0;
    int net2Move=0;
    for(int m=0;m<moveSet.size();m++)
    {
        MetaMove* move=moveSet[m];
        if(attemptDelMove(move,affectedVariables)==0)
        {
            successMove++;
            netScoreDelta=netScoreDelta+move->getScoreImprovement();
        }
    }
    for(map<int,INTINTMAP*>::iterator cIter=affectedVariables.begin();cIter!=affectedVariables.end();cIter++)
    {
        cIter->second->clear();
        delete cIter->second;
    }
    affectedVariables.clear();
    cout <<"Total successful moves " << successMove << " out of total "
    << moveSet.size() << " with net score improvement "
    << netScoreDelta
    << endl;
    return 0;
}*/

// affectedVars not needed, each target is added only once
int
MetaLearner::attemptMove(MetaMove* move)
{
    int specID=move->getConditionSetInd();
    SpeciesDataManager* sdm=speciesDataSet[specID];
    //vector<Variable*>& varSet=sdm->getVariableManager()->getVariableSet();
    //Variable* u=varSet[move->getSrcVertex()];
    //Variable* v=varSet[move->getTargetVertex()];
    int regi=move->getTFID();
    int targeti=move->getTargetID();
    FactorGraph* csGraph=sdm->getFactorGraph();
    SlimFactor* dFactor=csGraph->getFactorAt(move->getTargetVertex());
    dFactor->mergedMB.insert(move->getSrcVertex()); //Aug 23: dFactor->mergedMB[move->getSrcVertex()]=0;
    dFactor->mbScore=move->getTargetMBScore();
    dFactor->setMBWts(move->getSrcWeight());

    vector<int>* newEdgeStatus;//STRINTMAP* newEdgeStatus=NULL;
    char ogpair[20];
    sprintf(ogpair,"%d-%d",regi,targeti);
    string ogpairKey(ogpair);
    //cout <<"Attempting move " << ogpairKey << " in species " << speciesIDNameMap[specID] << endl;
    if(affectedOGPairs.find(ogpairKey)==affectedOGPairs.end())
    {
        newEdgeStatus=new vector<int>(speciesIDNameMap.size(),0);//new STRINTMAP;
        affectedOGPairs[ogpairKey]=newEdgeStatus;
    }
    else
    {
        newEdgeStatus=affectedOGPairs[ogpairKey];
    }
    //string& species=speciesIDNameMap[specID];
    (*newEdgeStatus)[specID]=1;
    return 0;
}

/*
int
MetaLearner::attemptMove(MetaMove* move,map<int,INTINTMAP*>& affectedVars)
{
    int specID=move->getConditionSetInd();
    SpeciesDataManager* sdm=speciesDataSet[speciesIDNameMap[specID]];
    VSET& varSet=sdm->getVariableManager()->getVariableSet();
    Variable* u=varSet[move->getSrcVertex()];
    Variable* v=varSet[move->getTargetVertex()];
    int regi=move->getTFID();
    int targeti=move->getTargetID();
    int uogno=ogr->getMappedOrthogroupID(u->getName().c_str(),speciesIDNameMap[specID].c_str()); //TF
    int vogno=ogr->getMappedOrthogroupID(v->getName().c_str(),speciesIDNameMap[specID].c_str()); //Target
    
    INTINTMAP* csVars=affectedVars[specID];
    //INTINTMAP* csVars=affectedVars[1];
    
    if(csVars->find(uogno)!=csVars->end())
     {
     cout <<"TF OG " << uogno  << " already used for "<< speciesIDNameMap[specID] << endl;
     return -1;
     }
    if(vogno==-1)
    {
        return 0;
    }
    if((vogno!=-1) && (csVars->find(vogno)!=csVars->end()))
    {
        cout <<"Target OG " << vogno  << " already used for "<< speciesIDNameMap[specID] << endl;
        return -1;
    }
    (*csVars)[vogno]=0;
    //(*csVars)[uogno]=0;
    FactorGraph* csGraph=sdm->getFactorGraph();
    SlimFactor* dFactor=csGraph->getFactorAt(move->getTargetVertex());
    dFactor->mergedMB[move->getSrcVertex()]=0;
    dFactor->mbScore=move->getTargetMBScore();
    dFactor->setMBWts(move->getSrcWeight());
    char ogpair[256];
    sprintf(ogpair,"%d-%d",uogno,vogno);
    string ogpairKey(ogpair);
    cout <<"Attempting move " << ogpairKey << " in species " << speciesIDNameMap[specID] << endl;
    if(edgeConditionMap.find(ogpairKey)==edgeConditionMap.end())
    {
        cout <<"No edge pair " << ogpairKey << endl;
        exit(0);
    } //shilu: edgeConditionMap contains all possible edges
    STRINTMAP* newEdgeStatus=NULL;
    if(affectedOGPairs.find({regi,targeti})==affectedOGPairs.end()) //if(affectedOGPairs.find(ogpairKey)==affectedOGPairs.end())
    {
        newEdgeStatus=new STRINTMAP;
        affectedOGPairs[{regi,targeti}]=newEdgeStatus; //affectedOGPairs[ogpairKey]=newEdgeStatus;
    }
    else
    {
        newEdgeStatus=affectedOGPairs[{regi,targeti}]; //newEdgeStatus=affectedOGPairs[ogpairKey];
    }
    string& species=speciesIDNameMap[specID];
    //(*newEdgeStatus)[species]=1;
    (*newEdgeStatus)[species]=(*newEdgeStatus)[species]+1;
    if((*newEdgeStatus)[species]>1)
    {
        cout <<"Adding edge to same family again for "<< species << " gene " << v->getName() << " TF " << u->getName() << endl;
    }
    return 0;
}*/


/*int
MetaLearner::attemptDelMove(MetaMove* move,map<int,INTINTMAP*>& affectedVars)
{
    int specID=move->getConditionSetInd();
    SpeciesDataManager* sdm=speciesDataSet[speciesIDNameMap[specID]];
    VSET& varSet=sdm->getVariableManager()->getVariableSet();
    Variable* u=varSet[move->getSrcVertex()];
    Variable* v=varSet[move->getTargetVertex()];
    int uogno=ogr->getMappedOrthogroupID(u->getName().c_str(),speciesIDNameMap[specID].c_str());
    int vogno=ogr->getMappedOrthogroupID(v->getName().c_str(),speciesIDNameMap[specID].c_str());
    
    INTINTMAP* csVars=affectedVars[1];
    if(csVars->find(uogno)!=csVars->end())
    {
        return -1;
    }
    if((vogno!=-1) && (csVars->find(vogno)!=csVars->end()))
    {
        return -1;
    }
    FactorGraph* csGraph=sdm->getFactorGraph();
    SlimFactor* dFactor=csGraph->getFactorAt(move->getTargetVertex());
    INTINTMAP_ITER vIter=dFactor->mergedMB.find(move->getSrcVertex());
    if(vIter==dFactor->mergedMB.end())
    {
        cout <<"No mb var " << move->getSrcVertex() << " in mb of " << move->getTargetVertex() << endl;
        return 0;
    }
    (*csVars)[uogno]=0;
    dFactor->mergedMB.erase(vIter->first);
    dFactor->mbScore=move->getTargetMBScore();
    dFactor->setMBWts(move->getSrcWeight());
    if(vogno=-1)
    {
        return 0;
    }
    (*csVars)[vogno]=0;
    char ogpair[256];
    sprintf(ogpair,"%d-%d",uogno,vogno);
    string ogpairKey(ogpair);
    if(edgeConditionMap.find(ogpairKey)==edgeConditionMap.end())
    {
        cout <<"No edge pair " << ogpairKey << endl;
        exit(0);
    }
    INTINTMAP* currEdgeStatus=edgeConditionMap[ogpairKey];
    STRINTMAP* newEdgeStatus=new STRINTMAP;
    affectedOGPairs[ogpairKey]=newEdgeStatus;
    for(INTINTMAP_ITER cIter=currEdgeStatus->begin();cIter!=currEdgeStatus->end();cIter++)
    {
        (*newEdgeStatus)[speciesIDNameMap[cIter->first]]=cIter->second;
    }
    string& species=speciesIDNameMap[specID];
    //(*newEdgeStatus)[species]=0;
    (*newEdgeStatus)[species]=(*newEdgeStatus)[species]-1;
    if((*newEdgeStatus)[species]>0)
    {
        cout <<"Found a duplicate for "<< species << " gene " << v->getName() << " ogkey " << ogpairKey << endl;
    }
    
    return 0;
}*/

int
MetaLearner::dumpAllGraphs(int currK,int foldid)
{
    cout <<"MetaLearner::dumpAllGraphs" << endl;
    char foldoutDirName[1024];
    char aFName[1024];
    //for(map<string,SpeciesDataManager*>::iterator eIter=speciesDataSet.begin();eIter!=speciesDataSet.end();eIter++)
    for(int i=0;i<speciesDataSet.size();i++){
        SpeciesDataManager* sdm=speciesDataSet[i];//eIter->second;
        const char* dirname=sdm->getOutputLoc();
        sprintf(foldoutDirName,"%s/fold%d",dirname,foldid);
        sprintf(aFName,"%s/var_mb_pw_k%d.txt",foldoutDirName,currK);
        ofstream oFile(aFName);
        FactorGraph* fg=sdm->getFactorGraph();
        vector<Variable*>& varSet=sdm->getVariableManager()->getVariableSet();
        vector<SlimFactor*>& factorSet=fg->getAllFactors();
        PotentialManager* potMgr=sdm->getPotentialManager();
        for(int aIter=0;aIter<factorSet.size();aIter++)
        {
            SlimFactor* sFactor=factorSet[aIter];
            potMgr->dumpVarMB_PairwiseFormat(sFactor, oFile, varSet);
        }
        //fg->dumpVarMB_PairwiseFormat(oFile,varSet);
        oFile.close();
        /*sprintf(aFName,"%s/net_ogspace_k%d.txt",foldoutDirName,currK);
        ofstream eFile(aFName);
        vector<SlimFactor*>& factorSet=fg->getAllFactors();
        for(int aIter=0;aIter<factorSet.size();aIter++)
        {
            SlimFactor* sFactor=factorSet[aIter];
            unordered_map<int,double>& mbWts=sFactor->mbWts;
            for(auto mIter=sFactor->mergedMB.begin();mIter!=sFactor->mergedMB.end();mIter++)
            {
                double wtval=mbWts[*mIter];
                int regulatorogid=ogr->getMappedOrthogroupID(varSet[*mIter]->getName().c_str(),speciesIDNameMap[i].c_str());
                int targetogid=ogr->getMappedOrthogroupID(varSet[aIter]->getName().c_str(),speciesIDNameMap[i].c_str());
                //eFile << varSet[mIter->first]->getName()<< "\t"
                //<< varSet[sFactor->vIds[0]]->getName() << "\t" << wtval
                eFile << "\t" << regulatorogid << "\t" << targetogid << "\t" << wtval<< endl;
            }
        }
        eFile.close();*/
    }
    return 0;
}

int
MetaLearner::showModelParameters(int foldid)
{
    char foldoutDirName[1024];
    char aFName[1024];
    //for(map<string,SpeciesDataManager*>::iterator eIter=speciesDataSet.begin();eIter!=speciesDataSet.end();eIter++)
    for(int i=0;i<speciesDataSet.size();i++){
        SpeciesDataManager* sdm=speciesDataSet[i];//eIter->second;
        const char* dirname=sdm->getOutputLoc();
        sprintf(foldoutDirName,"%s/fold%d",dirname,foldid);
        sprintf(aFName,"%s/modelparams.txt",foldoutDirName);
        ofstream oFile(aFName);
        FactorGraph* fg=sdm->getFactorGraph();
        PotentialManager* potMgr=sdm->getPotentialManager();
        vector<SlimFactor*>& slimFactorSet=fg->getAllFactors();
        vector<Variable*>& varSet=sdm->getVariableManager()->getVariableSet();
        for(auto sIter=slimFactorSet.begin();sIter!=slimFactorSet.end();sIter++)
        {
            SlimFactor* sFactor=*sIter;
            /*Potential* sPot=new Potential;
            sPot->setAssocVariable(varSet[sFactor->fId],Potential::FACTOR);
            for(auto mIter=sFactor->mergedMB.begin();mIter!=sFactor->mergedMB.end();mIter++)
            {
                Variable* aVar=varSet[*mIter];
                sPot->setAssocVariable(aVar,Potential::MARKOV_BNKT);
            }
            sPot->potZeroInit();
            //potMgr->populatePotential(sPot,false);
            potMgr->populatePotential_Eff(sPot,false);
            sPot->initMBCovMean();
            double mbcondvar=sPot->getCondVariance();
            double mbbias=sPot->getCondBias();
            unordered_map<int,double>& mbwt=sPot->getCondWeight();*/
            //The weight of the model is really a space filler.
            double mbcondvar=0;
            double mbbias=0;
            unordered_map<int,double> mbwt;
            potMgr->computePotentialMBCovMean(sFactor, mbcondvar, mbbias, mbwt);
            Variable* var=varSet[sFactor->fId];
            oFile<<"Var="<<var->getName()<<"\tWt=-1"<<"\tCondVar="<<mbcondvar << "\tCondBias="<<mbbias<<"\tCondWt=";
            for(auto dIter=mbwt.begin();dIter!=mbwt.end();dIter++)
            {
                if(dIter!=mbwt.begin())
                {
                    oFile <<",";
                }
                Variable* mbVar=varSet[dIter->first];
                oFile<< mbVar->getName()<<"="<<dIter->second;
            }
            oFile << endl;
            mbwt.clear();

        }
        oFile.close();
    }
    return 0;
}


/*
bool
MetaLearner::checkMBSize(int cid, int u, int currK)
{
    bool check=true;
    SpeciesDataManager* sdm=speciesDataSet[cid];
    vector<Variable*>& varSet=sdm->getVariableManager()->getVariableSet();
    Variable* sVar=varSet[u];
    FactorGraph* fg=sdm->getFactorGraph();
    SlimFactor* sFactor=fg->getFactorAt(u);
    //map<string,int>& regulatorSet=sdm->getRegulators();   //not needed?
    if(sFactor->mergedMB.size()>=currK)
    {
        check=false;
    }
    return check;
}*/




int
MetaLearner::generateData(int sampleCnt, int burnin,int foldid)
{
    map<int,ofstream*> newdataFiles;
    gsl_rng* rndgen=gsl_rng_alloc(gsl_rng_default);
    map<int,map<int,Potential*>*> potSet;
    //for(map<string,SpeciesDataManager*>::iterator gIter=speciesDataSet.begin();gIter!=speciesDataSet.end();gIter++)
    for(int i=0;i<speciesDataSet.size();i++){
        //int gid=speciesNameIDMap[gIter->first];
        FactorGraph* fg=speciesDataSet[i]->getFactorGraph();
        VariableManager* varMgr=speciesDataSet[i]->getVariableManager();
        vector<Variable*>& varSet=varMgr->getVariableSet();
        map<int,Potential*>* pSet=new map<int,Potential*>;
        potSet[i]=pSet;
        PotentialManager* potMgr=speciesDataSet[i]->getPotentialManager();
        for(map<string,int>::iterator vIter=subgraphVarSet.begin();vIter!=subgraphVarSet.end();vIter++)
        {
            int vId=vIter->second;
            SlimFactor* sFactor=fg->getFactorAt(vId);
            Potential* sPot=new Potential;
            (*pSet)[vId]=sPot;
            sPot->setAssocVariable(varSet[sFactor->fId],Potential::FACTOR);
            for(auto mIter=sFactor->mergedMB.begin();mIter!=sFactor->mergedMB.end();mIter++)
            {
                Variable* aVar=varSet[*mIter];
                sPot->setAssocVariable(aVar,Potential::MARKOV_BNKT);
            }
            sPot->potZeroInit();
            //potMgr->populatePotential(sPot,false);
            potMgr->populatePotential_Eff(sPot,false);
            sPot->initMBCovMean();
        }
    }
    //for(map<string,SpeciesDataManager*>::iterator sIter=speciesDataSet.begin();sIter!=speciesDataSet.end();sIter++)
    for(int i=0;i<speciesDataSet.size();i++){
        SpeciesDataManager* sdm=speciesDataSet[i];
        char fileName[1024];
        sprintf(fileName,"%s/newsamples_f%d.txt",sdm->getOutputLoc(),foldid);
        ofstream* oFile=new ofstream(fileName);
        newdataFiles[i]=oFile;  //speciesNameIDMap[sIter->first]
        //Write the model file
        sprintf(fileName,"%s/newmodel_f%d.txt",sdm->getOutputLoc(),foldid);
        ofstream mFile(fileName);
        mFile <<"NodeCnt\t"<< subgraphVarSet.size() << endl;
        mFile <<"ContinuousNodes";
        int nodeId=0;
        for(map<string,int>::iterator vIter=subgraphVarSet.begin();vIter!=subgraphVarSet.end();vIter++)
        {
            mFile <<"\t"<< nodeId;
            nodeId++;
        }
        mFile << endl;
        mFile <<"DiscreteNodes"<< endl;
        nodeId=0;
        for(map<string,int>::iterator vIter=subgraphVarSet.begin();vIter!=subgraphVarSet.end();vIter++)
        {
            mFile <<"NodeName="<< vIter->first<< "\tNodeID="<< nodeId<< "\tParents="
            << "\tChildren=" <<"\tValues=0,1,2" << endl;
            nodeId++;
        }
        mFile.close();
    }
    
    for(map<int,map<int,Potential*>*>::iterator psIter=potSet.begin();psIter!=potSet.end();psIter++)
    {
        int did=psIter->first;
        INTDBLMAP initialSample;
        INTDBLMAP saveSample;
        int iter=0;
        map<int,Potential*>* pSet=psIter->second;
        generateInitSample(initialSample,rndgen,pSet);
        
        while(iter<(sampleCnt+burnin))
        {
            int selectFirst=0;
            //We just need to align localmodelid with the outputfile to which the data is generated
            if(iter>=burnin)
            {
                if(iter==burnin)
                {
                    cout <<"Starting to write"<< endl;
                }
                writeEvidenceTo(newdataFiles[did],initialSample);
            }
            for(map<string,int>::iterator vIter=subgraphVarSet.begin(); vIter!=subgraphVarSet.end();vIter++)
            {
                saveSample[vIter->second]=initialSample[vIter->second];
                FactorGraph* fg=speciesDataSet[psIter->first]->getFactorGraph();
                SlimFactor* sFactor=fg->getFactorAt(vIter->second);
                Potential* sPot=(*pSet)[vIter->second];
                double sample_s=sPot->generateSample(initialSample,sFactor->fId,rndgen);
                int siter=0;
                while((sample_s<-160|| sample_s>160) && (siter<50))
                {
                    sample_s=sPot->generateSample(initialSample,sFactor->fId,rndgen);
                    siter++;
                }
                if(sample_s<-160)
                {
                    sample_s=-160;
                }
                else if(sample_s>160)
                {
                    sample_s=160;
                }
                if((isnan(sample_s))||(isinf(sample_s)))
                {
                    cout << "Found Nan/Inf for " << sFactor->fId <<"="<<sample_s  <<" Neighbors: ";
                    for(auto mbIter=sFactor->mergedMB.begin();mbIter!=sFactor->mergedMB.end();mbIter++)
                    {
                        cout << " "<< *mbIter <<"=" << initialSample[*mbIter];
                    }
                    cout << endl;
                }
                initialSample[sFactor->fId]=sample_s;
            }
            iter++;
        }
    }
    gsl_rng_free(rndgen);
    cout <<"Bottom clips " << endl;
    for(map<string,int>::iterator sIter=bottomClip.begin();sIter!=bottomClip.end();sIter++)
    {
        cout << sIter->first.c_str() << " " << sIter->second << endl;
    }
    cout <<"Top clips " << endl;
    for(map<string,int>::iterator sIter=topClip.begin();sIter!=topClip.end();sIter++)
    {
        cout << sIter->first.c_str() << " " << sIter->second << endl;
    }
    return 0;
}



int
MetaLearner::generateInitSample(INTDBLMAP& initialSample,gsl_rng* rndgen,map<int,Potential*>* pSet)
{
    double globalMean=0;
    globalVar=0;
    for(map<string,int>::iterator vIter=subgraphVarSet.begin();vIter!=subgraphVarSet.end();vIter++)
    {
        int vId=vIter->second;
        Potential* sPot=(*pSet)[vId];
        cout <<vId <<"\t" << sPot->getCondBias() <<"\t" << sPot->getCondVariance() << endl;
        globalMean=globalMean+sPot->getCondBias();
        globalVar=globalVar+sPot->getCondVariance();
    }
    double norm=(double)(subgraphVarSet.size()*speciesDataSet.size());
    globalMean=globalMean/norm;
    globalVar=globalVar/norm;
    for(map<string,int>::iterator vIter=subgraphVarSet.begin();vIter!=subgraphVarSet.end();vIter++)
    {
        double sampleval=gsl_ran_gaussian(rndgen,sqrt(globalVar));
        sampleval=sampleval+globalMean;
        initialSample[vIter->second]=sampleval;
    }
    return 0;
}


int
MetaLearner::writeEvidenceTo(ofstream* oFile,INTDBLMAP& sample)
{
    int nodeId=0;
    for(map<string,int>::iterator dIter=subgraphVarSet.begin();dIter!=subgraphVarSet.end();dIter++)
    {
        if(dIter!=subgraphVarSet.begin())
        {
            (*oFile)<<"\t";
        }
        double val=sample[dIter->second];
        if(val<-160)
        {
            val=-160;
            if(bottomClip.find(dIter->first)==bottomClip.end())
            {
                bottomClip[dIter->first]=1;
            }
            else
            {
                bottomClip[dIter->first]=bottomClip[dIter->first]+1;
            }
        }
        else if(val > 160)
        {
            val=160;
            if(topClip.find(dIter->first)==topClip.end())
            {
                topClip[dIter->first]=1;
            }
            else
            {
                topClip[dIter->first]=topClip[dIter->first]+1;
            }
        }
        (*oFile)<< nodeId << "=["<< exp(val) << "]";
        nodeId++;
    }
    (*oFile) << endl;
    return 0;
}
