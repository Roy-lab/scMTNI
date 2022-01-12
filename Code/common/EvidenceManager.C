#include <cassert>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
#include <math.h>
#include <string>
#include "Error.H"
#include "Variable.H"
#include "VariableManager.H"
#include "Evidence.H"
#include "EvidenceManager.H"
#include "Utils.H"

EvidenceManager::EvidenceManager()
{
	foldCnt=1;
	preRandomizeSplit=false;
	randseed=-1;
    testdataSize=0;
    traindataSize=0;
    dataMat = NULL;
}

EvidenceManager::~EvidenceManager()
{
    if (dataMat!=NULL)
    {
        delete dataMat;
    }
}

int
EvidenceManager::deleteData()
{
    delete dataMat;
    return 0;
}

int
EvidenceManager::setVariableManager(VariableManager* aPtr)
{
	vMgr=aPtr;
	return 0;
}

VariableManager*
EvidenceManager::getVariableManager()
{
    return vMgr;
}

Error::ErrorCode
EvidenceManager::loadEvidenceFromFile(const char* inFName)
{
	ifstream inFile(inFName);
	char buffer[80000];
	while(inFile.good())
	{
		inFile.getline(buffer,80000);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		else if(strchr(buffer,'#')!=NULL)
		{
			continue;
		}
		//All the evidences for each variable are stored in a map, indexed by the varId
		EMAP* evidMap=new EMAP;
		char* tok=strtok(buffer,"\t");
		//The toks take the form of varid and value
		while(tok!=NULL)
		{
			Evidence* evid;
			if(populateEvidence(&evid,tok)==-1)
			{
				cout <<"Error while populating evidence " << endl;
				return Error::DATAFILE_ERR;
			}
			(*evidMap)[evid->getAssocVariable()]=evid;
			tok=strtok(NULL,"\t");
		}
		//if(evidenceSet.size()<=20)
		//{
			evidenceSet.push_back(evidMap);
		//}
	}
	inFile.close();
	cout <<"Read " << evidenceSet.size() << " different datapoints " << endl;
	return Error::SUCCESS;
}

Error::ErrorCode
EvidenceManager::loadEvidenceFromTable(vector<string>& inputTable)
{
	/* Reads gene names and expression levels from a tab-separated file, where the first row is assumed (for now)
	   to be headers. In subsequent rows, the first column is a gene name and remaining columns are expression levels.
	   This method does what VariableManager::readVariables() and EvidenceManager::loadEvidenceFromFile_Continuous()
	   do for .model and .data files.
	*/

	/* For each sample, an EMAP (Evidence Map) is created that has an Evidence object per variable (gene) indexed by
	   variable ID. These EMAPs are saved in evidenceSet, a vector of EMAP*.
	*/
    
	// First, let's collect the values per sample. The table has them per gene...
	//vector<vector<double>> valuesPerSample;
	int nodeCount = inputTable.size();
	// How many samples do we have?
	vector<string> substrs = Utils::split(inputTable[0], '\t');
    //cout <<"Evidence " <<substrs[0] << endl;
	int sampleCount = substrs.size() - 1;  // not counting the gene name
    cout <<"EvidenceManager::loadEvidenceFromTable sampleCount=" << sampleCount << " nodeCount=" << nodeCount << endl;
	//valuesPerSample.reserve(sampleCount);
	/*for (int sample = 0; sample < sampleCount; sample++)
	{
		vector<double> sampleVec (nodeCount);
		valuesPerSample.push_back(sampleVec);
	}*/
    dataMat=new Matrix(nodeCount,sampleCount);
	int nodeNum = 0;
    vMgr = new VariableManager;
	for (auto line : inputTable)
	{
		substrs = Utils::split(line, '\t');
        // set variable name:
        string nodeName = substrs[0];
        vMgr->readVariablesFromData(nodeName,nodeNum);
		for (int sample = 0; sample < sampleCount; sample++)
		{
            double varVal =stod(substrs[sample + 1]);
            //valuesPerSample[sample][nodeNum] =varVal;
            dataMat->setValue(varVal,nodeNum,sample);
		}
		++nodeNum;
	}
    substrs.clear();
    /*
	// Now, let's process the data
	//double epsilon = numeric_limits<double>::min();  // so we don't try to take log of zero
	for (int sample = 0; sample < sampleCount; sample++)
	{
		EMAP* evidMap = new EMAP;
		for (int nodeNum = 0; nodeNum < nodeCount; nodeNum++)
		{
			Evidence* evidence = new Evidence;
			evidence->assocVariable(nodeNum);
			// double varVal = log(valuesPerSample[sample][nodeNum] + epsilon);
			double varVal = valuesPerSample[sample][nodeNum];
			evidence->setEvidVal(varVal);
			if (maxAbsValues.find(nodeNum) == maxAbsValues.end())
			{
				maxAbsValues[nodeNum] = fabs(varVal);
			}
			else
			{
				if (maxAbsValues[nodeNum] < fabs(varVal))
				{
					maxAbsValues[nodeNum] = fabs(varVal);
				}
			}
			(*evidMap)[nodeNum] = evidence;
            //cout << "data[" <<nodeNum << "," <<sample << "]=" <<dataMat->getValue(nodeNum,sample) << endl;
		}
		evidenceSet.push_back(evidMap);
	}
	cout <<"Read " << evidenceSet.size() << " different datapoints" << endl;*/
	return Error::SUCCESS;
}

Error::ErrorCode
EvidenceManager::loadEvidenceFromFile_Continuous(const char* inFName)
{
	ifstream inFile(inFName);
	char buffer[160000];
	int lineNo=0;
	while(inFile.good())
	{
		inFile.getline(buffer,160000);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		else if(strchr(buffer,'#')!=NULL)
		{
			continue;
		}
		if(lineNo>=500)
		{
			lineNo++;
		//	continue;
		}

		//All the evidences for each variable are stored in a map, indexed by the varId
		EMAP* evidMap=new EMAP;
		char* tok=strtok(buffer,"\t");
		//The toks take the form of varid and value
		while(tok!=NULL)
		{
			Evidence* evid;
			if(populateEvidence_Continuous(&evid,tok)==-1)
			{
				cout <<"Error while populating evidence " << endl;
				return Error::DATAFILE_ERR;
			}
			(*evidMap)[evid->getAssocVariable()]=evid;
			tok=strtok(NULL,"\t");
		}
		//if(evidenceSet.size()<=20)
		//{
			evidenceSet.push_back(evidMap);
		//}
	}
	inFile.close();
	cout <<"Read " << evidenceSet.size() << " different datapoints " << endl;
	//normalize();
	return Error::SUCCESS;
}

int
EvidenceManager::normalize()
{	
	for(int e=0;e<evidenceSet.size();e++)
	{
		EMAP* evidMap=evidenceSet[e];
		for(EMAP_ITER eIter=evidMap->begin();eIter!=evidMap->end();eIter++)
		{
			Evidence* evid=eIter->second;
			double aVal=evid->getEvidVal();
			double maxVal=maxAbsValues[eIter->first];
			aVal=aVal/maxVal;
			evid->setEvidVal(aVal);
		}
	}
	return 0;
}

//We create a matrix of randomized evidence, where each evidence has some value
//for the random variables. We populate the matrix one random variable at a time.
//We first generate a vector of permuted indices, in the range of 0 to the total
//number of evidences. Then we populate the part of the matrix associated with
//this variable by querying values from the original matrix in the order specified
//by the permuted indices

int
EvidenceManager::randomizeEvidence(gsl_rng* r)
{
	//First create all the evidence sets
	for(int i=0;i<evidenceSet.size();i++)
	{
		EMAP* evidMap=new EMAP;
		randEvidenceSet.push_back(evidMap);
	}
	//Populate variable wise
	vector<Variable*>& variableSet=vMgr->getVariableSet();
	int* randInds=new int[evidenceSet.size()];
	for(int vIter=0;vIter<variableSet.size();vIter++)
	{
		//generate a random vector of indices ranging from 0 to evidenceSet.size()-1
		populateRandIntegers(r,randInds,randEvidenceSet.size());	
		for(int i=0;i<evidenceSet.size();i++)
		{
			EMAP* evidMap=evidenceSet[randInds[i]];
			EMAP* randEvidMap=randEvidenceSet[i];
			Evidence* evid=(*evidMap)[vIter];
			(*randEvidMap)[vIter]=evid;
		}
	}
	return 0;
}

int 
EvidenceManager::getNumberOfEvidences()
{
	return evidenceSet.size();
}

EMAP* 
EvidenceManager::getEvidenceAt(int evId)
{
	if((evId>=evidenceSet.size()) && (evId<0))
	{
		return NULL;
	}
	return evidenceSet[evId];
}

// return data
Matrix*
EvidenceManager::getData()
{
    return dataMat;
}

EMAP* 
EvidenceManager::getRandomEvidenceAt(int evId)
{
	if((evId>=randEvidenceSet.size()) && (evId<0))
	{
		return NULL;
	}
	return randEvidenceSet[evId];
}


int
EvidenceManager::addEvidence(EMAP* evidSet)
{
	evidenceSet.push_back(evidSet);
	return 0;
}


int
EvidenceManager::addToEvidence(int eSetID, int vId, INTDBLMAP& evidData)
{
	EMAP* emap=evidenceSet[eSetID];
	Evidence* evid=new Evidence;
	evid->setData(evidData);
	(*emap)[vId]=evid;
	return 0;
}

int
EvidenceManager::dumpEvidenceSet(ostream& oFile)
{
	for(int i=0;i<evidenceSet.size();i++)
	{
		EMAP* evidMap=evidenceSet[i];
		for(EMAP_ITER eIter=evidMap->begin();eIter!=evidMap->end();eIter++)
		{
			Evidence* evid=eIter->second;
			if(eIter!=evidMap->begin())
			{
				oFile<<"\t";
			}
			evid->dumpEvidence(oFile);
		}
		oFile << endl;
	}
	return 0;
}

int
EvidenceManager::getMLSettings(ostream& oFile)
{
	for(int i=0;i<evidenceSet.size();i++)
	{
		EMAP* evidMap=evidenceSet[i];
		if(i==0)
		{
			for(EMAP_ITER eIter=evidMap->begin();eIter!=evidMap->end();eIter++)
			{
				if(eIter!=evidMap->begin())
				{
					oFile<<"\t";
				}
				oFile<< eIter->first;
			}
			oFile << endl;
		}

		for(EMAP_ITER eIter=evidMap->begin();eIter!=evidMap->end();eIter++)
		{
			if(eIter!=evidMap->begin())
			{
				oFile<<"\t";
			}
			Evidence* evid=eIter->second;
			oFile << evid->getMLVal();
		}
		oFile << endl;
	}
	return 0;
}


int 
EvidenceManager::setFoldCnt(int f)
{
	foldCnt=f;
	return 0;
}

int
EvidenceManager::generateValidationSet(const char* vFName, int vSetSize,gsl_rng* r)
{
	ifstream inFile(vFName);
	if(inFile.good())
	{
		char buffer[256];
		while(inFile.good())
		{
			inFile.getline(buffer,255);
			if(strlen(buffer)<=0)
			{
				continue;
			}
			int dId=atoi(buffer);
			validationIndex[dId]=0;
		}
		inFile.close();
	}
	else
	{
		populateRandIntegers(r,validationIndex,evidenceSet.size(),vSetSize);
		ofstream oFile(vFName);
		for(INTINTMAP_ITER vIter=validationIndex.begin();vIter!=validationIndex.end();vIter++)
		{
			oFile << vIter->first << endl;
		}
		oFile.close();
	}
	return 0;
}

int 
EvidenceManager::setPreRandomizeSplit()
{
	preRandomizeSplit=true;
	return 0;
}

int 
EvidenceManager::setPreRandomizeSplitSeed(int seed)
{
	randseed=seed;
	return 0;
}

/*
int 
EvidenceManager::splitData(int s)
{
	int testSetSize=(evidenceSet.size()-validationIndex.size())/foldCnt;
	int testStartIndex=s*testSetSize;
	int testEndIndex=(s+1)*testSetSize;
	if(s==foldCnt-1)
	{
		testEndIndex=evidenceSet.size()-validationIndex.size();
	}
	trainIndex.clear();
	testIndex.clear();
	int m=0;
	for(int i=0;i<evidenceSet.size();i++)
	{
		if(validationIndex.find(i)!=validationIndex.end())
		{
			continue;
		}
		if((m>=testStartIndex) && (m<testEndIndex))
		{
			testIndex[i]=0;
		}
		else
		{
			trainIndex[i]=0;
		}
		m++;
	}
	return 0;
}*/


int 
EvidenceManager::splitData(int s)
{
     
	int testSetSize=(dataMat->getColCnt()-validationIndex.size())/foldCnt;  //(evidenceSet.size()-validationIndex.size())/foldCnt
    if(foldCnt==1)
    {
        //testStartIndex=-1;
        //testEndIndex=-1;
        // only testdata, no training or validation data
        testdataSize=testSetSize;
        //cout <<"EvidenceManager::splitData preRandomizeSplit=" << preRandomizeSplit << " testdataSize=" << testdataSize << " dataMat->getColCnt()=" << dataMat->getColCnt()<<endl;
        return 0;
    }
	int testStartIndex=s*testSetSize;
	int testEndIndex=(s+1)*testSetSize;
	if(s==foldCnt-1)
	{
		testEndIndex=evidenceSet.size()-validationIndex.size();
	}
	trainIndex.clear();
	testIndex.clear();
	int m=0;
	int* randInds=NULL;
	if(preRandomizeSplit)
	{
		randInds=new int[evidenceSet.size()];
		//generate a random vector of indices ranging from 0 to evidenceSet.size()-1
		gsl_rng* r=gsl_rng_alloc(gsl_rng_default);
		if(randseed<0)
		{
			randseed = rand();
		}
		gsl_rng_set(r,randseed);
		populateRandIntegers(r,randInds,evidenceSet.size());	
		gsl_rng_free(r);
		cout <<"Random seed " << randseed << endl;
	}
	for(int i=0;i<evidenceSet.size();i++)
	{
		int eInd=i;
		if(randInds!=NULL)
		{
			eInd=randInds[i];
		}
		if(validationIndex.find(eInd)!=validationIndex.end())
		{
			continue;
		}
		if((m>=testStartIndex) && (m<testEndIndex))
		{
            testIndex.push_back(eInd); //testIndex[eInd]=0;
		}
		else
		{
            trainIndex.push_back(eInd); //trainIndex[eInd]=0;
		}
		m++;
	}
	if(preRandomizeSplit)
	{
		delete[] randInds;
	}
	return 0;
}


/*INTINTMAP&
EvidenceManager::getTrainingSet()
{
	return trainIndex;
}
INTINTMAP&
EvidenceManager::getTestSet()
{
	return testIndex;
}*/

vector<int>&
EvidenceManager::getTrainingSet()
{
    return trainIndex;
}

int
EvidenceManager::getTrainSetSize()
{
    return traindataSize;
}


vector<int>&
EvidenceManager::getTestSet()
{
    return testIndex;
}

int
EvidenceManager::getTestSetSize()
{
    return testdataSize;
}

INTINTMAP&
EvidenceManager::getValidationSet()
{	
	return validationIndex;
}

int 
EvidenceManager::populateEvidence(Evidence** evid,const char* evidStr)
{
	//first check for validity of evidStr
	if(strchr(evidStr,'=')==NULL)
	{
		return -1;
	}
	*evid=new Evidence;
	
	INTDBLMAP evidData;
	int currInd=0;
	int ttInd=0;
	int tokId=0;
	char tempTok[256];
	while(evidStr[currInd]!='\0')
	{
		if((evidStr[currInd]=='=') || 
		   (evidStr[currInd]==']') ||
		   (evidStr[currInd]==',')
		  )
		{
			tempTok[ttInd]='\0';
			ttInd=0;
			if(tokId==0)
			{
				//This is the variable
				int vId=atoi(tempTok);
				Variable* var=vMgr->getVariableAt(vId);
				var->initEvidence(evidData);
				(*evid)->assocVariable(vId);
			}
			else
			{
				char* pos=strchr(tempTok,'|');
				//Hard evidence
				if(pos==NULL)
				{
					int varVal=atoi(tempTok);
					if(evidData.find(varVal)==evidData.end())
					{
						cout <<"No value "<< varVal << " in the domain of  a variable" << endl;
						return -1;
					}
					evidData[varVal]=1.0;
					(*evid)->setType(Evidence::HARD);
				}
				else
				{
					*pos='\0';
					int varVal=atoi(tempTok);
					double varValProb=atof(pos+1);
					if(evidData.find(varVal)==evidData.end())
					{
						cout <<"No value "<< varVal << " in the domain of  a variable" << endl;
						return -1;
					}
					evidData[varVal]=varValProb;
					//Will be setting it multiple times but its ok for now.
					(*evid)->setType(Evidence::SOFT);
				}
			}
			tokId++;
		}
		else if(evidStr[currInd]!='[')
		{
			tempTok[ttInd]=evidStr[currInd];
			ttInd++;
		}
		currInd++;
	}
	(*evid)->setData(evidData);
	return 0;
}


int 
EvidenceManager::populateEvidence_Continuous(Evidence** evid,const char* evidStr)
{
	//first check for validity of evidStr
	if(strchr(evidStr,'=')==NULL)
	{
		return -1;
	}
	*evid=new Evidence;
	
	int currInd=0;
	int ttInd=0;
	int tokId=0;
	char tempTok[256];
	while(evidStr[currInd]!='\0')
	{
		if((evidStr[currInd]=='=') || 
		   (evidStr[currInd]==']') ||
		   (evidStr[currInd]==',')
		  )
		{
			tempTok[ttInd]='\0';
			ttInd=0;
			if(tokId==0)
			{
				//This is the variable
				int vId=atoi(tempTok);
				Variable* var=vMgr->getVariableAt(vId);
				(*evid)->assocVariable(vId);
			}
			else
			{
				char* pos=strchr(tempTok,'|');
				//Hard evidence
				if(pos==NULL)
				{
					double varVal=log(atof(tempTok));
					(*evid)->setEvidVal(varVal);
					int vId=(*evid)->getAssocVariable();
					if(maxAbsValues.find(vId)==maxAbsValues.end())
					{
						maxAbsValues[vId]=fabs(varVal);
					}
					else
					{
						double mVal=maxAbsValues[vId];
						if(mVal<fabs(varVal))
						{
							maxAbsValues[vId]=fabs(varVal);
						}
					}
				}
			}
			tokId++;
		}
		else if(evidStr[currInd]!='[')
		{
			tempTok[ttInd]=evidStr[currInd];
			ttInd++;
		}
		currInd++;
	}
	return 0;
}

int
EvidenceManager::dumpSummaryStat(ostream& oFile)
{
	//Need an evidence like object but to store the frequency over all evidences
	//rather than a single evidence
	map<int,INTDBLMAP*> summary;
	map<int,double> normFactors;
	for(int i=0;i<evidenceSet.size();i++)
	{
		EMAP* evidMap=evidenceSet[i];
		for(EMAP_ITER eIter=evidMap->begin();eIter!=evidMap->end();eIter++)
		{
			INTDBLMAP* evCnt=NULL;
			if(summary.find(eIter->first)==summary.end())
			{
				evCnt=new INTDBLMAP;
				summary[eIter->first]=evCnt;
			}
			else
			{
				evCnt=summary[eIter->first];
			}
			//Get data and add to evCnt
			INTDBLMAP& data=eIter->second->getData();
			for(INTDBLMAP_ITER idIter=data.begin();idIter!=data.end();idIter++)
			{
				if(evCnt->find(idIter->first)==evCnt->end())
				{
					(*evCnt)[idIter->first]=idIter->second;
				}
				else
				{
					(*evCnt)[idIter->first]=(*evCnt)[idIter->first]+idIter->second;
				}
				//Add the normalization factor for all the freq or exp. freq cnts
				if(normFactors.find(eIter->first)==normFactors.end())
				{
					normFactors[eIter->first]=idIter->second;
				}
				else
				{
					normFactors[eIter->first]=normFactors[eIter->first]+idIter->second;
				}
			}
		}
	}

	//Now iterate over the evidence summary, normalize and display values
	for(map<int,INTDBLMAP*>::iterator aIter=summary.begin();aIter!=summary.end();aIter++)
	{
		double normConst=normFactors[aIter->first];
		INTDBLMAP* evCnt=aIter->second;
		oFile <<"Distribution of "<< aIter->first;
		for(INTDBLMAP_ITER idIter=evCnt->begin();idIter!=evCnt->end();idIter++)
		{
			oFile << " " << idIter->first<<"=" << idIter->second/normConst;
		}
		oFile << endl;
	}
	return 0;
}


int 
EvidenceManager::populateRandIntegers(gsl_rng* r, int* randInds,int size)
{
	double step=1.0/(double)size;
	map<int,int> usedInit;
	for(int i=0;i<size;i++)
	{
		double rVal=gsl_ran_flat(r,0,1);
		int rind=(int)(rVal/step);
		while(usedInit.find(rind)!=usedInit.end())
		{
			rVal=gsl_ran_flat(r,0,1);
			rind=(int)(rVal/step);
		}
		usedInit[rind]=0;
		randInds[i]=rind;
	}
	usedInit.clear();
	return 0;
}


int 
EvidenceManager::populateRandIntegers(gsl_rng* r, INTINTMAP& randInds,int size, int subsetsize)
{
	double step=1.0/(double)size;
	for(int i=0;i<subsetsize;i++)
	{
		double rVal=gsl_ran_flat(r,0,1);
		int rind=(int)(rVal/step);
		while(randInds.find(rind)!=randInds.end())
		{
			rVal=gsl_ran_flat(r,0,1);
			rind=(int)(rVal/step);
		}
		randInds[rind]=0;
	}
	return 0;
}

