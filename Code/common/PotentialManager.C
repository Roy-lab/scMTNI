/*
 *  scMTNI: single-cell Multi-Task learning Network Inference
 *   Copyright 2022 Shilu Zhang (szhang256@wisc.edu) and  Sushmita Roy (sroy@biostat.wisc.edu)
 *
 *   Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, 
 *   including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, 
 *   subject to the following conditions:
 *
 *   The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 *
 *   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
 *   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *   */

#include <iostream>
#include <cstring>
#include <math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "CommonTypes.H"
#include "Error.H"
#include "Variable.H"
#include "Potential.H"
#include "Evidence.H"
#include "EvidenceManager.H"
#include "SlimFactor.H"
#include "PotentialManager.H"
#include <chrono>
#include "Utils.H"
#include "VariableManager.H"


using namespace std::chrono;

PotentialManager::PotentialManager()
{
	ludecomp=NULL;
	perm=NULL;
    data=NULL;
    meanMat=NULL;
    covMat=NULL;
    testdataSize=0;
}

PotentialManager::~PotentialManager()
{
	if(ludecomp!=NULL)
	{
		gsl_matrix_free(ludecomp);
	}
	if(perm!=NULL)
	{
		gsl_permutation_free(perm);
	}
	/*for(map<int,Potential*>::iterator pIter=potBuffer.begin();pIter!=potBuffer.end();pIter++)
	{
		delete pIter->second;
	}
	potBuffer.clear();*/
    if (data!=NULL)
    {
        delete data;
    }
		if (meanMat!=NULL)
    {
        delete meanMat;
    }
    if (covMat!=NULL)
    {
        delete covMat;
    }
}

int
PotentialManager::deleteData()
{
    if (data!=NULL)
    {
        delete data;
    }
    return 0;
}


int 
PotentialManager::setEvidenceManager(EvidenceManager* aPtr)
{
	evMgr=aPtr;
	return 0;
}

/*int
PotentialManager::setRestrictedNeighborSet(map<int,Variable*>& rSet)
{
	for(map<int,Variable*>::iterator vIter=rSet.begin();vIter!=rSet.end();vIter++)
	{
		restrictedNeighborSet[vIter->first]=vIter->second;
	}
	return 0;
}*/

int
PotentialManager::setOutputDir(const char* aDirName)
{
	strcpy(outputDir,aDirName);
	return 0;
}

//added from EvidenceManager
Error::ErrorCode
PotentialManager::loadEvidenceFromTable(vector<string>& inputTable)
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
    data=new Matrix(nodeCount,sampleCount);
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
            data->setValue(varVal,nodeNum,sample);
        }
        ++nodeNum;
    }
    testdataSize=sampleCount;
    cout <<"PotentialManager::loadEvidenceFromTable varCount=" << data->getRowCnt() << " sampleCount=" << data->getColCnt()  << endl;
    substrs.clear();
    return Error::SUCCESS;
}

Error::ErrorCode
PotentialManager::loadEvidenceFromTable(vector<string>& inputTable,unordered_set<string>& inputVariableNames)
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
    int nodeCount = inputVariableNames.size();
    // How many samples do we have?
    vector<string> substrs = Utils::split(inputTable[0], '\t');
    //cout <<"Evidence " <<substrs[0] << endl;
    int sampleCount = substrs.size() - 1;  // not counting the gene name
    data=new Matrix(nodeCount,sampleCount);
    int nodeNum = 0;
    vMgr = new VariableManager;
    for (auto line : inputTable)
    {
        substrs = Utils::split(line, '\t');
        // set variable name:
        string nodeName = substrs[0];
        if(inputVariableNames.find(nodeName)!=inputVariableNames.end())
        {
            vMgr->readVariablesFromData(nodeName,nodeNum);
            for (int sample = 0; sample < sampleCount; sample++)
            {
                double varVal =stod(substrs[sample + 1]);
                //valuesPerSample[sample][nodeNum] =varVal;
                data->setValue(varVal,nodeNum,sample);
            }
            ++nodeNum;
        }/*else{
            cout << nodeName << " is not in GeneList" << endl;
        }*/
        
    }
    testdataSize=sampleCount;
    cout <<"PotentialManager::loadEvidenceFromTable varCount=" << data->getRowCnt() << " sampleCount=" << data->getColCnt()  << endl;
    substrs.clear();
    return Error::SUCCESS;
}

VariableManager*
PotentialManager::getVariableManager()
{
    return vMgr;
}

/*
int
PotentialManager::initPooled()
{
	//INTINTMAP& trainEvidSet=evMgr->getTrainingSet();
	INTINTMAP& testEvidSet=evMgr->getTestSet();
	//estimateAllMeanCov(false,globalMean,globalCovar,testEvidSet,NULL,NULL);
	estimateAllMeanCov(false,globalMean,globalCovar,testEvidSet);
	ludecomp=gsl_matrix_alloc(MAXFACTORSIZE_ALLOC,MAXFACTORSIZE_ALLOC);
	perm=gsl_permutation_alloc(MAXFACTORSIZE_ALLOC);
	return 0;
}

int
PotentialManager::init()
{
	char mFName[1024];
	char sdFName[1024];
	sprintf(mFName,"%s/gauss_mean.txt",outputDir);
	sprintf(sdFName,"%s/gauss_std.txt",outputDir);
	ifstream inFile(mFName);
	if(inFile.good())
	{
		readAllMeanCov(mFName,sdFName);
	}
	else
	{
		INTINTMAP& trainEvidSet=evMgr->getTrainingSet();
		estimateAllMeanCov(false,globalMean,globalCovar,trainEvidSet,mFName,sdFName);
	}
	ludecomp=gsl_matrix_alloc(MAXFACTORSIZE_ALLOC,MAXFACTORSIZE_ALLOC);
	perm=gsl_permutation_alloc(MAXFACTORSIZE_ALLOC);
	return 0;
}
*/

int
PotentialManager::init(int f)
{
	char mFName[1024];
	char sdFName[1024];
	sprintf(mFName,"%s/gauss_mean_%d.txt",outputDir,f);
	sprintf(sdFName,"%s/gauss_std_%d.txt",outputDir,f);
	ifstream inFile(mFName);
	if(inFile.good())
	{
		readAllMeanCov(mFName,sdFName);
	}
	else
	{
		//INTINTMAP& trainEvidSet=evMgr->getTrainingSet();
		vector<int>& trainEvidSet=evMgr->getTestSet();  //INTINTMAP& trainEvidSet=evMgr->getTestSet();
		//estimateAllMeanCov(false,globalMean,globalCovar,trainEvidSet,mFName,sdFName);
		//estimateAllMeanCov(false,globalMean,globalCovar,trainEvidSet);
        estimateAllMeanCov_Eff(false,trainEvidSet);
	}
	ludecomp=gsl_matrix_alloc(MAXFACTORSIZE_ALLOC,MAXFACTORSIZE_ALLOC);
	perm=gsl_permutation_alloc(MAXFACTORSIZE_ALLOC);
	return 0;
}

int
PotentialManager::init()  // one fold only!
{
    //int trainEvidSetSize=evMgr->getTestSetSize();  //INTINTMAP& trainEvidSet=evMgr->getTestSet();
    //cout <<"PotentialManager::init() trainEvidSetSize=" << trainEvidSetSize << endl;
    //data is the data matrix which will have the variable by sample information
    //data=evMgr->getData();
    int varCnt=data->getRowCnt();//number of variables
    meanMat=new Matrix(varCnt,1);
    meanMat->setAllValues(0);
    covMat=new Matrix(varCnt,varCnt);
    covMat->setAllValues(-1);
    // estimate mean later
    //estimateAllMeanCov_Eff_onefold(false,trainEvidSetSize);
    //ludecomp=gsl_matrix_alloc(MAXFACTORSIZE_ALLOC,MAXFACTORSIZE_ALLOC);
    //perm=gsl_permutation_alloc(MAXFACTORSIZE_ALLOC);
    //cout << "PotentialManager::init one fold: sample=" << data->getColCnt() << " varibleCnt=" << varCnt << " MAXFACTORSIZE_ALLOC="<<MAXFACTORSIZE_ALLOC <<endl;
    return 0;
}


/*int
PotentialManager::initValidationSet(int validSetSize)
{
	char mFName[1024];
	char sdFName[1024];
	sprintf(mFName,"%s/gauss_mean_v%d.txt",outputDir,validSetSize);
	sprintf(sdFName,"%s/gauss_std_v%d.txt",outputDir,validSetSize);
	ifstream inFile(mFName);
	if(inFile.good())
	{
		readAllMeanCov(mFName,sdFName);
	}
	else
	{
		INTINTMAP& validationSet=evMgr->getValidationSet();
		estimateAllMeanCov(false,globalMean,globalCovar,validationSet,mFName,sdFName);
	}
	ludecomp=gsl_matrix_alloc(MAXFACTORSIZE_ALLOC,MAXFACTORSIZE_ALLOC);
	perm=gsl_permutation_alloc(MAXFACTORSIZE_ALLOC);
	return 0;
}


int
PotentialManager::initValidationSet(int validSetSize,int dId)
{
	char mFName[1024];
	char sdFName[1024];
	sprintf(mFName,"%s/gauss_mean_v%d.txt",outputDir,validSetSize);
	sprintf(sdFName,"%s/gauss_std_v%d.txt",outputDir,validSetSize);
	INTINTMAP& validationSet=evMgr->getValidationSet();
	estimateAllMeanCov(false,globalMean,globalCovar,validationSet,mFName,sdFName,dId);
	ludecomp=gsl_matrix_alloc(MAXFACTORSIZE_ALLOC,MAXFACTORSIZE_ALLOC);
	perm=gsl_permutation_alloc(MAXFACTORSIZE_ALLOC);
	return 0;
}*/

int
PotentialManager::reset()
{
    /*if (data!=NULL)
    {
        delete data;
    }*/
    if (meanMat!=NULL)
    {
        delete meanMat;
    }
    if (covMat!=NULL)
    {
        delete covMat;
    }
	/*globalMean.clear();
	for(map<int,INTDBLMAP*>::iterator cIter=globalCovar.begin();cIter!=globalCovar.end();cIter++)
	{
		cIter->second->clear();
		delete cIter->second;
	}
	globalCovar.clear();*/
	if(ludecomp!=NULL)
	{
		gsl_matrix_free(ludecomp);
	}
	if(perm!=NULL)
	{
		gsl_permutation_free(perm);
	}
	/*for(map<int,Potential*>::iterator pIter=potBuffer.begin();pIter!=potBuffer.end();pIter++)
	{
		delete pIter->second;
	}
	potBuffer.clear();*/
	return 0;
}

/*int
PotentialManager::initRandom()
{
	vector<int>& trainEvidSet=evMgr->getTestSet();
	estimateAllMeanCov(true,globalMean_Rand,globalCovar_Rand,trainEvidSet,NULL,NULL);
	return 0;
}

int
PotentialManager::estimateAllMeanCov(bool random, INTDBLMAP& gMean, map<int,INTDBLMAP*>& gCovar,INTINTMAP& trainEvidSet, const char* mFName, const char* sdFName,int leaveOutData)
{
	ofstream mFile;
	ofstream sdFile;

	if(!random)
	{
		if((mFName!=NULL) && (sdFName!=NULL))
		{
			mFile.open(mFName);
			sdFile.open(sdFName);
		}
	}

	int evidCnt=trainEvidSet.size();
	if(leaveOutData!=-1)
	{
		evidCnt=evidCnt-1;
	}
	//First get the mean and then the variance
	int dId=0;
	for(INTINTMAP_ITER eIter=trainEvidSet.begin();eIter!=trainEvidSet.end();eIter++)
	{
		if(dId==leaveOutData)
		{	
			dId++;
			continue;
		}
		EMAP* evidMap=NULL;
		if(random)
		{
			evidMap=evMgr->getRandomEvidenceAt(eIter->first);
		}
		else
		{
			evidMap=evMgr->getEvidenceAt(eIter->first);
		}
		for(EMAP_ITER vIter=evidMap->begin();vIter!=evidMap->end(); vIter++)
		{
			int vId=vIter->first;
			Evidence* evid=vIter->second;
			double val=evid->getEvidVal();
			if(gMean.find(vId)==gMean.end())
			{
				gMean[vId]=val;
			}
			else
			{
				gMean[vId]=gMean[vId]+val;
			}
		}
		dId++;	
	}
	//Now estimate the mean
	for(INTDBLMAP_ITER idIter=gMean.begin();idIter!=gMean.end();idIter++)
	{
		if(idIter->first==176)
		{
			//cout <<"Stop here: Variable " << idIter->first << " mean " << idIter->second << endl;
		}
		idIter->second=idIter->second/(double) evidCnt;
		if(!random)
		{
			if(mFile.good())
			{
				mFile<<idIter->first<<"\t" << idIter->second<< endl;
			}
		}
	}
	int covPair=0;
	//Now the variance
	for(INTINTMAP_ITER eIter=trainEvidSet.begin();eIter!=trainEvidSet.end();eIter++)
	{
		EMAP* evidMap=NULL;
		if(random)
		{
			evidMap=evMgr->getRandomEvidenceAt(eIter->first);
		}
		else
		{
			evidMap=evMgr->getEvidenceAt(eIter->first);
		}
		for(EMAP_ITER vIter=evidMap->begin();vIter!=evidMap->end(); vIter++)
		{
			int vId=vIter->first;
			Evidence* evid=vIter->second;
			double vval=evid->getEvidVal();
			double vmean=gMean[vId];
			INTDBLMAP* vcov=NULL;
			if(gCovar.find(vId)==gCovar.end())
			{
				vcov=new INTDBLMAP;
				gCovar[vId]=vcov;
			}
			else
			{
				vcov=gCovar[vId];
			}
			for(EMAP_ITER uIter=vIter;uIter!=evidMap->end();uIter++)
			{
				int uId=uIter->first;
				//Don't compute covariance of vId uId pairs that both are not in the restrictedNeighborSet, when
				//the restrictedNeighborSet is empty
				//if((!random) && (vId!=uId) && (restrictedNeighborSet.size()>0))
				//{
				//	if((restrictedNeighborSet.find(vId)==restrictedNeighborSet.end()) && (restrictedNeighborSet.find(uId)==restrictedNeighborSet.end()))
				//	{
				//		continue;
				//	}
				//}
				Evidence* evid1=uIter->second;
				double uval=evid1->getEvidVal();
				double umean=gMean[uId];
				double diffprod=(vval-vmean)*(uval-umean);
				INTDBLMAP* ucov=NULL;
				if(gCovar.find(uId)==gCovar.end())
				{
					ucov=new INTDBLMAP;
					gCovar[uId]=ucov;
				}
				else
				{
					ucov=gCovar[uId];
				}
				if(vcov->find(uId)==vcov->end())
				{
					covPair++;
					(*vcov)[uId]=diffprod;
				}
				else
				{
					(*vcov)[uId]=(*vcov)[uId]+diffprod;
				}
				if(uId!=vId)
				{
					if(ucov->find(vId)==ucov->end())
					{
						(*ucov)[vId]=diffprod;
					}
					else
					{
						(*ucov)[vId]=(*ucov)[vId]+diffprod;
					}
				}
			}
		}

	}
	cout <<"Total covariance pairs estimated " << covPair << endl;
	//Now estimate the variance
	for(map<int,INTDBLMAP*>::iterator idIter=gCovar.begin();idIter!=gCovar.end();idIter++)
	{
		INTDBLMAP* var=idIter->second;
		for(INTDBLMAP_ITER vIter=var->begin();vIter!=var->end();vIter++)
		{
			if(vIter->first==idIter->first)
			{
				//vIter->second=2*vIter->second/((double)(gCovar.size()-1));
				//vIter->second=2*vIter->second/((double)(evidCnt-1));
				//vIter->second=(0.001+vIter->second)/((double)(evidCnt-1));
				vIter->second=(vIter->second)/((double)(evidCnt-1));
				double variance=vIter->second;
				if(idIter->first==176)
				{
				//	cout <<"Stop here: Variable " << idIter->first << " variance " << idIter->second << endl;
				}
			}
			else
			{
				vIter->second=vIter->second/((double)(evidCnt-1));
				//vIter->second=vIter->second/((double)(gCovar.size()-1));
				//vIter->second=0;
			}
			if(!random)
			{
				if(sdFile.good())
				{
					sdFile<<idIter->first<<"\t" << vIter->first <<"\t" << vIter->second << endl;
				}
			}
		}
	}
	if(!random)
	{
		if(mFile.good())
		{
			mFile.close();
		}
		if(sdFile.good())
		{
			sdFile.close();
		}
	}	
	return 0;
}
*/
/*
int
PotentialManager::estimateAllMeanCov(bool random, INTDBLMAP& gMean, map<int,INTDBLMAP*>& gCovar,INTINTMAP& trainEvidSet)
{
	int evidCnt=trainEvidSet.size();
	//First get the mean and then the variance
	int dId=0;
	for(INTINTMAP_ITER eIter=trainEvidSet.begin();eIter!=trainEvidSet.end();eIter++)
	{
		EMAP* evidMap=NULL;
		if(random)
		{
			evidMap=evMgr->getRandomEvidenceAt(eIter->first);
		}
		else
		{
			evidMap=evMgr->getEvidenceAt(eIter->first);
		}
		for(EMAP_ITER vIter=evidMap->begin();vIter!=evidMap->end(); vIter++)
		{
			int vId=vIter->first;
			Evidence* evid=vIter->second;
			double val=evid->getEvidVal();
			if(gMean.find(vId)==gMean.end())
			{
				gMean[vId]=val;
			}
			else
			{
				gMean[vId]=gMean[vId]+val;
			}
		}
		dId++;	
	}
	//Now estimate the mean
	for(INTDBLMAP_ITER idIter=gMean.begin();idIter!=gMean.end();idIter++)
	{
		idIter->second=idIter->second/(double) evidCnt;
		INTDBLMAP* vcov=new INTDBLMAP;
		gCovar[idIter->first]=vcov;
	}
	return 0;
}*/

int
PotentialManager::estimateAllMeanCov_Eff(bool random, vector<int>& trainEvidSet)
{
    cout << "PotentialManager::estimateAllMeanCov_Eff start" << endl;
    int evidCnt=trainEvidSet.size();
    //First get the mean and then the variance
    ////data is the data matrix which will have the variable by sample information
    EMAP* evidMap=evMgr->getEvidenceAt(trainEvidSet[0]);//trainEvidSet.begin()->first);
    int varCnt=evidMap->size();
    if(data==NULL)
    {
        data=new Matrix(varCnt,trainEvidSet.size());
        meanMat=new Matrix(varCnt,1);
        meanMat->setAllValues(0);
        covMat=new Matrix(varCnt,varCnt);
        covMat->setAllValues(-1);
    }
    Matrix *data=evMgr->getData();
    cout << " evidCnt=" << evidCnt << " varCnt=" << varCnt ;
    for(int eIter=0;eIter<trainEvidSet.size();eIter++) //for(INTINTMAP_ITER eIter=trainEvidSet.begin();eIter!=trainEvidSet.end();eIter++)
    {
        EMAP* evidMap=NULL;
        if(random)
        {
            evidMap=evMgr->getRandomEvidenceAt(trainEvidSet[eIter]);
        }
        else
        {
            evidMap=evMgr->getEvidenceAt(trainEvidSet[eIter]);
        }
        //for(EMAP_ITER vIter=evidMap->begin();vIter!=evidMap->end(); vIter++)
        for(int vId=0;vId<evidMap->size();vId++)
        {
            //int vId=vIter->first;
            //Evidence* evid=vIter->second;
            Evidence* evid=(*evidMap)[vId];
            double val=evid->getEvidVal();
            data->setValue(val,vId,trainEvidSet[eIter]);
        }
    }
    //Done copying. Now we can go over the rows of data and get the means
    for(int i=0;i<varCnt;i++)
    {
        double s=0;
        for(int j=0;j<data->getColCnt();j++)
        {
            double val=data->getValue(i,j);
            s=s+val;
        }
        double sampleSize=(double) data->getColCnt();
        meanMat->setValue(s/sampleSize,i,0);
        cout << "i=" << i << " mean=" <<s/sampleSize << endl;
    }
    cout <<"PotentialManager::estimateAllMeanCov_Eff end" << endl;
    return 0;
}

/*
// shilu: work for data with only one fold
int
PotentialManager::estimateAllMeanCov_Eff_onefold(bool random, int trainEvidSetSize)
{
    int evidCnt=trainEvidSetSize; //trainEvidSet.size();
    //First get the mean and then the variance
    ////data is the data matrix which will have the variable by sample information
    //EMAP* evidMap=evMgr->getEvidenceAt(0); //trainEvidSet[0]);//trainEvidSet.begin()->first);
    data=evMgr->getData();
    int varCnt=data->getRowCnt();//evidMap->size();
    //data=new Matrix(varCnt,trainEvidSetSize);
    meanMat=new Matrix(varCnt,1);
    meanMat->setAllValues(0);
    covMat=new Matrix(varCnt,varCnt);
    covMat->setAllValues(-1);
    cout << "PotentialManager::estimateAllMeanCov_Eff one fold: evidCnt(sample)=" << data->getColCnt() << " varibleCnt=" << data->getRowCnt() << endl;

    //Now we can go over the rows of data and get the means
    //double sampleSize=(double) data->getColCnt();
    for(int i=0;i<varCnt;i++)
    {
        double mean=data->RowMean(i);
        meanMat->setValue(mean,i,0);
        for(int j=i;j<varCnt;j++)
        {
            estimateCovariance_Eff(i,j);
        }
    }
    //cout <<"PotentialManager::estimateAllMeanCov_Eff end" << endl;
    return 0;
}*/

/*int
PotentialManager::estimateCovariance(bool random,INTDBLMAP* vcov, int vId, int uId)
{
	if((uId==4 && vId==71) || (uId==71 && uId==4))
	{
		//cout << "Stop here " << endl;
	}
	vector<int>& trainEvidSet=evMgr->getTestSet();
	int evidCnt=trainEvidSet.size();
	INTDBLMAP* ucov=globalCovar[uId];
	for(int eIter=0;eIter<trainEvidSet.size();eIter++) //for(INTINTMAP_ITER eIter=trainEvidSet.begin();eIter!=trainEvidSet.end();eIter++)
	{
		EMAP* evidMap=NULL;
		if(random)
		{
			evidMap=evMgr->getRandomEvidenceAt(trainEvidSet[eIter]);
		}
		else
		{
			evidMap=evMgr->getEvidenceAt(trainEvidSet[eIter]);
		}
		Evidence* evid=(*evidMap)[vId];
		double vval=evid->getEvidVal();
		double vmean=globalMean[vId];
		Evidence* evid1=(*evidMap)[uId];
		double uval=evid1->getEvidVal();
		double umean=globalMean[uId];
		double diffprod=(vval-vmean)*(uval-umean);
		if(vcov->find(uId)==vcov->end())
		{
			(*vcov)[uId]=diffprod;
		}
		else
		{
			(*vcov)[uId]=(*vcov)[uId]+diffprod;
		}
		if(uId!=vId)
		{
			if(ucov->find(vId)==ucov->end())
			{
				(*ucov)[vId]=diffprod;
			}
			else
			{
				(*ucov)[vId]=(*ucov)[vId]+diffprod;
			}
		}
	}
	//cout <<"Total covariance pairs estimated " << covPair << endl;
	//Now estimate the variance
	if(uId==vId)
	{
		double ssd=(*vcov)[uId];
		(*vcov)[uId]=(0.001+ssd)/((double)(evidCnt-1));
        	//cout << "estimateCovariance: vId=" << vId << " uId=" << uId << " val=" << (*vcov)[uId] << endl;
	}
	else
	{
		double ssduv=(*ucov)[vId];
		(*ucov)[vId]=ssduv/((double)(evidCnt-1));
		(*vcov)[uId]=ssduv/((double)(evidCnt-1));
        	//cout << "estimateCovariance: vId=" << vId << " uId=" << uId << " covariance=" << (*ucov)[vId] << endl;
		//cout << "estimateCovariance: uId=" << uId << " vId=" << vId << " covariance=" << (*ucov)[vId] << endl;
    	}
	return 0;
}*/

int
PotentialManager::estimateCovariance_Eff(int uId, int vId)
{
    
    double vmean=meanMat->getValue(vId,0);
    double umean=meanMat->getValue(uId,0);

    /*double ssd=0;
    for(int i=0;i<data->getColCnt();i++)
    {
        double vval=data->getValue(vId,i);
        double uval=data->getValue(uId,i);
        double diffprod=(vval-vmean)*(uval-umean);
        ssd=ssd+diffprod;
    }*/
    
    double ssd=data->vectorMultiply(vId,vmean,uId,umean);
    //Now estimate the variance
    double var=ssd/((double)(data->getColCnt()-1));
    /*if(uId==vId)
    {
        var=(0.001+ssd)/((double)(data->getColCnt()-1));
    }
    else
    {
        var=ssd/((double)(data->getColCnt()-1));
    }*/
    covMat->setValue(var,uId,vId);
    covMat->setValue(var,vId,uId);
    return 0;
}


int
PotentialManager::readAllMeanCov(const char* mFName, const char* sdFName)
{

	ifstream mFile(mFName);
	ifstream sdFile(sdFName);
	char buffer[1024];
	while(mFile.good())
	{
		mFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		int vId;
		double mean=0;
		while(tok!=NULL)
		{	
			if(tokCnt==0)
			{
				vId=atoi(tok);	
			}
			else if(tokCnt==1)
			{
				mean=atof(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		//globalMean[vId]=mean;
        meanMat->setValue(mean,vId,0);
	}
	mFile.close();
	int lineNo=0;
	while(sdFile.good())
	{
		sdFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		int vId=0;
		int uId=0;
		double covariance=0;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				uId=atoi(tok);
			}
			else if(tokCnt==1)
			{
				vId=atoi(tok);
			}
			else if(tokCnt==2)
			{
				covariance=atof(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		/*INTDBLMAP* ucov=NULL;
		if(globalCovar.find(uId)==globalCovar.end())
		{
			ucov=new INTDBLMAP;
			globalCovar[uId]=ucov;
		}
		else
		{
			ucov=globalCovar[uId];
		}
		(*ucov)[vId]=covariance;*/
        covMat->setValue(covariance,uId,vId);
        covMat->setValue(covariance,vId,uId);
		if(uId==0 && vId==0)
		{
			cout << "Found uId=0 vId=0 covar="<< covariance << " at lineno " << lineNo  << endl;
		}
		lineNo++;
	}
	sdFile.close();
	return 0;
}

/*
Error::ErrorCode
PotentialManager::populatePotentialsSlimFactors(map<int,SlimFactor*>& factorSet,VSET& varSet)
{
	//The set of flags to keep status of the potentials that have been calculated
	map<int,bool> doneFlag;
	for(map<int,SlimFactor*>::iterator fIter=factorSet.begin();fIter!=factorSet.end();fIter++)
	{
		doneFlag[fIter->first]=false;
	}
	int popFId=0;
	for(map<int,SlimFactor*>::reverse_iterator rIter=factorSet.rbegin();rIter!=factorSet.rend();rIter++)
	{
		//If we have computed the potential for this flag move one
		if(doneFlag[rIter->first])
		{
			popFId++;
			continue;
		}
		SlimFactor* sFactor=rIter->second;
		if(sFactor->fId==176)
		{
			cout <<"Stop here " << endl;
		}
		//Otherwise create the potential
		Potential* aPotFunc=new Potential;
		for(int j=0;j<sFactor->vCnt;j++)
		{
			Variable* aVar=varSet[sFactor->vIds[j]];
			if(j==sFactor->vCnt-1)
			{
				aPotFunc->setAssocVariable(aVar,Potential::FACTOR);
			}
			else
			{
				aPotFunc->setAssocVariable(aVar,Potential::MARKOV_BNKT);
			}
		}
		aPotFunc->potZeroInit();
		//populatePotential(aPotFunc,false);
        populatePotential_Eff(aPotFunc,false);
		aPotFunc->calculateJointEntropy();
		sFactor->jointEntropy=aPotFunc->getJointEntropy();
		if(sFactor->jointEntropy<0)
		{
		//	sFactor->jointEntropy=0;
		//	cout <<"Negative entropy for " << sFactor->fId << endl;
		}
		doneFlag[rIter->first]=true;
		delete aPotFunc;
		if(popFId%100000==0)
		{
			cout <<"Done with " << factorSet.size()-popFId << " factors " << endl;
		}
		popFId++;
	}
	return Error::SUCCESS;
}*/

/*
int 
PotentialManager::estimateMarginalEntropies(map<int,SlimFactor*>& slimFactors,VSET& varSet,bool random)
{	
	for(map<int,SlimFactor*>::iterator aIter=slimFactors.begin();aIter!=slimFactors.end();aIter++)
	{
		SlimFactor* sFactor=aIter->second;
		if(sFactor->vCnt>1)
		{
			break;
		}
		Potential* aPotFunc=new Potential;
		Variable* aVar=varSet[sFactor->vIds[0]];
		aPotFunc->setAssocVariable(aVar,Potential::FACTOR);
		aPotFunc->potZeroInit();
		//populatePotential(aPotFunc,random);
        populatePotential_Eff(aPotFunc,random);
		aPotFunc->calculateJointEntropy();
		sFactor->jointEntropy=aPotFunc->getJointEntropy();
		delete aPotFunc;
	}
	return 0;
}*/

/*
//Estimate the random information for all factors of a particular size. Assume that
//randInfo is already allocated
Error::ErrorCode
PotentialManager::estimateRandomInfo(map<int,SlimFactor*>& factorSet, 
		VSET& varSet, vector<double>& randInfo, int fSize)
{
	int rInd=0;
	for(map<int,SlimFactor*>::iterator rIter=factorSet.begin();rIter!=factorSet.end();rIter++)
	{
		SlimFactor* sFactor=rIter->second;
		if(sFactor->vCnt<fSize)
		{
			continue;
		}
		else if(sFactor->vCnt>fSize)
		{
			break;
		}
		//Otherwise create a potential
		Potential* aPotFunc=new Potential;
		for(int j=0;j<sFactor->vCnt;j++)
		{
			Variable* aVar=varSet[sFactor->vIds[j]];
			if(j==sFactor->vCnt-1)
			{
				aPotFunc->setAssocVariable(aVar,Potential::FACTOR);
			}
			else
			{
				aPotFunc->setAssocVariable(aVar,Potential::MARKOV_BNKT);
			}
		}
		aPotFunc->potZeroInit();
		//populatePotential(aPotFunc,true);
        populatePotential_Eff(aPotFunc,true);
		aPotFunc->calculateJointEntropy();
		double rInfo=(-1)*aPotFunc->getJointEntropy();
		for(int j=0;j<sFactor->vCnt;j++)
		{
			SlimFactor* subFactor=factorSet[sFactor->vIds[j]];
			rInfo=rInfo+subFactor->jointEntropy;
		}
		randInfo.push_back(rInfo);
		rInd++;
		delete aPotFunc;
	}
	return Error::SUCCESS;

}*/

/*
int 
PotentialManager::populateFactor(map<int,SlimFactor*>& factorSet,VSET& varSet,SlimFactor* sFactor,bool random)
{
	Potential* aPotFunc=new Potential;
	//string fullConfStr;
	char confStr[CONSTR_LEN];
	for(int j=0;j<sFactor->vCnt;j++)
	{
		Variable* aVar=varSet[sFactor->vIds[j]];
		if(j==sFactor->vCnt-1)
		{
			
			aPotFunc->setAssocVariable(aVar,Potential::FACTOR);
		}
		else
		{
			aPotFunc->setAssocVariable(aVar,Potential::MARKOV_BNKT);
		}
		//sprintf(confStr,"-%d",sFactor->vIds[j]);
		//fullConfStr.append(confStr);
	}
	aPotFunc->potZeroInit();
	//populatePotential(aPotFunc,random);
    populatePotential_Eff(aPotFunc,random);
	aPotFunc->calculateJointEntropy();
	//jointEntropies[fullConfStr]=aPotFunc->getJointEntropy();
	double rInfo=(-1)*aPotFunc->getJointEntropy();
	for(int j=0;j<sFactor->vCnt;j++)
	{
		SlimFactor* subFactor=factorSet[sFactor->vIds[j]];
		rInfo=rInfo+subFactor->jointEntropy;
	}
	sFactor->mutualInfo=rInfo;
	
	delete aPotFunc;
	return 0;
}


int 
PotentialManager::populateFactor_Buffer(map<int,SlimFactor*>& factorSet,VSET& varSet,SlimFactor* sFactor,bool random)
{
	Potential* aPotFunc=NULL;
	if(potBuffer.find(sFactor->vCnt)==potBuffer.end())
	{
		aPotFunc=new Potential;
		aPotFunc->initMatrix(sFactor->vCnt);
		potBuffer[sFactor->vCnt]=aPotFunc;
	}
	else
	{
		aPotFunc=potBuffer[sFactor->vCnt];
		aPotFunc->resetVarSet();
	}
	//string fullConfStr;
	char confStr[CONSTR_LEN];
	for(int j=0;j<sFactor->vCnt;j++)
	{
		Variable* aVar=varSet[sFactor->vIds[j]];
		if(j==sFactor->vCnt-1)
		{
			
			aPotFunc->setAssocVariable(aVar,Potential::FACTOR);
		}
		else
		{
			aPotFunc->setAssocVariable(aVar,Potential::MARKOV_BNKT);
		}
		//sprintf(confStr,"-%d",sFactor->vIds[j]);
		//fullConfStr.append(confStr);
	}
	aPotFunc->potZeroInit_MeanOnly();
	//populatePotential(aPotFunc,random);
    populatePotential_Eff(aPotFunc,random);
	aPotFunc->calculateJointEntropy();
	//jointEntropies[fullConfStr]=aPotFunc->getJointEntropy();
	double rInfo=(-1)*aPotFunc->getJointEntropy();
	for(int j=0;j<sFactor->vCnt;j++)
	{
		SlimFactor* subFactor=factorSet[sFactor->vIds[j]];
		rInfo=rInfo+subFactor->jointEntropy;
	}
	sFactor->mutualInfo=rInfo;
	
	//delete aPotFunc;
	return 0;
}*/

/*
//update mean and covariance
int
PotentialManager::populatePotential(Potential* aPot, bool random)
{
    cout << "PotentialManager::populatePotential" << endl;
	VSET& potVars=aPot->getAssocVariables();
	for(VSET_ITER vIter=potVars.begin();vIter!=potVars.end(); vIter++)
	{
		double mean=0;
		double cov=0;
		INTDBLMAP* covar=NULL;
		if(random)
		{
			if(globalMean_Rand.find(vIter->first)==globalMean_Rand.end())
			{
				cout <<"No var with id " << vIter->first << endl;
				exit(0);
			}
			mean=globalMean_Rand[vIter->first];
			covar=globalCovar_Rand[vIter->first];
		}
		else
		{
			if(globalMean.find(vIter->first)==globalMean.end())
			{
				cout <<"No var with id " << vIter->first << endl;
				exit(0);
			}
			mean=globalMean[vIter->first];
			covar=globalCovar[vIter->first];
		}
		aPot->updateMean(vIter->first,mean);
		for(VSET_ITER uIter=vIter;uIter!=potVars.end();uIter++)
		{
			if(covar->find(uIter->first)==covar->end())
			{
				estimateCovariance(random,covar,vIter->first,uIter->first);
				//cout <<"No var " << uIter->first << " in covariance of " << vIter->first << endl;
				//exit(0);
			}
			double cval=(*covar)[uIter->first];
			aPot->updateCovariance(vIter->first,uIter->first,cval);
			aPot->updateCovariance(uIter->first,vIter->first,cval);
		}

	}
	//aPot->makeValidJPD();
	aPot->makeValidJPD(ludecomp,perm);
	return 0;
}*/

//precompute covMat and meanMat:
// compute mbcondVar[CondVariance], mbcondMean_Part[CondBias], mbcondMean_Vect[CondWeight]
double
PotentialManager::computePotentialMBCovMean(SlimFactor* sFactor, int& status)//unordered_map<int,double>& mbcondMean_Vect)
{
    int mbsize=sFactor->mergedMB.size();
    if(mbsize==0)
    {
        return 0;
    }
    vector <int> ParChildID(sFactor->mergedMB.begin(),sFactor->mergedMB.end());
    ParChildID.push_back(sFactor->fId); //the last one is target
    //extract covariance
    Matrix* covariance=new Matrix(ParChildID.size(),ParChildID.size());
    for(int i=0;i<ParChildID.size();i++) {
        //double mean=meanMat->getValue(vIter->first,0);
        //aPot->updateMean(vIter->first,mean);
        int vIter=ParChildID[i];  //index in covMat
        for(int j=i;j<ParChildID.size();j++)
        {
            int uIter=ParChildID[j];  //index in covMat
            double cval=covMat->getValue(vIter,uIter);
            covariance->setValue(cval,i,j);
            covariance->setValue(cval,j,i);
        }
    }
    //covariance->showMatrix(cout);
    int vId=sFactor->fId; //factorVariables.begin()->first;
    //Matrix* covariance=covMat->getSubMatrix(ParChildID);
    double determinant=covariance->detMatrix();
    if(determinant<=0) //<1e-50)
    {
        //cout <<"Negative/Zero Determinant " << determinant << " for Target=" << sFactor->fId <<endl;
        status=-1;
        return 0;
    }
    double mbcondVar=covMat->getValue(vId,vId); //covariance->getValue(vIdmId,vIdmId);
    //double mbcondMean_Part=meanMat->getValue(vId,0); //meanWrap[vId];
    
    Matrix* mbcov=covariance->getSubMatrix(0,0,mbsize,mbsize); //Matrix* mbcov=new Matrix(mbsize,mbsize);
    Matrix* mbmargvar=covariance->getSubMatrix(ParChildID.size()-1,0,1,mbsize); //Matrix* mbmargvar=new Matrix(1,mbsize); covariance(vId,j)
    double determinantP=mbcov->detMatrix();
    if(determinantP<=0) //<1e-50
    {
        //cout <<"Negative/Zero Determinant Parent " << determinantP << " for Target=" << sFactor->fId <<endl;
        status=-1;
        return 0;
    }
    //cout <<"mbcov:" <<endl;
    //mbcov->showMatrix(cout);
    //cout <<"mbmargvar:" <<endl;
    //mbmargvar->showMatrix(cout);
    Matrix* covInv=mbcov->invMatrix();
    Matrix* prod1=mbmargvar->multiplyMatrix(covInv);
    //cout << "mbcondMean_Vect get my conditional weight: " ;
    //unordered_map<int,double> mbcondMean_Vect;
    for(int i=0;i<ParChildID.size()-1;i++)  //start from parent only, last one is target
    {
        int varID=ParChildID[i];
        double aVal=prod1->getValue(0,i);
        double bVal=mbmargvar->getValue(0,i);
        mbcondVar=mbcondVar-(aVal*bVal);
        //mbcondMean_Vect[varID]=aVal;  // mbcondMean_Vect is conditional weight
        //double cVal= meanMat->getValue(varID,0); //meanWrap[aIter->second];
        //mbcondMean_Part=mbcondMean_Part-(cVal*aVal);
        //cout << " i=" << i<< " varID=" << ParChildID[i] <<" conditionalWeight[" << varID << "]=" << aVal <<endl;
    }
    if(mbcondVar<0){
        status=-1;
        return 0;
    }
    double jointll1 = computeLL_Tracetrick(ParChildID.size(),determinant);
    double jointll2 = computeLL_Tracetrick(mbsize,determinantP);
    double pll = jointll1 - jointll2;
    //cout <<"Target=" << sFactor->fId << " dataSize=" << testdataSize << " determinant=" << determinant << " dim=" << ParChildID.size() << " determinantParent=" << determinantP <<  " dimP=" <<mbsize << " var(target)=" <<covMat->getValue(vId,vId) <<" jointll1=" << jointll1 <<" jointll2=" << jointll2<<endl;
    if(isinf(pll)){
        status=-1;
        //cout <<"Target=" << sFactor->fId <<" determinant=" << determinant << " dim=" << ParChildID.size() << " determinantParent=" << determinantP <<  " dimP=" <<mbsize << " var(target)=" <<covMat->getValue(vId,vId) <<" mbcondVar=" << mbcondVar <<endl;
        return 0;
    }
    
    ParChildID.clear();
    delete covariance;
    delete mbcov;
    delete mbmargvar;
    delete covInv;
    delete prod1;
    return pll;
}

int
PotentialManager::dumpVarMB_PairwiseFormat(SlimFactor* sFactor, ofstream& oFile, VSET& varSet)
{
    //unordered_map<int,double>& mbWts;
    int mbsize=sFactor->mergedMB.size();
    if(mbsize==0)
    {
        return 0;
    }
    vector <int> ParChildID(sFactor->mergedMB.begin(),sFactor->mergedMB.end());
    ParChildID.push_back(sFactor->fId); //the last one is target
    //extract covariance
    Matrix* covariance=new Matrix(ParChildID.size(),ParChildID.size());
    //cout << "PotentialManager::computePotentialMBCovMean my covariance for target ID=" << sFactor->fId << endl;
    for(int i=0;i<ParChildID.size();i++) {
        int vIter=ParChildID[i];  //index in covMat
        for(int j=i;j<ParChildID.size();j++)
        {
            int uIter=ParChildID[j];  //index in covMat
            double cval=covMat->getValue(vIter,uIter);
            covariance->setValue(cval,i,j);
            covariance->setValue(cval,j,i);
        }
    }
    //covariance->showMatrix(cout);
    Matrix* mbcov=covariance->getSubMatrix(0,0,mbsize,mbsize); //Matrix* mbcov=new Matrix(mbsize,mbsize);
    Matrix* mbmargvar=covariance->getSubMatrix(ParChildID.size()-1,0,1,mbsize); //Matrix* mbmargvar=new Matrix(1,mbsize); covariance(vId,j)

    //cout <<"mbcov:" <<endl;
    //mbcov->showMatrix(cout);
    //cout <<"mbmargvar:" <<endl;
    //mbmargvar->showMatrix(cout);
    Matrix* covInv=mbcov->invMatrix();
    Matrix* prod1=mbmargvar->multiplyMatrix(covInv);
    //cout << "mbcondMean_Vect get my conditional weight: " ;
    //unordered_map<int,double> mbcondMean_Vect;
    for(int i=0;i<ParChildID.size()-1;i++)  //start from parent only, last one id=ParChildID.size()-1 is target
    {
        int varID=ParChildID[i];
        double aVal=prod1->getValue(0,i); // aVal is conditional weight
        //mbWts[varID]=aVal;
        oFile << varSet[varID]->getName()<< "\t"<< varSet[sFactor->fId]->getName() << "\t" << aVal << endl;
        //cout << " i=" << i<< " varID=" << ParChildID[i] <<" conditionalWeight[" << varID << "]=" << aVal <<endl;
    }
    //cout <<"determinant=" << determinant << " dim=" << ParChildID.size() << " determinantParent=" << determinantP <<  " dimP=" <<mbsize <<" jointll1=" << jointll1 << " jointll2=" <<jointll2 <<  endl;
    ParChildID.clear();
    delete covariance;
    delete mbcov;
    delete mbmargvar;
    delete covInv;
    delete prod1;
    return 0;
}

// compute mbcondVar[CondVariance], mbcondMean_Part[CondBias], mbcondMean_Vect[CondWeight]
// for MetaLearner::showModelParameters() final parameters
int
PotentialManager::computePotentialMBCovMean(SlimFactor* sFactor, double& mbcondVar, double& mbbias, unordered_map<int,double>& mbcondMean_Vect)
{
    int mbsize=sFactor->mergedMB.size();
    if(mbsize==0)
    {
        return 0;
    }
    vector <int> ParChildID(sFactor->mergedMB.begin(),sFactor->mergedMB.end());
    ParChildID.push_back(sFactor->fId); //the last one is target
    //extract covariance
    Matrix* covariance=new Matrix(ParChildID.size(),ParChildID.size());
    //cout << "PotentialManager::computePotentialMBCovMean my covariance for target ID=" << sFactor->fId << endl;
    for(int i=0;i<ParChildID.size();i++) {
        //double mean=meanMat->getValue(vIter->first,0);
        //aPot->updateMean(vIter->first,mean);
        int vIter=ParChildID[i];  //index in covMat
        for(int j=i;j<ParChildID.size();j++)
        {
            int uIter=ParChildID[j];  //index in covMat
            double cval=covMat->getValue(vIter,uIter);
            covariance->setValue(cval,i,j);
            covariance->setValue(cval,j,i);
        }
    }
    //covariance->showMatrix(cout);
    int vId=sFactor->fId; //factorVariables.begin()->first;
    //Matrix* covariance=covMat->getSubMatrix(ParChildID);
    mbcondVar=covMat->getValue(vId,vId); //covariance->getValue(vIdmId,vIdmId);
    mbbias=meanMat->getValue(vId,0); //meanWrap[vId];
    
    Matrix* mbcov=covariance->getSubMatrix(0,0,mbsize,mbsize); //Matrix* mbcov=new Matrix(mbsize,mbsize);
    Matrix* mbmargvar=covariance->getSubMatrix(ParChildID.size()-1,0,1,mbsize); //Matrix* mbmargvar=new Matrix(1,mbsize); covariance(vId,j)
    Matrix* covInv=mbcov->invMatrix();
    Matrix* prod1=mbmargvar->multiplyMatrix(covInv);
    //cout << "mbcondMean_Vect get my conditional weight: " ;
    for(int i=0;i<ParChildID.size()-1;i++)  //start from parent only, last one is target
    {
        int varID=ParChildID[i];
        double aVal=prod1->getValue(0,i);
        double bVal=mbmargvar->getValue(0,i);
        mbcondVar=mbcondVar-(aVal*bVal);
        mbcondMean_Vect[varID]=aVal;  // mbcondMean_Vect is conditional weight
        double cVal= meanMat->getValue(varID,0); //meanWrap[aIter->second];
        mbbias=mbbias-(cVal*aVal);
        //cout << " i=" << i<< " varID=" << ParChildID[i] <<" conditionalWeight[" << varID << "]=" << aVal <<endl;
    }
    ParChildID.clear();
    delete covariance;
    delete mbcov;
    delete mbmargvar;
    delete covInv;
    delete prod1;
    return 0;
}

double
PotentialManager::computeLL_Tracetrick(int dim, double determinant)  //sampleSize=testdataSize;
{
    double ll=testdataSize * (dim*log(2*PI)+log(determinant));
    double t=dim*(testdataSize-1);
    ll=(ll+t)*(-0.5);
    return ll;
}


int
PotentialManager::populatePotential_Eff(Potential* aPot, bool random)
{
    unordered_map<int,Variable*>& potVars=aPot->getAssocVariables();
    for(auto vIter=potVars.begin();vIter!=potVars.end(); vIter++)
    {
        double mean=meanMat->getValue(vIter->first,0);
        aPot->updateMean(vIter->first,mean);
        for(auto uIter=vIter;uIter!=potVars.end();uIter++)
        {
            double cval=covMat->getValue(vIter->first,uIter->first);
            if(cval==-1)
            {
                estimateCovariance_Eff(vIter->first,uIter->first);
                //cerr <<"No var " << uIter->first << " in covariance of " << vIter->first << endl;
                //exit(-1);
            }
            cval=covMat->getValue(vIter->first,uIter->first);
            aPot->updateCovariance(vIter->first,uIter->first,cval);
            aPot->updateCovariance(uIter->first,vIter->first,cval);
        }
    }

    //aPot->makeValidJPD();
    aPot->makeValidJPD(ludecomp,perm); //same as aPot->makeValidJPD();
    return 0;
}

double
PotentialManager::computeMeanVarPseudoLikelihood_onefold(int id) //SlimFactor* sFactor,VSET& varSet)
{
    // compute meanMat and Variance on covMat, PseudoLikelihood
    //int id=aVar->getID();  //id=sFactor->fId=aVar->getID();
    double vmean=data->RowMean(id);
    meanMat->setValue(vmean,id,0);
    double dataSetSize=data->getColCnt(); //evMgr->getTestSetSize();
    double ssd=data->vectorMultiply(id,vmean,id,vmean);
    //double variance=covMat->getValue(id,id);
    //add 1e-10 to avoid singularity issues:
    double variance=ssd/(dataSetSize-1.0)+1e-10;  //(0.001+ssd)/((double)(data->getColCnt()-1))
    covMat->setValue(variance,id,id);
    double pll=-0.5*ssd/variance-0.5*log(2.0*PI*variance)*dataSetSize;
    //cout << "PotentialManager::computeMeanVarPseudoLikelihood_onefold id="<<id<<" variance=" <<variance << " mean=" << vmean << " ssd=" << ssd << " dataSetSize=" << dataSetSize <<" pll=" << pll << endl;
    for(int j=id+1;j<data->getRowCnt();j++)
    {
        estimateCovariance_Eff(id,j);
    }
    return pll;
}

double
PotentialManager::getPseudoLikelihood(SlimFactor* sFactor,vector<Variable*>& varSet, bool train)
{
	Potential* aPotFunc=new Potential;
	Variable* aVar=varSet[sFactor->fId];
	aPotFunc->setAssocVariable(aVar,Potential::FACTOR);
	for(auto aIter=sFactor->mergedMB.begin();aIter!=sFactor->mergedMB.end();aIter++)  //INTINTMAP_ITER
	{
		Variable* aVar=varSet[*aIter];
		aPotFunc->setAssocVariable(aVar,Potential::MARKOV_BNKT);
	}
	aPotFunc->potZeroInit();
	//populatePotential(aPotFunc,false);
    populatePotential_Eff(aPotFunc,false);
	//This function creates a submatrix of the covariance matrix and inverts it
	aPotFunc->initMBCovMean();
	vector<int>* dataSet=NULL; //INTINTMAP* dataSet=NULL;
	if(train)
	{
		dataSet=&(evMgr->getTrainingSet());
	}
	else
	{
		dataSet=&(evMgr->getTestSet());
	}
	INTDBLMAP subData;
	double pll=0;
	int thresholded=0;
    for(int dIter=0;dIter<dataSet->size();dIter++) //for(INTINTMAP_ITER dIter=dataSet->begin();dIter!=dataSet->end();dIter++)
	{
		EMAP* evidMap=NULL;
		evidMap=evMgr->getEvidenceAt((*dataSet)[dIter]); //dIter->first
		Evidence* evid=(*evidMap)[sFactor->fId];
		double val=evid->getEvidVal();
		subData[sFactor->fId]=val;
		for(auto vIter=sFactor->mergedMB.begin();vIter!=sFactor->mergedMB.end(); vIter++) //INTINTMAP_ITER
		{
			int vId=*vIter;
			Evidence* evid=(*evidMap)[vId];
			double val=evid->getEvidVal();
			subData[vId]=val;
		}
		double cll=aPotFunc->getCondPotValueFor(subData);

		if(cll<1e-50)
		{
			cll=1e-50;
			thresholded++;
		}
		pll=pll+log(cll);
	}
	subData.clear();
	if(thresholded>0)
	{
	//	cout <<"Thresholded " << thresholded << " datapoints to 1e-50" << endl;
	}
	delete aPotFunc;
	return pll;
}

/*
double
PotentialManager::getPseudoLikelihood(SlimFactor* sFactor,VSET& varSet, bool train,int& status)
{
	status=0;
	Potential* aPotFunc=new Potential;
	Variable* aVar=varSet[sFactor->fId];
	aPotFunc->setAssocVariable(aVar,Potential::FACTOR);
	for(auto aIter=sFactor->mergedMB.begin();aIter!=sFactor->mergedMB.end();aIter++) //INTINTMAP_ITER
	{
		Variable* aVar=varSet[*aIter];
		aPotFunc->setAssocVariable(aVar,Potential::MARKOV_BNKT);
	}
	aPotFunc->potZeroInit();
	//populatePotential(aPotFunc,false);
    populatePotential_Eff(aPotFunc,false);
	//This function creates a submatrix of the covariance matrix and inverts it
	aPotFunc->initMBCovMean();
	if(aPotFunc->getCondVariance()<0)
	{
		status=-1;
		delete aPotFunc;
		return 0;
	}
	vector<int>* dataSet=NULL;
	if(train)
	{
		dataSet=&(evMgr->getTrainingSet());
	}
	else
	{
		dataSet=&(evMgr->getTestSet());
	}
	INTDBLMAP subData;
	double pll=0;
	int thresholded=0;
    for(int dIter=0;dIter<dataSet->size();dIter++) //for(INTINTMAP_ITER dIter=dataSet->begin();dIter!=dataSet->end();dIter++)
	{
		EMAP* evidMap=NULL;
		evidMap=evMgr->getEvidenceAt((*dataSet)[dIter]);
		Evidence* evid=(*evidMap)[sFactor->fId];
		double val=evid->getEvidVal();
		subData[sFactor->fId]=val;
		for(auto vIter=sFactor->mergedMB.begin();vIter!=sFactor->mergedMB.end(); vIter++) //INTINTMAP_ITER
		{
			int vId=*vIter;
			Evidence* evid=(*evidMap)[vId];
			double val=evid->getEvidVal();
			subData[vId]=val;
		}
		double cll=aPotFunc->getCondPotValueFor(subData);

		if(cll<1e-50)
		{
			cll=1e-50;
			thresholded++;
		}
		pll=pll+log(cll);
	}
	subData.clear();
	if(thresholded>0)
	{
	//	cout <<"Thresholded " << thresholded << " datapoints to 1e-50" << endl;
	}
	delete aPotFunc;
	return pll;
}


double
PotentialManager::getPseudoLikelihood(SlimFactor* sFactor,VSET& varSet, bool train,int& status, Potential* aPotFunc)
{
	status=0;
	Variable* aVar=varSet[sFactor->fId];
	aPotFunc->setAssocVariable(aVar,Potential::FACTOR);
	for(auto aIter=sFactor->mergedMB.begin();aIter!=sFactor->mergedMB.end();aIter++) //INTINTMAP_ITER
	{
		Variable* aVar=varSet[*aIter];
		aPotFunc->setAssocVariable(aVar,Potential::MARKOV_BNKT);
	}
	aPotFunc->potZeroInit();
	//populatePotential(aPotFunc,false);
    populatePotential_Eff(aPotFunc,false);
	//This function creates a submatrix of the covariance matrix and inverts it
	aPotFunc->initMBCovMean();
	if(aPotFunc->getCondVariance()<0)
	{
		status=-1;
		return 0;
	}
	vector<int>* dataSet=NULL;
	if(train)
	{
		dataSet=&(evMgr->getTrainingSet());
	}
	else
	{
		dataSet=&(evMgr->getTestSet());
	}
	INTDBLMAP subData;
	double pll=0;
	int thresholded=0;
	for(int dIter=0;dIter<dataSet->size();dIter++) //for(INTINTMAP_ITER dIter=dataSet->begin();dIter!=dataSet->end();dIter++)
	{
		EMAP* evidMap=NULL;
		evidMap=evMgr->getEvidenceAt((*dataSet)[dIter]);
		Evidence* evid=(*evidMap)[sFactor->fId];
		double val=evid->getEvidVal();
		subData[sFactor->fId]=val;
		for(auto vIter=sFactor->mergedMB.begin();vIter!=sFactor->mergedMB.end(); vIter++) //INTINTMAP_ITER
		{
			int vId=*vIter;
			Evidence* evid=(*evidMap)[vId];
			double val=evid->getEvidVal();
			subData[vId]=val;
		}
		double cll=aPotFunc->getCondPotValueFor(subData);

		if(cll<1e-50)
		{
			cll=1e-50;
			thresholded++;
		}
		pll=pll+log(cll);
	}
	subData.clear();
	return pll;
}
*/
/*
double 
PotentialManager::getGaussianLikelihood(map<int,SlimFactor*>& factorSet,VSET& varSet, bool train)
{
    cout << "PotentialManager::getGaussianLikelihood" << endl;
	Potential* aPotFunc=new Potential;
	for(map<int,SlimFactor*>::iterator fIter=factorSet.begin();fIter!=factorSet.end();fIter++)
	{
		Variable* aVar=varSet[fIter->first];
		if(fIter==factorSet.begin())
		{
			aPotFunc->setAssocVariable(aVar,Potential::FACTOR);
		}
		else
		{
			aPotFunc->setAssocVariable(aVar,Potential::MARKOV_BNKT);
		}
	}
	aPotFunc->potZeroInit();
	//Use the graph structure to update the particular elements of the covariance and mean of this potential
	for(map<int,SlimFactor*>::iterator fIter=factorSet.begin();fIter!=factorSet.end();fIter++)
	{
		SlimFactor* sFactor=fIter->second;

		double mean=0;
		double cov=0;
		INTDBLMAP* covar=NULL;
		if(globalMean.find(fIter->first)==globalMean.end())
		{
			cout <<"No var with id " << fIter->first << endl;
			exit(0);
		}
		mean=globalMean[fIter->first];
		covar=globalCovar[fIter->first];
		aPotFunc->updateMean(fIter->first,mean);
		double vval=(*covar)[fIter->first];
		aPotFunc->updateCovariance(fIter->first,fIter->first,vval);
		for(auto mIter=sFactor->mergedMB.begin();mIter!=sFactor->mergedMB.end();mIter++) //INTINTMAP_ITER
		{
			double cval=(*covar)[*mIter];
			aPotFunc->updateCovariance(fIter->first,*mIter,cval);
			aPotFunc->updateCovariance(*mIter,fIter->first,cval);
		}
	}

	aPotFunc->makeValidJPD();
	INTDBLMAP datapt;
	double gll=0;
	int thresholded=0;
	vector<int>* dataSet;
	if(train)
	{
		dataSet=(&evMgr->getTrainingSet());
	}
	else
	{
		dataSet=(&evMgr->getTestSet());
	}
	for(int dIter=0;dIter<dataSet->size();dIter++) //for(INTINTMAP_ITER dIter=dataSet->begin();dIter!=dataSet->end();dIter++)
	{
		EMAP* evidMap=NULL;
		evidMap=evMgr->getEvidenceAt((*dataSet)[dIter]);
		for(map<int,SlimFactor*>::iterator fIter=factorSet.begin();fIter!=factorSet.end();fIter++)
		{
			int vId=fIter->first;
			Evidence* evid=(*evidMap)[vId];
			double val=evid->getEvidVal();
			datapt[vId]=val;
		}
		double jll=aPotFunc->getJointPotValueFor(datapt);

		if(jll<1e-50)
		{
			jll=1e-50;
			thresholded++;
		}
		gll=gll+log(jll);
	}
	delete aPotFunc;
	return gll;
}

double
PotentialManager::getLikelihood(SlimFactor* sFactor,VSET& varSet)
{
	cout <<"Not implemented" << endl;
	return 0;
}


double
PotentialManager::getLikelihood(SlimFactor* sFactor,VSET& varSet,map<int,int>& visitedVertices )
{
	double dll=0;
	cout <<"Not implemented" << endl;
	return dll;
}


int
PotentialManager::estimateConditionalPotential(SlimFactor* sFactor,VSET& varSet,Potential** pot, STRDBLMAP& counts)
{
	cout <<"Not implemented" << endl;
	return 0;
}

int
PotentialManager::populatePotential(Potential* pot,STRDBLMAP& counts)
{
	cout <<"Not implemented" << endl;
	return 0;
}

int 
PotentialManager::estimateCanonicalPotential(SlimFactor* sFactor, VSET& variableSet,INTINTMAP& defInst,INTINTMAP& factorSubsets,map<int,SlimFactor*>& canonicalFactorSet)
{
	cout <<"Not implemented" << endl;
	return 0;
}


int 
PotentialManager::estimateCanonicalPotential_Abbeel(SlimFactor* sFactor, VSET& variableSet,INTINTMAP& defInst,INTINTMAP& factorSubsets,map<int,SlimFactor*>& canonicalFactorSet)
{
	cout <<"Not implemented" << endl;
	return 0;
}



int 
PotentialManager::estimateCanonicalPotential_Approximate(SlimFactor* sFactor, VSET& variableSet,INTINTMAP& defInst,INTINTMAP& factorSubsets,map<int,SlimFactor*>& canonicalFactorSet)
{
	cout <<"Not implemented" << endl;
	return 0;
}*/

/*
int
PotentialManager::resetPotFuncs()
{
	for(map<int,Potential*>::iterator pIter=potFuncs.begin();pIter!=potFuncs.end();pIter++)
	{
		delete pIter->second;
	}
	potFuncs.clear();
	return 0;
}



int 
PotentialManager::estimateCanonicalPotential_Joint(SlimFactor* sFactor, VSET& variableSet,INTINTMAP& defInst,INTINTMAP& factorSubsets,map<int,SlimFactor*>& canonicalFactorSet)
{
	cout <<"Not implemented" << endl;
	return 0;
}

Potential*
PotentialManager::getPotential(int fId)
{
	if(potFuncs.find(fId)==potFuncs.end())
	{
		return NULL;
	}
	return potFuncs[fId];
}


double 
PotentialManager::getConditionalEntropy(int vId,INTINTMAP& fVars,VSET& varSet)
{
	double condEntropy=0;
	//string fullConfStr;
	//string partConfStr;
	char confStr[CONSTR_LEN];

	double fullJointEntropy=0;

		Potential* potFunc=new Potential;
		for(INTINTMAP_ITER aIter=fVars.begin();aIter!=fVars.end();aIter++)
		{
			if(aIter==fVars.begin())
			{
				potFunc->setAssocVariable(varSet[aIter->first],Potential::FACTOR);
			}
			else
			{
				potFunc->setAssocVariable(varSet[aIter->first],Potential::MARKOV_BNKT);
			}
		}
		potFunc->potZeroInit();
		//populatePotential(potFunc,false);
        populatePotential_Eff(potFunc,false);
		potFunc->calculateJointEntropy();
		fullJointEntropy=potFunc->getJointEntropy();
		//jointEntropies[fullConfStr]=fullJointEntropy;
		delete potFunc;
	//}
	if(fVars.size()==1)
	{

		return fullJointEntropy;
	}
	double partJointEntropy=0;

		potFunc=new Potential;
		bool setFactorVar=false;
		for(INTINTMAP_ITER aIter=fVars.begin();aIter!=fVars.end();aIter++)
		{
			if(aIter->first==vId)
			{
				continue;
			}
			if(!setFactorVar)
			{
				setFactorVar=true;
				potFunc->setAssocVariable(varSet[aIter->first],Potential::FACTOR);
			}
			else
			{
				potFunc->setAssocVariable(varSet[aIter->first],Potential::MARKOV_BNKT);
			}
		}
		STRDBLMAP counts;
		potFunc->potZeroInit();
		//populatePotential(potFunc,false);
        populatePotential_Eff(potFunc,false);
		potFunc->calculateJointEntropy();
		partJointEntropy=potFunc->getJointEntropy();
		//jointEntropies[partConfStr]=partJointEntropy;
		delete potFunc;
	//}
	condEntropy=fullJointEntropy-partJointEntropy;
	//fullConfStr.clear();
	//partConfStr.clear();

	return condEntropy;
}


double 
PotentialManager::getSampleLikelihood(map<int,SlimFactor*>& factorSet, VSET& varSet, INTINTMAP* sample)
{
	double sampleLL=0;
	cout <<"Not implemented " <<endl;
	return sampleLL;
}


int 
PotentialManager::getVariableSample(INTINTMAP& jointConf,VSET& varSet,int vId,SlimFactor* sFactor, gsl_rng* r)
{
	Potential* pot=NULL;
	if(potFuncs.find(vId)==potFuncs.end())
	{
		pot=new Potential;
		STRDBLMAP counts;
		estimateConditionalPotential(sFactor,varSet,&pot,counts);
		potFuncs[vId]=pot;
	}
	else
	{
		pot=potFuncs[vId];
	}
	//int sample=pot->generateSample(jointConf,vId,r);
	int sample=-1;
	return sample;
}

int
PotentialManager::clearJointEntropies()
{
	jointEntropies.clear();
	return 0;
}

double
PotentialManager::estimateCanonicalValue(INTINTMAP& reqInst,INTINTMAP& defInst,INTINTMAP& allSubsets,map<int,SlimFactor*>& canFactors,Potential* condPot)
{
	double pVal=0;
	cout <<"Not implemented " << endl;
	return 0;
}


double
PotentialManager::estimateCanonicalValue_Joint(INTINTMAP& reqInst,INTINTMAP& defInst,INTINTMAP& allSubsets,map<int,SlimFactor*>& canFactors,Potential* jointPot)
{
	cout <<"Not implemented " << endl;
	return 0;
}*/
