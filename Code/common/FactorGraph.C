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

#include <cstring>
#include "Variable.H"

#include "Error.H"
#include "VariableManager.H"
#include "SlimFactor.H"
#include "Error.H"
#include "LatticeStructure.H"
//#include "FactorManager.H"
#include "FactorGraph.H"

FactorGraph::FactorGraph()
{
}

FactorGraph::~FactorGraph()
{
	for(auto fIter=factorSet.begin();fIter!=factorSet.end();fIter++)
	{
		delete *fIter;
	}
	factorSet.clear();
}


int 
FactorGraph::setFactor(SlimFactor* sFactor)
{
    //index same as sFactor->fId, start from 0
    factorSet.push_back(sFactor);// factorSet[sFactor->fId]=sFactor;
	return 0;
}

int 
FactorGraph::getFactorCnt()
{
	return factorSet.size();
}

SlimFactor* 
FactorGraph::getFactorAt(int fid)
{
	/*if(factorSet.find(fid)==factorSet.end())
	{
		return NULL;
	}*/
	return factorSet[fid];
}


vector<SlimFactor*>&
FactorGraph::getAllFactors()
{
	return factorSet;
}
/*
int 
FactorGraph::dumpFactors_ClusterFormat(const char* outputDir, int maxk, VSET& variableSet)
{
	char aFName[1024];
	sprintf(aFName,"%s/maximal_clusters_cw_k%d.txt",outputDir,maxk);
	ofstream oFile(aFName);
	for(map<int,SlimFactor*>::iterator aIter=factorSet.begin();aIter!=factorSet.end();aIter++)
	{
		if(aIter->second->vCnt==1)
		{
			continue;
		}
		aIter->second->showFactor(oFile,variableSet);
	}
	oFile.close();
	return 0;
}


int 
FactorGraph::dumpFactors_PairwiseFormat(const char* outputDir,int maxk,VSET& variableSet)
{
	char aFName[1024];
	sprintf(aFName,"%s/maximal_clusters_pw_k%d.txt",outputDir,maxk);
	ofstream oFile(aFName);
	for(map<int,SlimFactor*>::iterator aIter=factorSet.begin();aIter!=factorSet.end();aIter++)
	{
		if(showFactor_PW(aIter->second,oFile,variableSet)==-1)
		{
			oFile.close();
			return -1;;
		}
	}
	oFile.close();
	return 0;
}*/

/*
int 
FactorGraph::dumpVarMB_ClusterFormat(const char* outputDir,int maxk,VSET& variableSet,FactorManager* fMgr)
{
	char aFName[1024];
	sprintf(aFName,"%s/var_mb_cw_k%d.txt",outputDir,maxk);
	ofstream oFile(aFName);
	for(map<int,SlimFactor*>::iterator aIter=factorSet.begin();aIter!=factorSet.end();aIter++)
	{
		SlimFactor* sFactor=aIter->second;
		//Now show the Markov blankets of this factor
		for(INTDBLMAP_ITER idIter=sFactor->goodMBIDs.begin();idIter!=sFactor->goodMBIDs.end();idIter++)
		{
			//SlimFactor* mbFactor=fMgr->getFactorAt(idIter->first);
			SlimFactor* mbFactor=factorSet[aIter->first];
			oFile << variableSet[sFactor->vIds[0]]->getName();
			for(int i=0;i<mbFactor->vCnt;i++)
			{
				if(mbFactor->vIds[i]==sFactor->vIds[0])
				{
					continue;
				}
				oFile << "\t" << variableSet[mbFactor->vIds[i]]->getName();
			}
			oFile << endl;
		}
		oFile << endl;
	}
	
	oFile.close();
	return 0;
}*/

/*
int 
FactorGraph::dumpVarMB_PairwiseFormat(const char* outputDir,int maxk,VSET& variableSet,FactorManager* fMgr)
{
	char aFName[1024];
	sprintf(aFName,"%s/var_mb_pw_k%d.txt",outputDir,maxk);
	ofstream oFile(aFName);
	for(map<int,SlimFactor*>::iterator aIter=factorSet.begin();aIter!=factorSet.end();aIter++)
	{
		SlimFactor* sFactor=aIter->second;
		if(sFactor->vCnt>1)
		{
			break;
		}
		//Now show the Markov blankets of this factor
		/*for(INTDBLMAP_ITER idIter=sFactor->goodMBIDs.begin();idIter!=sFactor->goodMBIDs.end();idIter++)
		{
			SlimFactor* mbFactor=fMgr->getFactorAt(idIter->first);
			for(int i=0;i<mbFactor->vCnt;i++)
			{
				if(mbFactor->vIds[i]==sFactor->vIds[0])
				{
					continue;
				}
				oFile << variableSet[sFactor->vIds[0]]->getName()<< "\t" 
				<< variableSet[mbFactor->vIds[i]]->getName() << endl;
			}
		}*//*
		for(auto mIter=sFactor->mergedMB.begin();mIter!=sFactor->mergedMB.end();mIter++) //INTINTMAP_ITER
		{
			
			oFile << variableSet[sFactor->vIds[0]]->getName()<< "\t" 
			<< variableSet[*mIter]->getName() << endl;
		}
	}
	oFile.close();
	return 0;
}*/


int 
FactorGraph::dumpVarMB_PairwiseFormat(ofstream& oFile,VSET& variableSet)
{
	for(auto aIter=factorSet.begin();aIter!=factorSet.end();aIter++)
	{
		SlimFactor* sFactor=*aIter;//aIter->second
		/*if(sFactor->vCnt>1)
		{
			break;
		}*/
		unordered_map<int,double>& mbWts=sFactor->mbWts;
		for(auto mIter=sFactor->mergedMB.begin();mIter!=sFactor->mergedMB.end();mIter++)//INTINTMAP_ITER
		{
			/*if(mbWts.find(*mIter)==mbWts.end())
			{
				cout <<"No neighbor " << mIter->first << " for " << sFactor->fId << endl;
				exit(0);
			}*/
			double wtval=mbWts[*mIter];
            //sFactor->vIds[0] same as sFactor->fId
			oFile << variableSet[*mIter]->getName()<< "\t"
			<< variableSet[sFactor->fId]->getName() << "\t" << wtval
			<< endl;
		}
	}
	return 0;
}

/*
int
FactorGraph::readVarMB_PairwiseFormat(const char* aFName,VariableManager* vMgr,FactorManager* fMgr)
{
	ifstream inFile(aFName);
	if(!inFile.good())
	{
		cout <<"Bad file name " << aFName << endl;
		return -1;
	}
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		int vId=-1;
		int mbvId=-1;
		int tokCnt=0;
		char* tok=strtok(buffer,"\t");
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				vId=vMgr->getVarID(tok);
			}
			else if(tokCnt==1)
			{
				mbvId=vMgr->getVarID(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		SlimFactor* sFactor=factorSet[vId];
        sFactor->mergedMB.insert(mbvId); //sFactor->mergedMB[mbvId]=0;
	}
	inFile.close();
	return 0;
}*/

/*
int 
FactorGraph::dumpCandidateVarMB_PairwiseFormat(const char* outputDir,int maxk,VSET& variableSet,FactorManager* fMgr)
{
	char aFName[1024];
	sprintf(aFName,"%s/candidatevar_mb_pw_k%d.txt",outputDir,maxk);
	ofstream oFile(aFName);
	for(map<int,SlimFactor*>::iterator aIter=factorSet.begin();aIter!=factorSet.end();aIter++)
	{
		SlimFactor* sFactor=aIter->second;
		if(sFactor->vCnt>1)
		{
			break;
		}
		//Now show the Markov blankets of this factor
		for(INTDBLMAP_ITER idIter=sFactor->candidateNeighbours.begin();idIter!=sFactor->candidateNeighbours.end();idIter++)
		{
			oFile << variableSet[sFactor->vIds[0]]->getName()<< "\t" 
			<< variableSet[idIter->first]->getName() << "\t" << idIter->second << "\t" << sFactor->marginalEntropy-idIter->second << endl;
		}
	}
	oFile.close();
	return 0;
}*/

/*
//Here we just check if the mutual consistency check holds
bool 
FactorGraph::isConsistent()
{
	bool consistent=true;
	for(map<int,SlimFactor*>::iterator fIter=factorSet.begin();fIter!=factorSet.end();fIter++)
	{
		SlimFactor* sFactor=fIter->second;
		for(auto vIter=sFactor->mergedMB.begin();vIter!=sFactor->mergedMB.end();vIter++) //INTINTMAP_ITER
		{
			SlimFactor* mFactor=factorSet[*vIter];
			if(mFactor->mergedMB.find(fIter->first)==mFactor->mergedMB.end())
			{
				consistent=false;
				cout <<"Inconsistency. factor: "<< sFactor->fId << " mbvar: " << *vIter <<endl;
			}
		}
		if(!consistent)
		{
			break;
		}
	}
	return consistent;
}

int 
FactorGraph::showFactor_PW(SlimFactor* sFactor,ostream& oFile, VSET& variableSet)
{
	map<string,int> shownPairs;
	//First if this factor is of size 2, then there might be a pairwise interaction
	//we should capture
	for(int i=0;i<sFactor->vCnt;i++)
	{
		int uid=sFactor->vIds[i];
		for(int j=i+1;j<sFactor->vCnt;j++)
		{
			int vid=sFactor->vIds[j];
			string fKey;
			char key[256];
			if(uid<vid)
			{
				sprintf(key,"-%d-%d",uid,vid);
			}
			else 
			{
				sprintf(key,"-%d-%d",vid,uid);
			}
			fKey.append(key);
			if(shownPairs.find(fKey)==shownPairs.end())
			{
				oFile << variableSet[sFactor->vIds[0]]->getName() <<"\tunknown\t"
				<< variableSet[sFactor->vIds[1]]->getName() << endl;
				shownPairs[fKey]=0;
			}
		}
	}
	//Now use every variable in the factor with every variable in the Markov blanket and create a pair
	for(int i=0;i<sFactor->vCnt;i++)
	{
		int vid=sFactor->vIds[i];
		for(auto vIter=sFactor->mergedMB.begin();vIter!=sFactor->mergedMB.end();vIter++) //INTINTMAP_ITER
		{
			string fKey;
			char key[256];
			int uid=*vIter;
			if(uid<vid)
			{
				sprintf(key,"-%d-%d",uid,vid);
			}
			else 
			{
				sprintf(key,"-%d-%d",vid,uid);
			}
			fKey.append(key);
			if(shownPairs.find(fKey)!=shownPairs.end())
			{
				continue;
			}
			oFile << variableSet[uid]->getName() <<"\tunknown\t"
			<< variableSet[vid]->getName() << endl;
			shownPairs[fKey]=0;
		}
	}
	return 0;
}*/
