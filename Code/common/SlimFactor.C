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

#include <math.h>
#include "Variable.H"
#include "Potential.H"
#include "SlimFactor.H"

SlimFactor::SlimFactor()
{
	mbScore=-1;
	//mutualInfo=0;
	//refCnt=0;
}

SlimFactor::SlimFactor(int fSize)
{
	//vIds=new int[fSize];
	//vCnt=fSize;
	//secondPId=-1;
	//mutualInfo=0;
	//jointEntropy=0;
	//confidence=0;
	fId=-1;
	mbScore=-1;
	//refCnt=0;
}

SlimFactor::~SlimFactor()
{
	//delete[] vIds;
	//goodMBIDs.clear();
	mergedMB.clear();
	//tabulist.clear();
	//candidateNeighbours.clear();
	//candidateNeighbours_vect.clear();
	//mbSubsetStartInd.clear();
	mbWts.clear();
}

/*
//This returns the maximal subsets of this factor. That is
//all subsets with vCnt-1 factors
int 
SlimFactor::generateMaximalSubsets(int** subsets)
{
	for(int i=0;i<vCnt;i++)
	{
		int varId=0;
		for(int j=0;j<vCnt;j++)
		{
			if(i!=j)
			{
				subsets[i][varId]=vIds[j];
				varId++;
			}
		}
	}

	return 0;
}

//Right now I am doing a linear search to find the variables that are not in this factor
//Need to update this later
int
SlimFactor::getSetDiff(SlimFactor* aFactor, int* diffVIds, int& diffSize)
{
	diffSize=0;
	for(int i=0;i<aFactor->vCnt;i++)
	{
		if(!isMemberVariable(aFactor->vIds[i]))
		{
			diffVIds[diffSize]=aFactor->vIds[i];
			diffSize++;
		}
	}
	return 0;
}

bool
SlimFactor::isMemberVariable(int v)
{
	//Do a linear search of v in vIds. This is ok since our clusters are small when this function is called
	int startIndex=0;
	int endIndex=vCnt-1;
	if((vIds[startIndex]> v) || (vIds[endIndex]<v))
	{
		return false;
	}
	bool hit=false;
	for(int i=0;i<=endIndex;i++)
	{
		if(vIds[i]==v)
		{
			hit=true;
		}
	}
	return hit;
}*/

/*
int 
SlimFactor::genMBSubsets(int ssSize)
{
	if(mbSubsets.size()==0)
	{
		mbSubsetStartInd[1]=0;
		//Initialize
		for(auto vIter=mergedMB.begin();vIter!=mergedMB.end();vIter++)
		{
			INTINTMAP* sset=new INTINTMAP;
			mbSubsets.push_back(sset);
			(*sset)[vIter->first]=0;
		}
		mbSubsetStartInd[ssSize]=mbSubsets.size();
	}
	//At each call we grow our subsets by one element
	int startInd=mbSubsetStartInd[ssSize-1];
	int endInd=mbSubsetStartInd[ssSize];
	for(int i=startInd;i<endInd;i++)
	{
		INTINTMAP* oldset=mbSubsets[i];
		int lastVid=oldset->rbegin()->first;
		auto vIter=mergedMB.begin();
		while((vIter->first<=lastVid) && (vIter!=mergedMB.end()))
		{
			vIter++;
		}
		//Use the remaining variables from vIter to create new subsets
		for(;vIter!=mergedMB.end();vIter++)
		{
			INTINTMAP* newset=new INTINTMAP;
			for(INTINTMAP_ITER uIter=oldset->begin();uIter!=oldset->end();uIter++)
			{
				(*newset)[uIter->first]=0;
			}
			(*newset)[vIter->first]=0;
			mbSubsets.push_back(newset);
		}
	}
	mbSubsetStartInd[ssSize+1]=mbSubsets.size();
	return 0;
}*/

/*
bool
SlimFactor::allEntriesInsignificant()
{
	bool insignificant=true;
	STRDBLMAP& canonicalVals=canonicalParams->getJointPotTable();
	STRDBLMAP_ITER aIter=canonicalVals.begin();
	while((aIter!=canonicalVals.end()) && (insignificant))
	{
		if(aIter->second!=1.0)
		{
			insignificant=false;
		}
		aIter++;
	}
	return insignificant;
}*/

/*
int
SlimFactor::thresholdToOne(double threshold)
{
	int onedEntries=0;
	STRDBLMAP& canonicalVals=canonicalParams->getJointPotTable();
	for(STRDBLMAP_ITER aIter=canonicalVals.begin();aIter!=canonicalVals.end();aIter++)
	{
		if(fabs(log(aIter->second))<=threshold)
		{
			aIter->second=1.0;
			onedEntries++;
		}
	}
	if(onedEntries>0)
	{
		//cout <<"Oned "<< onedEntries << " of a total of " << canonicalVals.size() << endl;
	}
	return 0;
}*/

/*
int 
SlimFactor::showFactor(ostream& oFile, VSET& variableSet, bool newLine)
{
	for(int i=0;i<vCnt;i++)
	{
		if(i>0)
		{
			oFile << "-";
		}
		oFile << variableSet[vIds[i]]->getName();
	}
	if(newLine)
	{
		oFile << endl;
	}
	return 0;
}*/

int
SlimFactor::setMBWts(unordered_map<int,double>& wts)
{
	for(auto wIter=wts.begin();wIter!=wts.end();wIter++)
	{
		mbWts[wIter->first]=wIter->second;
	}
	return 0;
}
