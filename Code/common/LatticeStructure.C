#include <iostream>
#include "LatticeStructure.H"


LatticeStructure::LatticeStructure()
{
}

LatticeStructure::~LatticeStructure()
{
	for(MATRIX_ITER mIter=mSubset.begin();mIter!=mSubset.end();mIter++)
	{
		mIter->second->clear();
		delete mIter->second;
	}
	mSubset.clear();
	for(MATRIX_ITER mIter=mSuperset.begin();mIter!=mSuperset.end();mIter++)
	{
		mIter->second->clear();
		delete mIter->second;
	}
	mSuperset.clear();
}

//Add supersetId as an entry to the row corresponding to subsetId
int 
LatticeStructure::addSubset(int subsetId,int supersetId)
{
	INTINTMAP* superSets=NULL;
	if(mSubset.find(subsetId)==mSubset.end())
	{
		superSets=new INTINTMAP;
		mSubset[subsetId]=superSets;
	}
	else
	{
		superSets=mSubset[subsetId];
	}
	(*superSets)[supersetId]=0;
	//factor with ID supersetId will be the superset for all 
	//subsets of the factor with ID subsetId
	/*if(mSuperset.find(subsetId)!=mSuperset.end())
	{
		INTINTMAP* subSets=mSuperset[subsetId];
		for(INTINTMAP_ITER iIter=subSets->begin();iIter!=subSets->end();iIter++)
		{
			int subsubSetId=iIter->first;
			addSubset(subsubSetId,supersetId);
		}
	}*/
	return 0;
}

//Add subsetId as an entry to the row corresponding to supersetId
int 
LatticeStructure::addSuperset(int supersetId,int subsetId)
{
	INTINTMAP* subSets=NULL;
	if(mSuperset.find(supersetId)==mSuperset.end())
	{
		subSets=new INTINTMAP;
		mSuperset[supersetId]=subSets;
	}
	else
	{
		subSets=mSuperset[supersetId];
	}
	(*subSets)[subsetId]=0;
	//factor with ID supersetId will be the superset for all 
	//subsets of the factor with ID subsetId
	/*if(mSuperset.find(subsetId)!=mSuperset.end())
	{
		INTINTMAP* subSets=mSuperset[subsetId];
		for(INTINTMAP_ITER iIter=subSets->begin();iIter!=subSets->end();iIter++)
		{
			int subsubsetId=iIter->first;
			addSuperset(supersetId,subsubsetId);
		}
	}*/
	return 0;
}

INTINTMAP* 
LatticeStructure::getSubsets(int fId)
{
	if(mSuperset.find(fId)==mSuperset.end())
	{
		return NULL;
	}
	return mSuperset[fId];
}

INTINTMAP* 
LatticeStructure::getSupersets(int fId)
{
	if(mSubset.find(fId)==mSubset.end())
	{
		return NULL;
	}
	return mSubset[fId];
}

int
LatticeStructure::getSupersets(int fId, INTINTMAP& supersets, int level)
{
	if(level==0)
	{
		INTINTMAP* sSet=getSupersets(fId);
		if(sSet==NULL)
		{
			return 0;
		}
		for(INTINTMAP_ITER fIter=sSet->begin();fIter!=sSet->end();fIter++)
		{
			supersets[fIter->first]=0;
		}
	}
	else
	{
		int currlevel=1;
		getSuperSets(fId,supersets,currlevel,level);
	}
	return 0;
}

int
LatticeStructure::getSuperSets(int fId,INTINTMAP& supersets,int currlevel,int level)
{
	if(currlevel>level)
	{
		return 0;
	}
	INTINTMAP* immeSuperset=getSupersets(fId);
	if(immeSuperset==NULL)
	{
		return 0;
	}
	for(INTINTMAP_ITER ssIter=immeSuperset->begin();ssIter!=immeSuperset->end();ssIter++)
	{
		if(currlevel==level)
		{
			supersets[ssIter->first]=0;
		}
		else
		{
			getSuperSets(ssIter->first,supersets,currlevel+1,level);
		}
	}
	return 0;
}

int
LatticeStructure::getAllSubsets(int fId,INTINTMAP& ssIds)
{
	INTINTMAP* immeSubset=getSubsets(fId);
	if(immeSubset==NULL)
	{
		return 0;
	}
	for(INTINTMAP_ITER ssIter=immeSubset->begin();ssIter!=immeSubset->end();ssIter++)
	{
		ssIds[ssIter->first]=0;
		getAllSubsets(ssIter->first,ssIds);
	}
	return 0;
}

int
LatticeStructure::deleteFromLattice(int fId)
{
	//We need to get all factors that are supersets of this factor and delete them from there rows in the mSuperset
	//Finally we want to delete this factor from mSubset.
	MATRIX_ITER sIter=mSubset.find(fId);
	if(sIter!=mSubset.end())
	{
		INTINTMAP* supersets=sIter->second;
		/*for(INTINTMAP_ITER tIter=supersets->begin();tIter!=supersets->end();tIter++)
		{
			if(mSuperset.find(tIter->first)==mSuperset.end())
			{
				cout <<"No such factor id in mSuperset matrix" << endl;
				return -1;
			}
			INTINTMAP* aSuperSet=mSuperset[tIter->first];
			//Now erase the factor corresponding to fId from here
			INTINTMAP_ITER eIter=aSuperSet->find(fId);
			if(eIter==aSuperSet->end())
			{
				cout <<"Inconsistent lattice structure " << endl;
				return -1;
			}
			aSuperSet->erase(eIter);
		}*/
		sIter->second->clear();
		delete sIter->second;
		mSubset.erase(sIter);
	}

	//We need to get all factors that are subsets of this factor, delete them from corresponding rows of mSubset
	sIter=mSuperset.find(fId);
	if(sIter!=mSuperset.end())
	{
		INTINTMAP* subsets=sIter->second;
		for(INTINTMAP_ITER tIter=subsets->begin();tIter!=subsets->end();tIter++)
		{
			if(mSubset.find(tIter->first)==mSubset.end())
			{
				cout <<"No such factor id in mSubset matrix" << endl;
				return -1;
			}
			INTINTMAP* aSubSet=mSubset[tIter->first];
			//Now erase the factor corresponding to fId from here
			INTINTMAP_ITER eIter=aSubSet->find(fId);
			if(eIter==aSubSet->end())
			{
				cout <<"Inconsistent lattice structure " << endl;
				return -1;
			}
			aSubSet->erase(eIter);
		}
		sIter->second->clear();
		delete sIter->second;
		mSuperset.erase(sIter);
	}
	
	return 0;
}

int
LatticeStructure::clear()
{
	for(MATRIX_ITER sIter=mSubset.begin();sIter!=mSubset.end();sIter++)
	{
		sIter->second->clear();
		delete sIter->second;
	}
	mSubset.clear();
	for(MATRIX_ITER mIter=mSuperset.begin();mIter!=mSuperset.end();mIter++)
	{
		mIter->second->clear();
		delete mIter->second;
	}
	mSuperset.clear();
	return 0;
}
