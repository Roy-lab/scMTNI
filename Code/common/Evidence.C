#include "Evidence.H"

Evidence::Evidence()
{
	hiddenFlag=false;
	hardValue=-1;
}

Evidence::~Evidence()
{
}

int
Evidence::assocVariable(int id)
{
	varID=id;
	return 0;
}

int
Evidence::getAssocVariable()
{
	return varID;
}

//Set the likelihood vector for this variable
int 
Evidence::setData(INTDBLMAP& data)
{
	for(INTDBLMAP_ITER idIter=data.begin();idIter!=data.end();idIter++)
	{
		evidData[idIter->first]=idIter->second;
		if(eType==Evidence::HARD)
		{
			if(idIter->second>0)
			{
				hardValue=idIter->first;
			}
		}
	}
	return 0;
}

int
Evidence::addToData(int valId, double val)
{
	evidData[valId]=val;
	return 0;
}

INTDBLMAP& 
Evidence::getData()
{
	return evidData;
}

int 
Evidence::setType(Evidence::EvidenceType aType)
{
	eType=aType;
	return 0;
}

Evidence::EvidenceType
Evidence::getType()
{
	return eType;
}

int
Evidence::dumpEvidence(ostream& oFile)
{
	oFile << varID <<"=[";
	for(INTDBLMAP_ITER idIter=evidData.begin();idIter!=evidData.end();idIter++)
	{
		if(idIter!=evidData.begin())
		{
			oFile<<",";
		}
		oFile << idIter->first<<"|"<<idIter->second;
		
	}
	oFile<<"]";
	return 0;
}


int 
Evidence::makeHidden()
{
	hiddenFlag=true;
	return 0;
}

bool 
Evidence::isHidden()
{
	return hiddenFlag;
}

int
Evidence::getMLVal()
{
	int maxVal=-1;
	double maxProb=0.0;
	for(INTDBLMAP_ITER aIter=evidData.begin();aIter!=evidData.end();aIter++)
	{
		if(aIter->second>maxProb)
		{
			maxProb=aIter->second;
			maxVal=aIter->first;
		}
	}
	return maxVal;
}

int
Evidence::getHardEvidVal()
{
	return hardValue;
}


int 
Evidence::setEvidVal(double aval)
{
	contValue=aval;
	return 0;
}

double 
Evidence::getEvidVal()
{
	return contValue;
}
