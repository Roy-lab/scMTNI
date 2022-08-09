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
