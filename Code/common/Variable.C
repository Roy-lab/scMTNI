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

#include "Variable.H"

Variable::Variable()
{
}

Variable::~Variable()
{
}

int 
Variable::setName(const char* aStr)
{
	name.append(aStr);
	return 0;
}

const string& 
Variable::getName()
{
	return name;
}

int
Variable::setID(int aId)
{
	vId=aId;
	return 0;
}

int
Variable::getID()
{
	return vId;
}

int 
Variable::setNextValue(int aVal)
{
	valueSet.push_back(aVal);
	return 0;
}

int 
Variable::setValues(INTVECT& aValSet)
{
	for(int i=0;i<aValSet.size();i++)
	{
		valueSet.push_back(aValSet[i]);
	}
	return 0;
}

INTVECT& 
Variable::getValues()
{
	return valueSet;
}

int 
Variable::getValueCnt()
{
	return valueSet.size();
}


bool 
Variable::isValidValue(int aVal)
{
	bool found=false;
	int i=0;
	while((i<valueSet.size()) && (!found))
	{
		if(aVal==valueSet[i])
		{
			found=true;
		}
		i++;
	}
	return found;
}


int
Variable::initEvidence(INTDBLMAP& evidData)
{
	for(int i=0;i<valueSet.size();i++)
	{
		evidData[valueSet[i]]=0.0;
	}
	return 0;
}

