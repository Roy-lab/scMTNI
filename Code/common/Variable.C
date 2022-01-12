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

