#include "Error.H"

const char* Error::errorString[]={"Success", "Variable Schema", "Datafile", "Variable", 
				"SampleCnt", "Factor Creation",  "MarkovBlanket", 
				"CorruptEvidence", "Null Potential", "Missing Var Assignment",
				"Missing Configuration", "Unknown"};

Error::Error()
{
}

Error::~Error()
{
}

const char*
Error::getErrorString(int eCode)
{
	if(eCode >= Error::UNKNOWN)
	{
		return errorString[Error::UNKNOWN];
	}
	return errorString[eCode];
}
