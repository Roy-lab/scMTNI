#include <cstring>
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>
#include "Error.H"
#include "Utils.H"
#include "Variable.H"
#include "VariableManager.H"


VariableManager::VariableManager()
{
}

VariableManager::~VariableManager()
{
}

Error::ErrorCode
VariableManager::readVariablesFromTable(vector<string>& inputTable)
{
	int nodeCount = inputTable.size();
	int nodeID = 0;
	for (auto line : inputTable)
	{
		vector<string> substrs = Utils::split(line, '\t');
		string nodeName = substrs[0];
		Variable* var = new Variable;
		var->setID(nodeID);
		var->setName(nodeName.c_str());
        variableSet.push_back(var);//variableSet[nodeID]=var;
		varNameIDMap[nodeName] = nodeID;
		++nodeID;
        
	}
	cout << "Read information about " << nodeCount << " variables" << endl;
	return Error::SUCCESS;
}

Error::ErrorCode
VariableManager::readVariablesFromData(string nodeName,int nodeID)
{
    Variable* var = new Variable;
    var->setID(nodeID);
    var->setName(nodeName.c_str());
    variableSet.push_back(var);
    varNameIDMap[nodeName] = nodeID;
    return Error::SUCCESS;
}

/*
//Reads the schema of the variables

Error::ErrorCode
VariableManager::readVariables(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[40000];
	int nodeCnt=0;
	while(inFile.good())
	{
		inFile.getline(buffer,40000);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		if(strstr(buffer,"NodeCnt")!=NULL)
		{
			char* pos=strchr(buffer,'\t');
			if(pos==NULL)
			{
				cout <<"Bad file format" << endl;
				return Error::VARSCHEMA_ERR;
			}
			nodeCnt=atoi(pos+1);
		}
		else if((strstr(buffer,"ContinuousNodes")!=NULL)||
			(strstr(buffer,"DiscreteNodes")!=NULL))	
		{
			char* tok=strtok(buffer,"\t");
			int tokCnt=0;

			while(tok!=NULL)
			{
				if(tokCnt>0)
				{
					int nid=atoi(tok);
					Variable* var=new Variable;
					var->setID(nid);
					variableSet[nid]=var;
				}
				tokCnt++;
				tok=strtok(NULL,"\t");
			}
		}
		else if(strstr(buffer,"NodeName")!=NULL)
		{
			char* tok=strtok(buffer,"\t");
			int tokCnt=0;
			int currNodeId;
			Variable* var=NULL;
			char varName[256];
			while(tok!=NULL)
			{
				//Node name
				if(tokCnt==0)
				{
					char* pos=strchr(tok,'=');
					if(pos==NULL)
					{
						cout << "Bad format in nodename " << endl;
						return Error::VARSCHEMA_ERR;
					}
					strcpy(varName,pos+1);
				}
				//Node id
				else if(tokCnt==1)
				{
					char* pos=strchr(tok,'=');
					if(pos==NULL)
					{
						cout <<"Bad format" << endl;
						return Error::VARSCHEMA_ERR;
					}
					currNodeId=atoi(pos+1);
					var=variableSet[currNodeId];
					var->setName(varName);
					string varKey(varName);
					varNameIDMap[varKey]=currNodeId;
				}
				//Parent or child ids fields are 2 and 3, which we must ignore for structure learning
				//value fields is tokCnt=4
				/*else if(tokCnt==4)
				{
					char* pos=strchr(tok,'=');
					if(pos==NULL)
					{
						cout <<"Bad format" << endl;
						return Error::VARSCHEMA_ERR;
					}
					int index=1;
					char tempId[5];
					int tempIndex=0;
					while(*(pos+index)!='\0')
					{
						if(*(pos+index)==',')
						{
							tempId[tempIndex]='\0';
							int valField=atoi(tempId);
							var->setNextValue(valField);
							tempIndex=0;
						}
						else
						{
							tempId[tempIndex]=*(pos+index);
							tempIndex++;
						}
						index++;
					}
					if(tempIndex>0)
					{
						tempId[tempIndex]='\0';
						int value=atoi(tempId);
						var->setNextValue(value);
					}
				}*//*
				tok=strtok(NULL,"\t");
				tokCnt++;
			}
		}
	}
	inFile.close();
	cout <<"Read information about " << variableSet.size() << " variables " << endl;
	return Error::SUCCESS;
}*/

int
VariableManager::getVarID(const char* varName)
{
	string varKey(varName);
        /* comment by shilu
	if(varNameIDMap.find(varKey)==varNameIDMap.end())
	{
		return -1;
	}*/
	int vId=varNameIDMap[varKey];
	return vId;
}

bool 
VariableManager::isValid(int varID,int varVal)
{
	Variable* rVar=variableSet[varID];
	return rVar->isValidValue(varVal);
}

vector<Variable*>&
VariableManager::getVariableSet()
{
	return variableSet;
}


Variable* 
VariableManager::getVariableAt(int vId)
{
	/*if(variableSet.find(vId)==variableSet.end())
	{
		cout << "Illegal variable id " << vId << endl;
		return NULL;
	}*/
	return variableSet[vId];
}
