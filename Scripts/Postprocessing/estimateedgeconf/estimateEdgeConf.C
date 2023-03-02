#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <string.h>
#include <stdlib.h>

using namespace std;

map<string,double> edgeSet;

int readNetwork(const char*);
int readNetworkSet(const char*);
int estimateConfidence();
int filterHighConfidenceEdges(double,const char*);
int showConfidenceAllEdges(const char*);
int networkCnt=0;

int
main(int argc, const char** argv)
{
	if(argc!=5)
	{
		cout <<"Usage: estEdgeConf filenamelist confidence outputfile [filterededges|alledges]" << endl;
		return 0;
	}
	readNetworkSet(argv[1]);
	estimateConfidence();
	if(strcmp(argv[4],"filterededges")==0)
	{
		filterHighConfidenceEdges(atof(argv[2]),argv[3]);
	}
	else
	{
		showConfidenceAllEdges(argv[3]);
	}
	return 0;
}

int readNetworkSet(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		readNetwork(buffer);
		networkCnt=networkCnt+1;
	}
	inFile.close();
	return 0;
}

int
readNetwork(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[1024];
	map<string,int> localNetwork;
	while(inFile.good())
	{
		inFile.getline(buffer,1024);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		string gene1;
		string gene2;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				gene1.append(tok);
			}
			else if(tokCnt==1)
			{
				gene2.append(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		string key;
	//	if(gene1.compare(gene2)<=0)
	//	{
			key.append(gene1);
			key.append("\t");
			key.append(gene2);
	//	}
	//	else
	//	{
	//		key.append(gene2);
	//		key.append("\t");
	//		key.append(gene1);
	//	}
		localNetwork[key]=0;
	}
	for(map<string,int>::iterator eIter=localNetwork.begin();eIter!=localNetwork.end();eIter++)
	{
		if(edgeSet.find(eIter->first)==edgeSet.end())
		{
			edgeSet[eIter->first]=1;
		}
		else
		{
			edgeSet[eIter->first]=edgeSet[eIter->first]+1;
		}
	}
	localNetwork.clear();
	inFile.close();
	return 0;
}

int
estimateConfidence()
{
	for(map<string,double>::iterator eIter=edgeSet.begin();eIter!=edgeSet.end();eIter++)
	{
		eIter->second=eIter->second/((double)networkCnt);
	}
	return 0;
}

int
filterHighConfidenceEdges(double confThresh, const char* aFName)
{
	ofstream oFile(aFName);
	for(map<string,double>::iterator eIter=edgeSet.begin();eIter!=edgeSet.end();eIter++)
	{
		if(eIter->second<confThresh)
		{
			continue;
		}
		oFile << eIter->first << endl;
	}
	oFile.close();
	return 0;
}

int 
showConfidenceAllEdges(const char* aFSuff)
{
	char aFName[1024];
	sprintf(aFName,"%sdist.txt",aFSuff);
	ofstream oFile(aFName);
	sprintf(aFName,"%salledge.txt",aFSuff);
	ofstream eFile(aFName);
	
	map<int,double> edgeConfDist;
	double minConf=1/(double)networkCnt;
	double maxConf=1.0;
	double binSize=(maxConf-minConf)/20.0;
	for(map<string,double>::iterator eIter=edgeSet.begin();eIter!=edgeSet.end();eIter++)
	{
		double econf=eIter->second;
		int bin=(int)((econf-minConf)/binSize);
		if(bin>20)
		{
			cout<<" stop here bin=" <<  bin << " confidence " << econf << " for edge " << eIter->first << endl;
		}
		if(edgeConfDist.find(bin)==edgeConfDist.end())
		{
			edgeConfDist[bin]=1;
		}
		else
		{
			edgeConfDist[bin]=edgeConfDist[bin]+1;
		}
		eFile << eIter->first<< "\t"<< econf << endl;
	}
	double cdf=0;
	for(map<int,double>::iterator aIter=edgeConfDist.begin();aIter!=edgeConfDist.end();aIter++)
	{
		double pdf=aIter->second/((double)edgeSet.size());
		cdf=cdf+pdf;
		oFile << minConf+(binSize*aIter->first) <<"\t"<< pdf<< "\t" << cdf<< endl;
	}
	oFile.close();
	eFile.close();
	return 0;
}
