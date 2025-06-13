#include <iostream>
#include <fstream>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "Framework.H"
int sortfunc_less(const void* first, const void* second);
int sortfunc_more(const void* first, const void* second);
int* sortingind=NULL;
double* sortedvals=NULL;

Framework::Framework()
{
}

Framework::~Framework()
{
}

int 
Framework::readNetwork(const char* aFName)
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
		if(strstr(buffer,"YAR009C")!=NULL)
		{	
			cout <<"Found TF "<<  buffer<< endl;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		string tf;
		string tgt;
		double regwt=0;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				tf.append(tok);
			}
			else if(tokCnt==1)
			{
				tgt.append(tok);	
			}
			else if(tokCnt==2)
			{
				regwt=atof(tok);
				regwt=fabs(regwt);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		map<string,double>* regprog=NULL;
		if(nw.find(tf)==nw.end())
		{
			regprog=new map<string,double>;
			nw[tf]=regprog;
		}
		else 
		{
			regprog=nw[tf];
		}
		(*regprog)[tgt]=regwt;
	}
	cout <<"Read "<< nw.size() << " TFs" << endl;
	inFile.close();
	return 0;
}

int 
Framework::convertToPercentile(const char* aFName, bool less)
{
	ofstream oFile(aFName);
	for(map<string,map<string,double>*>::iterator gIter=nw.begin();gIter!=nw.end();gIter++)
	{
		vector<string> sortedEdge;
		sort(gIter->second,sortedEdge,less);
		map<string,double>* edgeWts=gIter->second;
		double currScore=(*edgeWts)[sortedEdge[0]];
		map<string,double> edgeEqual;
		for(int i=0;i<sortedEdge.size();i++)
		{
			double newScore=(*edgeWts)[sortedEdge[i]];
			bool update=false;
			if(less)
			{
				if(newScore<currScore)
				{
					update=true;
				}
			}
			else
			{
				if(newScore>currScore)
				{
					update=true;
				}
			}
			if(update)
			{
				double lessScores=(double)(sortedEdge.size()-(i));
				double percentilerank=(double)edgeEqual.size();
				if(edgeEqual.size()>1)
				{
					percentilerank=0.5*percentilerank;
				}
				percentilerank=(lessScores+percentilerank)/((double)sortedEdge.size());
				for(map<string,double>::iterator sIter=edgeEqual.begin();sIter!=edgeEqual.end();sIter++)
				{
					oFile <<gIter->first <<"\t"<< sIter->first <<"\t" << percentilerank << endl;
				}
				edgeEqual.clear();
				currScore=newScore;
				edgeEqual[sortedEdge[i]]=0;
			}
			else
			{	
				edgeEqual[sortedEdge[i]]=0;
			}
		}	
		double percentilerank=(double)edgeEqual.size();
		if(edgeEqual.size()>1)
		{
			percentilerank=0.5*percentilerank;
		}
		percentilerank=percentilerank/((double)sortedEdge.size());
		for(map<string,double>::iterator sIter=edgeEqual.begin();sIter!=edgeEqual.end();sIter++)
		{
			oFile <<gIter->first <<"\t"<< sIter->first <<"\t" << percentilerank << endl;
		}
	}
	oFile.close();
	return 0;
}

int
Framework::sort(map<string,double>* toSort, vector<string>& sorted, bool less )
{
	sortedvals=new double[toSort->size()];
	sortingind=new int[toSort->size()];
	int t=0;
	map<int,string> edgeIndex;
	for(map<string,double>::iterator eIter=toSort->begin();eIter!=toSort->end();eIter++)
	{
		sortedvals[t]=eIter->second;
		sortingind[t]=t;
		edgeIndex[t]=eIter->first;
		t++;
	}
	if(less)
	{
		qsort(sortingind,toSort->size(),sizeof(int),&sortfunc_less);
	}
	else
	{
		qsort(sortingind,toSort->size(),sizeof(int),&sortfunc_more);
	}
	for(int e=0;e<toSort->size();e++)
	{
		int eid=sortingind[e];
		sorted.push_back(edgeIndex[eid]);
	}
	delete[] sortingind;	
	delete[] sortedvals;
	return 0;
}

int 
sortfunc_less(const void* first, const void* second)
{
	int ind1=*((int*)first);	
	int ind2=*((int*)second);
	double val1=sortedvals[ind1];
	double val2=sortedvals[ind2];
	int compstat=0;
	if(val1>val2)
	{
		compstat=-1;
	}
	else if(val1<val2)
	{
		compstat=1;
	}
	return compstat;
}

int 
sortfunc_more(const void* first, const void* second)
{
	int ind1=*((int*)first);	
	int ind2=*((int*)second);
	double val1=sortedvals[ind1];
	double val2=sortedvals[ind2];
	int compstat=0;
	if(val1<val2)
	{
		compstat=-1;
	}
	else if(val1>val2)
	{
		compstat=1;
	}
	return compstat;
}


int
main(int argc, const char** argv)
{
	if(argc!=4)
	{
		cout <<"Usage: convertToDream in.txt out.txt scoreorder[incr|decr]" << endl;
		return 0;
	}
	Framework fw;
	fw.readNetwork(argv[1]);
	if(strcmp(argv[3],"incr")==0)
	{
		fw.convertToPercentile(argv[2],false);
	}
	else
	{
		fw.convertToPercentile(argv[2],true);
	}
	return 0;
}
