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

#include <fstream>
#include <iostream>
#include <cstring>
#include <queue>
#include <math.h>

#include <sys/timeb.h>
#include <sys/time.h>
#include <time.h>

#include "Vertex.H"
#include "Graph.H"


Graph::Graph()
{
}

Graph::~Graph()
{
}

//The format of the file is with an edge information per line
int 
Graph::makeGraph(const char* fName)
{
	ifstream inFile(fName);
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)==0)
		{
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		Vertex* v=NULL;
		Vertex* u=NULL;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				string vKey(tok);
				if(vertexList.find(vKey)==vertexList.end())
				{
					v=new Vertex;
					v->setName(tok);
					vertexList[vKey]=v;
				}
				else
				{
					v=vertexList[vKey];
				}
			}
			else if(tokCnt==1)
			{
				string uKey(tok);
				if(vertexList.find(uKey)==vertexList.end())
				{
					u=new Vertex;
					u->setName(tok);
					vertexList[uKey]=u;
				}
				else
				{
					u=vertexList[uKey];
				}
			}
			tokCnt++;
			tok=strtok(NULL,"\t");
		}
		if(v!=NULL && u!=NULL)
		{
			if(v->getOutDegree()<2 && u->getOutDegree()<2)
			{
				v->setOutNeighbour(u);
				u->setOutNeighbour(v);
			}
		}
	}
	inFile.close();
	cout <<"Number of nodes: "<< vertexList.size() <<   endl;
	return 0;
}

//The format of the file is with an edge information per line
int 
Graph::makeGraph(const char* fName,int nCnt)
{
	ifstream inFile(fName);
	char buffer[1024];
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)==0)
		{
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		Vertex* v=NULL;
		Vertex* u=NULL;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				string vKey(tok);
				if(vertexList.find(vKey)==vertexList.end())
				{
					v=new Vertex;
					v->setName(tok);
					vertexList[vKey]=v;
				}
				else
				{
					v=vertexList[vKey];
				}
			}
			else if(tokCnt==1)
			{
				string uKey(tok);
				if(vertexList.find(uKey)==vertexList.end())
				{
					u=new Vertex;
					u->setName(tok);
					vertexList[uKey]=u;
				}
				else
				{
					u=vertexList[uKey];
				}
			}
			tokCnt++;
			tok=strtok(NULL,"\t");
		}
		if(v!=NULL && u!=NULL)
		{
			if(v->getOutDegree()<nCnt && u->getOutDegree()<nCnt)
			{
				v->setOutNeighbour(u);
				u->setOutNeighbour(v);
			}
		}
	}
	inFile.close();
	cout <<"Number of nodes: "<< vertexList.size() <<   endl;
	return 0;
}

int
Graph::getVertexCnt()
{
	return vertexList.size();
}

int
Graph::dumpGraph()
{
	ofstream oFile("graph.sif");
	map<string,int> doneVertices;
	for(map<string,Vertex*>::iterator aIter=vertexList.begin();aIter!=vertexList.end();aIter++)
	{
		NINFO_MAP& neighbrs=aIter->second->getImmediateNeighbours();
		for(NINFO_MAP_ITER nmIter=neighbrs.begin();nmIter!=neighbrs.end();nmIter++)
		{
			if(doneVertices.find(nmIter->first)==doneVertices.end())
			{
				oFile << aIter->first.c_str() <<" interacts " << nmIter->first.c_str() << endl;
			}
		}
		doneVertices[aIter->first]=0;
	}
	oFile.close();
	return 0;
}

//This is a bad way of doing stuff. 
int 
Graph::obtainConnectivity()
{
	map<string,Vertex*>::iterator aIter;
	/*struct timeb begintime;
	struct timeb endtime;*/
	struct timeval begintime;
	struct timeval endtime;
	gettimeofday(&begintime,NULL);
	for(aIter=vertexList.begin();aIter!=vertexList.end();aIter++)
	{
		Vertex* v=aIter->second;
		v->findReachableNodes();
	}
	gettimeofday(&endtime,NULL);

	cout << "Time elapsed" << endtime.tv_sec-begintime.tv_sec<< " seconds and " << endtime.tv_usec-begintime.tv_usec << " micro secs" << endl;

	return 0;
}

int 
Graph::showConnectivity()
{
	for(map<string,Vertex*>::iterator aIter=vertexList.begin(); aIter!=vertexList.end();aIter++)
	{
		Vertex* v=aIter->second;
		v->showReachability(cout);
	}
	return 0;
}

//For a graph to be connected every vertex must be reachable to every other vertex
bool 
Graph::isConnected()
{
	bool isReachable=true;
	map<string,Vertex*>::iterator aIter=vertexList.begin();
	while(aIter!=vertexList.end() && isReachable)
	{
		Vertex* u=aIter->second;
		map<string,Vertex*>::iterator bIter=vertexList.begin();
		while(bIter!=vertexList.end() && isReachable)
		{
			if(aIter!=bIter)
			{
				Vertex* v=bIter->second;
				if(!u->isReachable(v))
				{
					isReachable=false;
				}
			}
			bIter++;
		}
		aIter++;
	}
	return isReachable;
}

int 
Graph::getDegreeDist()
{
	//true for indegree
	getDegreeDist(true,indegreeDist);
	getDegreeDist(false,outdegreeDist);
	return 0;
}

int 
Graph::showDegreeDist(const char* fName)
{
	ofstream oFile(fName);
	oFile <<"InDegree\tNodeCnt"<< endl;
	for(map<int,int>::iterator aIter=indegreeDist.begin();aIter!=indegreeDist.end();aIter++)
	{
		oFile <<aIter->first <<"\t" << aIter->second << endl;
	}
	oFile <<"OutDegree\tNodeCnt" << endl;
	for(map<int,int>::iterator aIter=outdegreeDist.begin();aIter!=outdegreeDist.end();aIter++)
	{
		oFile <<aIter->first <<"\t" << aIter->second << endl;
	}
	
	oFile.close();
}

Vertex*
Graph::getVertex(const char* vname)
{
	Vertex* v=NULL;
	string key(vname);
	map<string,Vertex*>::iterator aIter=vertexList.find(key);
	if(aIter==vertexList.end())
	{
		cout <<"Warning! Null pointer returned for " << vname << endl;
	}
	else
	{
		v=aIter->second;
	}
	return v;
}


int
Graph::getDegreeDist(bool degreeType,map<int,int>& degreeDist)
{
	for(map<string,Vertex*>::iterator aIter=vertexList.begin();aIter!=vertexList.end();aIter++)
	{
		Vertex* v=aIter->second;
		int d=0;
		if(degreeType)
		{
			d=v->getInDegree();
		}
		else
		{
			d=v->getOutDegree();
		}
		if(degreeDist.find(d)==degreeDist.end())
		{
			degreeDist[d]=1;
		}
		else
		{
			degreeDist[d]=degreeDist[d]+1;
		}
	}
	return 0;
}
