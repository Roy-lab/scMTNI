#include <queue>
#include <iostream>
#include <cstring>
#include <math.h>
#include <stdlib.h>
#include "Graph.H"
#include "Vertex.H"
Vertex::Vertex()
{
	betweenness=0;
}

Vertex::~Vertex()
{
}

int
Vertex::setName(const char* aName)
{
	strcpy(vname,aName);
	return 0;
}

const char*
Vertex::getName()
{
	return vname;
}

int 
Vertex::setInNeighbour(Vertex* v)
{
	string nKey(v->getName());
	inNeighbours[nKey]=v;
	//reachableNeighbours[nKey]=1;
	return 0;
}

int
Vertex::setOutNeighbour(Vertex* v)
{
	string nKey(v->getName());
	outNeighbours[nKey]=v;
	return 0;
}

//Do a Breadth first search
int 
Vertex::findReachableNodes()
{
	//lazy bum me is going to STL queue of vertices
	queue<Vertex*> reachableQ;
	//Initialize reachableQ with in neighbours
	for(NINFO_MAP_ITER nmIter=outNeighbours.begin();nmIter!=outNeighbours.end();nmIter++)
	{
		reachableQ.push(nmIter->second);
		reachableNeighbours[nmIter->first]=1;
		nmIter->second->addReachableFromNode(vname,1);
	}

	while(!reachableQ.empty())
	{
		Vertex* v=reachableQ.front();
		//get the neighbours of aNode
		NINFO_MAP& vNeigh=v->getImmediateNeighbours();
		//This node, u is reachable to all the neighbours of v by 1 + dist(u,v)
		string vkey(v->getName());
		if(reachableNeighbours.find(vkey)==reachableNeighbours.end())
		{
			cout <<"Expected a node in reachable list, did not find " << endl;
			exit(0);
		}
		int dist_uv=reachableNeighbours[vkey];
		for(NINFO_MAP_ITER nmIter=vNeigh.begin();nmIter!=vNeigh.end();nmIter++)
		{
			if(strcmp(nmIter->first.c_str(),vname)==0)
			{
				continue;
			}
			//Check for cycles
			if(reachableNeighbours.find(nmIter->first)==reachableNeighbours.end())
			{
				//Everytime the node v acts as an intermediary, increment the betweenness count
				v->incrBetweenness(vname,nmIter->first.c_str());
				reachableNeighbours[nmIter->first]=dist_uv+1;
				reachableQ.push(nmIter->second);
				shortestPathCnts[nmIter->first]=1;
				nmIter->second->addReachableFromNode(vname,dist_uv+1);
			}
			else 
			{
				//If we have encountered this vertex before and the length is equal
				//to the length of the previously encounterd path then increment betweenness
				int currDist=reachableNeighbours[nmIter->first];
				if(currDist==(dist_uv+1))
				{
					v->incrBetweenness(vname,nmIter->first.c_str());
					shortestPathCnts[nmIter->first]=shortestPathCnts[nmIter->first]+1;
				}
				if(currDist>(dist_uv+1))
				{
					cout <<"Error! A shorter path found later on!! BFS search is faulty" << endl;
					exit(0);
				}
			}
		}
		reachableQ.pop();
	}
	return 0;
}


NINFO_MAP& 
Vertex::getImmediateNeighbours()
{
	return outNeighbours;
}

NDIST_MAP&
Vertex::getReachableNeighbours()
{
	return reachableNeighbours;
}

bool
Vertex::isReachable(Vertex* v)
{
	string vKey(v->getName());
	if(reachableNeighbours.find(vKey)==reachableNeighbours.end())
	{
		return false;
	}
	return true;
}

bool
Vertex::isReachable(string& vKey)
{
	if(reachableNeighbours.find(vKey)==reachableNeighbours.end())
	{
		return false;
	}
	return true;
}

int
Vertex::getPathLength(Vertex* v)
{
	int pathLen=-1;
	string vKey(v->getName());
	if(reachableNeighbours.find(vKey)!=reachableNeighbours.end())
	{
		pathLen=reachableNeighbours[vKey];
	}
	return pathLen;
}

int
Vertex::getInDegree()
{
	return inNeighbours.size();
}

int
Vertex::getOutDegree()
{
	return outNeighbours.size();
}

int 
Vertex::setGraph(Graph* aPtr)
{
	gPtr=aPtr;
	return 0;
}

int
Vertex::getReachableOutDegree()
{
	return reachableNeighbours.size();
}


int
Vertex::getReachableInDegree()
{
	return reachableFromNode.size();
}

int 
Vertex::setQueueStatus(Vertex::QueueStatus aStat)
{
	qStatus=aStat;
	return 0;
}

Vertex::QueueStatus 
Vertex::getQueueStatus()
{
	return qStatus;
}

int 
Vertex::addReachableNeighbour(Vertex* v,int dist)
{
	//Then for all my current reachable neighbours make them reachable from v using their current distance from this
	//node and this node's distance from the neighbours
	for(NDIST_MAP_ITER nmIter=reachableNeighbours.begin();nmIter!=reachableNeighbours.end();nmIter++)
	{
		Vertex* neighbour=gPtr->getVertex(nmIter->first.c_str());
		int currDist=nmIter->second;
		neighbour->setReachableNeighbour(v,dist+currDist);
		v->setReachableNeighbour(neighbour,dist+currDist);
	}
	setReachableNeighbour(v,dist);
	//v->setReachableNeighbour(this,dist);
	return 0;
}

int
Vertex::addReachableFromNode(const char* regName, int dist)
{
	string regNameKey(regName);
	reachableFromNode[regName]=dist;
	return 0;
}


int
Vertex::setReachableNeighbour(Vertex* v,int dist)
{
	//Check for self-loops
	if(strcmp(v->getName(),vname)==0)
	{
		return 0;
	}
	string key(v->getName());
	if(reachableNeighbours.find(key)==reachableNeighbours.end())
	{
		reachableNeighbours[key]=dist;
	}
	else
	{
		int origDist=reachableNeighbours[key];
		if(dist<origDist)
		{
			reachableNeighbours[key]=dist;
		}
	}
	return 0;
}

int 
Vertex::showReachability(ostream& oFile)
{
	oFile << vname;
	for(NDIST_MAP_ITER nmIter=reachableNeighbours.begin();nmIter!=reachableNeighbours.end();nmIter++)
	{
		oFile <<" "<< nmIter->first.c_str() <<"(" << nmIter->second << ")";
	}
	oFile<<endl;
	return 0;
}

bool
Vertex::isOutNeighbour(Vertex* putNeighbr)
{
	string aKey(putNeighbr->getName());
	bool found=false;
	if(outNeighbours.find(aKey)!=outNeighbours.end())
	{
		found=true;
	}
	return found;
}


bool
Vertex::isInNeighbour(Vertex* putNeighbr)
{
	string aKey(putNeighbr->getName());
	bool found=false;
	if(inNeighbours.find(aKey)!=inNeighbours.end())
	{
		found=true;
	}
	return found;
}


int
Vertex::incrBetweenness(const char* src, const char* dest)
{
	string sKey;
	sKey.append(src);
	sKey.append("#");
	sKey.append(dest);
	if(betweennessCnt.find(sKey)==betweennessCnt.end())
	{
		betweennessCnt[sKey]=1;
	}
	else
	{
		betweennessCnt[sKey]=betweennessCnt[sKey]+1;
	}
	return 0;
}


int
Vertex::getShortestPathCnt(string& destKey)
{
	if(shortestPathCnts.find(destKey)==shortestPathCnts.end())
	{
		cout <<"No shortest path between " << vname << " and " << destKey.c_str() << endl;
		exit(0);
	}
	return shortestPathCnts[destKey];
}

int
Vertex::computeAvgPathLength()
{
	avgPathLen=0;
	if(reachableNeighbours.size()==0)
	{
		return 0;
	}
	for(NDIST_MAP_ITER nmIter=reachableNeighbours.begin();nmIter!=reachableNeighbours.end();nmIter++)
	{
		avgPathLen=avgPathLen+nmIter->second;
	}
	avgPathLen=avgPathLen/(double) reachableNeighbours.size();
	return 0;
}


int
Vertex::computeSdPathLength()
{
	sdPathLen=0;
	if(reachableNeighbours.size()==0)
	{
		return 0;
	}
	for(NDIST_MAP_ITER nmIter=reachableNeighbours.begin();nmIter!=reachableNeighbours.end();nmIter++)
	{
		double diff=(double)nmIter->second-avgPathLen;
		sdPathLen=sdPathLen+(diff*diff);
	}
	sdPathLen=sqrt(sdPathLen/((double) (reachableNeighbours.size()-1)));
	return 0;
}

int
Vertex::computePureReachableNodes(map<string,Vertex*>& vertexList)
{
	pureTargets=0;
	for(NDIST_MAP_ITER nmIter=reachableNeighbours.begin();nmIter!=reachableNeighbours.end();nmIter++)
	{
		Vertex* tVertex=vertexList[nmIter->first];
		if(tVertex->getOutDegree()>0)
		{
			continue;
		}
		pureTargets++;
	}
	return 0;
}

int
Vertex::computeAvgTargetIndegree()
{
	avgTargetIndegree=0;
	avgTargetReachableIndegree=0;
	for(NINFO_MAP_ITER nmIter=outNeighbours.begin();nmIter!=outNeighbours.end();nmIter++)
	{
		avgTargetIndegree=avgTargetIndegree+(double)nmIter->second->getInDegree();
		avgTargetReachableIndegree=avgTargetReachableIndegree+(double)nmIter->second->getReachableInDegree();
	}
	avgTargetIndegree=avgTargetIndegree/((double)outNeighbours.size());
	avgTargetReachableIndegree=avgTargetReachableIndegree/((double)outNeighbours.size());
	return 0;

}


int
Vertex::computeAvgTFOutdegree()
{
	avgTFOutdegree=0;
	avgTFReachableOutdegree=0;
	for(NINFO_MAP_ITER nmIter=inNeighbours.begin();nmIter!=inNeighbours.end();nmIter++)
	{
		avgTFOutdegree=avgTFOutdegree+(double)nmIter->second->getOutDegree();
		avgTFReachableOutdegree=avgTFReachableOutdegree+(double)nmIter->second->getReachableOutDegree();
	}
	avgTFOutdegree=avgTFOutdegree/((double)inNeighbours.size());
	avgTFReachableOutdegree=avgTFReachableOutdegree/((double)inNeighbours.size());
	return 0;

}



double 
Vertex::getAvgCommonTargets()
{
	return avgCommonTargets;
}

double 
Vertex::getSdCommonTargets()
{
	return sdCommonTargets;
}

double
Vertex::getAvgPathLength()
{
	return avgPathLen;
}

double
Vertex::getSdPathLength()
{
	return sdPathLen;
}


double 
Vertex::getAvgTargetIndegree()
{
	return avgTargetIndegree;
}

double 
Vertex::getAvgTargetReachableIndegree()
{
	return avgTargetReachableIndegree;
}

int
Vertex::getPureReachableNodes()
{
	return pureTargets;
}

double 
Vertex::getAvgTFOutdegree()
{
	return avgTFOutdegree;
}

double 
Vertex::getAvgTFReachableOutdegree()
{
	return avgTFReachableOutdegree;
}

