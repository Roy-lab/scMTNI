#include <queue>
#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include "Graph.H"
#include "Vertex.H"
Vertex::Vertex()
{
	betweenness=0;
	vOrder=-1;
	predecessor=NULL;
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


//Do a Breadth first search but preserve the topological sort order
int 
Vertex::findReachableNodes_Ordered(int maxSubGraphSize)
{
	//lazy bum me is going to STL queue of vertices
	queue<Vertex*> reachableQ;
	int currSubGraphSize=1;
	//Initialize reachableQ with in neighbours
	for(NINFO_MAP_ITER nmIter=outNeighbours.begin();nmIter!=outNeighbours.end();nmIter++)
	{
		if(nmIter->second->getPredecessor()!=NULL)
		{
			continue;
		}
		if(nmIter->second->getOrder()>vOrder)
		{
			reachableQ.push(nmIter->second);
			reachableNeighbours[nmIter->first]=1;
			if(nmIter->second->setPredecessor(this)==-1)
			{
				return -1;
			}
		}
	}

	while((!reachableQ.empty()) && (currSubGraphSize<maxSubGraphSize))
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
			if(nmIter->second->getPredecessor()!=NULL)
			{
				continue;
			}
			int neighbourOrder=nmIter->second->getOrder();
			if( (neighbourOrder<=v->getOrder()) || (neighbourOrder <=vOrder))
			{
				continue;
			}
			reachableNeighbours[nmIter->first]=dist_uv+1;
			reachableQ.push(nmIter->second);
			if(nmIter->second->setPredecessor(v)==-1)
			{
				return -1;
			}
		}
		reachableQ.pop();
		currSubGraphSize++;
	}
	return 0;
}


//Do a neighbourhood walk 
int 
Vertex::doNeighbourhoodWalk(int maxRadius)
{
	//lazy bum me is going to STL queue of vertices
	queue<Vertex*> reachableQ;
	int currRadius=1;
	//Initialize reachableQ with in neighbours
	string myname(vname);
	map<string,STRVECT*> parentPaths;
	for(NINFO_MAP_ITER nmIter=outNeighbours.begin();nmIter!=outNeighbours.end();nmIter++)
	{
		if(nmIter->second==this)
		{
			continue;
		}
		reachableQ.push(nmIter->second);
		//Make a path
		NDIST_MAP* path=new NDIST_MAP;
		vector<string>* ord_path=new vector<string>;
		ord_path->push_back(myname);
		ord_path->push_back(nmIter->first);
		(*path)[myname]=0;
		(*path)[nmIter->first]=0;
		string pkey;
		generateKey(path,pkey);
		pathSet[pkey]=path;
		pathSet_ordered[pkey]=ord_path;
		//Save the endpoints of the path
		STRVECT* sVect=new STRVECT;
		sVect->push_back(nmIter->first);
		parentPaths[pkey]=sVect;
	}

	while((parentPaths.size()>0) && (currRadius<maxRadius))
	{
		map<string,STRVECT*> tempParents;
		for(map<string,STRVECT*>::iterator pIter=parentPaths.begin();pIter!=parentPaths.end();pIter++)
		{
			//Get the path and extend it
			NDIST_MAP* ppath=pathSet[pIter->first];
			vector<string>* ord_path=pathSet_ordered[pIter->first];
			STRVECT* endPoints=pIter->second;
			for(int i=0;i<endPoints->size();i++)
			{
				Vertex* v=gPtr->getVertex((*endPoints)[i].c_str());
				//get the neighbours of aNode
				NINFO_MAP& vNeigh=v->getImmediateNeighbours();
				for(NINFO_MAP_ITER nmIter=vNeigh.begin();nmIter!=vNeigh.end();nmIter++)
				{
					if(ppath->find(nmIter->first)!=ppath->end())
					{
						//Found a cycle, but make sure that nmIter->first is not
						//just the previous vertex v
						if((strcmp(nmIter->first.c_str(),myname.c_str())==0) &&(ppath->size()>2))
						{
							string cycleKey;
							generateKey(ppath,cycleKey);
							cycleSet[cycleKey]=ord_path;
						}
						continue;
					}
					//Make  a new path
					NDIST_MAP* newpath=new NDIST_MAP;
					vector<string>* ord_newpath=new vector<string>;

					for(int j=0;j<ord_path->size();j++)
					{
						string& existName=(*ord_path)[j];
						(*newpath)[existName]=0;
						ord_newpath->push_back(existName);
					}
					(*newpath)[nmIter->first]=0;
					ord_newpath->push_back(nmIter->first);
					string pkey;
					generateKey(newpath,pkey);
					//generateKey(ord_newpath,pkey);
					if(pathSet.find(pkey)==pathSet.end())
					{
						pathSet[pkey]=newpath;
						pathSet_ordered[pkey]=ord_newpath;
						STRVECT* sVect=new STRVECT;
						sVect->push_back(nmIter->first);
						tempParents[pkey]=sVect;
					}
					else
					{	
						//Update only the endpoints
						if(tempParents.find(pkey)!=tempParents.end())
						{
							STRVECT* sVect=tempParents[pkey];
							sVect->push_back(nmIter->first);
						}
						else
						{
							cout <<"Path " << pkey.c_str() << " proposed, present in pathSet, but absent in tempParents" << endl;
						}
						delete newpath;
						delete ord_newpath;
					}
				}
			}
		}
		parentPaths.clear();
		for(map<string,STRVECT*>::iterator aIter=tempParents.begin();aIter!=tempParents.end();aIter++)
		{
			parentPaths[aIter->first]=aIter->second;
		}
		tempParents.clear();
		currRadius++;
	}
	return 0;
}

map<string,NDIST_MAP*>&
Vertex::getPathSet()
{	
	return pathSet;
}

map<string,vector<string>*>&
Vertex::getPathSet_Ordered()
{
	return pathSet_ordered;
}

map<string,vector<string>*>&
Vertex::getCycleSet()
{
	return cycleSet;
}

int 
Vertex::setOrder(int anOrder)
{
	vOrder=anOrder;
	return 0;
}

int 
Vertex::getOrder()
{
	return vOrder;
}

int 
Vertex::setPredecessor(Vertex* vPtr)
{
	if(predecessor!=NULL)
	{
		cout <<"Predecessor pointer for " << vname << " already set to " << predecessor->getName() << endl;
		cout <<"Failed to set to " << vPtr->getName() << endl;
		return -1;
	}
	predecessor=vPtr;
	return 0;
}

int
Vertex::resetPredecessor()
{
	predecessor=NULL;
	return 0;
}
Vertex* 
Vertex::getPredecessor()
{
	return predecessor;
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


int 
Vertex::showReachability_Ordered(ostream& oFile)
{
	for(NDIST_MAP_ITER nmIter=reachableNeighbours.begin();nmIter!=reachableNeighbours.end();nmIter++)
	{
		Vertex* reachableNode=gPtr->getVertex(nmIter->first.c_str());
		vector<Vertex*>	path;
		Vertex* currNode=reachableNode;
		Vertex* predNode=currNode->getPredecessor();
		while(predNode!=this)
		{
			path.push_back(predNode);
			currNode=predNode;
			predNode=currNode->getPredecessor();
		}
		oFile << vname;
		for(int i=0;i<path.size() ;i++)
		{
			oFile <<"\t"<< path[i]->getName();
		}
		oFile <<"\t"<< nmIter->first.c_str() << endl;
	}
	return 0;
}


int 
Vertex::showReachability_Ordered_All(ostream& oFile)
{
	int currDist=1;
	int shownCnt=1;
	while(shownCnt>0)
	{
		shownCnt=0;
		for(NDIST_MAP_ITER nmIter=reachableNeighbours.begin();nmIter!=reachableNeighbours.end();nmIter++)
		{
			if(nmIter->second!=currDist)
			{
				continue;
			}
			Vertex* reachableNode=gPtr->getVertex(nmIter->first.c_str());
			vector<Vertex*>	path;
			Vertex* currNode=reachableNode;
			Vertex* predNode=currNode->getPredecessor();
			while(predNode!=this)
			{
				path.push_back(predNode);
				currNode=predNode;
				predNode=currNode->getPredecessor();
			}
			oFile << vname;
			for(int i=0;i<path.size() ;i++)
			{
				oFile <<"-"<< path[i]->getName();
			}
			oFile <<"-"<< nmIter->first.c_str() << endl;
			shownCnt++;
		}
		currDist++;
	}
	return 0;
}


int 
Vertex::showReachability_Ordered_All(ostream& oFile1, ostream& oFile2)
{
	int currDist=1;
	int shownCnt=1;
	while(shownCnt>0)
	{
		shownCnt=0;
		for(NDIST_MAP_ITER nmIter=reachableNeighbours.begin();nmIter!=reachableNeighbours.end();nmIter++)
		{
			if(nmIter->second!=currDist)
			{
				continue;
			}
			Vertex* reachableNode=gPtr->getVertex(nmIter->first.c_str());
			vector<Vertex*>	path;
			path.push_back(this);
			Vertex* currNode=reachableNode;
			Vertex* predNode=currNode->getPredecessor();
			while(predNode!=this)
			{
				path.push_back(predNode);
				currNode=predNode;
				predNode=currNode->getPredecessor();
			}
			path.push_back(reachableNode);
			for(int i=0;i<path.size() ;i++)
			{
				if(i>0)
				{
					oFile1 <<"#";
				}
				oFile1 << path[i]->getName();
			}
			oFile1 << endl;
			bool startDisp=false;
			//Get the edge set associated with these vertices
			for(int i=0;i<path.size();i++)
			{
				Vertex* u=path[i];
				for(int j=i+1;j<path.size();j++)
				{
					Vertex* v=path[j];
					if(u->isOutNeighbour(v))
					{
						if(!startDisp)
						{
							startDisp=true;
						}
						else
						{
							oFile2<<"\t";
						}
						oFile2<< u->getName()<<"#"<< v->getName();
					}
				}
			}
			oFile2<< endl;
			shownCnt++;
		}
		currDist++;
	}
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
Vertex::computeCloseness(map<string,Vertex*>& allVertices)
{
	/*closeness=(double) reachableNeighbours.size();
	if(closeness==0)
	{
		return 0;
	}
	double totalPathLen=0;
	//This definition is taken from Crucitti
	for(NDIST_MAP_ITER aIter=reachableNeighbours.begin();aIter!=reachableNeighbours.end();aIter++)
	{
		double pathLen=aIter->second;
		totalPathLen=totalPathLen+pathLen;
	}
	closeness=closeness/totalPathLen;
	double deflen=(double)allVertices.size()+1;
	double disconnPathLen=(allVertices.size()-reachableNeighbours.size()-1)*deflen;*/
	closeness=(double) allVertices.size();
	double totalPathLen=0;
	double defLen=(double) (allVertices.size()+1);
	for(map<string,Vertex*>::iterator aIter=allVertices.begin();aIter!=allVertices.end();aIter++)
	{
		if(reachableNeighbours.find(aIter->first)==reachableNeighbours.end())
		{
			totalPathLen=totalPathLen+defLen;
		}
		else
		{
			totalPathLen=totalPathLen+reachableNeighbours[aIter->first];
		}
	}
	closeness=closeness/totalPathLen;	
	//closeness=closeness_Conn/(totalPathLen + disconnPathLen);
	return 0;
}

int
Vertex::computeBetweenness(map<string,Vertex*>& allVertices, int totalPathsinGraph)
{
	//Here we look at the betweennessCnt map. The betweenness is the sum of between coefficient
	//over all short paths on which this node occurs
	betweenness=0;
	bridgeness=0;
	int totalPaths=0;
	for(map<string,int>::iterator aIter=betweennessCnt.begin();aIter!=betweennessCnt.end();aIter++)
	{
		//Get the first and second vertices from the key
		int pos=aIter->first.find('#');
		string firstName=aIter->first.substr(0,pos);
		string secondName=aIter->first.substr(pos+1,aIter->first.length()-1);
		if(allVertices.find(firstName)==allVertices.end())
		{
			cout << "No vertex by name " << firstName.c_str() << endl;
			exit(0);
		}
		if(allVertices.find(secondName)==allVertices.end())
		{
			cout << "No vertex by name " << secondName.c_str() << endl;
			exit(0);
		}
		Vertex* srcVertex=allVertices[firstName];
		//Get the number of shortestpaths between srcVertex and destVertex
		int spathCnt=srcVertex->getShortestPathCnt(secondName);
		double bCoeff=(double)aIter->second/(double) spathCnt;
		betweenness=betweenness+bCoeff;
		totalPaths=totalPaths+aIter->second;
		bridgeness=bridgeness+aIter->second;
	}
	if(betweenness>0)
	{
		betweennessNorm=betweenness/totalPathsinGraph;
	}
	return 0;
}

int
Vertex::computeOverlappingTargets(map<string,Vertex*>& tfNodes)
{
	avgCommonTargets=0;
	int downstreamTF=0;
	for(map<string,Vertex*>::iterator aIter=tfNodes.begin();aIter!=tfNodes.end();aIter++)
	{
		if(aIter->second==this)
		{
			continue;
		}
		Vertex* othertf=aIter->second;
		int overlapCnt=0;
		//Tf must not be downstream. Otherwise obviously we will have lots of shared nodes
		if(reachableNeighbours.find(aIter->first)!=reachableNeighbours.end())
		{
			downstreamTF++;
			continue;
		}
		
		//Of my reachable nodes how many are reachable by this tf
		for(NDIST_MAP_ITER nmIter=reachableNeighbours.begin();nmIter!=reachableNeighbours.end();nmIter++)
		{
			if(othertf->isReachable((string&)nmIter->first))
			{
				overlapCnt++;
			}
		}
		//Of my reachable nodes how many are reachable by this tf in one step
		/*for(NINFO_MAP_ITER nmIter=outNeighbours.begin();nmIter!=outNeighbours.end();nmIter++)
		{
			if(othertf->getPathLength(nmIter->second)==1)
			{
				overlapCnt++;
			}
		}*/
		if(overlapCnt>0)
		{
			overlappingTargets[aIter->first]=overlapCnt;
			avgCommonTargets=avgCommonTargets+overlapCnt;
		}
	}
	avgCommonTargets=avgCommonTargets/(double)tfNodes.size();
	sdCommonTargets=0;
	for(map<string,Vertex*>::iterator aIter=tfNodes.begin();aIter!=tfNodes.end();aIter++)
	{
		int overlapCnt=0;
		if(overlappingTargets.find(aIter->first)!=overlappingTargets.end())
		{
			overlapCnt=overlappingTargets[aIter->first];
		}
		double diff=overlapCnt-avgCommonTargets;
		sdCommonTargets=sdCommonTargets+(diff*diff);
	}
	sdCommonTargets=sqrt(sdCommonTargets/((double) (tfNodes.size()-1)));
	avgCommonTargets=avgCommonTargets/((double) outNeighbours.size());
	return 0;
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



int
Vertex::computeRegPathRatio()
{
	//regPathRatio=((double)inNeighbours.size())/((double) reachableFromNode.size());
	regPathRatio=((double) reachableFromNode.size())/((double) inNeighbours.size());
	return 0;
}

int
Vertex::computeRegulatibility()
{
	regulatibility=0;
	for(map<string,int>::iterator aIter=reachableFromNode.begin();aIter!=reachableFromNode.end();aIter++)
	{
		regulatibility=regulatibility+(double) aIter->second;
	}
	regulatibility=regulatibility/((double) reachableFromNode.size());
	return 0;
}


double 
Vertex::getCloseness()
{
	return closeness;
}

double 
Vertex::getBetweenness()
{
	return betweenness;
}

double
Vertex::getBetweennessNorm()
{
	return betweennessNorm;
}

double
Vertex::getBridgeness()
{
	return bridgeness;
}

int
Vertex::getCoRegulators()
{
	return overlappingTargets.size();
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

double 
Vertex::getRegPathRatio()
{
	return regPathRatio;
}

double 
Vertex::getRegulatibility()
{
	return regulatibility;
}

int
Vertex::generateKey(NDIST_MAP* path, string& key)
{
	for(NDIST_MAP_ITER pIter=path->begin();pIter!=path->end();pIter++)
	{
		if(pIter!=path->begin())
		{
			key.append("-");
		}
		key.append(pIter->first.c_str());
	}
	return 0;
}


int
Vertex::generateKey(vector<string>* path, string& key)
{
	for(int i=0;i<path->size();i++)
	{
		if(i!=0)
		{
			key.append("-");
		}
		key.append((*path)[i].c_str());
	}
	return 0;
}
