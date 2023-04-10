/*This class is a node of a graph*/
#ifndef _VERTEX
#define _VERTEX

#include <map>
#include <string>
#include <vector>
#include <iostream>

using namespace std;
//This map is used to store stuff like the distance between the current
//vertex and all its neighbours
class Vertex;
typedef map<string,Vertex*> NINFO_MAP;
typedef map<string,Vertex*>::iterator NINFO_MAP_ITER;

typedef vector<string> STRVECT;

//Distance to different nodes
typedef map<string,int> NDIST_MAP;
typedef map<string,int>::iterator NDIST_MAP_ITER;
class Graph;
class Vertex
{
	public:
		Vertex();
		~Vertex();
		int setName(const char*);
		const char* getName();
		//Set by reading the interaction graph
		int setInNeighbour(Vertex*);
		int setOutNeighbour(Vertex*);
		//Do BFS starting from this node
		int findReachableNodes();
		//Do BFS but preserved the topological sort order
		int findReachableNodes_Ordered(int);
		int doNeighbourhoodWalk(int);
		map<string,NDIST_MAP*>& getPathSet();
		map<string,vector<string>*>& getPathSet_Ordered();
		map<string,vector<string>*>& getCycleSet();

		int setOrder(int);
		int getOrder();
		int setPredecessor(Vertex*);
		int resetPredecessor();
		Vertex* getPredecessor();

		
		NINFO_MAP& getImmediateNeighbours();
		NDIST_MAP& getReachableNeighbours();
		bool isReachable(Vertex*);
		bool isReachable(string&);
		int getPathLength(Vertex*);
		int getInDegree();	
		int getOutDegree();	
		
		int getReachableOutDegree();
		int getReachableInDegree();
		typedef enum
		{
			//Not been on the queue yet
			NO_QUEUE,
			//Has been pushed on the queue
			ON_QUEUE,
			//Has been popped off the queue
			OFF_QUEUE
		}QueueStatus;
		int setGraph(Graph*);
		int setQueueStatus(Vertex::QueueStatus);
		Vertex::QueueStatus getQueueStatus();
		int addReachableNeighbour(Vertex*,int);
		int addReachableFromNode(const char*,int);
		int setReachableNeighbour(Vertex*,int);
		int showReachability(ostream&);
		int showReachability_Ordered(ostream&);
		int showReachability_Ordered_All(ostream&);
		int showReachability_Ordered_All(ostream&,ostream&);

		//Clique related stuff
		bool isInNeighbour(Vertex*);
		bool isOutNeighbour(Vertex*);

		int incrBetweenness(const char*, const char*);
		int getShortestPathCnt(string&);
		
		int computeCloseness(map<string,Vertex*>&);
		int computeBetweenness(map<string,Vertex*>&,int);
		int computeOverlappingTargets(map<string,Vertex*>&);
		int computeAvgPathLength();
		int computeSdPathLength();
		int computeAvgTargetIndegree();
		int computePureReachableNodes(map<string,Vertex*>&);
		int computeAvgTFOutdegree();
		
		int computeRegPathRatio();
		int computeRegulatibility();

		double getCloseness();
		double getBetweenness();
		double getBetweennessNorm();
		double getBridgeness();
		int getCoRegulators();

		double getAvgCommonTargets();
		double getSdCommonTargets();
		
		double getAvgPathLength();
		double getSdPathLength();
		double getAvgTargetIndegree();
		double getAvgTargetReachableIndegree();
		int getPureReachableNodes();
		double getAvgTFOutdegree();
		double getAvgTFReachableOutdegree();
		
		double getRegPathRatio();
		double getRegulatibility();
		
	private:
		int generateKey(NDIST_MAP*,string&);
		int generateKey(vector<string>*,string&);
		//These are the in neighours
		NINFO_MAP inNeighbours;
		NINFO_MAP outNeighbours;
		NDIST_MAP reachableNeighbours;
		//The number of shortestpaths that exist between this node and other nodes
		NDIST_MAP shortestPathCnts;
		char vname[256];
		QueueStatus qStatus;
		Graph* gPtr;

		//Get closeness to vertices to all vertices
		double closeness;
		double betweenness;
		double betweennessNorm;
		
		//This is the number of shortest paths that go through this node
		double bridgeness;
		double regPathRatio;
		double regulatibility;
		double avgPathLen;
		double sdPathLen;
		double avgTargetIndegree;
		double avgTargetReachableIndegree;
		double avgTFOutdegree;
		double avgTFReachableOutdegree;
		int pureTargets;
		double avgCommonTargets;
		double sdCommonTargets;
		//The key here is a pair of vertices
		//and the value is the number of times this vertex is on the shortest path between
		//these two vertices
		map<string,int> betweennessCnt;
		map<string,int> overlappingTargets;
		map<string,int> reachableFromNode;

		int vOrder;
		Vertex* predecessor;

		map<string,NDIST_MAP*> pathSet;
		map<string,vector<string>*> pathSet_ordered;
		map<string,vector<string>*> cycleSet;
};
#endif
