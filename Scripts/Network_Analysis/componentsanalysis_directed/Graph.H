#ifndef _GRAPH
#define _GRAPH
#include <map>
#include <string>
#include <vector>
using namespace std;

class Vertex;
class Clique;
class Graph
{
	public:
		Graph();
		~Graph();
		int makeGraph(const char*);
		int dumpGraph();
		map<string,Vertex*>& getVertexList();
		int showNodeProperties(const char*);
		int obtainConnectivity();
		int showConnectivity();
		int showPathLenFrequency(const char*);

		int obtainConnectivity_Efficient();
		bool isConnected();

		int getConnectedComponents(const char*);
		int getConnectedComponents_Linewise(const char*);
		int getConnectedComponents(vector<map<string,Vertex*>*>&);
		
		//Compute average distance between a pair of nodes that are connected
		double computeAveragePathLength();
		int computeAverageBetweenness();
		int computeAveragePureReachabilityDiff();
		int getTotalPathCnt();
		//Measure the distance between graphs
		int getPathLength(const char*, const char*,int&);
		int computeGraphDistance(Graph*);
		int computeKLDivergence(Graph*);
		int computePathLenFrequency();
		int smoothPathLenDist(Graph* );
		int computeSymmetricKLDivergence(Graph*);
		double getKLDivergence(Graph*);
		int getDegreeDist();
		
		int showDegreeDist(const char*);
		int getNodesWithDegreeAtleast(int,const char*,const char*);

		Vertex* getVertex(const char*);

		int getAssortativity();
		int showEdgeDegrees(const char*);

		int computeDifferentNodeProperties();

		//Find the number of reachable nodes that are common between
		//different transcription factors. We want this number to be low
		//in order for the network to be modular
		int computeTargetOverlap();
		//Here we want to see if study the importance of a target node.
		//We want to study the number of paths that terminate in a target node
		//and how this changes with the indegree of the node. The number of regulatory
		//paths is atleast equal to the indegree of a node. We can look at the ratio
		//of the indegree to the number of paths terminating in the node.
		//A number close to 1 means that the number of regulatory paths is close to
		//the number regulators. But a number smaller than 1 means that this node
		//can be tweaked by many more regulators other than its immediate nodes.
		int computeAvgRegPathRatio();
	
		int doTopologicalSort();
		int genSubGraphs(int maxSubGSize, ofstream& oFile);
		int genSubGraphs(int maxSubGSize, ofstream& oFile1, ofstream& oFile2);
		int genNeighbourhood(int maxRadius, ofstream& oFile1, ofstream& oFile2);
		int genNeighbourhood_AllAtaTime(int maxRadius, ofstream& oFile1, ofstream& oFile2);

		int findAllCliques(ofstream& oFile1,ofstream& oFile2);
		int findMaximalCliques();
		int showCliques(ofstream& oFile1, ofstream & oFile2);
		int showCliques(ofstream& oFile);

		int genCycles(int,ofstream&,ofstream&);
		string& getEdgeHeader();
		vector<string>* getEdgeAttribs(string &);
	private:
		int resetPredecessor();
		int getDegreeDist(bool,map<int,int>&);
		int dfsVisit(Vertex*,int);

		map<int,int> pathLenFreqCnt;
		map<int,double> pathLenDist;
		map<string,Vertex*> vertexList;
		map<string,Vertex*> untouchedVertexMap;
		map<int,int> indegreeDist;
		map<int,int> outdegreeDist;

		double harmonicMeanLength;
		double outIncorrCoeff;
		double reachoutIncorrCoeff;
		double outOutcorrCoeff;
		int totalPaths;
		double avgPathLength;
		double avgBetweenness;
		double avgBetweennessNorm;
		double avgBridgeness;
		double sdPathLength;
		double avgRegPathRatio;
		double sdRegPathRatio;
		double avgRegulatibility;
		double sdRegulatibility;
		double avgPureReachabilityDiff;

		//Stuff related to the topological sort
		map<string,int> visitFlag;
		vector<Vertex*> topSorted;
		vector<Vertex*> topSorted_R;
		vector<Clique*> cliqueSet;

		//Store edge attributes
		map<string,vector<string>*> edgeAttribs; //The string in each vector just stores all the attributes in one string
		string edgeHeader;
};

#endif
