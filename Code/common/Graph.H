#ifndef _GRAPH
#define _GRAPH
#include <map>
#include <string>
#include <vector>

using namespace std;

class Vertex;
class Graph
{
	public:
		Graph();
		~Graph();
		int makeGraph(const char*);
		int makeGraph(const char*,int);
		int getVertexCnt();
		int dumpGraph();
		int obtainConnectivity();
		int showConnectivity();

		bool isConnected();
		int getDegreeDist();
		int showDegreeDist(const char*);

		Vertex* getVertex(const char*);

	private:
		int getDegreeDist(bool,map<int,int>&);

		map<string,Vertex*> vertexList;
		map<string,Vertex*> untouchedVertexMap;
		map<int,int> indegreeDist;
		map<int,int> outdegreeDist;
		
};

#endif
