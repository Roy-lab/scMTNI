/*Stores a clique of vertices. This can be considered as a compound node. We need to inherit from*/

#ifndef _CLIQUE
#define _CLIQUE
#include <map>
#include <fstream>
using namespace std;

class Vertex;
class Clique
{
	public:
		Clique();
		~Clique();
		int addMember(Vertex*);
		bool inClique(Vertex*);
		int showClique(ofstream&);
		map<string,Vertex*>& getMembers();
	private:
		map<string,Vertex*> cliqueMembers;
};

#endif
