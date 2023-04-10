#include <iostream>

using namespace std;

#include "Vertex.H"
#include "Clique.H"


Clique::Clique()
{
}

Clique::~Clique()
{
}

int 
Clique::addMember(Vertex* anode)
{
	string akey(anode->getName());
	cliqueMembers[akey]=anode;
	return 0;
}


//Check if the node n is in this clique, i.e. it is connected to all the nodes in this
//clique
bool
Clique::inClique(Vertex* n)
{
	bool yes=true;
	map<string,Vertex*>::iterator nodeIter=cliqueMembers.begin();
	while(nodeIter!=cliqueMembers.end() && yes)
	{
		Vertex* member=nodeIter->second;
		if((!member->isOutNeighbour(n)) || (!n->isOutNeighbour(member)))
		{
			yes=false;
		}
		else
		{
			nodeIter++;
		}
	}
	return yes;
}

int
Clique::showClique(ofstream& oFile)
{
	for(map<string,Vertex*>::iterator cIter=cliqueMembers.begin();cIter!=cliqueMembers.end();cIter++)
	{
		if(cIter!=cliqueMembers.begin())
		{
			oFile <<"#";
		}
		oFile << cIter->first.c_str();
	}
	oFile << endl;
	return 0;
}


map<string,Vertex*>& 
Clique::getMembers()
{
	return cliqueMembers;
}
