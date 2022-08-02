#ifndef _ORTHOGROUP_
#define _ORTHOGROUP_
#include <string>
#include <map>
#include <vector>
using namespace std;

class GeneMap;
class MappedOrthogroup
{
	public:
	MappedOrthogroup();
	~MappedOrthogroup();
	//int setMembers(map<string,string>&);
    int addMembers(string&);
	int setID(int);
	int getID();
	//int incrCnt();
	//int getCnt();

	vector<string>& getOrthoMembers();
    //map<string,GeneMap*>& getOrthoMembers();
	//GeneMap* getSpeciesHits(const char* specName);
	//STRINTMAP* getSpeciesHitsForGene(const char* srcSpecName, const char* targetSpecName, const char* geneName);
	//map<int,map<string,string>*>& getGeneSets();
	private:
    //map<string,string> orthoMembers;
    vector<string> orthoMembers;
    //map<string,GeneMap*> orthoMembers;
	//map<string,int> scerHits;
	//map<int,map<string,string>*> geneSets;
	//map<string,int> specList;
	int oid;
	//int cnts;
};
#endif
