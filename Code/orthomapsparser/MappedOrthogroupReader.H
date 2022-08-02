#ifndef _ORTHOGROUP_READER
#define _ORTHOGROUP_READER
#include <map>
using namespace std;

class MappedOrthogroup;
class MappedOrthogroupReader
{
	public:
	MappedOrthogroupReader();
	~MappedOrthogroupReader();
	int readFile(const char*);
	int readSpeciesMapping(const char*);
	int getMappedOrthogroupID(const char*, const char*);
	MappedOrthogroup* getMappedOrthogroup(const char*, const char*);
	map<int,MappedOrthogroup*>& getMappedOrthogroups();
	STRINTMAP* getOrtholog(const char* srcSpecName, const char* geneName, const char* targetSpecName);

	private:
	int addMembers(char*,MappedOrthogroup*);
	int generateGeneOrthoMap();
	map<int,MappedOrthogroup*> orthogroupSet;
	map<string,STRINTMAP*> geneOrthoMap;
	map<int,string> speciesIDNameMap;
	bool persistentDebug;
};
#endif
