/*Stores for every species, a list of genes. For each gene it stores its orthology relationships in several species*/
#ifndef _GENEMAP_
#define _GENEMAP_

#include <map>
#include <string>

using namespace std;
typedef map<string,int> STRINTMAP;
typedef map<string,int>::iterator STRINTMAP_ITER;
class GeneMap
{
	public:
	GeneMap();
	~GeneMap();
	int addPair(const string& srcGene, const string& targetSpecies, const string& targetGene);
    int addPair(const string& srcGene, const string& targetGene);
	map<string,int>* getHits(const char* srcGene, const char* targetSpecies);
	map<string,map<string,STRINTMAP*>*>& getGeneSet();

	private:
	map<string,map<string,STRINTMAP*>*> geneSet;
    map<string,STRINTMAP*> geneSet0;
};
#endif
