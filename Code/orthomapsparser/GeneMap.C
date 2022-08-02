#include <iostream>
#include "GeneMap.H"

GeneMap::GeneMap()
{
}

GeneMap::~GeneMap()
{
}
//remove targetSpecies by shilu
int
GeneMap::addPair(const string& srcGene, const string& targetGene)
{
    STRINTMAP* hitsForGene=NULL;
    if(geneSet0.find(srcGene)==geneSet0.end())
    {
        hitsForGene=new STRINTMAP;
        geneSet0[srcGene]=hitsForGene;
    }
    else
    {
        hitsForGene=geneSet0[srcGene];
    }
    (*hitsForGene)[targetGene]=0;
    return 0;
}
/*
 //geneSet[srcGene][targetSpecies]=targetGene
 int
 GeneMap::addPair(const string& srcGene, const string& targetSpecies, const string& targetGene)
 {
     map<string,STRINTMAP*>* hitsForGene=NULL;
     if(geneSet.find(srcGene)==geneSet.end())
     {
         hitsForGene=new map<string,STRINTMAP*>;
         geneSet[srcGene]=hitsForGene;
     }
     else
     {
         hitsForGene=geneSet[srcGene];
     }
     if(targetSpecies.length()==0)
     {
         return 0;
     }
     STRINTMAP* hitsInSpecies=NULL;
     if(hitsForGene->find(targetSpecies)==hitsForGene->end())
     {
         hitsInSpecies=new STRINTMAP;
         (*hitsForGene)[targetSpecies]=hitsInSpecies;
     }
     else
     {
         hitsInSpecies=(*hitsForGene)[targetSpecies];
     }
     (*hitsInSpecies)[targetGene]=0;
     return 0;
 }
 */

//geneSet[srcGene][targetSpecies][targetGene]=0
int
GeneMap::addPair(const string& srcGene, const string& targetSpecies, const string& targetGene)
{
    cout <<"GeneMap::addPair srcGene=" << srcGene << " targetSpecies=" << targetSpecies << " targetGene=" << targetGene << endl;
	map<string,STRINTMAP*>* hitsForGene=NULL;
	if(geneSet.find(srcGene)==geneSet.end())
	{
		hitsForGene=new map<string,STRINTMAP*>;
		geneSet[srcGene]=hitsForGene;
	}
	else
	{
		hitsForGene=geneSet[srcGene];
	}
	if(targetSpecies.length()==0)
	{
		return 0;
	}
	STRINTMAP* hitsInSpecies=NULL;
	if(hitsForGene->find(targetSpecies)==hitsForGene->end())
	{
		hitsInSpecies=new STRINTMAP;
		(*hitsForGene)[targetSpecies]=hitsInSpecies;	
	}
	else
	{
		hitsInSpecies=(*hitsForGene)[targetSpecies];	
	}
	(*hitsInSpecies)[targetGene]=0;
	return 0;
}


STRINTMAP* 
GeneMap::getHits(const char* srcGene, const char* targetSpecies)
{
	string sKey(srcGene);
	if(geneSet.find(sKey)==geneSet.end())
	{
		cout <<"Weird! No  gene " << sKey << " in its orthogroups"<< endl;
		return NULL;
	}
	map<string,STRINTMAP*>* hits=geneSet[sKey];
	string tKey(targetSpecies);
	if(hits->find(tKey)==hits->end())
	{
		return NULL;
	} 
	return (*hits)[tKey];
}

map<string,map<string,STRINTMAP*>*>&
GeneMap::getGeneSet()
{
	return geneSet;
}
