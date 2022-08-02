#include <iostream>
#include "GeneMap.H"
#include "MappedOrthogroup.H"


MappedOrthogroup::MappedOrthogroup()
{
	//cnts=0;
}

MappedOrthogroup::~MappedOrthogroup()
{
}

/*
int
MappedOrthogroup::setMembers(map<string,string>& speciesGeneMap)
{
    map<string,string>* geneSet=new map<string,string>;
    int setid=geneSets.size();
    geneSets[setid]=geneSet;
    for(map<string,string>::iterator sIter=speciesGeneMap.begin();sIter!=speciesGeneMap.end();sIter++)
    {
        (*geneSet)[sIter->first]=sIter->second;
        GeneMap* speciesHit_s=NULL;
        if(orthoMembers.find(sIter->first)==orthoMembers.end())
        {
            speciesHit_s=new GeneMap;
            orthoMembers[sIter->first]=speciesHit_s;
        }
        else
        {
            speciesHit_s=orthoMembers[sIter->first];
        }
        if(speciesGeneMap.size()==1)
        {
            speciesHit_s->addPair(sIter->second,"","");
        }
        map<string,string>::iterator tIter=sIter;
        tIter++;
        for(;tIter!=speciesGeneMap.end();tIter++)
        {
            GeneMap* speciesHit_t=NULL;
            speciesHit_s->addPair(sIter->second,tIter->first,tIter->second);
            if(orthoMembers.find(tIter->first)==orthoMembers.end())
            {
                speciesHit_t=new GeneMap;
                orthoMembers[tIter->first]=speciesHit_t;
            }
            else
            {
                speciesHit_t=orthoMembers[tIter->first];
            }
            speciesHit_t->addPair(tIter->second,sIter->first,sIter->second);
        }
        cout << endl;
    }
    return 0;
}
*/
/*
int 
MappedOrthogroup::setMembers(map<string,string>& speciesGeneMap)
{
	map<string,string>* geneSet=new map<string,string>;
    //cout <<"MappedOrthogroup::setMembers speciesGeneMap.size=" << speciesGeneMap.size() << endl;
	for(map<string,string>::iterator sIter=speciesGeneMap.begin();sIter!=speciesGeneMap.end();sIter++)
	{
        orthoMembers[sIter->first]=sIter->second;
        //cout <<"speciesGeneMap["<<sIter->first <<"]=" <<sIter->second << endl;
	}
    //orthoMembers first is species, second is gene
    for(auto &mem: orthoMembers0)
    {
        cout <<"orthoMembers0: " <<mem << endl;
    }
	return 0;
}*/

int
MappedOrthogroup::addMembers(string &geneName)
{
    orthoMembers.push_back(geneName);
    //orthoMembers first is species, second is gene
    return 0;
}

int 
MappedOrthogroup::setID(int aid)
{
	oid=aid;
	return 0;
}

int 
MappedOrthogroup::getID()
{
	return oid;
}

vector<string>&
MappedOrthogroup::getOrthoMembers()
{
    return orthoMembers;
}

/*
int 
MappedOrthogroup::incrCnt()
{
	cnts=cnts+1;
	return 0;
}

int 
MappedOrthogroup::getCnt()
{
	return cnts;
}*/

/*
map<string,GeneMap*>&
MappedOrthogroup::getOrthoMembers()
{
    return orthoMembers;
}
 
GeneMap* 
MappedOrthogroup::getSpeciesHits(const char* specName)
{
	string key(specName);
	if(orthoMembers.find(key)==orthoMembers.end())
	{
		return NULL;
	}
	return orthoMembers[key];
}
 
STRINTMAP* 
MappedOrthogroup::getSpeciesHitsForGene(const char* srcSpecName, const char* targetSpecName, const char* geneName)
{
	string key(srcSpecName);
	if(orthoMembers.find(key)==orthoMembers.end())
	{
		cout <<"No hits for species  " << srcSpecName << endl;
		return NULL;
	}
	GeneMap* geneMap=orthoMembers[key];
	STRINTMAP* hits=geneMap->getHits(geneName,targetSpecName);
	return hits;
}

map<int,map<string,string>*>& 
MappedOrthogroup::getGeneSets()
{
	return geneSets;
}
*/
