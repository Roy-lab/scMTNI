#include <map>
#include <string>
#include <iostream>
#include <fstream>
using namespace std;
#include "GeneMap.H" 
#include "MappedOrthogroup.H"
#include "MappedOrthogroupReader.H"

int
getOrthologousNetwork(MappedOrthogroupReader* mor, const char*,const char*, const char*);

int
main(int argc, const char** argv)
{
	if(argc!=6)
	{
		cout <<"Usage: testOGR orthogroupmapping specidmap network srcname targetname " << endl;
		return 0;
	}
	MappedOrthogroupReader mor;
	mor.readSpeciesMapping(argv[1]);
	mor.readFile(argv[2]);
	getOrthologousNetwork(&mor,argv[3],argv[4],argv[5]);
	return 0;
}

int
getOrthologousNetwork(MappedOrthogroupReader* mor,const char* aFName, const char* srcSpecName, const char* targetSpecName)
{
	ifstream inFile(aFName);
	char buffer[1024];
	int total=0;
	int hit=0;
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		string tfname;
		string targetname;
		while(tok!=NULL)
		{
			if(tokCnt==0)
			{
				tfname.append(tok);
			}
			else
			{
				targetname.append(tok);
			}
			tok=strtok(NULL,"\t");
			tokCnt++;
		}
		if(strcmp(tfname.c_str(),"MATA1")==0)
		{
			continue;
		}
		total++;
		STRINTMAP* tforthos=mor->getOrtholog(srcSpecName,tfname.c_str(),targetSpecName);
		STRINTMAP* tgtorthos=mor->getOrtholog(srcSpecName,targetname.c_str(),targetSpecName);
		if((tforthos==NULL)||(tgtorthos==NULL))
		{
			continue;
		}
		hit++;
		for(STRINTMAP_ITER tfIter=tforthos->begin();tfIter!=tforthos->end();tfIter++)
		{
			for(STRINTMAP_ITER tgtIter=tgtorthos->begin();tgtIter!=tgtorthos->end();tgtIter++)
			{
			//	cout << tfIter->first <<"\t" <<tgtIter->first << endl;
			}
		}
	}
	inFile.close();
	cout <<"Found " << hit << " out of " << total << " interactions " << endl;
	return 0;
}
