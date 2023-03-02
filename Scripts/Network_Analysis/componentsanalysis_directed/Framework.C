#include <iostream>
#include <stdlib.h>
using namespace std;
#include "Vertex.H"
#include "Clique.H"
#include "Graph.H"
#include "Framework.H"

Framework::Framework()
{
}

Framework::~Framework()
{
}

int
Framework::init(const char* inFName)
{
	g.makeGraph(inFName);
	return 0;
}

int
Framework::start(const char* outFName, int hubthreshold)
{
	g.obtainConnectivity();
	if(g.isConnected())
	{
		cout <<"Graph is connected" << endl;
	}
	else
	{
		cout <<"Graph is not connected" << endl;
	}
	vector<map<string,Vertex*>*> components;
	g.getConnectedComponents(components);
	char allCompFName[1024];
	sprintf(allCompFName,"%s_allcomps.txt",outFName);
	ofstream oFile(allCompFName);
	//oFile <<"<Root>" << endl;
	int gcid=-1;
	int maxsize=0;
	//Show only giancomponent
	for(int c=0;c<components.size();c++)
	{
		//char compFName[1024];
		//sprintf(compFName,"%s/comp%d.sif",outFName,c);
		//ofstream cFile(compFName);
		map<string,Vertex*>* acomponent=components[c];
		if(maxsize<acomponent->size())
		{
			maxsize=acomponent->size();
			gcid=c;
		}
		for(map<string,Vertex*>::iterator mIter=acomponent->begin();mIter!=acomponent->end();mIter++)
		{
			if(mIter!=acomponent->begin())
			{
				oFile <<"\t";
			}
			oFile << mIter->first;
		} 
		oFile<<  endl;
	}
	char compFName[1024];
	sprintf(compFName,"%s_giantcomponent.txt",outFName);
	ofstream cFile(compFName);
	cFile<<g.getEdgeHeader() <<endl;
	map<string,Vertex*>* acomponent=components[gcid];
		for(map<string,Vertex*>::iterator vIter=acomponent->begin();vIter!=acomponent->end();vIter++)
		{
			map<string,Vertex*>& neighborSet=vIter->second->getImmediateNeighbours();
			//show only those neighbors associated with at least hubthreshold edges
			if(neighborSet.size()<hubthreshold)
			{
				continue;
			}
			for(map<string,Vertex*>::iterator nIter=neighborSet.begin();nIter!=neighborSet.end();nIter++)
			{
				char buff[1024];
				sprintf(buff,"%s\t%s",vIter->first.c_str(),nIter->first.c_str());
				string eKey;
				eKey.append(buff);
				vector<string>* attribs=g.getEdgeAttribs(eKey);
				if(attribs!=NULL)
				{
					for(int i=0;i<attribs->size();i++)
					{
						cFile << eKey.c_str()<< "\t"<< (*attribs)[i].c_str() << endl;
					}
				}
				sprintf(buff,"%s\t%s",nIter->first.c_str(),vIter->first.c_str());
				eKey.append(buff);
				attribs=g.getEdgeAttribs(eKey);
				if(attribs!=NULL)
				{
					for(int i=0;i<attribs->size();i++)
					{
						cFile << eKey.c_str()<< "\t"<< (*attribs)[i].c_str() << endl;
					}
				}
				
			}
			
		} 
		cFile.close();
		//oFile <<"<Samples></Samples>" << endl;
		//oFile << "</Module>" << endl;
	cout <<"Giant component is "<< gcid << " with " << maxsize << " nodes" << endl;
		
	oFile.close();
	return 0;
}

int
main(int argc, const char** argv)
{
	if(argc!=4)
	{
		cout <<"Usage: getConnectedComponents infile outfile hubnodefilter" << endl;
		return 0;
	}
	Framework fw;
	fw.init(argv[1]);
	//fw.start(argv[2]);
	fw.start(argv[2],atoi(argv[3]));
	return 0;
}
