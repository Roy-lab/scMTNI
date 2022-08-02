#include <fstream>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "SpeciesDistance.H"

SpeciesDistance::SpeciesDistance()
{
	root=NULL;
	proot=0.5;
}

SpeciesDistance::~SpeciesDistance()
{
    for(auto pIter=speciesSet.begin();pIter!=speciesSet.end();pIter++)
    {
        delete pIter->second;
    }
    speciesSet.clear();
    binEdgeKeyProbMap.clear();
    speciesNameIDMap.clear();
}

int
SpeciesDistance::setProot(double aVal)
{
    proot=aVal;
    return 0;
}


int
SpeciesDistance::setspeciesNameIDMap(unordered_map<string,int>& spNameIDMap){
    speciesNameIDMap=spNameIDMap;
}

int 
SpeciesDistance::readSpeciesTree(const char* aFName)
{
	ifstream inFile(aFName);
	char buffer[1024];
	//double p_gain = pgain;
	//double p_maintain_edge = pmaintain;
	while(inFile.good())
	{
		inFile.getline(buffer,1023);
		if(strlen(buffer)<=0)
		{
			continue;
		}
		if(strchr(buffer,'#')!=NULL)
		{
			continue;
		}
		char* tok=strtok(buffer,"\t");
		int tokCnt=0;
		string childSpeciesName;
		string parentSpeciesName;
		//string childtype;
		double p_gain=0;
		double p_maintain_edge=0;
        while(tok!=NULL)
        {
            if(tokCnt==0)
            {
                childSpeciesName.append(tok);
            }
            else if(tokCnt==1)
            {
                parentSpeciesName.append(tok);
            }
            else if(tokCnt==2)
            {
                p_gain=atof(tok);
            }
            else if(tokCnt==3)
            {
                p_maintain_edge=1-atof(tok);
            }
            tok=strtok(NULL,"\t");
            tokCnt++;
        }
		SpeciesDistance::Species* childSpecies=NULL;
		SpeciesDistance::Species* parentSpecies=NULL;
		if(speciesSet.find(childSpeciesName)==speciesSet.end())
		{
			childSpecies=new SpeciesDistance::Species;
			childSpecies->name.append(childSpeciesName.c_str());
			childSpecies->parent=NULL;
			speciesSet[childSpeciesName]=childSpecies;
		}
		else
		{
			childSpecies=speciesSet[childSpeciesName];
		}
		if(speciesSet.find(parentSpeciesName)==speciesSet.end())
		{
			parentSpecies=new SpeciesDistance::Species;
			parentSpecies->name.append(parentSpeciesName.c_str());
			parentSpecies->parent=NULL;
			speciesSet[parentSpeciesName]=parentSpecies;
		}
		else
		{
			parentSpecies=speciesSet[parentSpeciesName];
		}
		childSpecies->parent=parentSpecies;
		if(parentSpecies->parent==NULL)
		{
			root=parentSpecies;
		}
		//Probabilities are for the child node
		childSpecies->p_maintain_edge=p_maintain_edge;
		childSpecies->p_gain=p_gain;
        parentSpecies->children.push_back(childSpecies);
        cout << "("<<childSpeciesName << "|" << parentSpeciesName << ") p_maintain_edge=" << p_maintain_edge << " p_gain=" << p_gain << endl;
	}
        cout <<"speciesSet.size()=" <<speciesSet.size() << " Root is " << root->name << endl;
	inFile.close();
	return 0;
}


//The probability that an edge is maintained in a child given that the edge is present in the ancestor
double 
SpeciesDistance::getProbGainGain(string& spName)
{
	/*if(speciesSet.find(spName)==speciesSet.end())
	{
		cout <<"No species with name " << spName.c_str() << endl;
		exit(0);
	}*/
	Species* species=speciesSet[spName];
	double gain_gain=species->p_maintain_edge;
	return gain_gain;
}

//The probability that an edge is gained in a child given that the edge is absent in the ancestor
double 
SpeciesDistance::getProbGainLoss(string& spName)
{
	/*if(speciesSet.find(spName)==speciesSet.end())
	{
		cout <<"No species with name " << spName.c_str() << endl;
		exit(0);
	}*/
	Species* species=speciesSet[spName];
	double gain_loss=species->p_gain;
	return gain_loss;
}

//The probability than an edge is lost in a child given that the edge was present in the ancestor
double 
SpeciesDistance::getProbLossGain(string& spName)
{
	/*if(speciesSet.find(spName)==speciesSet.end())
	{
		cout <<"No species with name " << spName.c_str() << endl;
		exit(0);
	}*/
	Species* species=speciesSet[spName];
	double loss_gain=1-species->p_maintain_edge;
	return loss_gain;
}

//The  probability that an edge is lost in a child given the edge was absent in the ancestor
double 
SpeciesDistance::getProbLossLoss(string& spName)
{
	/*if(speciesSet.find(spName)==speciesSet.end())
	{
		cout <<"No species with name " << spName.c_str() << endl;
		exit(0);
	}*/
	Species* species=speciesSet[spName];
	double loss_loss=1-species->p_gain;
	return loss_loss;
}

// new getEdgeStatusProb by Shilu
double
SpeciesDistance::getEdgeStatusProb(vector<int>& edgeStatus)
{
    /*cout << "Edge status: ";
    for(map<string,int>::iterator eIter=edgeStatus.begin();eIter!=edgeStatus.end();eIter++)
    {
        cout <<eIter->first <<"=" << eIter->second <<"| " ;
    }*/
    double edgeprior=0;
    int binkey=0;
    int binpos=0;
    for(int eIter=0;eIter<edgeStatus.size();eIter++)
    {
        if(edgeStatus[eIter]!=0)
        {
            int val=(int)pow(2,binpos);
            binkey=binkey+val;
        }
        binpos=binpos+1;
    }
    if(binEdgeKeyProbMap.find(binkey)!=binEdgeKeyProbMap.end())
    {
        return binEdgeKeyProbMap[binkey];
    }

    /*cout << "binkey=" << binkey << " [Edge status: ";
    for(int eIter=0;eIter<edgeStatus.size();eIter++)
    {
        cout <<edgeStatus[eIter] <<" ";
    }
    cout << "] " ;*/
    //Start with the root
    //changed 6/4 to test setting the prior probability of an edge in the root node = .2
    double score;
    if(edgeStatus[0]==0)
    {
        score = 1-proot;
    }
    else
    {
        score = proot;
    }
    edgeprior=score;
    //cout << "P(G)=P(" << root->name <<"=" <<edgeStatus[0] <<")[" << score <<"]";
    for(int i=0;i<root->children.size();i++)
    {
        double childrenScore=getSubTree(edgeStatus[0],root->children[i],edgeStatus);
        edgeprior=edgeprior*childrenScore;
        //cout << " * P(" << root->children[i]->name <<"=" <<edgeStatus[speciesNameIDMap[root->children[i]->name]] <<"|" << root->name <<"=" <<edgeStatus[0] << ")[" <<  childrenScore << "]";
    }
    //cout << "=" << edgeprior << endl;
    binEdgeKeyProbMap[binkey]=edgeprior;
    return edgeprior;
}


// new getSubTree function by Shilu
double
SpeciesDistance::getSubTree(bool parentEdge, Species* child, vector<int>& edgeStatus)
{
    double score=0;
    int status=edgeStatus[speciesNameIDMap[child->name]];

    if(child->children.empty())
    {
        //This is a leaf node
        if(parentEdge)
        {
            if(status==0)
            {
                score=getProbLossGain(child->name);
            }
            else
            {
                score=getProbGainGain(child->name);
            }
        }
        else
        {
            if(status==0)
            {
                score=getProbLossLoss(child->name);
            }
            else
            {
                score=getProbGainLoss(child->name);
            }
        }
        //cout << "&&&&Reach LeafNode! P(" << child->name << "=" << status <<"|" <<child->parent->name << "=" <<  parentEdge << ")= " << score <<endl;
    }
    else
    {
        //This is a parent node

        // its own score:
        if(edgeStatus[speciesNameIDMap[child->parent->name]])
        {
            if(status==0)
            {
                score=getProbLossGain(child->name);
            }
            else
            {
                score=getProbGainGain(child->name);
            }
        }
        else
        {
            if(status==0)
            {
                score=getProbLossLoss(child->name);
            }
            else
            {
                score=getProbGainLoss(child->name);
            }
        }
        //cout << "Current " << child->name << " score=" << score << endl;
        //cout << " * P(" << child->name <<"=" <<edgeStatus[child->name] <<"|" << child->parent->name <<"=" <<edgeStatus[child->parent->name] << ")" ;
        // children score
        for(int i=0;i<child->children.size();i++)
        {
            // edgeStatus[child->name] status for current node, i.e. parent edge status for children
            double childrenScore=getSubTree(edgeStatus[speciesNameIDMap[child->name]],child->children[i],edgeStatus);
            score=score*childrenScore;
            //cout << " * P(" << child->children[i]->name <<"=" <<edgeStatus[speciesNameIDMap[child->children[i]->name]] <<"|" << child->name <<"=" <<edgeStatus[speciesNameIDMap[child->name]] << ")" ;
        }
    }
    return score;
}

SpeciesDistance::Species*
SpeciesDistance::getRoot()
{
	return root;
}
