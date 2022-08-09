/*
 *  scMTNI: single-cell Multi-Task learning Network Inference
 *   Copyright 2022 Shilu Zhang (szhang256@wisc.edu) and  Sushmita Roy (sroy@biostat.wisc.edu)
 *
 *   Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, 
 *   including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, 
 *   subject to the following conditions:
 *
 *   The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
 *
 *   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
 *   IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 *   */

#include <iostream>
#include <math.h>
#include "Evidence.H"

#include "Variable.H"
#include "Potential.H"

#include "gsl/gsl_randist.h"

Potential::Potential()
{
	mean=NULL;
	covariance=NULL;
	inverse=NULL;
	matrixInd=0;
	mbcondVar=0;
	mbcondMean_Part=0;
}

Potential::~Potential()
{
	varSet.clear();
	factorVariables.clear();
	markovBlnktVariables.clear();
	if(covariance!=NULL)
	{
		delete covariance;
	}
	if(inverse!=NULL)
	{
		delete inverse;
	}
	if(mean!=NULL)
	{
		delete mean;
	}
	meanWrap.clear();
	vIDMatIndMap.clear();
	matIndvIDMap.clear();
	mbcondMean_Vect.clear();
}

int
Potential::resetVarSet()
{
	matrixInd=0;
	varSet.clear();
	factorVariables.clear();
	markovBlnktVariables.clear();
	meanWrap.clear();
	vIDMatIndMap.clear();
	matIndvIDMap.clear();
	return 0;
}


int
Potential::initMatrix(int k)
{
	covariance=new Matrix(k,k);
	mean=new Matrix(k,1);
	return 0;
}

int 
Potential::setAssocVariable(Variable* var,Potential::VariableRole vRole)
{
    varSet[var->getID()]=var;
	vIDMatIndMap[var->getID()]=matrixInd;
	matIndvIDMap[matrixInd]=var->getID();
	matrixInd++;
	switch(vRole)
	{
		case Potential::FACTOR:
		{
			factorVariables[var->getID()]=0;
			break;
		}
		case Potential::MARKOV_BNKT:
		{
			markovBlnktVariables[var->getID()]=0;
			break;
		}
	}
	return 0;
}

unordered_map<int,Variable*>&
Potential::getAssocVariables()
{
	return varSet;
}

//The goal of this function is to simply change the variable on which we must condition
//on when computing the conditional potential
int
Potential::updateFactorVariable(int newfVar)
{
	if(factorVariables.size()!=1)
	{
		cout << "Too many factor variables " << endl;
		return -1;
	}
	int currfVar=factorVariables.begin()->first;
	if(currfVar==newfVar)
	{
		return 0;
	}
	Variable* fVar=varSet[currfVar];
	INTINTMAP_ITER vIter=markovBlnktVariables.find(newfVar);
	if(vIter==markovBlnktVariables.end())
	{
		cout <<"Did not find variable " << newfVar <<" in the mb set " << endl;
		return -1;
	}
	markovBlnktVariables.erase(vIter);
	markovBlnktVariables[currfVar]=0;
	factorVariables.clear();
	factorVariables[newfVar]=0;
	return 1;
}

// initialize meanWrap[INTDBLMAP], covariance[matrix], mean[1d matrix]
int
Potential::potZeroInit()
{
	for(auto vIter=varSet.begin();vIter!=varSet.end();vIter++)
	{
		meanWrap[vIter->first]=0;
	}
	int row=varSet.size();
	covariance=new Matrix(row,row);
	covariance->setAllValues(0);
	mean=new Matrix(row,1);
	mean->setAllValues(0);
	return 0;
}


int
Potential::potZeroInit_MeanOnly()
{
	for(auto vIter=varSet.begin();vIter!=varSet.end();vIter++)
	{
		meanWrap[vIter->first]=0;
	}
	return 0;
}

//Initialize the mean and variance using user specified values
int
Potential::potCustomInit(double mVal, double varVal)
{
	covariance=new Matrix(varSet.size(),varSet.size());
	mean=new Matrix(varSet.size(),1);
	for(auto vIter=varSet.begin();vIter!=varSet.end();vIter++)
	{
        meanWrap[vIter->first]=mVal;
        int ind=vIDMatIndMap[vIter->first];
        mean->setValue(mVal,ind,0);
	}
	for(auto vIter=varSet.begin();vIter!=varSet.end();vIter++)
	{
		int i=vIDMatIndMap[vIter->first];
		for(auto uIter=varSet.begin();uIter!=varSet.end();uIter++)
		{
			int j=vIDMatIndMap[uIter->first];
			covariance->setValue(varVal,i,j);
		}
	}
	return 0;
}

//Get the joint prob value for a particular configuration
double 
Potential::getJointPotValueFor(INTDBLMAP& varConf)
{
	string aKey;
	double pVal=0;
	Matrix* valMat=new Matrix(varSet.size(),1);
	for(INTDBLMAP_ITER idIter=varConf.begin();idIter!=varConf.end();idIter++)
	{
		int i=vIDMatIndMap[idIter->first];
		valMat->setValue(idIter->second,i,0);
	}
	Matrix* meanDiff=valMat->subtractMatrix(mean);
	Matrix* diffT=meanDiff->transMatrix();
	Matrix* p1=diffT->multiplyMatrix(inverse);
	Matrix* p2=p1->multiplyMatrix(meanDiff);
	double prod=p2->getValue(0,0);
	pVal=exp(-0.5*prod);
	pVal=pVal/normFactor;
	delete meanDiff;
	delete diffT;
	delete p1;
	delete p2;
	return pVal;
}


double 
Potential::getJointPotValueFor(STRINTMAP& varConf)
{
	cout <<"Not implemented " << endl;
	return 0;
}


double 
Potential::getJointPotValueForConf(string& varConf)
{
	cout <<"Not implemented " << endl;
	return 0;
}


int 
Potential::updateMean(int vID,double mVal)
{
	int mID=vIDMatIndMap[vID];
	meanWrap[vID]=mVal;
	mean->setValue(mVal,mID,0);
	return 0;
}

int 
Potential::updateCovariance(int vID,int uID,double sVal)
{
	int mID=vIDMatIndMap[vID];
	int nID=vIDMatIndMap[uID];
	covariance->setValue(sVal,mID,nID);
	return 0;
}


int
Potential::dumpPotential(ostream& oFile)
{
	oFile <<"Mean"<< endl;
	if(meanWrap.size()>0)
	{
		for(INTDBLMAP_ITER sIter=meanWrap.begin();sIter!=meanWrap.end();sIter++)
		{
			oFile <<varSet[sIter->first]->getName() << "\t" << sIter->second << endl;
		}
	}
	return 0;
}

int
Potential::showVariables(ostream& oFile)
{
	for(auto vIter=varSet.begin();vIter!=varSet.end();vIter++)
	{
		oFile <<" " << vIter->second->getName();
	}
	oFile << endl;
	return 0;
}

//The mean vector is ok. Now we need to get the invariance of the covariance and determinant.
//Also compute the normalization factor
int
Potential::makeValidJPD()
{
	inverse=covariance->invMatrix();
	//inverse->showMatrix(cout);
	determinant=covariance->detMatrix();
    //cout << "Potential::makeValidJPD() determinant=" << determinant <<endl;
	if(determinant <0)
	{
		cout <<"Negative Determinant " << determinant << endl;
	}
	double n=((double)varSet.size())/2.0;
	normFactor=pow(2*PI,n)*determinant;
	normFactor=sqrt(normFactor);
	return 0;
}


int
Potential::makeValidJPD(gsl_matrix* ludecomp, gsl_permutation* p)
{
        //cout << "Potential::makeValidJPD ludecomp=" << ludecomp->size1 << " col=" << ludecomp->size2 << " psize=" << p->size << " covariance=" << endl;
	if(inverse!=NULL)
	{
		delete inverse;
	}
        //covariance->showMatrix();
	inverse=covariance->invMatrix(ludecomp,p);
	determinant=covariance->detMatrix(ludecomp,p);
    //cout << "Potential::makeValidJPD(ludecomp,p) determinant=" << determinant << " ludecomp.size=" << ludecomp->size1 << " p.size=" << p->size <<endl;
        //cout << " determinant=" << determinant << endl;
	if(determinant <0)
	{
		cout <<"Negative Determinant " << determinant << endl;
		covariance->showMatrix();
		for(INTINTMAP_ITER vIter=matIndvIDMap.begin();vIter!=matIndvIDMap.end();vIter++)
		{
			int vId=vIter->second;
			cout << vIter->first << " " << vId << " "<< varSet[vId]->getName() << endl;
		}
	}
	double n=((double)varSet.size())/2.0;
	normFactor=pow(2*PI,n)*determinant;
	normFactor=sqrt(normFactor);
	return 0;
}

//This is the conditional entropy and is calculated as H(FVar|MarkovBlnkVar)
int
Potential::calculateEntropy()
{
	double entropy=0;
	double jtEntropy=0;
	//Need to create a smaller version of the covariance matrix discarding the row and columns for the
	//factorVariable
	Matrix* neighbourCov=new Matrix(markovBlnktVariables.size(),markovBlnktVariables.size());
	int fmind=vIDMatIndMap[factorVariables.begin()->first];
	for(INTINTMAP_ITER vIter=markovBlnktVariables.begin();vIter!=markovBlnktVariables.end();vIter++)
	{
		int i=vIDMatIndMap[vIter->first];
		for(INTINTMAP_ITER uIter=markovBlnktVariables.begin();uIter!=markovBlnktVariables.end();uIter++)
		{
			int j=vIDMatIndMap[uIter->first];
			double cv=covariance->getValue(i,j);
			if((i==fmind) || (j==fmind))
			{
				cout <<"Factor and Markov Blanket variables are the same!!" << endl;
				return -1;
			}
			if(i>fmind)
			{
				i--;
			}
			if(j>fmind)
			{
				j--;
			}
			neighbourCov->setValue(cv,i,j);
		}
	}
	double nDet=neighbourCov->detMatrix();
	double commFact=1+log(2*PI);
	double p=((double)markovBlnktVariables.size());
	double mbEntropy=0.5*((p*commFact) + log(nDet));
	double n=((double)varSet.size());
	jointEntropy=0.5*((n*commFact) + log(determinant));
	potEntropy=jointEntropy-mbEntropy;
	return 0;
}


int
Potential::calculateJointEntropy()
{
	double commFact=1+log(2*PI);
	double n=((double)varSet.size());
	jointEntropy=0.5*((n*commFact) + log(determinant));
	return 0;
}


double
Potential::getEntropy()
{
	return potEntropy;
}

double
Potential::getJointEntropy()
{
	return jointEntropy;
}


double
Potential::generateSample(INTDBLMAP& jointConf, int vId,gsl_rng* r)
{
	if(jointConf.find(factorVariables.begin()->first)==jointConf.end())
	{
		cout <<"Fatal error! No variable assignment for " << factorVariables.begin()->first << endl;
		exit(0);
	}
	double newmean=0;
	for(auto aIter=mbcondMean_Vect.begin();aIter!=mbcondMean_Vect.end();aIter++)
	{
		if(jointConf.find(aIter->first)==jointConf.end())
		{
			cout <<"Fatal error! No variable assignment for " << aIter->first << endl;
			exit(0);
		}
		double aval=jointConf[aIter->first];
		newmean=newmean+(aval*aIter->second);
	}
	newmean=newmean+mbcondMean_Part;
	double x=gsl_ran_gaussian(r,sqrt(mbcondVar));
	x=x+newmean;
	return x;
}


int 
Potential::copyMe(Potential** apot)
{
	*apot=new Potential;
	for(INTINTMAP_ITER vIter=matIndvIDMap.begin();vIter!=matIndvIDMap.end();vIter++)
	{
		Variable* v=varSet[vIter->second];
		if(factorVariables.find(vIter->second)!=factorVariables.end())
		{
			(*apot)->setAssocVariable(v,Potential::FACTOR);
		}
		else
		{
			(*apot)->setAssocVariable(v,Potential::MARKOV_BNKT);
		}
	}
	
	for(INTDBLMAP_ITER sdIter=meanWrap.begin();sdIter!=meanWrap.end();sdIter++)
	{
		(*apot)->updateMean(sdIter->first,sdIter->second);
	}
	for(INTDBLMAP_ITER uIter=meanWrap.begin();uIter!=meanWrap.end();uIter++)
	{
		int i=vIDMatIndMap[uIter->first];
		for(INTDBLMAP_ITER vIter=meanWrap.begin();vIter!=meanWrap.end();vIter++)
		{
			int j=vIDMatIndMap[vIter->first];
			double cv=covariance->getValue(i,j);
			(*apot)->updateCovariance(i,j,cv);
		}
	}
	return 0;
}

// compute mbcondVar[CondVariance], mbcondMean_Part[CondBias], mbcondMean_Vect[CondWeight]
int 
Potential::initMBCovMean()
{
    int vId=factorVariables.begin()->first;
    int vIdmId=vIDMatIndMap[vId];
    mbcondVar=covariance->getValue(vIdmId,vIdmId);
    mbcondMean_Part=meanWrap[vId];
    //cout << "My target ID is vId=" << vId << " vIdmId=" << vIdmId << " mbcondVar=" << mbcondVar << " mbcondMean_Part="  << mbcondMean_Part <<" meanWrap.size=" << meanWrap.size() <<" covariance="  << endl;
    //covariance->showMatrix(cout);
    if(markovBlnktVariables.size()==0)
    {
        return 0;
    }
    Matrix* mbcov=new Matrix(markovBlnktVariables.size(),markovBlnktVariables.size());
    Matrix* mbmargvar=new Matrix(1,markovBlnktVariables.size());
    INTINTMAP localMatIDMap;
    for(INTDBLMAP_ITER uIter=meanWrap.begin();uIter!=meanWrap.end();uIter++)
    {
        int i=vIDMatIndMap[uIter->first];
        int inew=i;
        if(i>vIdmId)
        {
            inew--;
        }
        for(INTDBLMAP_ITER vIter=meanWrap.begin();vIter!=meanWrap.end();vIter++)
        {
            if(vIter->first==vId)
            {
                continue;
            }
            int j=vIDMatIndMap[vIter->first];
            double cv=covariance->getValue(i,j);
            if(j>vIdmId)
            {
                j--;
            }
            if(uIter->first==vId)
            {
                mbmargvar->setValue(cv,0,j);
            }
            else
            {
                mbcov->setValue(cv,inew,j);
            }
        }
        if(uIter->first!=vId)
        {
            localMatIDMap[inew]=uIter->first;
        }
    }
    //cout <<"mbcov:" <<endl;
    //mbcov->showMatrix(cout);
    //cout <<"mbmargvar:" <<endl;
    //mbmargvar->showMatrix(cout);
    Matrix* covInv=mbcov->invMatrix();
    Matrix* prod1=mbmargvar->multiplyMatrix(covInv);
    //cout << "Potential::initMBCovMean get my conditional weight: " ;
    for(INTINTMAP_ITER aIter=localMatIDMap.begin();aIter!=localMatIDMap.end();aIter++)
    {
        double aVal=prod1->getValue(0,aIter->first);
        double bVal=mbmargvar->getValue(0,aIter->first);
        mbcondVar=mbcondVar-(aVal*bVal);
        mbcondMean_Vect[aIter->second]=aVal;
        double cVal=meanWrap[aIter->second];
        mbcondMean_Part=mbcondMean_Part-(cVal*aVal);
        //cout << " aIter->first=" << aIter->first<<" mbcondMean_Vect[" << aIter->second << "]=" << aVal << " ";
    }
    //cout <<endl;
    localMatIDMap.clear();
    if(mbcondVar<1e-10)
    {
        //mbcondVar=1e-10;
    }
    delete mbcov;
    delete mbmargvar;
    delete covInv;
    delete prod1;
    
    return 0;
}

// shilu added: return log(pval)
double
Potential::getCondPotValueFor(double x)
{
    double norm=0.5*log(2*PI*mbcondVar);
    double dev=-1.0*(x-mbcondMean_Part)*(x-mbcondMean_Part)/(2*mbcondVar);  //newmean=mbcondMean_Part
    double pval=dev-norm;
    return pval;
}

double 
Potential::getCondPotValueFor(INTDBLMAP& assignment)
{
	/*if(assignment.find(factorVariables.begin()->first)==assignment.end())
	{
		cout <<"Fatal error! No variable assignment for " << factorVariables.begin()->first << endl;
		exit(0);
	}*/
	double newmean=0;
    //mbcondMean_Vect.size=0
	for(auto aIter=mbcondMean_Vect.begin();aIter!=mbcondMean_Vect.end();aIter++)
	{
		/*if(assignment.find(aIter->first)==assignment.end())
		{
			cout <<"Fatal error! No variable assignment for " << aIter->first << endl;
			exit(0);
		}*/
		double aval=assignment[aIter->first];
		newmean=newmean+(aval*aIter->second);
	}
    //cout << "mbcondMean_Vect.size="<<mbcondMean_Vect.size() << " mbcondMean_Part="<< mbcondMean_Part << endl;
	newmean=newmean+mbcondMean_Part;
	//double normsq=2*PI*mbcondVar;
	double norm=sqrt(2*PI*mbcondVar);
	//norm=sqrt(normsq);
	double x=assignment[factorVariables.begin()->first];
	double dev=(x-newmean)*(x-newmean);
	dev=dev/(2*mbcondVar);
	double eval=exp(-1.0*dev);
	double pval=eval/norm;
	return pval;

}



double 
Potential::getCondPotValueFor(map<int,Evidence*>* evidMap)
{
	if(evidMap->find(factorVariables.begin()->first)==evidMap->end())
	{
		cout <<"Fatal error! No variable assignment for " << factorVariables.begin()->first << endl;
		exit(0);
	}
	double newmean=0;
	for(auto aIter=mbcondMean_Vect.begin();aIter!=mbcondMean_Vect.end();aIter++)
	{
		if(evidMap->find(aIter->first)==evidMap->end())
		{
			cout <<"Fatal error! No variable assignment for " << aIter->first << endl;
			exit(0);
		}
		Evidence* evid=(*evidMap)[aIter->first];
		double aval=evid->getEvidVal();
		newmean=newmean+(aval*aIter->second);
	}
	newmean=newmean+mbcondMean_Part;
	double normsq=2*PI*mbcondVar;
	double norm=sqrt(2*PI*mbcondVar);
	norm=sqrt(normsq);
	Evidence* fevid=(*evidMap)[factorVariables.begin()->first];
	double x=fevid->getEvidVal();
	double dev=(x-newmean)*(x-newmean);
	dev=dev/(2*mbcondVar);
	double eval=exp(-1.0*dev);
	double pval=eval/norm;
	return pval;
}

// added by vperiyasamy
double
Potential::computeLL_Tracetrick(int sampleSize)
{
	//double ll=0;
	int dim=covariance->getRowCnt();
	//double constant=sampleSize*dim*log(2*PI);
	double ll=((double) sampleSize) * (dim*log(2*PI)+log(determinant));
	/*Matrix sos=Matrix(vIDMatIndMap.size(),vIDMatIndMap.size());
	for(int i=0;i<dim;i++)
	{
		for(int j=0;j<dim;j++)
		{
			sos.setValue(covariance->getValue(i,j)*(sampleSize-1),i,j);
		}
	}
	//sos.subtractScalar(0.001);
	Matrix* m=inverse->multiplyMatrix(&sos);
	double t=0;
	for(int i=0;i<dim;i++)
	{
		t=t+m->getValue(i,i);
	}  *///commented by shilu
    double t=dim*(sampleSize-1);
	ll=(ll+t)*(-0.5);  //(ll+t+constant)*(-0.5);
    //cout << "computeLL_Tracetrick dim=" << dim << " sampleSize=" << sampleSize<<" PI=" << PI <<" determinant=" << determinant <<" log(determinant)=" << log(determinant) <<" ll=" << ll <<endl;
	//delete m;
	return ll;
}


double
Potential::getCondVariance()
{
	return mbcondVar;
}

double
Potential::getCondBias()
{
	return mbcondMean_Part;
}

unordered_map<int,double>&
Potential::getCondWeight()
{
	return mbcondMean_Vect;
}

double
Potential::getDeterm()
{
    return determinant;
}
