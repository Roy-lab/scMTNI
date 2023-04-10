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
#include "MetaMove.H"
MetaMove::MetaMove()
{
}

MetaMove::~MetaMove()
{
	//conditionSet.clear();
	srcWeight.clear();
}

int 
MetaMove::setScoreImprovement(double aVal)
{
	scoreDelta=aVal;
	return 0;
}

int
MetaMove::setSrcMBScore(double aScore)
{
	mbscore=aScore;
	return 0;
}

int
MetaMove::setTargetMBScore(double aScore)
{
	targetMBScore=aScore;
	return 0;
}


int 
MetaMove::setSrcVertex(int vid)
{
	src=vid;
	return 0;
}

int
MetaMove::setTargetVertex(int vid)
{	
	target=vid;
	return 0;
}

int
MetaMove::setTFID(int vid)
{
    regi=vid;
    return 0;
}

int
MetaMove::setTargetID(int vid)
{
    targeti=vid;
    return 0;
}

/*int
MetaMove::setConditionSet(INTINTMAP& vSet)
{
	for(INTINTMAP_ITER vIter=vSet.begin();vIter!=vSet.end();vIter++)
	{
		conditionSet[vIter->first]=vIter->second;
	}

	return 0;
}*/

int
MetaMove::setConditionSetInd(int aind)
{
	conditionSetInd=aind;
	return 0;
}

int 
MetaMove::setSrcWeight(unordered_map<int,double>& awt)
{
	for(auto wIter=awt.begin();wIter!=awt.end();wIter++)
	{
		srcWeight[wIter->first]=wIter->second;
	}
	return 0;
}

int 
MetaMove::getSrcVertex()
{
	return src;
}

int
MetaMove::getTargetVertex()
{
	return target;
}

int
MetaMove::getTFID()
{
    return regi;
}

int
MetaMove::getTargetID()
{
    return targeti;
}

/*INTINTMAP&
MetaMove::getConditionSet()
{
	return conditionSet;
}*/

int
MetaMove::getConditionSetInd()
{	
	return conditionSetInd;
}

double 
MetaMove::getSrcMBScore()
{
	return mbscore;
}

double
MetaMove::getTargetMBScore()
{
	return targetMBScore;
}


double 
MetaMove::getScoreImprovement()
{
	return scoreDelta;
}

unordered_map<int,double>&
MetaMove::getSrcWeight()
{
	return srcWeight;
}

