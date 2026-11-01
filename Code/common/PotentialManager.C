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

#include <cassert>
#include <iostream>
#include <cstring>
#include <math.h>
#include <fstream>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "CommonTypes.H"
#include "Error.H"
#include "Variable.H"
#include "SlimFactor.H"
#include "PotentialManager.H"
#include <chrono>
#include "Utils.H"
#include "VariableManager.H"

using namespace std::chrono;

PotentialManager::PotentialManager()
{
    meanMat = NULL;
    covMat = NULL;
    sampleCount = 0;
}

PotentialManager::~PotentialManager()
{
	if (meanMat != NULL) {
        delete meanMat;
    }
    if (covMat != NULL) {
        delete covMat;
    }
}

void
PotentialManager::loadEvidenceFromTable(string fileName, const vector<string>& variableList)
{
    // Reads gene names and expression levels from a tab-separated file, where the first row is assumed (for now)
    // to be headers. In subsequent rows, the first column is a gene name and remaining columns are expression levels.

    if (fileName.find('.') == std::string::npos) {
        fileName = fileName + ".table";
        std::cout << "no .table in input expression filename, add .table suffix:" << fileName << endl;
    }

    cout << "readEvidenceTable:" << fileName << endl;
    ifstream inFile(fileName);
    assert(inFile);

    bool headerLine = true;
    vector<string> inputTable;
    string inputLine;

    while (getline(inFile, inputLine)) {
        if (inputLine.empty() || inputLine.find('#') == 0 || headerLine) {
            headerLine = false;
            continue;
        }
        inputTable.push_back(inputLine);
    }

    int varCount = variableList.size();

    // Count the columns and subtract 1 (for the gene name column) to get sample count.
    vector<string> substrs = Utils::split(inputTable[0], '\t');
    sampleCount = substrs.size() - 1;

    vMgr = new VariableManager;
    Matrix data(varCount, sampleCount);

    // Load sample data into matrix, and variables into VariableManager

    for (auto line : inputTable) {
        substrs = Utils::split(line, '\t');
        string varName = substrs[0];

        vector<string>::const_iterator iter = find(variableList.begin(), variableList.end(), varName);

        // Variable list restricts which genes we include, so if it isn't present, skip it.
        if (iter == variableList.end()) {
            continue;
        }

        int varIndex = distance(variableList.begin(), iter);
        vMgr->addVariable(varName, varIndex);

        for (int sample = 0; sample < sampleCount; sample++) {
            double varVal = stod(substrs[sample + 1]);
            data.setValue(varVal, varIndex, sample);
        }
    }

    // Store means and covariances.

    meanMat = new Matrix(varCount, 1);
    meanMat->setAllValues(0);

    covMat = new Matrix(varCount, varCount);
    covMat->setAllValues(-1);

    for (int i = 0; i < varCount; i++) {
        double vMean = data.RowMean(i);
        meanMat->setValue(vMean, i, 0);

        double ssd = data.vectorMultiply(i, vMean, i, vMean);

        //add 1e-10 to avoid singularity issues:
        double variance = ssd / (sampleCount - 1.0) + 1e-10;
        covMat->setValue(variance, i, i);

        for(int j = i + 1; j < varCount; j++) {
            double uMean = meanMat->getValue(j, 0);
            double ssd = data.vectorMultiply(i, vMean, j, uMean);
            double var = ssd / (sampleCount - 1.0);
            covMat->setValue(var, j, i);
            covMat->setValue(var, i, j);
        }
    }

    cout << "PotentialManager::loadEvidenceFromTable varCount=" << varCount << " sampleCount=" << sampleCount << endl;
}

VariableManager*
PotentialManager::getVariableManager()
{
    return vMgr;
}

double
PotentialManager::computeConditionalLL(SlimFactor* sFactor, int& status)
{
    int mbsize = sFactor->mergedMB.size();
    if (mbsize == 0) {
        return 0;
    }

    vector<int> mbVars(sFactor->mergedMB.begin(), sFactor->mergedMB.end());
    Matrix* covariance = createMBCovarianceMatrix(sFactor, mbVars);

    double determinant = covariance->detMatrix();
    if (determinant <= 0) {
        status = -1;
        return 0;
    }

    int vId = sFactor->fId;
    double mbcondVar = covMat->getValue(vId, vId);
    
    Matrix* mbcov = covariance->getSubMatrix(0, 0, mbsize, mbsize);
    Matrix* mbmargvar = covariance->getSubMatrix(mbsize, 0, 1, mbsize);

    double determinantP = mbcov->detMatrix();
    if (determinantP <= 0) {
        status = -1;
        return 0;
    }

    Matrix* covInv = mbcov->invMatrix();
    Matrix* prod1 = mbmargvar->multiplyMatrix(covInv);

    for (int i = 0; i < mbVars.size(); i++) {
        int varID = mbVars[i];
        double aVal = prod1->getValue(0, i);
        double bVal = mbmargvar->getValue(0, i);
        mbcondVar = mbcondVar - (aVal * bVal);
    }

    if (mbcondVar < 0) {
        status=-1;
        return 0;
    }

    double jointll1 = computeJointLL(mbsize + 1, determinant);
    double jointll2 = computeJointLL(mbsize, determinantP);
    double pll = jointll1 - jointll2;

    if (isinf(pll)) {
        status=-1;
        return 0;
    }

    delete covariance;
    delete mbcov;
    delete mbmargvar;
    delete covInv;
    delete prod1;

    return pll;
}

void
PotentialManager::dumpVarMB(SlimFactor* sFactor, ofstream& oFile)
{
    int mbsize = sFactor->mergedMB.size();
    if (mbsize == 0) {
        return;
    }

    vector<int> mbVars(sFactor->mergedMB.begin(), sFactor->mergedMB.end());

    Matrix* covariance = createMBCovarianceMatrix(sFactor, mbVars);
    Matrix* mbcov = covariance->getSubMatrix(0, 0, mbsize, mbsize);
    Matrix* mbmargvar = covariance->getSubMatrix(mbsize, 0, 1, mbsize);
    Matrix* covInv = mbcov->invMatrix();
    Matrix* prod1 = mbmargvar->multiplyMatrix(covInv);

    Variable* target = vMgr->getVariable(sFactor->fId);

    for(int i = 0; i < mbVars.size(); i++) {
        Variable* parent = vMgr->getVariable(mbVars[i]);
        double conditionalWeight = prod1->getValue(0, i);
        oFile << parent->getName() << "\t"<< target->getName() << "\t" << conditionalWeight << endl;
    }

    delete covariance;
    delete mbcov;
    delete mbmargvar;
    delete covInv;
    delete prod1;
}

void
PotentialManager::computePotentialMBCovMean(SlimFactor* sFactor, double& mbcondVar, double& mbbias, unordered_map<int,double>& mbcondMean_Vect)
{
    int mbsize = sFactor->mergedMB.size();
    if (mbsize == 0) {
        return;
    }

    vector<int> mbVars(sFactor->mergedMB.begin(), sFactor->mergedMB.end());

    Matrix* covariance = createMBCovarianceMatrix(sFactor, mbVars);

    int vId = sFactor->fId;
    mbcondVar = covMat->getValue(vId, vId);
    mbbias = meanMat->getValue(vId, 0);

    Matrix* mbcov = covariance->getSubMatrix(0, 0, mbsize, mbsize);
    Matrix* mbmargvar = covariance->getSubMatrix(mbsize, 0, 1, mbsize);
    Matrix* covInv = mbcov->invMatrix();
    Matrix* prod1 = mbmargvar->multiplyMatrix(covInv);

    for(int i = 0; i < mbVars.size(); i++) {
        int varID = mbVars[i];
        double conditionalWeight = prod1->getValue(0, i);
        double bVal = mbmargvar->getValue(0, i);
        mbcondVar = mbcondVar - (conditionalWeight * bVal);
        mbcondMean_Vect[varID] = conditionalWeight;
        double cVal = meanMat->getValue(varID, 0);
        mbbias = mbbias - (cVal * conditionalWeight);
    }

    delete covariance;
    delete mbcov;
    delete mbmargvar;
    delete covInv;
    delete prod1;
}

double
PotentialManager::computeJointLL(int dim, double determinant)
{
    double ll = sampleCount * (dim * log(2 * M_PI) + log(determinant));
    double t = dim * (sampleCount - 1);
    ll = -0.5 * (ll + t);
    return ll;
}

double
PotentialManager::computeUnivariateLL(int id)
{
    double mean = meanMat->getValue(id, 0);
    double variance = covMat->getValue(id, id);
    double pll = -0.5 * (sampleCount - 1.0) - 0.5 * log(2.0 * M_PI * variance) * sampleCount;
    return pll;
}

Matrix*
PotentialManager::createMBCovarianceMatrix(SlimFactor* factor, vector<int>& mbVars)
{
    int parentCount = mbVars.size();
    int varCount = parentCount + 1;
    int childID = factor->fId;
    int childIndex = varCount - 1;
    Matrix* covariance = new Matrix(varCount, varCount);
    for(int i = 0; i < parentCount; i++) {
        int vID = mbVars[i];
        for(int j = i; j < parentCount; j++) {
            int uID = mbVars[j];
            double cval = covMat->getValue(vID, uID);
            covariance->setValue(cval, i, j);
            covariance->setValue(cval, j, i);
        }
        double cval = covMat->getValue(vID, childID);
        covariance->setValue(cval, i, childIndex);
        covariance->setValue(cval, childIndex, i);
    }
    double cval = covMat->getValue(childID, childID);
    covariance->setValue(cval, childIndex, childIndex);
    return covariance;
}
