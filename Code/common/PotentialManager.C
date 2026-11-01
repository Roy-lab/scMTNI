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

void
PotentialManager::computeLLs(int factorID, vector<int>& existingParents, vector<int>&candidateParents, unordered_map<int, double>&scores) {

	if (existingParents.size() == 0) {
		computeSingleParentLLs(factorID, sampleCount, candidateParents, scores);
		return;
	}

	// In order to compute likelihood, we need the variance of the factor conditioned upon existing parents and the candidate parent.
	// We use the following framing
	// A : Factor
	// B : existing parents
	// C : candidate parent
	// Var(A | B, C) = Var(A | B) - Cov(AC | B) * Var(C | B)^-1 * Cov(CA | B)
	// Var (A | B) = Var(A) - Cov(AB) * Var(B)^-1 * Cov(BA)
	// Var (C | B) = Var(C) - Cov(CB) * Var(B)^-1 * Cov(BC)
	// Cov (AC | B) = Cov(AC) - Cov(AB) * Var(B)^-1 * Cov(BC)

	// This allows us to precompute Var(A | B), Var(B)^-1, and Cov(AB) * Var(B)^-1, because they don't rely on C. Since we only add a single
	// candidate parent at a time, the per candidate work doesn't require any matrix inverts, only scalar and vector multiplication.

	int parentCount = existingParents.size();

	// Var(A)
	double factorVariance = covMat->getValue(factorID, factorID);

	// Cov(AB)
	gsl_vector* existingParentMarginalVariances = gsl_vector_alloc(parentCount);

	// Cov(BB)
	gsl_matrix* existingParentCovariances = gsl_matrix_alloc(parentCount, parentCount);

	for (int i = 0; i < parentCount; i++) {
		int varAID = existingParents[i];
		double marginalCovariance = covMat->getValue(factorID, varAID);
		gsl_vector_set(existingParentMarginalVariances, i, marginalCovariance);

		for (int j = i; j < parentCount; j++) {
			int varBID = existingParents[j];
			double covariance = covMat->getValue(varAID, varBID);
			gsl_matrix_set(existingParentCovariances, i, j, covariance);
			gsl_matrix_set(existingParentCovariances, j, i, covariance);
		}
	}

	gsl_permutation* permutation = gsl_permutation_alloc(parentCount);

	int signum=0;
	gsl_linalg_LU_decomp(existingParentCovariances, permutation, &signum);

	// Var(B)^-1
	gsl_matrix* parentCovInverse = gsl_matrix_alloc(parentCount, parentCount);
	gsl_linalg_LU_invert(existingParentCovariances, permutation, parentCovInverse);

	gsl_matrix_free(existingParentCovariances);
	gsl_permutation_free(permutation);

	// Cov(AB) * Var(B)^-1
	gsl_vector* prod = gsl_vector_alloc(parentCount);
	gsl_blas_dgemv(CblasTrans, 1, parentCovInverse, existingParentMarginalVariances, 0, prod);

	// Cov(AB) * Var(B)^-1 * Cov(BA)
	double dot = 0.0;
	gsl_blas_ddot(prod, existingParentMarginalVariances, &dot);

	// Var(A|B)
	double existingConditionalVariance = factorVariance - dot;

	gsl_vector_free(existingParentMarginalVariances);

	for (int i = 0; i < candidateParents.size(); i++) {
		int candidateID = candidateParents[i];

		// Var(C)
		double candidateVariance = covMat->getValue(candidateID, candidateID);

		// Cov(BC)
		gsl_vector* candidateMarginalVariances = gsl_vector_alloc(parentCount);

		for (int i = 0; i < parentCount; i++) {
			int varAID = existingParents[i];
			double marginalCovariance = covMat->getValue(candidateID, varAID);
			gsl_vector_set(candidateMarginalVariances, i, marginalCovariance);
		}

		// Cov(BC) * Var(B)^-1
		gsl_vector* candidateProd = gsl_vector_alloc(parentCount);
		gsl_blas_dgemv(CblasTrans, 1, parentCovInverse, candidateMarginalVariances, 0, candidateProd);

		// Cov(BC) * Var(B)^-1 * Cov(CB)
		double dot = 0.0;
		gsl_blas_ddot(candidateProd, candidateMarginalVariances, &dot);

		// Var(C|B)
		double candidateConditionalVariance = candidateVariance - dot;

		if (candidateConditionalVariance < 1e-10) {
			continue;
		}

		// Cov(AC)
		double candidateFactorCovariance = covMat->getValue(candidateID, factorID);

		// Cov(AB) * Var(B)^-1 * Cov(BC)
		dot = 0.0;
		gsl_blas_ddot(prod, candidateMarginalVariances, &dot);

		gsl_vector_free(candidateProd);
		gsl_vector_free(candidateMarginalVariances);

		// Cov(AC | B)
		double factorAndCandidateConditionalCov = candidateFactorCovariance - dot;

		// Var(A | B, C)
		double finalVariance = existingConditionalVariance - factorAndCandidateConditionalCov * factorAndCandidateConditionalCov / candidateConditionalVariance;

		if (finalVariance < 1e-5) {
			finalVariance = 1e-5;
		}

		if(isnan(finalVariance) || isinf(finalVariance)) {
			continue;
		}

		scores[candidateID] = -0.5 * ((sampleCount - 1) + sampleCount * log(2 * M_PI) + sampleCount * log(finalVariance));
	}

	gsl_matrix_free(parentCovInverse);
	gsl_vector_free(prod);
}

// Computes a gaussian log likelihood of factor id conditioned upon a single parent, for each candidate parent.
// The single parent case is broken out into a separate function because it can be done with fast scalar arithmetic.
void
PotentialManager::computeSingleParentLLs(int factorID, int sampleSize, vector<int>&candidateParents, unordered_map<int, double>&scores) {

	double factorVariance = covMat->getValue(factorID, factorID);

	for (int i = 0; i < candidateParents.size(); i++) {
		int candidateID = candidateParents[i];

		double candidateVariance = covMat->getValue(candidateID, candidateID);

		if (candidateVariance < 1e-10) {
			continue;
		}

		double factorCandidateCov = covMat->getValue(factorID, candidateID);
		double finalVariance = factorVariance - factorCandidateCov * factorCandidateCov / candidateVariance;

		if (finalVariance < 1e-5) {
			finalVariance = 1e-5;
		}

		if(isnan(finalVariance) || isinf(finalVariance)) {
			continue;
		}

		scores[candidateID] = -0.5 * ((sampleSize - 1) + sampleSize * log(2 * M_PI) + sampleSize * log(finalVariance));
	}
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
