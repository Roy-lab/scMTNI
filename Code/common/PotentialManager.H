/* The goal of this class is to take all the MarkovBlankets and the factors and
 *  the data and create all potentials involving valid pairings MarkovBlankets 
 *  and factors. We create all putative potentials upfront because we want to
 *  iterate over the data once and update counts of configurations in one go*/

#ifndef _POTENTIAL_MANAGER
#define _POTENTIAL_MANAGER

#include <vector>
#include <unordered_map>
#include "Matrix.H"

using namespace std;

class SlimFactor;
class VariableManager;

class PotentialManager
{
    public:
        PotentialManager();
        ~PotentialManager();

        void loadEvidenceFromTable(string fileName, const vector<string>& variableList);
        VariableManager* getVariableManager();
        double computeUnivariateLL(int);
        double computeConditionalLL(SlimFactor* sFactor, int&);
        void computePotentialMBCovMean(SlimFactor* sFactor, double& mbcondvar, double& mbbias, unordered_map<int,double>& mbcondMean_Vect);
        void dumpVarMB(SlimFactor* sFactor, ofstream& oFile);

    private:
        double computeJointLL(int dim, double determinant);
        Matrix* createMBCovarianceMatrix(SlimFactor* factor, vector<int>& mbVars);

        Matrix* meanMat;
        Matrix* covMat;
        int sampleCount;
        VariableManager* vMgr;
};
#endif
