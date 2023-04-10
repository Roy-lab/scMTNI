#ifndef _DISTANCE
#define _DISTANCE

#include<vector>

using namespace std;

class Distance
{
	public:
	Distance();
	~Distance();
	double computeSymmKLDivergence(double m1, double v1, double m2, double v2);
	double computeKLDivergence(double m1, double v1, double m2, double v2);
	double computeZstat(double m1,double v1,double m2, double v2,int sampleCnt);
	double computeCC(vector<double>& v1, vector<double>& v2);
	
	
};
#endif
