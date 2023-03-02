#include <math.h>
#include "Distance.H"


Distance::Distance()
{
}

Distance::~Distance()
{
}

double
Distance::computeSymmKLDivergence(double m1, double v1, double m2, double v2)
{
	double kl1=computeKLDivergence(m1,v1,m2,v2);
	double kl2=computeKLDivergence(m2,v2,m1,v1);
	double klsim=(kl1+kl2)/2;
	return klsim;
}

double 
Distance::computeKLDivergence(double m1, double v1, double m2, double v2)
{
	double t1=log(v2/v1);
	double t2=((m1-m2)*(m1-m2))/(v2);
	double t3=v1/v2;
	double kld=(t1+t2+t3-1)/2;
	return kld;
}

double 
Distance::computeZstat(double m1,double v1,double m2, double v2,int sampleCnt)
{
	double temp1=m1-m2;
	double temp2=sqrt((v1+v2)/sampleCnt);
	double zstat=temp1/temp2;
	return zstat;
}

double
Distance::computeCC(vector<double>& v1, vector<double>& v2)
{
	double cc=0;
	double m1=0;
	for(int i=0;i<v1.size();i++)
	{
		m1=m1+v1[i];
	}
	m1=m1/v1.size();

	double m2=0;
	for(int i=0;i<v2.size();i++)
	{
		m2=m2+v2[i];
	}
	m2=m2/v2.size();
	
	double xx=0;
	double yy=0;
	double xy=0;
	double oppRel=0;
	for(int i=0;i<v1.size();i++)
	{
		double diff1=v1[i]-m1;
		xx=xx + (diff1 * diff1);
		double diff2=v2[i]-m2;
		yy=yy+ (diff2 * diff2);
		xy=xy+(diff1*diff2);  	
		if( (diff1*diff2) < 0)
		{
			oppRel++;
		}
	}
	cc=sqrt((xy*xy)/(xx*yy));
	double threshold=v1.size()/2.0;

	if(oppRel > threshold)
	{
		cc=cc*(-1);
	}
	return cc;
}
