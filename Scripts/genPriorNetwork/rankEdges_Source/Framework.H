#ifndef _FRAMEWORK_
#define _FRAMEWORK_

#include <map>
#include <string>
#include <vector>
using namespace std;
class Framework
{
	public:
	Framework();
	~Framework();
	int readNetwork(const char*);
	int convertToPercentile(const char*,bool);
	private:
	int sort(map<string,double>*,vector<string>&,bool);
	map<string,map<string,double>*> nw;
	
};

#endif
