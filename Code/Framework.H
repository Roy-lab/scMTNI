/*The class which brings it all together*/
#ifndef _FRAMEWORK
#define _FRAMEWORK

class MetaLearner;
class SpeciesDistance;

class Framework
{
	public:
		Framework();
		~Framework();
		Error::ErrorCode init(int, char**);
		int start();

	private:

		MetaLearner metaLearner;
		SpeciesDistance spData;
		double epsThreshold;
		char actionStr[256];
};
#endif
