/*The class which brings it all together*/
#ifndef _FRAMEWORK
#define _FRAMEWORK

//All the managers need forward declarations
class MetaLearner;
class MappedOrthogroupReader;
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
		MappedOrthogroupReader mor;
		SpeciesDistance spData;
		double epsThreshold;
		char actionStr[256];
		int cvCnt;
};
#endif
