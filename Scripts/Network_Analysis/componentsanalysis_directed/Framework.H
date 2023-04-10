#ifndef _FRAMEWORK_
#define _FRAMEWORK_
class Framework
{
	public:
		Framework();
		~Framework();
		int init(const char*);
		int start(const char*,int);
	private:
		Graph g;

};
#endif
