#ifndef _ERROR
#define _ERROR
class Error
{
	public:
		Error();
		~Error();
		typedef enum
		{
			SUCCESS,
			VARSCHEMA_ERR,
			DATAFILE_ERR,
			VARIABLE_ERR,
			SAMPLECNT_ERR,
			FACTOR_CREAT_ERR,
			MB_SEARCH_ERR,
			EVID_ERR,
			NULL_POT_ERR,
			VAR_ASSIGN_ERR,
			POT_QUERY_ERR,
			UNKNOWN
		}ErrorCode;
		
		static const char* getErrorString(int);
	private:
		static const char* errorString[];
};
#endif
