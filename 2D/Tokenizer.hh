#ifndef __TOKENIZER_H__
#define __TOKENIZER_H__

#include <stdio.h>
#include <string>
#include <stack>

class Tokenizer {
public:

	// static version of the tokenizer:
	static char * read_token( FILE * fp);
	static long read_long( FILE * fp, char * comment);
	static double read_double( FILE * fp, char * comment);
	static char * read_string( FILE * fp, char * comment);
	static int is_double( const char * str);
	static int is_double( const std::string & str);
	static double to_double( const char * str);
	static double to_double( const std::string & str);
	static int is_long( const char * str);
	static long to_long( const char * str);

	// dynamic version (with push & pop)
public:
	Tokenizer( FILE * fp);
	~Tokenizer();
	void push( const std::string & s);
	bool error( void);
	std::string read_token( void);
	long read_long( const char * comment);
	double read_double( const char * comment);
	std::string read_string( const char * comment);

private:
	FILE * fptr;
	std::stack<std::string> token_list;
	bool error_flag;
};



#endif
