#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <ctype.h>
#include "Tokenizer.hh"
#include "die.hh"

char * Tokenizer::read_token( FILE * fp) {
	static char buff[ 4096];
	// skip all white spaces and all comments
	int c; int in_comment = 0;
	while( 1) {
		c = fgetc( fp);
		if( c == EOF) return NULL;
		if( in_comment) {
			if( c == '\n') in_comment = 0;
			continue;
		}
		if( c == '#') { in_comment = 1; continue; }
		if( isspace(c)) continue;
		break;
	}

	// if we are reading a string - handle that separately
	if( c == '"') {
		long n = 0;
		while( 1) {
			c = fgetc( fp);
			if( c == '"') break;
			if( c == EOF) break;
			if( long( n) >= long( sizeof( buff))) break;
			buff[n] = c;
			n ++;
		}
		if( c == EOF) return NULL;
		if( long(n) >= long(sizeof( buff))) return NULL;
		buff[n] = '\0';
		return buff;
	}

	// at this point we have skipped all white spaces and comments
	long n = 0;
	int first_char = c;
	while( 1) {
		// if c is EOF we are done
		if( c == EOF) break;
		// if c is a blank, we are done
		if( isspace( c)) break;
		// if the buffer is full we are done
		if( long(n) >= long( sizeof( buff) -1)) break;
		// put c in the buffer
		buff[n] = c;
		n ++;
		// get the next character
		c = fgetc( fp);
		// if initial character  was a special character
		// (non-alphanum), we are done
		if( ! isalnum( first_char) &&
		    first_char != '-' &&
		    first_char != '+' &&
		    first_char != '.') break;
		// if the next character is a special one, we are done
		if( ! isalnum( c) &&
		    c != '-' &&
		    c != '+' &&
		    c != '.' &&
		    c != '_') break;
	}
	ungetc( c, fp);
	buff[n] = '\0';
//	fprintf( stderr, "*** token = '%s'\n", buff);
	if( c == EOF) return NULL;
	return buff;
}

long Tokenizer::read_long( FILE * fp, char * comment) {
	char * token = read_token( fp);
	if( token == NULL) die( "got NULL in read_long(%s)", comment);
	long n;
	int res = sscanf( token, "%ld", & n);
	if( res != 1) die( "read_long(%s)", comment);
	return n;
}

double Tokenizer::read_double( FILE * fp, char * comment) {
	char * token = read_token( fp);
	if( token == NULL) die( "got NULL in read_double(%s)", comment);
	double n;
	int res = sscanf( token, "%lf", & n);
	if( res != 1) die( "read_double(%s) - '%s' not a double", comment,
			   token);
	return n;
}

char * Tokenizer::read_string( FILE * fp, char * comment) {
	char * token = read_token( fp);
	if( token == NULL) die( "got NULL in read_string(%s)", comment);
	return token;
}

int Tokenizer::is_double( const char * str)
{
	double d;
	if( 1 == sscanf( str, "%lf", & d)) return 1; else return 0;
}

int Tokenizer::is_double( const std::string & str)
{
	return is_double( str.c_str());
}

double Tokenizer::to_double( const char * str)
{
	if( str == NULL) die( "to_double( NULL)");
	double d;
	if( 1 != sscanf( str, "%lf", & d))
		die( "Tokenizer::to_double('%s')\n", str);
	return d;
}

double Tokenizer::to_double( const std::string & str)
{
	return to_double( str.c_str());
}

int Tokenizer::is_long( const char * str)
{
	long n;
	if( 1 == sscanf( str, "%ld", & n)) return 1; else return 0;
}

long Tokenizer::to_long( const char * str)
{
	if( str == NULL) die( "to_long( NULL)");
	long n;
	if( 1 != sscanf( str, "%ld", & n))
		die( "Tokenizer::to_long('%s')\n", str);
	return n;
}


// ======================================================================
// The dynamic version of Tokenizer
// ======================================================================

Tokenizer::Tokenizer( FILE * fp)
{
	fptr = fp;
	error_flag = 0;
}

Tokenizer::~Tokenizer( )
{
}

bool Tokenizer::error( void) {
	return error_flag;
}

void Tokenizer::push( const std::string & s)
{
	token_list.push( s);
}

long Tokenizer::read_long( const char * comment)
{
	std::string token = read_token();
	if( error()) die( "got NULL in read_long(%s)", comment);
	long n;
	int res = sscanf( token.c_str(), "%ld", & n);
	if( res != 1) die( "read_long(%s)", comment);
	return n;
}

double Tokenizer::read_double( const char * comment)
{
	std::string token = read_token();
	if( error()) die( "got NULL in read_double(%s)", comment);
	double n;
	int res = sscanf( token.c_str(), "%lf", & n);
	if( res != 1) die( "read_double(%s) - '%s' not a double", comment,
			   token.c_str());
	return n;
}

std::string Tokenizer::read_string( const char * comment)
{
	std::string token = read_token();
	if( error()) die( "got NULL in read_string(%s)", comment);
	return token;
}

std::string Tokenizer::read_token( void)
{
	if( token_list.empty()) {
		char * str = read_token( fptr);
		if( str == NULL) {
			error_flag = 1;
			return std::string( "<<<ERROR>>>");
		}
		return std::string( str);
	} else {
		std::string res = token_list.top();
		token_list.pop();
		return res;
	}
}
