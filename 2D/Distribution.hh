#ifndef __DISTRIBUTION_HH__
#define __DISTRIBUTION_HH__

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "Tokenizer.hh"
#include "die.hh"

class Distribution { 

public:

	enum DistType { Poisson,
			Normal,
			Uniform,
			Constant,
			Exponential,
			Undefined };

	DistType get_dist_type( void) {
		return type;
	}

	Distribution() {
		type = Undefined;
	}

	double get_random_number( void) {
		assert( type != Undefined);
		
		if( type == Uniform)
			return get_uniform_rnd( n[0], n[1]);
		else if( type == Exponential)
			return get_exponential_rnd( n[0]);
		else if( type == Normal)
			return get_normal_rnd( n[0], n[1]);
		else if( type == Poisson)
			return get_poisson_rnd( n[0]);
		else if( type == Constant)
			return n[0];
		else
			assert( 0);
		// just so that the compiler does not complain
		return 0;
	}

	int read_from_file( FILE * fp) {
		// get the first token
		char * token = Tokenizer::read_string( fp, "Distr::token1");
		if( token == NULL) die( "Distribution: no data.");
//		fprintf( stderr, "Distribution: token1 = '%s'\n", token);
		// if the first token is a number, consider this to be
		// a constant distribution
		if( Tokenizer::is_double( token)) {
			type = Constant;
			n[0] = Tokenizer::to_double( token);
			return 0;
		}
		
		// the first token is not a number - try one of the
		// distributions
		if( strcasecmp( token, "constant") == 0)
			type = Constant;
		else if( strcasecmp( token, "poisson") == 0)
			type = Poisson;
		else if( strcasecmp( token, "normal") == 0)
			type = Normal;
		else if( strcasecmp( token, "exponential") == 0)
			type = Exponential;
		else if( strcasecmp( token, "uniform") == 0)
			type = Uniform;
		else die( "Distribution::wrong distribution type: '%s'.",
			  token);
		
		// skip the bracket
		token = Tokenizer::read_string( fp, "Distr::L-brack");
		if( token == NULL || strcmp( token, "(") != 0)
			die( "Distribution: left bracket missing.");
		// read bunch of numbers separated by commas, depending on
		// the distribution type
		long n_params = 0;
		if( type == Exponential) n_params = 1;
		else if( type == Uniform) n_params = 2;
		else if( type == Poisson) n_params = 1;
		else if( type == Normal) n_params = 2;
		else assert( 0);
		long count = 0;
		while( 1) {
			// read a number
			n[count] = Tokenizer::read_double( fp, "n[]");
			count ++;
			if( count == n_params) break;
			// skip a comma
			char * token = Tokenizer::read_string( fp, "comma");
			if( token == NULL || strcmp( token, ",") != 0)
				die( "Distribution: comma missing.");
		}

		// skip the right bracket
		token = Tokenizer::read_string( fp, "right bracket");
		if( token == NULL || strcmp( token, ")") != 0)
			die( "Distribution: right bracket missing '%s'.",
			     token);
		
		return 0;
	}

	void print( FILE * fp) {
		if( type == Exponential)
			fprintf( fp, "Exponential(%f)", n[0]);
		else if( type == Normal)
			fprintf( fp, "Normal(%f,%f)", n[0], n[1]);
		else if( type == Poisson)
			fprintf( fp, "Poisson(%f)", n[0]);
		else if( type == Uniform)
			fprintf( fp, "Uniform(%f,%f)", n[0], n[1]);
		else if( type == Constant)
			fprintf( fp, "Constant(%f)", n[0]);
		else fprintf( fp, "UnknownDistro()");
	}

	static double get_exponential_rnd( double avg) {
		double dum;
		do
			dum = drand48();
		while( dum == 0.0);
		return -log( dum);
	}

	static double get_uniform_rnd( double n1, double n2) {
		return drand48() * (n2-n1) + n1;
	}

        static double get_normal_rnd( double mean, double var);

        static double get_poisson_rnd( double xm);

private:

	double n[10];		// holds parameters

	DistType type;		// type of distribution
};

#endif
