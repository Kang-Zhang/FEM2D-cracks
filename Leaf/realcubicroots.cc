#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "xmath.h"
#include "realcubicroots.h"

// extern int ns_debug_flag;

/*
static void print_in_hex( double val) {
	fprintf( stderr, "\t");
	unsigned char * ptr = (unsigned char *) & val;
	for( long i = 0 ; i < long( sizeof( val)) ; i ++) {
		unsigned int c = ptr[i];
		fprintf( stderr, "%x%x", c/16, c%16);
	}
	fprintf( stderr, "\n");
}
*/

void
get_cubic_roots(
	double a1, double a2, double a3,
	double & l1, double & l2, double & l3)
{
/*
	if( ns_debug_flag) {
		fprintf( stderr, "\n\nGET_CUBIC_ROOTS\n");
		fprintf( stderr, "  a1=%.30f\n", a1);
		print_in_hex( a1);
		fprintf( stderr, "  a2=%.30f\n", a2);
		print_in_hex( a2);
		fprintf( stderr, "  a3=%.30f\n", a3);
		print_in_hex( a3);
		
		double three = 3.0;
		double one = 1.0;
		double one_third = one / three;
		double one_third_sq = one_third * one_third;
		fprintf( stderr, "  (1/3)^2 = %.30f\n", one_third_sq);
	}
*/

        double Q = (a1 * a1 - 3 * a2) / 9.0;
        double R = (2 * a1 * a1 * a1 - 9 * a1 * a2 + 27 * a3) / 54.0;
        double Qcubed = Q * Q * Q;

/*
	if( ns_debug_flag) {
		fprintf( stderr, "  Q=%.30f\n", Q);
		fprintf( stderr, "  R=%.30f\n", R);
		fprintf( stderr, "  Qcubed=%.30f\n", Qcubed);
	}
*/

	double ratio = R / my_sqrt(Qcubed);
	if( isnan( ratio)) {
		if( Q*R < 0) ratio = -1;
		else if( Q*R > 0) ratio = 1;
		else ratio = 0;
	} else if( ratio < -1) ratio = -1;
	else if( ratio > 1) ratio = 1;

/*
	if( ns_debug_flag) {
		fprintf( stderr, "  ratio=%.20f\n", ratio);
	}
*/

	double theta = my_acos( ratio);

/*
	if( ns_debug_flag) {
		fprintf( stderr, "  theta=%.20f\n", theta);
	}
*/

	assert( ! isnan( theta));
	double sqrtQ = my_sqrt(Q);
	l1 = -2 * sqrtQ * cos(theta / 3) - a1 / 3;
	l2 = -2 * sqrtQ * cos((theta + 2 * M_PI) / 3) - a1 / 3;
	l3 = -2 * sqrtQ * cos((theta + 4 * M_PI) / 3) - a1 / 3;

/*
	if( ns_debug_flag) {
		fprintf( stderr, "  l1=%.20f\n", l1);
		fprintf( stderr, "  l2=%.20f\n", l2);
		fprintf( stderr, "  l3=%.20f\n", l3);
	}
*/

}
