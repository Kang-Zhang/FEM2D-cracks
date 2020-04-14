#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

double my_sqrt( double d)
{
	if( d < 0) {
		fprintf( stderr, "Warning: my_sqrt(%.20f)\n", d);
		d = 0;
	}
//	assert( d >= 0);
	return sqrt( d);
}

double my_acos( double d)
{
	if( d < -1 || d > 1) {
		fprintf( stderr, "Warning: my_acos(%.20f)\n", d);
		if( d < -1) d = -1;
		if( d > 1) d = 1;
	}
//	assert( -1.0 <= d);
//	assert( d <= 1.0);
	return acos( d);
}

