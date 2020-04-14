#ifndef __LINFUNC_HH__
#define __LINFUNC_HH__

double linfunc( double x1, double v1, double x2, double v2, double x) {
	if( x < x1) return v1;
	if( x > x2) return v2;
	if( x1 == x2) return v1;
	return ((x-x1)/(x2-x1)) * (v2-v1) + v1;
}


#endif
