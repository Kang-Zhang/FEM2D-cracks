#ifndef __STRESSTENSOR_HH__
#define __STRESSTENSOR_HH__

#include "Vector3d.hh"

class StressTensor {
public:
	double sxx, syy, szz, syz, sxz, sxy;

	StressTensor( void);

	StressTensor( double p_sxx, double p_syy, double p_szz,
		      double p_syz, double p_sxz, double p_sxy);
	
	void get_eigenvalues( double & a1, double & a2, double & a3);

	long get_eigenvectors( double a1, Vector3d & v1, Vector3d & v2);

	double eigensystem_err( double a, Vector3d & v1);

	void set( double p_sxx, double p_syy, double p_szz,
		  double p_syz, double p_sxz, double p_sxy);

	void add( StressTensor & s);
	void multiply_by_scalar( double a);
	void print( FILE * fp, char * str, char * fmtstr = NULL);
	double diff_sq( StressTensor & s);
	void make_eigen_tensor( Vector3d & v);
	char * to_str( char * fmt_str = NULL);
};

// overloaded operators - to make our lives simpler
// ----------------------------------------------------------------------
//
// multiplication - multiply a tensor by a vector on the right

Vector3d operator * ( StressTensor & s, Vector3d & v);

#endif
