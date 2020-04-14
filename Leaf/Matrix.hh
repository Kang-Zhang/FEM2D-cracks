#ifndef __MATRIX_HH__
#define __MATRIX_HH__

#include "Vector.hh"
#include <math.h>

class Matrix {

private:

public:
	
	long n_rows;
	long n_cols;
	double ** data;

	Matrix();

	Matrix( long p_rows , long p_cols, double * vals = NULL);
	Matrix( const Matrix & M);
	~Matrix();

	void reset( void);

	void copy( const Matrix & A);
	void resize( long new_n_rows, long new_n_cols);
	void add( long r, long c, double val);
	void add( Matrix & m);
	void sub( Matrix & m);
	void set( long r, long c, double val);
	void transpose( Matrix & res) const;
	void multiply_by_scalar( double val);
	double determinant( void);
	double get( long r, long c) const;
	void inverse( Matrix & res);
	void multiply( const Matrix & m, Matrix & res) const;
	void mult( const Vector_double & a, Vector_double & b) const;
	void print( char * str, char * fmtstr = NULL) const;
	void mathematica_print( FILE * fp, char * name, char * fmtstr = NULL);
	int is_symmetric( void) const;
	long get_bandwidth( void) const;
	int is_zero( double d) const { return fabs( d) < 1e-10; }
	double operator() (long r, long c) const;
	double & operator() (long r, long c);
	Matrix & operator= ( const Matrix & m);
};

Matrix inverse( const Matrix & A);
Matrix pseudo_inverse( const Matrix & A);
void svd( const Matrix & A, Matrix & U, Vector_double & S, Matrix & V);
Matrix transpose( const Matrix & A);
Matrix operator* ( const Matrix & m1, const Matrix & m2);

#endif
