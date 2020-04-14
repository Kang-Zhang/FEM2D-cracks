#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
//#include <malloc.h>
#include <math.h>
#include "pinv.hh"
#include "die.hh"
#include "Matrix.hh"
#include "gauss_jordan.hh"

void Matrix::reset( void)
// ----------------------------------------------------------------------
// set all elements of the matrix to 0
// ----------------------------------------------------------------------
{
	/*
	for( long r = 0 ; r < n_rows ; r ++)
		for( long c = 0 ; c < n_cols ; c ++)
			data[r][c] = 0.0;
	*/

	long n = n_rows * n_cols;
	double * r = data[0];
	for( long i = 0 ; i < n ; i ++) {
		* r = 0.0;
		r ++;
	}
}

void Matrix::set( long r, long c, double val)
// ----------------------------------------------------------------------
// set the given element of the matrix to the specified value
// ----------------------------------------------------------------------
{
	assert( r >= 0 );
	assert( c >= 0 );
	assert( r < n_rows);
	assert( c < n_cols);
	assert( ! isnan( val));

	data[r][c] = val;
}

void Matrix::add( long r, long c, double val)
// ----------------------------------------------------------------------
// add the specified value to the given element of the matrix
// ----------------------------------------------------------------------
{
	assert( r >= 0 );
	assert( c >= 0 );
	assert( r < n_rows);
	assert( c < n_cols);
	assert( ! isnan( val));

	data[r][c] += val;
}

Matrix::Matrix( long p_rows, long p_cols, double * vals)
// ----------------------------------------------------------------------
// create a matrix p_rows x p_cols
// ----------------------------------------------------------------------
{
	n_rows = p_rows;
	n_cols = p_cols;

//	data = (double **) malloc( sizeof( double *) * n_rows);
	data = new double *[n_rows];
	assert( data != NULL);

	data[0] = new double[n_cols * n_rows];
	for( long r = 1 ; r < n_rows ; r ++)
		data[r] = data[0] + r*n_cols;

	// if initial values are present, use these to initialize the matrix,
	// otherwise set it to 0's
	if( vals == NULL)
		reset();
	else {
		long ind = 0;
		for( long r = 0 ; r < n_rows ; r ++)
			for( long c = 0 ; c < n_cols ; c ++)
				data[r][c] = vals[ind++];
	}
}

Matrix::Matrix( void)
{
	data = NULL;
	n_rows = -1;
	n_cols = -1;
}

Matrix::Matrix( const Matrix & M)
// ----------------------------------------------------------------------
// create a copy of matrix M
// ----------------------------------------------------------------------
{
//	fprintf( stderr, "Matrix:: copy constructor\n");

	n_rows = M.n_rows;
	n_cols = M.n_cols;

	data = new double *[n_rows];
	assert( data != NULL);

	data[0] = new double[n_cols * n_rows];
	for( long r = 1 ; r < n_rows ; r ++)
		data[r] = data[0] + r*n_cols;
	for( long r = 0 ; r < n_rows ; r ++)
		for( long c = 0 ; c < n_cols ; c ++)
			data[r][c] = M.data[r][c];
}

void Matrix::copy( const Matrix & A)
{
	assert( n_rows >= A.n_rows);
	assert( n_cols >= A.n_cols);

	for( long r = 0 ; r < A.n_rows ; r ++)
		for( long c = 0 ; c < A.n_cols ; c ++)
			data[r][c] = A.data[r][c];
}

Matrix::~Matrix()
{
	if( data != NULL) {
		delete [] data[0];
		delete [] data;
	}
}

void Matrix::transpose( Matrix & res) const
// ----------------------------------------------------------------------
// calculate a transpose of the current matrix, put the result in the
// provided space 'res'
// ----------------------------------------------------------------------
{
	assert( res.n_cols == n_rows);
	assert( res.n_rows == n_cols);

	for( long r = 0 ; r < n_rows ; r ++)
		for( long c = 0 ; c < n_cols ; c ++)
			res.data[ c][ r] = data[r][c];
}

void Matrix::multiply_by_scalar( double val)
// ----------------------------------------------------------------------
// multiply the whole matrix by a given scalar
// ----------------------------------------------------------------------
{
	for( long r = 0 ; r < n_rows ; r ++)
		for( long c = 0 ; c < n_cols ; c ++)
			data[r][c] *= val;
}

double Matrix::get( long r, long c) const
// ----------------------------------------------------------------------
// return the value of the given element in the matrix
// ----------------------------------------------------------------------
{
	assert( r >= 0 );
	assert( c >= 0 );
	assert( r < n_rows);
	assert( c < n_cols);

	return data[r][c];
}

double Matrix::operator() (long r, long c) const
{
	assert( r >= 0 && r < n_rows);
	assert( c >= 0 && c < n_cols);
	return data[r][c];
}

double & Matrix::operator() (long r, long c)
{
	assert( r >= 0 && r < n_rows);
	assert( c >= 0 && c < n_cols);
	return data[r][c];
}

Matrix & Matrix::operator= ( const Matrix & M)
{
//	fprintf( stderr, "Matrix:: assign operator\n");

	if( data != NULL) {
		delete [] data[0];
		delete [] data;
	}
	
	n_rows = M.n_rows;
	n_cols = M.n_cols;

	data = new double *[n_rows];
	assert( data != NULL);

	data[0] = new double[n_cols * n_rows];
	for( long r = 1 ; r < n_rows ; r ++)
		data[r] = data[0] + r*n_cols;
	for( long r = 0 ; r < n_rows ; r ++)
		for( long c = 0 ; c < n_cols ; c ++)
			data[r][c] = M.data[r][c];
	return * this;
}

Matrix operator* ( const Matrix & m1, const Matrix & m2)
{
	// make sure the matrices have proper dimensions
	assert( m1.n_cols == m2.n_rows);

	// create a result
	Matrix res( m1.n_rows , m2.n_cols );
	m1.multiply( m2, res);
	return res;
}

Matrix inverse( const Matrix & A)
{
	Matrix B = A;
	gauss_jordan( B);
	return B;
}

Matrix pseudo_inverse( const Matrix & A)
{
	fprintf( stderr, "Matrix::pseudo_inverse...\n");
	Matrix res = A;
	pinvT( res.data, res.n_rows, res.n_cols);
	return transpose(res);
}

/*
Matrix pseudo_inverse( const Matrix & _A)
{
	fprintf( stderr, "Matrix::pseudo_inverse...\n");
	// name the dimensions
	long m = _A.n_rows; long n = _A.n_cols;
//	if( m < n) m = n;
	Matrix A(m,n);
	fprintf( stderr, "1\n");
	A.copy( _A);
	fprintf( stderr, "2\n");

	// get the SVD of A
	Matrix U, V; Vector_double S(n);
	svd( A, U, S, V);
	// get the max. singular value
	double max = S(0);
	for( long i = 1 ; i < n ; i ++)
		if( max < S(i)) max = S(i);
	double eps = 1e-6/max;
	// construct the result
	Matrix res( A.n_cols, A.n_rows);
	for( long i = 0 ; i < A.n_cols ; i ++) {
		for( long j = 0 ; j < A.n_rows ; j ++) {
			res(j,i) = 0.0;
			for( long k = 0 ; k < A.n_cols ; k ++) {
				if( S(k) > eps)
					res(j,i) += V(i,k) * U(j,k) * (1/S(k));
			}
		}
	}
	return res;
}
*/

/*
void svd( const Matrix & A, Matrix & U, Vector_double & S, Matrix & V)
{
	fprintf( stderr, "STD...\n");
	A.print( "A=");
	Array2D<double> Ajama( A.n_rows, A.n_cols);
	// populate Ajama;
	for( long i = 0 ; i < A.n_rows ; i ++)
		for( long j = 0 ; j < A.n_cols ; j ++)
			Ajama[i][j] = A(i,j);

	cerr << "A before = " << Ajama << "\n";
	JAMA::SVD<double> jsvd(Ajama);
	cerr << "A after = " << Ajama << "\n";
	Array2D<double> Uj, Vj; Array1D<double> Sj;
	jsvd.getU( Uj);
	jsvd.getSingularValues( Sj);
	jsvd.getV( Vj);
	// populate the results
	U = Matrix( Uj.dim1(), Uj.dim2());
	for( long i = 0 ; i < U.n_rows ; i ++)
		for( long j = 0 ; j < U.n_cols ; j ++)
			U(i,j) = Uj[i][j];
	S = Vector_double( Sj.dim());
	for( long i = 0 ; i < S.size ; i ++)
		S(i) = Sj[i];
	V = Matrix( Vj.dim1(), Vj.dim2());
	for( long i = 0 ; i < V.n_rows ; i ++)
		for( long j = 0 ; j < V.n_cols ; j ++)
			V(i,j) = Vj[i][j];
	U.print( "U=");
	S.print( "S=");
	V.print( "V=");
}
*/

Matrix transpose( const Matrix & A)
{
	Matrix B( A.n_cols, A.n_rows);
	A.transpose( B);
	return B;
}

Matrix transpose( const Matrix & A);

void Matrix::multiply( const Matrix & m, Matrix & res) const
// ----------------------------------------------------------------------
// multiply the current matrix by 'm' and store the result in 'res'
// ----------------------------------------------------------------------
{
	// make sure current matrix, m and res have compatible dimensions
	assert( n_cols == m.n_rows);
	assert( n_rows == res.n_rows);
	assert( m.n_cols == res.n_cols);

	// perform the multiplication
	for( long r = 0 ; r < res.n_rows ; r ++)
		for( long c = 0 ; c < res.n_cols ; c ++)
		{
			res.data[r][c]=0;
			for( long k = 0 ; k < n_cols ; k ++)
				res.data[r][c] += data[r][k]*m.data[k][c];
		}
}

void Matrix::mult( const Vector_double & a, Vector_double & b) const
// ----------------------------------------------------------------------
// Multiply the matrix from the right by 'a' and put the result into b
// ----------------------------------------------------------------------
{
	assert( n_cols == a.size);
	assert( n_rows == b.size);
	for( long row = 0 ; row < n_rows ; row ++) {
		b.array[row] = 0;
		for( long col = 0 ; col < n_cols ; col ++) {
			b.array[row] += data[row][col] * a.array[col];
		}
	}
}

void Matrix::add( Matrix & m)
// ----------------------------------------------------------------------
// adds matrix m to itself
// ----------------------------------------------------------------------
{
	assert( n_rows == m.n_rows);
	assert( n_cols == m.n_cols);

	for( long r = 0 ; r < n_rows ; r ++)
		for( long c = 0 ; c < n_cols ; c ++)
			data[r][c] += m.data[r][c];
}

void Matrix::sub( Matrix & m)
// ----------------------------------------------------------------------
// subtract matrix m from itself
// ----------------------------------------------------------------------
{
	assert( n_rows == m.n_rows);
	assert( n_cols == m.n_cols);

	for( long r = 0 ; r < n_rows ; r ++)
		for( long c = 0 ; c < n_cols ; c ++)
			data[r][c] -= m.data[r][c];
}

void Matrix::inverse( Matrix & res)
// ----------------------------------------------------------------------
// calculates the inverse of a matrix and puts the result in res
// ----------------------------------------------------------------------
{
	assert( n_cols == n_rows);
	assert( res.n_cols == n_cols);
	assert( res.n_rows == n_rows);
	assert( n_cols > 0 || n_cols < 4);

	if( n_rows == 3) {
		double a = data[0][0];
		double b = data[0][1];
		double c = data[0][2];
		double d = data[1][0];
		double e = data[1][1];
		double f = data[1][2];
		double g = data[2][0];
		double h = data[2][1];
		double i = data[2][2];
		double det = -c*e*g+b*f*g+c*d*h-a*f*h-b*d*i+a*e*i;
		if( det == 0) {
			fprintf( stderr,
				 "Matrix::inverse() - non-invertible\n");
			print( "    matrix = ");
			assert( det != 0.0);
		}
		res.data[0][0] = (-f*h + e*i)/det;
		res.data[0][1] = ( c*h - b*i)/det;
		res.data[0][2] = (-c*e + b*f)/det;
		res.data[1][0] = ( f*g - d*i)/det;
		res.data[1][1] = (-c*g + a*i)/det;
		res.data[1][2] = ( c*d - a*f)/det;
		res.data[2][0] = (-e*g + d*h)/det;
		res.data[2][1] = ( b*g - a*h)/det;
		res.data[2][2] = (-b*d + a*e)/det;
		return;
	}

	if( n_rows == 2) {
		double a = data[0][0];
		double b = data[0][1];
		double c = data[1][0];
		double d = data[1][1];
		double det = -b*c+a*d;
		res.data[0][0] = d/det;
		res.data[0][1] = -b/det;
		res.data[1][0] = -c/det;
		res.data[1][1] = a/det;
		return;
	}

	if( n_rows == 1) {
		double a = data[0][0];
		res.data[0][0] = 1/a;
		return;
	}

	assert( 0);
}

double Matrix::determinant( void)
// ----------------------------------------------------------------------
// calculates the determinant of the matrix
// ----------------------------------------------------------------------
{
	assert( n_cols == n_rows);
	assert( n_cols > 0 || n_cols < 4);

	if( n_rows == 3) {
		double a = data[0][0];
		double b = data[0][1];
		double c = data[0][2];
		double d = data[1][0];
		double e = data[1][1];
		double f = data[1][2];
		double g = data[2][0];
		double h = data[2][1];
		double i = data[2][2];
		return -c*e*g+b*f*g+c*d*h-a*f*h-b*d*i+a*e*i;
	}

	if( n_rows == 2) {
		double a = data[0][0];
		double b = data[0][1];
		double c = data[1][0];
		double d = data[1][1];
		return -b*c+a*d;
	}

	if( n_rows == 1) {
		double a = data[0][0];
		return a;
	}

	assert( 0);
	// just so that the compiler does not complain:
	return 0;
}

int Matrix::is_symmetric( void) const
// ----------------------------------------------------------------------
// returns 1 if the matrix is symmetric, 0 otherwise
// ----------------------------------------------------------------------
{
	double max_diff = 0.0;

	// only a square matrix can ever be symmetric
	if( n_cols != n_rows) return 0;

	for( long r = 0 ; r < n_rows ; r ++)
		for( long c = r + 1 ; c < n_cols ; c ++) {
			double diff = fabs(data[r][c] - data[c][r]);
			if( diff > max_diff) max_diff = diff;
		}
	return is_zero( max_diff);
}

long Matrix::get_bandwidth( void) const
// ----------------------------------------------------------------------
// calculates and returns the bandwidth of the matrix
//
// returns: -1 if the matrix is not square
//           0 if the matrix is all zeros
//           n the bandwidth of the matrix
// ----------------------------------------------------------------------
{
	// first - make sure the matrix is square
	if( n_cols != n_rows) return -1;

	for( long b = n_cols ; b >= 1 ; b --)
		for( long r = 0 ; r < n_cols - b + 1; r ++) {
			long c = b + r - 1;
			if( ! is_zero( get(r,c)) || ! is_zero( get(c,r)))
				return b;
		}
	return 0;
}

void Matrix::print( char * str, char * fmtstr) const
{
	if( fmtstr == NULL) fmtstr = " %8.4f";
	long len = strlen( str);

	// print the matrix row by row
	for( long r = 0 ; r < n_rows ; r ++) {
		if( r == n_rows / 2)
			printf( "%s", str);
		else
			for( long i = 0 ; i < len ; i ++)
				printf( " ");
		if( r == 0) printf( " /");
		else if( r == n_rows-1) printf( " \\");
		else printf( "| ");
		for( long c =0 ; c < n_cols ; c ++)
			printf( fmtstr, data[r][c]);
		if( r == 0) printf( " \\");
		else if( r == n_rows-1) printf( " /");
		else printf( "  |");
		printf( "\n");
	}
}

void Matrix::mathematica_print( FILE * fp, char * name, char * fmtstr)
{
	if( fmtstr == NULL) fmtstr = " %8.4f";

	fprintf( fp, "%s=\n", name);
	fprintf( fp, "{\n");
	for( long i = 0 ; i < n_rows ; i ++)
	{
		fprintf( fp, "  {");
		fprintf( fp, fmtstr, data[i][0]);
		for( long j = 1 ; j < n_cols ; j ++) {
			fprintf( fp, ", ");
			fprintf( fp, fmtstr, data[i][j]);
		}
		if( i < n_rows - 1)
			fprintf( fp, " },\n");
		else
			fprintf( fp, " }\n");
	}
	fprintf( fp, "}\n");
}


