#include <assert.h>
#include <algorithm>
#include "die.h"
#include "gauss_jordan.hh"

int gauss_jordan( Matrix & A, Matrix & B)
// ======================================================================
// Linear equation solution by Gauss-Jordan elimination
//
// a[1..n][1..n] is the input matrix. b[1..n][1..m] is input
// containing the m right-hand side vectors. On output, a is replaced by
// its matrix inverse, and b is replaced by the corresponding set of
// solution vectors.
//
// Returns:
//   - A is replaced by its inverse
//   - B is the solution
//   - result of the function: 0 = success
//                            -1 = singular matrix encountered
//
// This was copied & modified from the 'numerical recepies online'.
// ======================================================================
{
	// Make sure the sizes are compatible.
	assert( A.n_rows == A.n_cols);
	assert( A.n_cols == B.n_rows);

	// Name the dimensions of the matrices for convenience.
	long n = A.n_rows;
	long m = B.n_cols;

	// The integer arrays ipiv, indxr and indxc are used for
	// bookkeeping on the pivoting.
	long indxc[n];
	long indxr[n];
	long ipiv[n];
	for( long i = 0 ; i < n ; i ++) ipiv[i] = 0;
	// This is the main loop over the columns to be reduced.
	long irow = -1, icol = -1;
	for( long i = 0 ; i < n ; i ++) {
		double big = 0.0;
		// This is the outer loop of the search for a pivot
		// element.
		for( long j = 0 ; j < n ; j ++) {
			if( ipiv[j] == 1) continue;
			for( long k = 0 ; k < n ; k ++) {
				if( ipiv[k] == 0) {
					double fa = fabs(A(j,k));
					if( fa >= big) {
						big = fa;
						irow = j;
						icol = k;
					}
				}
				else if( ipiv[k] > 1)
					return -1;
//					die( "gaussj: Singular matrix [1].");
			}
		}
		ipiv[ icol] ++;

		// We now have the pivot element, so we interchange rows,
		// if needed, to put the pivot element on the diagonal. The
		// columns are not physically interchanged, only relabeled:
		// indxc[i], the column of the ith pivot element, is the
		// ith column that is reduced, while indxr[i] is the row in
		// which that pivot element was originally located. If
		// indxr[i] 6= indxc[i] there is an implied column
		// interchange. With this form of bookkeeping, the solution
		// b's will end up in the correct order, and the inverse
		// matrix will be scrambled by columns.
		if( irow != icol) {
			for( long j = 0 ; j < n ; j ++)
				std::swap( A(irow,j), A(icol,j));
			for( long j = 0 ; j < m ; j ++)
				std::swap( B(irow,j), B(icol,j));
		}

		// We are now ready to divide the pivot row by the pivot
		// element, located at irow and icol.
		indxr[i] = irow;
		indxc[i] = icol;
		if( A(icol,icol) == 0.0)
			return -1;
//			die( "gaussj: Singular matrix [2].");
		double pivinv = 1/A(icol,icol);
		A(icol,icol) = 1.0;
		for( long j = 0 ; j < n ; j ++) A(icol,j) *= pivinv;
		for( long j = 0 ; j < m ; j ++) B(icol,j) *= pivinv;

		// Next, we reduce the rows, except for the pivot one, of
		// course.
		for( long j = 0 ; j < n ; j ++) {
			if( j != icol) {
				double dum = A(j,icol);
				A(j,icol) = 0;
				for( long k = 0 ; k < n ; k ++)
					A(j,k) -= A(icol,k) * dum;
				for( long k = 0 ; k < m ; k ++)
					B(j,k) -= B(icol,k) * dum;
			}
		}
	}
	// This is the end of the main loop over columns of the
	// reduction. It only remains to unscramble the solution in view of
	// the column interchanges. We do this by interchanging pairs of
	// columns in the reverse order that the permutation was built up.
	for( long j = n-1; j >= 0; j--) {
		if( indxr[j] != indxc[j])
			for( long k = 0; k < n; k++)
				std::swap( A(k,indxr[j]), A(k,indxc[j]));
	}
	
	// we are done
	return 0;
}

int gauss_jordan( Matrix & A)
// ======================================================================
// Linear equation solution by Gauss-Jordan elimination
//
// a[1..n][1..n] is the input matrix.On output, a is replaced by
// its matrix inverse.
//
// Returns:
//   - A is replaced by its inverse
//   - result of the function: 0 = success
//                            -1 = singular matrix encountered
//
// This was copied & modified from the 'numerical recepies online'.
// ======================================================================
{
	// Make sure the sizes are compatible.
	assert( A.n_rows == A.n_cols);

	// Name the dimensions of the matrices for convenience.
	long n = A.n_rows;

	// The integer arrays ipiv, indxr and indxc are used for
	// bookkeeping on the pivoting.
	long indxc[n];
	long indxr[n];
	long ipiv[n];
	for( long i = 0 ; i < n ; i ++) ipiv[i] = 0;
	// This is the main loop over the columns to be reduced.
	long irow = -1, icol = -1;
	for( long i = 0 ; i < n ; i ++) {
		double big = 0.0;
		// This is the outer loop of the search for a pivot
		// element.
		for( long j = 0 ; j < n ; j ++) {
			if( ipiv[j] == 1) continue;
			for( long k = 0 ; k < n ; k ++) {
				if( ipiv[k] == 0) {
					double fa = fabs(A(j,k));
					if( fa >= big) {
						big = fa;
						irow = j;
						icol = k;
					}
				}
				else if( ipiv[k] > 1)
					return -1;
//					die( "gaussj: Singular matrix [1].");
			}
		}
		ipiv[ icol] ++;

		// We now have the pivot element, so we interchange rows,
		// if needed, to put the pivot element on the diagonal. The
		// columns are not physically interchanged, only relabeled:
		// indxc[i], the column of the ith pivot element, is the
		// ith column that is reduced, while indxr[i] is the row in
		// which that pivot element was originally located. If
		// indxr[i] 6= indxc[i] there is an implied column
		// interchange. With this form of bookkeeping, the solution
		// b's will end up in the correct order, and the inverse
		// matrix will be scrambled by columns.
		if( irow != icol) {
			for( long j = 0 ; j < n ; j ++)
				std::swap( A(irow,j), A(icol,j));
		}

		// We are now ready to divide the pivot row by the pivot
		// element, located at irow and icol.
		indxr[i] = irow;
		indxc[i] = icol;
		if( A(icol,icol) == 0.0)
			return -1;
//			die( "gaussj: Singular matrix [2].");
		double pivinv = 1/A(icol,icol);
		A(icol,icol) = 1.0;
		for( long j = 0 ; j < n ; j ++) A(icol,j) *= pivinv;

		// Next, we reduce the rows, except for the pivot one, of
		// course.
		for( long j = 0 ; j < n ; j ++) {
			if( j != icol) {
				double dum = A(j,icol);
				A(j,icol) = 0;
				for( long k = 0 ; k < n ; k ++)
					A(j,k) -= A(icol,k) * dum;
			}
		}
	}
	// This is the end of the main loop over columns of the
	// reduction. It only remains to unscramble the solution in view of
	// the column interchanges. We do this by interchanging pairs of
	// columns in the reverse order that the permutation was built up.
	for( long j = n-1; j >= 0; j--) {
		if( indxr[j] != indxc[j])
			for( long k = 0; k < n; k++)
				std::swap( A(k,indxr[j]), A(k,indxc[j]));
	}
	
	// we are done
	return 0;
}
