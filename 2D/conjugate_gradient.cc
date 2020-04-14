#include <assert.h>
#include <math.h>
#include <time.h>
#include "conjugate_gradient.hh"

int conjugate_gradient_diagpre(
	const SparseMatrix & A,
	const Vector_double & B,
	Vector_double & X,
	double max_err)
// ----------------------------------------------------------------------
// uses Gauss-Seidel iterative method to solve a system of linear equations
//     [A].[X] = [B]
// for [X], where [X] is the initial guess
// ----------------------------------------------------------------------
{
	fprintf( stderr, "Preconditioned CG!!!!\n");
	max_err /= 10;

	// benchmarking:
	clock_t start_time = clock();

	long n = A.n_rows;

	assert( n == A.n_cols);
	assert( n == B.size);
	assert( B.size == X.size);

	// Calculate a preconditioner M - an inverser of the diagonal elements
	// of A
	Vector_double M( n);
	for( long i = 0 ; i < n ; i ++) {
		double val = A.get(i,i);
		if( val < max_err/100) val = 1.0;
		M.array[i] = 1/val;
	}

// initialize:
//
//      g0 = M*A.x0 - M*b     d0 = -g0
//
// iterate while (g[k]^t . g[k]) is large:
//
//      alpha0 = g0^T . g0 / (d0^T . M*A . d0)
//
//      x1 = x0 + alpha0.d0
//
//      g1 = g0 + a0 . M*A . d0
//
//      beta0 = g1^T . g1 / (g0^T . g0)
//
//      d1 = - g1 + beta0 . d0
//
//      x0 = x1
//      g0 = g1
//      d0 = d1


	// setup x0 = X
	Vector_double x0( X);
	// initialize g0 = A.X - B
	Vector_double g0(n);
	A.mult( X, g0);
	// M*A.X
	for( long i = 0 ; i < n ; i ++) g0.array[i] *= M.array[i];
	// MB = M*B
	Vector_double MB( B);
	for( long i = 0 ; i < n ; i ++) MB.array[i] *= M.array[i];
	g0.sub( MB);
	// initialize d0 = -g0;
	Vector_double d0( g0);
	d0.multiply_by_scalar(-1);

	// construct some of the vectors used inside the loops,
	// so that we don't have to construct them all the time
	Vector_double Ad0(n);
	Vector_double x1(n);	
	Vector_double g1(n);
	Vector_double d1(n);

	// enter the iteration loop
	long n_iterations = 0;
	double gprod0 = 0;
	double bnorm = B.norm();
	double err;
	if( bnorm < max_err) bnorm = 1.0;
//	fprintf( stderr, "\t- CG bnorm = %.20f\n", bnorm);
	while( 1) {
		// calculate gprod0 = g0.g0 = norm(g0)
		gprod0 = g0.norm();
		if( isnan( gprod0)) {
			fprintf( stderr, "CG: ERROR - gprod0 is nan!!!\n");
			fprintf( stderr, "    iterations = %ld\n",
				 n_iterations);
			// here is some debug stuff - to help us in
			// debugging
			//
			// check the coefficient matrix, maybe one of the
			// entries is nan
			for( long i = 0 ; i < A.n_rows ; i ++) {
			  SparseMatrix::Row * r = A.rows[i];
			  for( long j = 0 ; j < r-> n_entries ; j ++) {
			    if( ! isnan( r-> entries[j].val)) continue;
			    fprintf( stderr, "    A[%ld,%ld] is nan\n",
				     i, r-> entries[j].col);
			  }
			}
			// check the B matrix as well
			for( long i = 0 ; i < B.size ; i ++) {
				if( ! isnan( B(i))) continue;
				fprintf( stderr, "    B[%ld] is nan\n",
					 i);
			}
			assert( 0);
		}
		// exit the loop if gprod0/norm(B) is small enough or the
		// number of iterations is very large
//              err = gprod0 / bnorm;
// 		if( err < max_err || n_iterations > A.n_rows * 2) {
// 			if( n_iterations > 0) break;
// 			break;
// 		}
		// exit loop if the max. of the error squared is small enough
// 		double g0_max_sq = g0.max_sq();
// 		err = g0_max_sq / bnorm;
		err = fabs( g0.max_abs());
		if( err < max_err ||
		    n_iterations > 2*A.n_rows)
		{
			break;
		}

		n_iterations ++;

		if( n_iterations > A.n_rows / 2)
			fprintf( stderr, "CG: gprod0/bnorm = %.30f\n",
				 gprod0/bnorm);

		// calculate alpha0
		//      alpha0 = g0^T . g0 / (d0^T . M * A . d0)
		// --------------------------------------------------
		// calculate Ad0 = A.d0
		A.mult( d0, Ad0);
		// Ad0 = M*A.d0
		for( long i = 0 ; i < n ; i ++) Ad0.array[i] *= M.array[i];
		// figure out d0.Ad0
		double d0TAd0 = d0.dot_product( Ad0);
		// figure out alpha0
		double alpha0 = gprod0 / d0TAd0;

		// calculate x1
		//      x1 = x0 + alpha0.d0
		x1 = d0;
		x1.multiply_by_scalar(alpha0);
		x1.add(x0);

		// calculate g1
		//      g1 = g0 + a0 . A . d0
		//         = g0 + a0 . Ad0
		// --------------------------------------------------
		g1 = Ad0;
		g1.multiply_by_scalar(alpha0);
		g1.add(g0);

		// calculate beta0
		//      beta0 = g1^T . g1 / (g0^T . g0)
		//            = norm(g1) / gprod0
		// --------------------------------------------------
		double beta0 = g1.norm() / gprod0;

		// calculate d1
		//      d1 = - g1 + beta0 . d0
		// --------------------------------------------------
		d1 = d0;
		d1.multiply_by_scalar( beta0);
		d1.sub( g1);

		// x0 <- x1, g0 <- g1, d0 <- d1
		x0 = x1;
		g0 = g1;
		d0 = d1;
	}

	// put the result into X
	X.copy( x0);
	// benchmarking:
	clock_t end_time = clock();
// 	fprintf( stderr,
// 		 "\t- CG() n_it=%ld gprod=%.20f in %.3fs (%.3fs/it)\n",
// 		 n_iterations, gprod0/bnorm,
// 		 double( end_time - start_time) / CLOCKS_PER_SEC,
// 		 (double(end_time-start_time)/CLOCKS_PER_SEC)/n_iterations);
	fprintf( stderr,
		 "\t- CG() n_it=%ld err=%.20f in %.3fs (%.3fs/it)\n",
		 n_iterations, err,
		 double( end_time - start_time) / CLOCKS_PER_SEC,
		 (double(end_time-start_time)/CLOCKS_PER_SEC)/n_iterations);
	return 0;
}

int conjugate_gradient(
	const SparseMatrix & A,
	const Vector_double & B,
	Vector_double & X,
	double max_err)
// ----------------------------------------------------------------------
// uses Gauss-Seidel iterative method to solve a system of linear equations
//     [A].[X] = [B]
// for [X], where [X] is the initial guess
// ----------------------------------------------------------------------
{
	max_err /= 10;

	// benchmarking:
	clock_t start_time = clock();

	// TBD:
	// Precondition matrix A & B

	long n = A.n_rows;

	assert( n == A.n_cols);
	assert( n == B.size);
	assert( B.size == X.size);

// initialize:
//
//      g0 = A.x0 - b     d0 = -g0
//
// iterate while (g[k]^t . g[k]) is large:
//
//      alpha0 = g0^T . g0 / (d0^T . A . d0)
//
//      x1 = x0 + alpha0.d0
//
//      g1 = g0 + a0 . A . d0
//
//      beta0 = g1^T . g1 / (g0^T . g0)
//
//      d1 = - g1 + beta0 . d0
//
//      x0 = x1
//      g0 = g1
//      d0 = d1


	// setup x0 = X
	Vector_double x0( X);
	// initialize g0 = A.X - B
	Vector_double g0(n);
	A.mult( X, g0);
	g0.sub( B);
	// initialize d0 = -g0;
	Vector_double d0( g0);
	d0.multiply_by_scalar(-1);

	// construct some of the vectors used inside the loops,
	// so that we don't have to construct them all the time
	Vector_double Ad0(n);
	Vector_double x1(n);	
	Vector_double g1(n);
	Vector_double d1(n);

	// enter the iteration loop
	long n_iterations = 0;
	double gprod0 = 0;
	double bnorm = B.norm();
	double err;
	if( bnorm < max_err) bnorm = 1.0;
//	fprintf( stderr, "\t- CG bnorm = %.20f\n", bnorm);
	while( 1) {
		// calculate gprod0 = g0.g0 = norm(g0)
		gprod0 = g0.norm();
		if( isnan( gprod0)) {
			fprintf( stderr, "CG: ERROR - gprod0 is nan!!!\n");
			fprintf( stderr, "    iterations = %ld\n",
				 n_iterations);
			// here is some debug stuff - to help us in
			// debugging
			//
			// check the coefficient matrix, maybe one of the
			// entries is nan
			for( long i = 0 ; i < A.n_rows ; i ++) {
			  SparseMatrix::Row * r = A.rows[i];
			  for( long j = 0 ; j < r-> n_entries ; j ++) {
			    if( ! isnan( r-> entries[j].val)) continue;
			    fprintf( stderr, "    A[%ld,%ld] is nan\n",
				     i, r-> entries[j].col);
			  }
			}
			// check the B matrix as well
			for( long i = 0 ; i < B.size ; i ++) {
				if( ! isnan( B(i))) continue;
				fprintf( stderr, "    B[%ld] is nan\n",
					 i);
			}
			assert( 0);
		}
		// exit the loop if gprod0/norm(B) is small enough or the
		// number of iterations is very large
//              err = gprod0 / bnorm;
// 		if( err < max_err || n_iterations > A.n_rows * 2) {
// 			if( n_iterations > 0) break;
// 			break;
// 		}
		// exit loop if the max. of the error squared is small enough
// 		double g0_max_sq = g0.max_sq();
// 		err = g0_max_sq / bnorm;
		err = fabs( g0.max_abs());
		if( err < max_err ||
		    n_iterations > 2*A.n_rows)
		{
			break;
		}

		n_iterations ++;

		if( n_iterations > A.n_rows)
			fprintf( stderr, "CG: gprod0/bnorm = %.30f\n",
				 gprod0/bnorm);

		// calculate alpha0
		//      alpha0 = g0^T . g0 / (d0^T . A . d0)
		// --------------------------------------------------
		// calculate Ad0 = A.d0
		A.mult( d0, Ad0);
		// figure out d0.Ad0
		double d0TAd0 = d0.dot_product( Ad0);
		// figure out alpha0
		double alpha0 = gprod0 / d0TAd0;

		// calculate x1
		//      x1 = x0 + alpha0.d0
		x1 = d0;
		x1.multiply_by_scalar(alpha0);
		x1.add(x0);

		// calculate g1 = the error of the new solution
		//      g1 = g0 + a0 . A . d0
		//         = g0 + a0 . Ad0
		// --------------------------------------------------
		g1 = Ad0;
		g1.multiply_by_scalar(alpha0);
		g1.add(g0);

		// calculate beta0
		//      beta0 = g1^T . g1 / (g0^T . g0)
		//            = norm(g1) / gprod0
		// --------------------------------------------------
		double beta0 = g1.norm() / gprod0;

		// calculate d1
		//      d1 = - g1 + beta0 . d0
		// --------------------------------------------------
		d1 = d0;
		d1.multiply_by_scalar( beta0);
		d1.sub( g1);

		// x0 <- x1, g0 <- g1, d0 <- d1
		x0 = x1;
		g0 = g1;
		d0 = d1;
	}

	// put the result into X
	X.copy( x0);
	// benchmarking:
	clock_t end_time = clock();
// 	fprintf( stderr,
// 		 "\t- CG() n_it=%ld gprod=%.20f in %.3fs (%.3fs/it)\n",
// 		 n_iterations, gprod0/bnorm,
// 		 double( end_time - start_time) / CLOCKS_PER_SEC,
// 		 (double(end_time-start_time)/CLOCKS_PER_SEC)/n_iterations);
	fprintf( stderr,
		 "\t- CG() n_it=%ld err=%.20f in %.3fs (%.3fs/it)\n",
		 n_iterations, err,
		 double( end_time - start_time) / CLOCKS_PER_SEC,
		 (double(end_time-start_time)/CLOCKS_PER_SEC)/n_iterations);
	return 0;
}

