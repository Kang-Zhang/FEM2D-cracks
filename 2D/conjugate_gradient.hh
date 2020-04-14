#ifndef __CONJUGATE_GRADIENT_HH__
#define __CONJUGATE_GRADIENT_HH__

#include <stdio.h>
#include <stdlib.h>
#include "Matrix.hh"
#include "SparseMatrix.hh"

int conjugate_gradient_diagpre(
	const SparseMatrix & A,
	const Vector_double & B,
	Vector_double & X,
	double max_err);

int conjugate_gradient(
	const SparseMatrix & A,
	const Vector_double & B,
	Vector_double & X,
	double max_err);

#endif

