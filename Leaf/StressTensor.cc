#include <math.h>
#include <assert.h>
#include <math.h>
#include "xmath.h"
#include "StressTensor.hh"
#include "nullspace.h"
#include "realcubicroots.h"

StressTensor::StressTensor(
	double p_sxx, double p_syy, double p_szz,
	double p_syz, double p_sxz, double p_sxy)
{
	sxx = p_sxx;
	syy = p_syy;
	szz = p_szz;
	syz = p_syz;
	sxz = p_sxz;
	sxy = p_sxy;
}

StressTensor::StressTensor( void)
{
	sxx = 0.0;
	syy = 0.0;
	szz = 0.0;
	syz = 0.0;
	sxz = 0.0;
	sxy = 0.0;
}

void
StressTensor::set(
	double p_sxx, double p_syy, double p_szz,
	double p_syz, double p_sxz, double p_sxy)
{
	sxx = p_sxx;
	syy = p_syy;
	szz = p_szz;
	syz = p_syz;
	sxz = p_sxz;
	sxy = p_sxy;
}

void
StressTensor::add( StressTensor & s)
{
	sxx += s.sxx;
	syy += s.syy;
	szz += s.szz;
	syz += s.syz;
	sxz += s.sxz;
	sxy += s.sxy;
}

void
StressTensor::multiply_by_scalar( double a)
{
	sxx *= a;
	syy *= a;
	szz *= a;
	syz *= a;
	sxz *= a;
	sxy *= a;
}

#ifdef DONT_COMPILE
void 
StressTensor::get_eigenvalues_old(
	double & l1, double & l2, double & l3)
// ----------------------------------------------------------------------
// calculate 3 eigenvalues from the stress tensor, return them in sorted
// order ( l1 < l2 < l3)
// ----------------------------------------------------------------------
{
	// the determinant of the matrix:
	//     stress - a.Identity
	// is a cubic polynomial of 'a'. Setting it to zero and solving
	// for 'a' yields 3 values --> eigenvalues of stress
	//
	// the determinant (according to  mathematica) is:
	// 
	//    Power(a,3) + Power(sxz,2)*syy - 2*sxy*sxz*syz + sxx*Power(syz,2)
	//    + Power(a,2)*(-sxx - syy - szz) + Power(sxy,2)*szz -
	//    sxx*syy*szz + a*(-Power(sxy,2) - Power(sxz,2) + sxx*syy
	//    - Power(syz,2) + (sxx + syy)*szz);
	// 
	// we extract the coefficients c1, c2, c3, so that
	// det = a^3 + c1*a^2 + c2*a + c3

	assert( ! isnan( sxx));
	assert( ! isnan( syy));
	assert( ! isnan( szz));
	assert( ! isnan( syz));
	assert( ! isnan( sxz));
	assert( ! isnan( sxy));

//	fprintf( stderr, "s.. = %.10f %.10f %.10f %.10f %.10f %.10f\n",
//		 sxx, syy, szz, syz, sxz, sxy);

	double c1 = -sxx - syy - szz;
	double c2 = -sxy*sxy - sxz*sxz + sxx*syy - syz*syz + (sxx + syy)*szz;
	double c3 = sxz*sxz*syy - 2*sxy*sxz*syz + sxx*syz*syz + 
		+ sxy*sxy*szz -	sxx*syy*szz;

//	fprintf( stderr, "c123= %.10f %.10f %.10f\n",
//		 c1, c2, c3);

	// solve for roots
	double Q = (c1 * c1 - 3 * c2) / 9;
	double R = (2 * c1 * c1 * c1 - 9 * c1 * c2 + 27 * c3) / 54;
	double Qcubed = Q * Q * Q;
	double d = Qcubed - R * R;

	// three real roots - this should always be the case for
	// symmetric matrices
	if (d >= 0) {
		double cosval = R / my_sqrt(Qcubed);
		if( isnan( cosval)) cosval = 0;
		else if( cosval < -1) cosval = -1;
		else if( cosval > 1) cosval = 1;
		double theta = my_acos( cosval);
		assert( ! isnan( theta));
		double sqrtQ = my_sqrt(Q);
		l1 = -2 * sqrtQ * cos( theta           / 3) - c1 / 3;
		l2 = -2 * sqrtQ * cos((theta + 2 * M_PI) / 3) - c1 / 3;
		l3 = -2 * sqrtQ * cos((theta + 4 * M_PI) / 3) - c1 / 3;
	} else {
		// this should never happen - one root only
		double e = pow(my_sqrt(-d) + fabs(R), 1. / 3.);
		if (R > 0)
			e = -e;
		l1 = (e + Q / e) - c1 / 3.;
		l2 = l1;
		l3 = l1;
		fprintf( stderr, "****** one root only (%f)**********\n",
			 l1);
		fprintf( stderr, "c123= %.10f %.10f %.10f\n",
			 c1, c2, c3);
	}

	// debug:
	assert( ! isnan( l1));
	assert( ! isnan( l2));
	assert( ! isnan( l3));
	// test solutions
	double err1 = pow(l1,3) + c1*pow(l1,2) + c2*l1 + c3;
	double err2 = pow(l2,3) + c1*pow(l2,2) + c2*l2 + c3;
	double err3 = pow(l3,3) + c1*pow(l3,2) + c2*l3 + c3;
	fprintf( stderr, "err = %f %f %f\n", err1, err2, err3);

	// sort l1, l2 and l3
	if( l1 > l2) { double t = l1; l1 = l2; l2 = t; }
	if( l2 > l3) { double t = l2; l2 = l3; l3 = t; }
	if( l1 > l2) { double t = l1; l1 = l2; l2 = t; }

	return;
}
#endif

void 
StressTensor::get_eigenvalues(
	double & l1, double & l2, double & l3)
// ----------------------------------------------------------------------
// calculate 3 eigenvalues from the stress tensor, return them in sorted
// order ( l1 < l2 < l3)
// ----------------------------------------------------------------------
{
	// the determinant of the matrix:
	//     stress - a.Identity
	// is a cubic polynomial of 'a'. Setting it to zero and solving
	// for 'a' yields 3 values --> eigenvalues of stress
	//
	// the determinant (according to  mathematica) is:
	// 
	//    Power(a,3) + Power(sxz,2)*syy - 2*sxy*sxz*syz + sxx*Power(syz,2)
	//    + Power(a,2)*(-sxx - syy - szz) + Power(sxy,2)*szz -
	//    sxx*syy*szz + a*(-Power(sxy,2) - Power(sxz,2) + sxx*syy
	//    - Power(syz,2) + (sxx + syy)*szz);
	// 
	// we extract the coefficients c1, c2, c3, so that
	// det = a^3 + c1*a^2 + c2*a + c3

	assert( ! isnan( sxx));
	assert( ! isnan( syy));
	assert( ! isnan( szz));
	assert( ! isnan( syz));
	assert( ! isnan( sxz));
	assert( ! isnan( sxy));

//	fprintf( stderr, "s.. = %.10f %.10f %.10f %.10f %.10f %.10f\n",
//		 sxx, syy, szz, syz, sxz, sxy);

	double c1 = -sxx - syy - szz;
	double c2 = -sxy*sxy - sxz*sxz + sxx*syy - syz*syz + (sxx + syy)*szz;
	double c3 = sxz*sxz*syy - 2*sxy*sxz*syz + sxx*syz*syz + 
		+ sxy*sxy*szz -	sxx*syy*szz;

	// get cubic roots
	get_cubic_roots( c1, c2, c3, l1, l2, l3);

//	fprintf( stderr, "roots = %f %f %f\n", l1, l2, l3);

	// debug:
	assert( ! isnan( l1));
	assert( ! isnan( l2));
	assert( ! isnan( l3));

	// test solutions
//	double err1 = pow(l1,3) + c1*pow(l1,2) + c2*l1 + c3;
//	double err2 = pow(l2,3) + c1*pow(l2,2) + c2*l2 + c3;
//	double err3 = pow(l3,3) + c1*pow(l3,2) + c2*l3 + c3;
//	fprintf( stderr, "err = %f %f %f\n", err1, err2, err3);

	// sort l1, l2 and l3
	if( l1 > l2) { double t = l1; l1 = l2; l2 = t; }
	if( l2 > l3) { double t = l2; l2 = l3; l3 = t; }
	if( l1 > l2) { double t = l1; l1 = l2; l2 = t; }

	return;
}

void StressTensor::print( FILE * fp, char * str, char * fmtstr)
{
	if( fmtstr == NULL) fmtstr = " %.4f";

	fprintf( fp, "%s: (", str);
	fprintf( fp, fmtstr, sxx); fprintf( fp, ",");
	fprintf( fp, fmtstr, syy); fprintf( fp, ",");
	fprintf( fp, fmtstr, szz); fprintf( fp, ",");
	fprintf( fp, fmtstr, syz); fprintf( fp, ",");
	fprintf( fp, fmtstr, sxz); fprintf( fp, ",");
	fprintf( fp, fmtstr, sxy); fprintf( fp, ")\n");
}

long
StressTensor::get_eigenvectors(
	double a1,
	Vector3d & v1, Vector3d & v2)
// ----------------------------------------------------------------------
// given an eigenvalue 'a1', find the corresponding eigenvectors
// ----------------------------------------------------------------------
{
	// method:
	//
	// - prepare a temp. matrix, equal to the tensor - a1*Identity
	// - find a null space of this temp. matrix
	//
	double tmp1[9];
	tmp1[0] = sxx - a1;
	tmp1[1] = sxy;
	tmp1[2] = sxz;
	tmp1[3] = sxy;
	tmp1[4] = syy - a1;
	tmp1[5] = syz;
	tmp1[6] = sxz;
	tmp1[7] = syz;
	tmp1[8] = szz - a1;

// TBD: this is just testing right now... scaling the matrix by
//      dividing it by its largest element
/*
	double max = fabs(tmp1[0]);
	for( long i = 1 ; i < 9 ; i ++)
		if( max < fabs(tmp1[i])) max = tmp1[i];
	if( max > 1e-6)
		for( long i = 0 ; i < 9 ; i ++)
			tmp1[i] /= max;
*/


	// create room for results
	double tmp2[9];
	// call the nullspace routine
	long nv = NullSpace( tmp1, tmp2, 1e-6, 3);
	
	// extract results (the first two rows of tmp2)
	v1.x = tmp2[0];
	v1.y = tmp2[1];
	v1.z = tmp2[2];
	v2.x = tmp2[3];
	v2.y = tmp2[4];
	v2.z = tmp2[5];

	// return number of eigenvectors
	return nv;
}

double
StressTensor::eigensystem_err(
	double a, Vector3d & v)
// ----------------------------------------------------------------------
// Reports an error of eigenvalue - eigenvector pair
//    - calculates the vector r = [s]*v - v*a
//    - reports its length squared
// ----------------------------------------------------------------------
{
//	printf( "eigensystem_err( %f, (%f,%f,%f))\n",
//		a, v.x, v.y, v.z);
	double rx = -(a*v.x) + sxx*v.x + sxy*v.y + sxz*v.z;
	double ry = sxy*v.x - a*v.y + syy*v.y + syz*v.z;
	double rz = sxz*v.x + syz*v.y + (-a + szz)*v.z;
	double res = rx*rx + ry*ry + rz*rz;
//	printf( "%.60f\n", res);

	return res;
}

double
StressTensor::diff_sq(
	StressTensor & s)
// ----------------------------------------------------------------------
// Calculates the difference between itself and another tensor, squares
// them and then sums them up
// ----------------------------------------------------------------------
{
	double sum =
		(sxx - s.sxx) * (sxx - s.sxx) +
		(syy - s.syy) * (syy - s.syy) +
		(szz - s.szz) * (szz - s.szz) +
		(syz - s.syz) * (syz - s.syz) +
		(sxz - s.sxz) * (sxz - s.sxz) +
		(sxy - s.sxy) * (sxy - s.sxy);
	return sum;
		
}

void
StressTensor::make_eigen_tensor(
	Vector3d & v)
// ----------------------------------------------------------------------
// makes a tensor from vector v by calculating [v]*transpose[v] / |v|
// ----------------------------------------------------------------------
{
	double l = v.length();
	if( l < 1e-10) {
		sxx = 0;
		syy = 0;
		szz = 0;
		sxy = 0;
		sxz = 0;
		syz = 0;
	} else {
		sxx = v.x * v.x / l;
		syy = v.y * v.y / l;
		szz = v.z * v.z / l;
		sxy = v.x * v.y / l;
		sxz = v.x * v.z / l;
		syz = v.y * v.z / l;
	}
}

char * StressTensor::to_str( char * fmtstr)
{
	if( fmtstr == NULL) fmtstr = "%.20f";
	static char res[ 4096];
	char buff[ 4096];
	sprintf( buff, "<%s,%s,%s,%s,%s,%s>",
		 fmtstr, fmtstr, fmtstr, fmtstr, fmtstr, fmtstr);
	sprintf( res, buff, sxx, syy, szz, syz, sxz, sxy);
	return res;
}

Vector3d operator * ( StressTensor & s, Vector3d & v)
// ----------------------------------------------------------------------
// calculate stress tensor and vector product
// ----------------------------------------------------------------------
{
	Vector3d res(
		s.sxx * v.x + s.sxy * v.y + s.sxz * v.z,
		s.sxy * v.x + s.syy * v.y + s.syz * v.z,
		s.sxz * v.x + s.syz * v.y + s.szz * v.z);
	return res;
}
