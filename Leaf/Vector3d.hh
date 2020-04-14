#ifndef __VECTOR3D_HH__
#define __VECTOR3D_HH__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

class Vector3d {

public:

	double x, y, z;

	Vector3d() {
		x = 0.0; y = 0.0; z = 0.0;
	}

	Vector3d( double px, double py, double pz) {
		x = px; y = py; z = pz;
	}

	void assign( Vector3d & v) {
		x = v.x; y = v.y; z = v.z;
	}

	void print( FILE * fp) {
		fprintf( fp, "<%f,%f,%f>", x, y, z);
	}

	void cross_product( Vector3d & v) {
		double tx = y*v.z - z*v.y;
		double ty = z*v.x - x*v.z;
		double tz = x*v.y - y*v.x;
		x = tx;
		y = ty;
		z = tz;
	}

	void scale( double a) {
		x *= a;
		y *= a;
		z *= a;
	}

	double scalar_product( Vector3d & v) {
		return x*v.x + y*v.y + z*v.z;
	}

	double cos_product( Vector3d & v) {
		return scalar_product( v) / ( length() * v.length());
	}
	
	double length( void) {
		return sqrt( x*x + y*y + z*z);
	}
	
	void normalize( void) {
		double l = sqrt( x*x + y*y + z*z);
		if( fabs( l) > 1e-10) {
			x /= l; y /= l; z /= l;
		}
	}
	
	void add( Vector3d & v) {
		x += v.x;
		y += v.y;
		z += v.z;
	}

	void sub( Vector3d & v) {
		x -= v.x;
		y -= v.y;
		z -= v.z;
	}

	void set( double px, double py, double pz) {
		x = px; y = py; z = pz;
	}

	void get( double & px, double & py, double & pz) {
		px = x; py = y; pz = z;
	}

	static double angle( Vector3d v1, Vector3d v2) {
		v1.normalize();
		v2.normalize();
		double prod = v1.scalar_product( v2);
		return fabs( acos( prod)) * 180 / M_PI;
	}

	Vector3d project_onto( const Vector3d & v) {
		double bx = this-> x;
		double by = this-> y;
		double bz = this-> z;
		double ax = v.x;
		double ay = v.y;
		double az = v.z;
		double AB = ax*bx + ay*by + az*bz;
		double AA = ax*ax + ay*ay + az*az;
		Vector3d res( ax * AB / AA, ay * AB / AA, az * AB / AA);
		return res;
	}

	char * to_str( char * fmtstr = NULL) {
		if( fmtstr == NULL) fmtstr = "%f";
		static char res[ 4096];
		char buff[ 4096];
		sprintf( buff, "<%s,%s,%s>", fmtstr, fmtstr, fmtstr);
		sprintf( res, buff, x, y, z);
		return res;
	}

	Vector3d & operator= ( const Vector3d & v) {
		x = v.x;
		y = v.y;
		z = v.z;
		return * this;
	}

	int is_degenerate( void) const {
		return isnan(x) || isnan(y) || isnan(z);
	}

	int is_zero( void) const {
		return (x == 0) && (y == 0) && ( z == 0);
	}

	static double get_angle( const Vector3d & v1, const Vector3d & v2);
};


// expression templates
// ----------------------------------------------------------------------
// Vector3d cross_product( Vector3d & v1, Vector3d & v2);
Vector3d cross_product( const Vector3d & v1, const Vector3d & v2);
// vector addition
Vector3d operator+( const Vector3d & v1, const Vector3d & v2);
// vector subtraction
Vector3d operator-( const Vector3d & v1, const Vector3d & v2);
// scalar product
double operator*( const Vector3d & v1, const Vector3d & v2);
// multiplication by a scalar
Vector3d operator*( const Vector3d & v1, double a);
Vector3d operator*( double a, const Vector3d & v1);
// division by a scalar
Vector3d operator/( const Vector3d & v1, double a);
// project v1 onto v2
Vector3d project( const Vector3d & v1, const Vector3d & v2);

#endif
