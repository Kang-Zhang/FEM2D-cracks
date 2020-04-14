#include <assert.h>
#include "Vector3d.hh"

/*
Vector3d cross_product( Vector3d & v1, Vector3d & v2)
{
	return Vector3d(
		-(v2.y*v1.z) + v1.y*v2.z,
		v2.x*v1.z - v1.x*v2.z,
		-(v2.x*v1.y) + v1.x*v2.y);
}
*/

Vector3d cross_product( const Vector3d & v1, const Vector3d & v2)
{
	return Vector3d(
		-(v2.y*v1.z) + v1.y*v2.z,
		v2.x*v1.z - v1.x*v2.z,
		-(v2.x*v1.y) + v1.x*v2.y);
}

Vector3d operator+( const Vector3d & v1, const Vector3d & v2)
{
	return Vector3d( v1.x + v2.x,
			 v1.y + v2.y,
			 v1.z + v2.z);
}

Vector3d operator-( const Vector3d & v1, const Vector3d & v2)
{
	return Vector3d( v1.x - v2.x,
			 v1.y - v2.y,
			 v1.z - v2.z);
}

Vector3d operator*( const Vector3d & v1, double a)
{
	return Vector3d( v1.x * a, v1.y * a, v1.z * a);
}

Vector3d operator*( double a, const Vector3d & v1)
{
	return Vector3d( v1.x * a, v1.y * a, v1.z * a);
}

Vector3d operator/( const Vector3d & v1, double a)
{
	return Vector3d( v1.x / a, v1.y / a, v1.z / a);
}

Vector3d project( const Vector3d & v1, const Vector3d & v2)
{
		double bx = v1.x;
		double by = v1.y;
		double bz = v1.z;
		double ax = v2.x;
		double ay = v2.y;
		double az = v2.z;
		double AB = ax*bx + ay*by + az*bz;
		double AA = ax*ax + ay*ay + az*az;

		Vector3d res( ax * AB / AA, ay * AB / AA, az * AB / AA);
 		if( isnan( res.x) || isnan( res.y) || isnan( res.z)) {
 			res.x = 0;
 			res.y = 0;
 			res.z = 0;
 		}

		return res;
}

// scalar product
double operator*( const Vector3d & v1, const Vector3d & v2)
{
	return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

double Vector3d::get_angle( const Vector3d & v1, const Vector3d & v2)
{
	Vector3d u1 = v1; u1.normalize();
	Vector3d u2 = v2; u2.normalize();
	return fabs( acos( u1*u2)) * 180 / M_PI;
}
