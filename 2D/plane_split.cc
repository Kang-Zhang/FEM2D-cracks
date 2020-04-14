#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "plane_split.hh"
#include "Vector3d.hh"

double split_edge_by_plane_raw(
	const Vector3d & norm,
	const Vector3d & p0,
	const Vector3d & p1,
	const Vector3d & p2)
// ---------------------------------------------------------------------------
// Given a plane defined by a normal <n> and a point on the plane P0,
// determine where the plane intersects line segment P1-P2
//
// The result:
//
//   - the isoparametric coordinate of the intersection (only valid if
//     between 0 and 1
// ---------------------------------------------------------------------------
{
	// calculate denominator
	double d = norm * (p2-p1);

	// calculate top
	double top = norm * (p0-p1);

	// calculate intersection point
	double res = top / d;

	return res;
}

double split_edge_by_plane(
	double nx, double ny, double nz, // normal of the plane
	double x0, double y0, double z0, // P0
	double x1, double y1, double z1, // P1
	double x2, double y2, double z2) // P2
// ---------------------------------------------------------------------------
// Given a plane defined by a normal <n> and a point on the plane P0,
// determine where the plane intersects line segment P1-P2
//
// The result:
//
//   - the sign of the result will determine where P1 lies with respect
//     to the plane: positive sign means above, negative sign means below,
//                   zero means right on the plane
//   
//   - the absolute value: if bigger than 1.0, the plane does not intersect
//     the edge. If between 0 and 1, the result indicates isoparametric
//     coordinate of the intersection, 0 being at P1 and 1 at P2
//
// ---------------------------------------------------------------------------
{
	// calculate v = p2 - p1;
	double vx = x2 - x1;
	double vy = y2 - y1;
	double vz = z2 - z1;

	// calculate the sign 'sign':
	// -------------------------------------------------------
	// if P1 is above, sign = 1, otherwise sign = -1
	double sign;
	if( nx * (x1-x0) + ny * (y1-y0) + nz * (z1-z0) > 0)
		sign = 1;
	else
		sign = -1;

	// calculate the intersection 'r'
	// -------------------------------------------------------
	double r;

	// calculate denominator
	double d = nx * vx + ny * vy + nz * vz;

	// there will be intersection only if d is not zero
	if( fabs( d) > 1e-6) {
		// calculate top
		double top = nx*(x0 - x1) + ny*(y0 - y1) + nz*(z0 - z1);

		// calculate intersection point
		double res = top / d;

		if( res >= 0 && res <= 1)
			r = res;
		else
			r = 1000;
	} else {
		r = 1000;
	}

	// return result
	return sign * r;
}

double
get_angle(
	double x1, double y1, double z1, // point P1
	double x2, double y2, double z2, // point P2
	double x3, double y3, double z3  // point P3
	)
// ----------------------------------------------------------------------
// - calculate the angle between lines (P1-P2) and (P1-P3) in degrees
// - the angle is the absolute value!
// ----------------------------------------------------------------------
{
	double vx = x2 - x1;
	double vy = y2 - y1;
	double vz = z2 - z1;
	double dv = sqrt( vx * vx + vy * vy + vz * vz);
	vx /= dv; vy /= dv; vz /= dv;
	double ux = x3 - x1;
	double uy = y3 - y1;
	double uz = z3 - z1;
	double du = sqrt( ux * ux + uy * uy + uz * uz);
	ux /= du; uy /= du; uz /= du;

	double cosalpha = vx * ux + vy * uy + vz * uz;
	if( cosalpha < -1) cosalpha = -1;
	if( cosalpha > 1) cosalpha = 1;

	double alpha = fabs( acos( cosalpha) * 180.0 / M_PI);

	return alpha;
}

int where_is_triangle_wrt_plane_angle(
	double nx, double ny, double nz, // plane normal N
	double x0, double y0, double z0, // P0
	double x1, double y1, double z1, // P1
	double x2, double y2, double z2, // P2
	double min_angle,	         // min. split angle
	double & s                       // where to return the resulting
				         // coordinate
	)
// ----------------------------------------------------------------------
// calculates whether and where a plane intersects an edge
// - plane is defined by a normal <nx,ny,nz> and a point <x0,y0,z0>
// - edge is defined by <x1,y1,z1> and <x2,y2,z2>
// - if the intersection exists, but the resulting angle at P0 with
//   (P1 or P2) and intersection point is less than a min. angle then
//   we return no intersection
// - in 's' we return 2 results:
//    - The sign of 's' determines:
//          - sign(s) < 0 : P1 is below the plane
//          - sign(s) > 0 : P1 is above the plane
//          - sign(s) = 0 : P1 is on the plane
//    - the absolute value of 's' determines the iso-parametric coordinate
//      of the intersection (provided it exists). fabs(s) is always between
//      0 and 1. '0' denotes intersection at P1, while '1' represents
//      intersection at P2.
// - return value:
//      1 = the whole triangle is above the plane, or whatever is below
//          the plane forms an angle smaller than min_angle
//     -1 = the while triangle is below the plane, or whatever is above
//          the plane forms an angle smaller than min_angle
//      0 = the triangle is split by the plane, and both angles formed
//          by the split are bigger than min_angle
// - remark:
//      - if both angles formed by the split are smaller than min_angle,
//        then the result is 'above=1' if bigger angle formed above, and
//        'below=-1' otherwise
// ----------------------------------------------------------------------
{
	// find the isoparametric coordinate on the edge p1-p2 where it is
	// intersected by the stress plane
	s = split_edge_by_plane(
		nx, ny, nz,
		x0, y0, z0,
		x1, y1, z1,
		x2, y2, z2);

	fprintf( stderr, "slit edge returned = %.20f\n", s);

	// if the whole edge is above or below, declare result immediately
// 	if( fabs(s) > 1)
// 		if( s > 0) return 1;
// 		else if( s < 0) return -1;
	if( s > 1) return 1;
	if( s < -1) return -1;

	// calculate where the intersection point would be
	double ix = x1 + (x2-x1) * fabs(s);
	double iy = y1 + (y2-y1) * fabs(s);
	double iz = z1 + (z2-z1) * fabs(s);

	// calculate the two angles at node 0 we would get if we split
	// the triangle
	double alpha1 = get_angle( x0, y0, z0,
				   ix, iy, iz,
				   x1, y1, z1);
	double alpha2 = get_angle( x0, y0, z0,
				   ix, iy, iz,
				   x2, y2, z2);

	fprintf( stderr, "alpha1 = %f\n", alpha1);
	fprintf( stderr, "alpha2 = %f\n", alpha2);

	// if both angles are smaller than min_angle, return verdict based
	// on which angle is bigger: above or below the plane
	if( alpha1 < min_angle && alpha2 < min_angle) {
		if( alpha1 >= alpha2)
			// if P1's angle is larger, return position of P1
			return (s < 0) ? -1 : 1;
		else
			// if P2's angle is larger, return position of P2
			return (s < 0) ? 1 : -1;
	}

	fprintf( stderr, "aaa 1\n");

	// if alpha1 is smaller than min_angle, then return position of P2
	if( alpha1 < min_angle)
		return (s < 0) ? 1 : -1;

	fprintf( stderr, "aaa 2\n");

	// if alpha2 is smaller than min_angle, then return position of P1
	if( alpha2 < min_angle)
		return (s < 0) ? -1 : 1;

	fprintf( stderr, "aaa 3\n");

	// both angles are bigger than min_angle, so return 'intersection'
	return 0;
}

int where_is_triangle_wrt_plane_length(
	double nx, double ny, double nz, // plane normal
	double x0, double y0, double z0, // vertex through which plane passes
	double x1, double y1, double z1, // second vertex
	double x2, double y2, double z2, // second vertex
	double min_length,	         // min. split edge length
	double & s                       // returns the isoparametric
					 // coordinate of the split,
					 // positive if P1 is above and
					 // negative is P1 is below
	)
{
	// find the isoparametric coordinate on the edge p1-p2 where it is
	// intersected by the stress plane
	s = split_edge_by_plane(
		nx, ny, nz,
		x0, y0, z0,
		x1, y1, z1,
		x2, y2, z2);

	// if the whole edge is above or below, declare result immediately
	if( fabs(s) > 1)
		if( s > 0) return 1;
		else if( s < 0) return -1;

	// calculate the lengths of the split edge
	double l = sqrt(
		(x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
	// calculate lengths of the split angles
	double l1 = fabs( s) * l;
	double l2 = (1-fabs(s)) * l;
	
	// if the min. length is not satisfied, declare the whole element e
	// either below or above, depending on which length is shorter and
	// which of p1 and p2 is above and below
	if( l1 < min_length || l2 < min_length) {
		if( ( l1 >= l2 && s >= 0) ||
		    ( l2 >= l1 && s <= 0))
		{
			return 1;
		} else {
			return -1;
		}
	} else {
		return 0;
	}
}

#ifdef SELFTEST

int main( int, char **)
{
	fprintf( stderr, "Test1: %f\n",
		 split_edge_by_plane_raw(
			 Vector3d( -1, 1, 0),
			 Vector3d( 0, 0, 0),
			 Vector3d( 0, 1, 0),
			 Vector3d( 2, 0, 0))
		);

	fprintf( stderr, "Test1: %f\n",
		 split_edge_by_plane_raw(
			 Vector3d( 0, 0, 2),
			 Vector3d( 0, 0, 1),
			 Vector3d( 1, 1, 2),
			 Vector3d( -1, -1, -1)
			 )
		);
/*
	printf( "Test1: %f\n",
		split_edge_by_plane(
			0, 1, 0,
			0, 2, 0,
			3, 0, 0,
			5, 3, 0));
	printf( "Test1: %f\n",
		split_edge_by_plane(
			0, 1, 0,
			0, 2, 0,
			5, 3, 0,
			3, 0, 0));
	printf( "Test1: %f\n",
		split_edge_by_plane(
			0, 1, 0,
			0, 2, 0,
			5, 3, 0,
			1, 5, 0));
	printf( "Test1: %f\n",
		split_edge_by_plane(
			0, 1, 0,
			0, 2, 0,
			1, 5, 0,
			5, 3, 0));

	int n = 1;
	{
		double s;
		int where = where_is_triangle_wrt_plane_angle(
			-1, 1, 0,
			1, 1, 0,
			4, 1, 0,
			4, 2, 0,
			0,
			s);
		printf( "Triangle%d: %d, %f\n", n++, where, s);
	}

	{
		double s;
		int where = where_is_triangle_wrt_plane_angle(
			1, -1, 0,
			1, 1, 0,
			4, 1, 0,
			4, 2, 0,
			0,
			s);
		printf( "Triangle%d: %d, %f\n", n++, where, s);
	}
	{
		double s;
		int where = where_is_triangle_wrt_plane_angle(
			-1, 3, 0,
			1, 1, 0,
			4, 1, 0,
			4, 2, 0,
			0,
			s);
		printf( "Triangle%d: %d, %f\n", n++, where, s);
	}
	{
		double s;
		int where = where_is_triangle_wrt_plane_angle(
			0, 1, 0,
			1, 1, 0,
			4, 1, 0,
			4, 2, 0,
			0,
			s);
		printf( "Triangle%d: %d, %f\n", n++, where, s);
	}
	{
		double s;
		int where = where_is_triangle_wrt_plane_angle(
			-1, 4, 0,
			1, 1, 0,
			4, 1, 0,
			4, 2, 0,
			0,
			s);
		printf( "Triangle%d: %d, %f\n", n++, where, s);
	}
	{
		double s;
		int where = where_is_triangle_wrt_plane_angle(
			-1, 400, 0,
			1, 1, 0,
			4, 1, 0,
			4, 2, 0,
			0,
			s);
		printf( "Triangle%d: %d, %f\n", n++, where, s);
	}
	{
		double s;
		int where = where_is_triangle_wrt_plane_angle(
			-1, 400, 0,
			1, 1, 0,
			4, 1, 0,
			4, 2, 0,
			10,
			s);
		printf( "Triangle%d: %d, %f\n", n++, where, s);
	}
	{
		double s;
		int where = where_is_triangle_wrt_plane_angle(
			1, -400, 0,
			1, 1, 0,
			4, 1, 0,
			4, 2, 0,
			10,
			s);
		printf( "Triangle%d: %d, %f\n", n++, where, s);
	}
	{
		double s;
		int where = where_is_triangle_wrt_plane_angle(
			-1, 3.1, 0,
			1, 1, 0,
			4, 1, 0,
			4, 2, 0,
			10,
			s);
		printf( "Triangle%d: %d, %f\n", n++, where, s);
	}
	{
		double s;
		int where = where_is_triangle_wrt_plane_angle(
			1, -3.1, 0,
			1, 1, 0,
			4, 1, 0,
			4, 2, 0,
			10,
			s);
		printf( "Triangle%d: %d, %f\n", n++, where, s);
	}
*/
}

#endif
