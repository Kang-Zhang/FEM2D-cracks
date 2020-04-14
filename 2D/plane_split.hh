#ifndef __PLANE_SPLIT_HH__
#define __PLANE_SPLIT_HH__

#include "Vector3d.hh"

double split_edge_by_plane_raw(
	const Vector3d & norm,
	const Vector3d & p0,
	const Vector3d & p1,
	const Vector3d & p2);

double
get_angle(
	double x1, double y1, double z1, // point P1
	double x2, double y2, double z2, // point P2
	double x3, double y3, double z3  // point P3
	);

double split_edge_by_plane(
	double nx, double ny, double nz,
	double x0, double y0, double z0,
	double x1, double y1, double z1,
	double x2, double y2, double z2);

int where_is_triangle_wrt_plane_angle(
	double nx, double ny, double nz, // plane normal
	double x0, double y0, double z0, // vertex through which plane passes
	double x1, double y1, double z1, // second vertex
	double x2, double y2, double z2, // second vertex
	double min_angle,	         // min. split angle
	double & s                       // returns the isoparametric
					 // coordinate of the split,
					 // positive if P1 is above and
					 // negative is P1 is below
	);
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
	);
#endif
