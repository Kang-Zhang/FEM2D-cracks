#ifndef __GL_DRAW_TETRA__
#define __GL_DRAW_TETRA__

#include <math.h>
#include "glDrawTriangle.hh"

void glDrawTetra( double cx, double cy, double cz, double d)
// ----------------------------------------------------------------------
// draw a tetra centered at cx, cy, cx with edge size d
// ----------------------------------------------------------------------
{
	Vector3d p1( cx + d*(-0.5), cy + d*(-0.288675), cz + d*(-0.204124));
	Vector3d p2( cx + d*( 0.5), cy + d*(-0.288675), cz + d*(-0.204124));
	Vector3d p3( cx + d*(   0), cy + d*( 0.57735 ), cz + d*(-0.204124));
	Vector3d p4( cx + d*(   0), cy + d*( 0.0     ), cz + d*( 0.612372));

	glDrawTriangle( p2, p1, p3);
	glDrawTriangle( p1, p2, p4);
	glDrawTriangle( p2, p3, p4);
	glDrawTriangle( p3, p1, p4);
}

#endif
