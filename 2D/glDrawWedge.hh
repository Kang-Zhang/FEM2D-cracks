#ifndef __GL_DRAW_WEDGE__
#define __GL_DRAW_WEDGE__

#include "glDrawTriangle.hh"

void glDrawQuad( const Vector3d & p1,
		 const Vector3d & p2,
		 const Vector3d & p3,
		 const Vector3d & p4,
		 const long ns = 4)
{
	// prepare an array of interpolated values
	Vector3d a[ns+1][ns+1];
	double dt = 1.0/ns;
	double ds = 1.0/ns;
	for( long i = 0 ; i <= ns ; i ++) {
		for( long j = 0 ; j <= ns ; j ++) {
			double t = dt * i;
			double s = ds * j;
			a[i][j] =
				p1 * (1-t)*(1-s) +
				p2 * (  t)*(1-s) +
				p3 * (  t)*(  s) +
				p4 * (1-t)*(  s);
		}
	}

	// draw the sub-quads
	for( long i = 0 ; i < ns ; i ++) {
		for( long j = 0 ; j < ns ; j ++) {
			Vector3d pp1 = a[i][j];
			Vector3d pp2 = a[i+1][j];
			Vector3d pp3 = a[i+1][j+1];
			Vector3d pp4 = a[i][j+1];

			Vector3d c = (pp1+pp2+pp3+pp4) * (1.0/4);
			
			glDrawTriangle( pp1, pp2, c);
			glDrawTriangle( pp2, pp3, c);
			glDrawTriangle( pp3, pp4, c);
			glDrawTriangle( pp4, pp1, c);
		}
	}
}

void glDrawWedge( Vector3d (&p)[6])
// ----------------------------------------------------------------------
// draw a cube centered at cx, cy, cx with edge size d
// ----------------------------------------------------------------------
{
	glDrawTriangle( p[0], p[1], p[2]); // top face
	glDrawTriangle( p[5], p[4], p[3]); // bottom face
	glDrawQuad(p[0], p[1], p[4], p[3], 1);
	glDrawQuad(p[1], p[2], p[5], p[4], 1);
	glDrawQuad(p[2], p[0], p[3], p[5], 1);
// 	glDrawTriangle( p[3], p[4], p[0]); // front face
// 	glDrawTriangle( p[0], p[4], p[1]);
// 	glDrawTriangle( p[4], p[5], p[1]); // right face
// 	glDrawTriangle( p[1], p[5], p[2]);
// 	glDrawTriangle( p[5], p[3], p[0]); // left face
// 	glDrawTriangle( p[5], p[0], p[2]);
}

#endif
