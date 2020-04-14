#ifndef __GL_DRAW_TRIANGLE__
#define __GL_DRAW_TRIANGLE__

#include "Vector3d.hh"

void glDrawTriangle( Vector3d & p1, Vector3d & p2, Vector3d & p3)
// ----------------------------------------------------------------------
// draw a cube centered at cx, cy, cx with edge size d
// ----------------------------------------------------------------------
{
	// figure out the normal
	Vector3d n = cross_product( p2-p1, p3-p1);
	n.normalize();
	glNormal3d( n.x, n.y, n.z);
	glBegin( GL_TRIANGLES);
	glVertex3d( p1.x, p1.y, p1.z);
	glVertex3d( p2.x, p2.y, p2.z);
	glVertex3d( p3.x, p3.y, p3.z);
	glEnd();
}

#endif
