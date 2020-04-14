#include <stdio.h>
#include <stdlib.h>
#include "Vector3d.hh"
#include "projected_area.hh"

double
projected_area(
	Vector3d & pn,
	TriangleList & tlist)
// ----------------------------------------------------------------------
// calculates the area of a closed manifold projected onto a plane defined
// by the normal pn
// ----------------------------------------------------------------------
{
	double area = 0.0;
	for( long i = 0 ; i < tlist.n_triangles ; i ++) {
		TriangleList::Triangle & t = tlist.triangles[i];
		// calculate the surface normal of this triangle
		Vector3d tn = cross_product(t.p[1]-t.p[0],t.p[2]-t.p[0]);
		// if the scalar product of tn and pn is negative, this
		// triangle is facing away, so ignor it
		if( tn*pn < 0) continue;
		// otherwise add the area of this triangle multiplied
		// by the cosine of tn,pn to the total area
		double cosangle = tn*pn / (tn.length()*pn.length());
		area += (tn.length()/2.0) * cosangle;
	}
	return area;
}
