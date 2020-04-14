#ifndef __PROJECTED_AREA_HH__
#define __PROJECTED_AREA_HH__

#include <assert.h>
#include "Vector3d.hh"

class Plane {
public:
	Vector3d norm;
	Vector3d pt;

	void project_onto( Vector3d & p) const {
		Vector3d v = p-pt;
		// project v onto normal
		Vector3d pv = project( v, norm);
		// subtract the projected vector from p to get result
		p = v - pv + pt;
	}
};

class TriangleList {
public:
	class Triangle {
	public:
		Vector3d p[3];
	};
	long n_triangles;
	long n_alloc;
	Triangle * triangles;

	TriangleList( long p_n_alloc) {
		n_triangles = 0;
		n_alloc = p_n_alloc;
		triangles = new Triangle[ n_alloc];
	}

	~TriangleList() {
		delete [] triangles;
	}

	void add( Vector3d & p1, Vector3d & p2, Vector3d & p3) {
		assert( n_triangles < n_alloc);
		triangles[n_triangles].p[0] = p1;
		triangles[n_triangles].p[1] = p2;
		triangles[n_triangles].p[2] = p3;
		n_triangles ++;
	}
		
};

double
projected_area(
	Vector3d & pn,		// area normal
	TriangleList & tlist);

#endif
