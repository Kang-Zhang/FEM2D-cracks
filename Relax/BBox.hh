#ifndef __BBOX_HH__
#define __BBOX_HH__

#include <cmath>
#include <vector>
#include "Vector.hpp"

using VectLib::V2D;

class BBox {
public:
    BBox() {
    }
    ~BBox() {
    }
    const std::vector <V2D> & pts () { return _pts; }
    void add (const V2D & p) {
	_pts.push_back (p);
	// record min/max
	if (_pts.size() == 1) {
	    minx = maxx = p.x();
	    miny = maxy = p.y();
	} else {
	    if( p.x() < minx) minx = p.x();
	    if( p.x() > maxx) maxx = p.x();
	    if( p.y() < miny) miny = p.y();
	    if( p.y() > maxy) maxy = p.y();
	}
    }
/*
    bool is_inside (const V2D & p) {
	// point must be on the left-hand side of all edges
	// forming the bounding box - only then is it on the inside
	for (size_t i = 0 ; i < _pts.size() ; i ++) {
	    V2D & p1 = _pts[i];
	    V2D & p2 = _pts[(i+1) % _pts.size()];
	    V2D v1 = p2 - p1;
	    V2D v2 = p - p1;
	    if (v1.x() * v2.y() - v2.x() * v1.y() < 0) return false;
	}
	return true;
    }
*/
    double dist_to_boundary (const V2D & p) {
	double min = -1;
	for (size_t i = 0 ; i < _pts.size() ; i ++) {
	    V2D & p1 = _pts[i];
	    V2D & p2 = _pts[(i+1) % _pts.size()];
	    double d = VectLib::dist_to_line_seg (p, p1, p2);
//	    std::cerr << "\t" << d << "\n";
	    if (min < 0 || d < min) min = d;
/*
	    V2D v12 = p2 - p1;
	    double u = (p1-p) * v12 / (v12 * v12);
	    if (u < 0);
*/
	}
	return min;
    }
    bool between (double x, double x1, double x2) {
	if (x1 <= x && x <= x2) return true;
	if (x2 <= x && x <= x1) return true;
	return false;
    }
    bool is_inside (const V2D & p) {
//	std::cerr << "is_inside (" << p.x() << "," << p.y() << ")\n";
	// determines whether a point is inside a polygon or not.
	// This is achieved by counting the number of intersections of a horizontal
	// ray starting at point p continuing x+ with the edges of the polygon
	bool res = false;
	for (size_t i = 0 ; i < _pts.size() ; i ++) {
//	    std::cerr << "\t" << i << " ";
	    V2D & p1 = _pts[i];
	    V2D & p2 = _pts[(i+1)%_pts.size()];
	    // special case - point directly located on a horizontal edge
	    if (p1.y() == p2.y() && p1.y() == p.y() && between (p.x(), p1.x(), p2.x())) {
//		std::cerr << "special case applies\n";
		return true;
	    }
	    // if edge completely above or below, no intersection
	    if (std::min(p1.y(),p2.y()) >= p.y()) {
//		std::cerr << "completely above\n";
		continue;
	    }
	    if (std::max(p1.y(),p2.y()) < p.y()) {
//		std::cerr << "completely below.\n";
		continue;
	    }
	    // the full ray intersects somewhere at point x:
	    double x = ((p.y() - p1.y()) * (p2.x() - p1.x())) / (p2.y()-p1.y()) + p1.x();
//	    std::cerr << "x = " << x << "\n";
	    if (x >= p.x()) res = ! res;
	}
	return res;
    }

    V2D snap_inside (const V2D & p) {
	// if the point is outside of a bounding box, it finds
	// the closest point on the boundary of the bounding box:
	//    - consider each boundary edge
	//    - if point is on the wrong side of the edge
	//         - find iso coordinate of a point on the edge
	//           closest to this point
	//         - if the iso coordinate is outside of the edge,
	//           snap the coordinate to one of the endpoints
	//         - move the point to this iso coordinate
	if (is_inside (p)) return p;
	V2D res = _pts[0];
	for (size_t i = 0 ; i < _pts.size() ; i ++) {
	    double x1 = _pts[i].x(); double y1 = _pts[i].y();
	    double x2 = _pts[(i+1) % _pts.size()].x(); double y2 = _pts[(i+1) % _pts.size()].y();
	    double x3 = p.x(); double y3 = p.y();
	    // if the point is inside of this edge, don't
	    // do anything
	    double v1x = x2 - x1;
	    double v1y = y2 - y1;
	    double v2x = x3 - x1;
	    double v2y = y3 - y1;
	    if( v1x * v2y  - v2x * v1y >= 0) continue;
	    // otherwise snap...
	    double den = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
	    double nom = x1*x1+x2*x3-x1*(x2+x3)+(y1-y2)*(y1-y3);
	    double t = 0;
	    if( den > 1e-10) t = nom / den;
	    if( t < 0) t = 0;
	    if( t > 1) t = 1;
	    V2D tmp (t * (x2-x1) + x1, t * (y2-y1) + y1);
	    V2D tmpdiff = tmp - p;
	    V2D resdiff = res - p;
	    if (tmpdiff * tmpdiff < resdiff * resdiff) res = tmp;
	}
	return res;
    }

    V2D random_point() {
	return _pts[0];
	while( 1) {
	    V2D p (drand48() * (maxx-minx) + minx,
		   drand48() * (maxy-miny) + miny);
	    if (is_inside (p)) return p;
	}
    }
    double get_area( void) {
	// return the are of the polygon
	double area = 0;
	for (size_t i = 0 ; i < _pts.size() ; i ++) {
	    double x1 = _pts[i].x(); double y1 = _pts[i].y();
	    double x2 = _pts[(i+1) % _pts.size()].x(); double y2 = _pts[(i+1) % _pts.size()].y();
	    area += (x2-x1)*(y2+y1)/2;
	}
	return fabs( area);
    }
    
private:
    std::vector <V2D> _pts;
public:
    double minx, miny;	// min/max values
    double maxx, maxy;	//      ||
};

#endif
    
