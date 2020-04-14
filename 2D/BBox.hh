#ifndef __BBOX_HH__
#define __BBOX_HH__

#include <math.h>

class BBox {
public:
	BBox() {
		n = 0;
		x = NULL;
		y = NULL;
	}
	~BBox() {
		if( x != NULL) free( x);
		if( y != NULL) free( y);
	}
	void add( double p_x, double p_y) {
		x = (double *) realloc( x, sizeof( double) * (n+1));
		y = (double *) realloc( y, sizeof( double) * (n+1));
		assert( x != NULL && y != NULL);
		x[ n] = p_x;
		y[ n] = p_y;
		// record min/max
		if( n == 0) {
			minx = maxx = p_x;
			miny = maxy = p_y;
		} else {
			if( p_x < minx) minx = p_x;
			if( p_x > maxx) maxx = p_x;
			if( p_y < miny) miny = p_y;
			if( p_y > maxy) maxy = p_y;
		}
		n ++;
	}
	void reset( void) {
		n = 0;
		if( x != NULL) {
			free( x); x = NULL;
			free( y); y = NULL;
		}
	}
	void print( FILE * fp) {
		fprintf( fp, "%ld", n);
		for( long i = 0 ; i < n ; i ++)
			fprintf( fp, " %f %f", x[i], y[i]);
	}
	int is_inside( double px, double py) {
		// point must be on the left-hand side of all edges
		// forming the bounding box - only then is it on the inside
		for( long i = 0 ; i < n ; i ++) {
			double x1 = x[i]; double x2 = x[(i+1) % n];
			double y1 = y[i]; double y2 = y[(i+1) % n];
			double v1x = x2 - x1;
			double v1y = y2 - y1;
			double v2x = px - x1;
			double v2y = py - y1;
			if( v1x * v2y  - v2x * v1y < 0) return 0;
		}
		return 1;
	}
	void snap_inside( double & px, double & py) {
		// if the point is outside of a bounding box, it finds
		// the closest point on the boundary of the bounding box:
		//    - consider each boundary edge
		//    - if point is on the wrong side of the edge
		//         - find iso coordinate of a point on the edge
		//           closest to this point
		//         - if the iso coordinate is outside of the edge,
		//           snap the coordinate to one of the endpoints
		//         - move the point to this iso coordinate
		for( long i = 0 ; i < n ; i ++) {
/*			double x1 = x[i]; double x2 = x[(i+1) % n];
			double y1 = y[i]; double y2 = y[(i+1) % n];
			double v1x = x2 - x1;
			double v1y = y2 - y1;
			double v2x = px - x1;
			double v2y = py - y1;
			// if the point is inside of this edge, don't
			// do anything
			if( v1x * v2y  - v2x * v1y >= 0) continue;
			// oterwise calculate the projection
			double v12 = v1x * v2x + v1y * v2y;
			double v11 = v1x * v1x + v1y * v1y;
			double k = 0;
			if( v11 > 1e-10) k = v12 / v11;
			double nvx = k * v1x;
			double nvy = k * v1y;
			px = x1 + nvx;
			py = y1 + nvy;
*/
			double x1 = x[i]; double y1 = y[i];
			double x2 = x[(i+1) % n]; double y2 = y[(i+1) % n];
			double x3 = px; double y3 = py;
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
			px = t * (x2-x1) + x1;
			py = t * (y2-y1) + y1;
		}
		
	}
	void random_point( double & x, double & y) {
		while( 1) {
			x = drand48() * (maxx-minx) + minx;
			y = drand48() * (maxy-miny) + miny;
			if( is_inside( x, y))
				return;
		}
	}
	double get_area( void) {
		// return the are of the polygon
		double area = 0;
		for( long i = 0 ; i < n ; i ++) {
			double x1 = x[i]; double y1 = y[i];
			double x2 = x[(i+1) % n]; double y2 = y[(i+1) % n];
			area += (x2-x1)*(y2+y1)/2;
		}
		return fabs( area);
	}

private:
	long n;			// number of points
	double * x, * y;	// arrays of points
public:
	double minx, miny;	// min/max values
	double maxx, maxy;	//      ||
};

#endif
