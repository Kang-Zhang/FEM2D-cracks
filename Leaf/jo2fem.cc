#include <cstdio>
#include <cstdlib>
#include <vector>
#include <utility>
#include <GL/gl.h>
#include <OGLUIall.h>
#include <cstring>
#include <stdarg.h>
#include "Vector3d.hh"
#include "Tokenizer.hh"
#include "die.h"
#include "gauss_jordan.hh"
#include "glDrawCube.h"
#include "Materials.hh"
#include "TriangulatorDelaunay2D.hh"
#include "StressTensor.hh"

// global variables
typedef StressTensor GrowthTensor;
typedef TriangulatorDelaunay2D::Triangle Triangle;
typedef std::vector<Triangle> TList;
static double curr_time = 0;
static long time_line_height = 20;
static int debug = 0;
bool draw_grid = true;
bool draw_lines = true;
bool draw_cubes = false;
bool draw_numbers = false;
bool draw_autozoom = true;
bool draw_antialias = true;
bool draw_tensors = false;
bool draw_growth_rates = false;

class Point {
public:
	double x,y;
	bool is_new;

	Point() { x = y = 0; }
	Point( double px, double py) { x = px; y = py; }
	int load( Tokenizer & tok);
	bool member_of( const std::vector<Point> pts) const;
	void print( FILE * fp) const;
};
typedef std::vector<Point> PList;

int Point::load( Tokenizer & tok)
{
	// read in the id & the coordinates
	(void) tok.read_double( "id");
	x = tok.read_double( "x");
	y = tok.read_double( "y");
	is_new = 0;

	return tok.error();
}

void Point::print( FILE * fp) const
{
	fprintf( fp, "%10f %10f\n", x, y);
}

class Measurement {
public:
	// time of the measurement
	double time;
	// the set of points
	PList pts;

	// parser for the measurement
	int load( Tokenizer & tok);

	// output method
	void print( FILE * fp);

	// request a triangulation of this measurement
	TList triangulation( void);
};

TList Measurement::triangulation( void)
{
	TriangulatorDelaunay2D triangulator;
	for( size_t i = 0 ; i < pts.size() ; i ++)
		triangulator.add_point( pts[i].x, pts[i].y);
	return triangulator.triangulate();
}

int Measurement::load( Tokenizer & tok)
{
	// get the time of the measurement
	time = tok.read_double( "time");

	// get the number of points
	long n_pts = tok.read_long( "num. of points");

	// read in the points
	for( long i = 0 ; i < n_pts ; i ++) {
		Point p;
		if( debug)
			fprintf( stderr, "\t\t- loading point %ld\n", i+1);
		p.load( tok);
		pts.push_back( p);
	}

	return tok.error();
}

void Measurement::print( FILE * fp)
{
	fprintf( fp, "%f # time of this measurement\n", time);
	fprintf( fp, "%ld # number of points for this measurement\n",
		 long( pts.size()));
	
	for( size_t i = 0 ; i < pts.size() ; i ++)
		pts[i].print( fp);
}

class Data {
public:
	std::vector<Measurement> mes;
	TList tlist;

	int load( Tokenizer & tok);
	void print( FILE * fp);
	void center( void);
	void extrapolate( void);
	void extrapolate2( void);
	void extrapolate3( void);
	void get_coord( size_t ind, double t,
			double & x, double & y, int & is_new);
	void get_triangle_pos (Triangle & t,
			       double time,
			       Vector3d (& pos)[3]);
	Vector3d get_triangle_center (Triangle & t, double time);
	double get_triangle_growth_rate (Triangle & t, double time, double dt);
	GrowthTensor get_growth_tensor (Triangle & t, double time, double dt);

	double min_growth_rate;
	double max_growth_rate;

	void calculate_growth_rates ();
	void triangulate ();
};

void Data::get_coord( size_t ind, double t,
		      double & x, double & y, int & is_new)
{
	assert( ind >= 0 || ind < mes.size());
	// adjust t (in case it is out of bounds)
	t = std::max<double> (t , 0.0);
	t = std::min<double> (t , double (mes.size ()));
	// get the measurements on the left and right
	long li = long(t);
	long ri = li+1;
	if( ri >= long(mes.size())) ri = mes.size()-1;
	// get the coordinates of both points
	double x1 = mes[li].pts[ind].x;
	double y1 = mes[li].pts[ind].y;
	double x2 = mes[ri].pts[ind].x;
	double y2 = mes[ri].pts[ind].y;
	// interpolate
	double tf = t - li;
	x = (1-tf) * x1 + (tf) * x2;
	y = (1-tf) * y1 + (tf) * y2;
	is_new = mes[li].pts[ind].is_new;
}

void Data::get_triangle_pos (Triangle & t,
			     double time,
			     Vector3d (& pos)[3])
{
	int is_new; // dummy
	get_coord( t.p1, time, pos[0].x, pos[0].y, is_new);
	get_coord( t.p2, time, pos[1].x, pos[1].y, is_new);
	get_coord( t.p3, time, pos[2].x, pos[2].y, is_new);
}

Vector3d Data::get_triangle_center (Triangle & t, double time)
{
	Vector3d pos[3];
	get_triangle_pos (t, time, pos);
	Vector3d res = (pos[0] + pos[1] + pos[2]) / 3.0;
	return res;
}

double Data::get_triangle_growth_rate (Triangle & t, double time, double dt)
{
	double t0 = time;
	double t1 = t0 + dt;
	if( t1 > mes.size() - 1)
	{
		t1 = mes.size() - 1;
		t0 = t1 - dt;
	}
	Vector3d pos0[3]; get_triangle_pos (t, t0, pos0);
	Vector3d pos1[3]; get_triangle_pos (t, t1, pos1);

	double area0 = cross_product(pos0[1]-pos0[0],pos0[2]-pos0[0])
		.length()/2;
	double area1 = cross_product(pos1[1]-pos1[0],pos1[2]-pos1[0])
		.length()/2;
	double rate = (area1-area0)/area1/dt;
	return rate;
}

void Data::triangulate ()
{
	// triangulate the last result
	tlist = mes[mes.size () - 1].triangulation();
}

void Data::calculate_growth_rates ()
{
	double dt = 0.05;
	// initial min/max
	min_growth_rate = max_growth_rate =
		get_triangle_growth_rate( tlist[0], 0, dt);
	for (double time = 0.0 ; time < mes.size() - 1 ; time += dt)
	{
		for( size_t i = 0 ; i < tlist.size() ; i ++) {
			double r = get_triangle_growth_rate (
				tlist[i], time, dt);
			max_growth_rate = std::max<double> (max_growth_rate,r);
			min_growth_rate = std::min<double> (min_growth_rate,r);
		}
	}
	max_growth_rate = std::min<double> (max_growth_rate, 3.0);
	min_growth_rate = std::max<double> (min_growth_rate, -1.0);
	min_growth_rate = std::min<double> (min_growth_rate, -0.1);
}

GrowthTensor Data::get_growth_tensor (Triangle & t, double time, double dt)
{
	int is_new;
	double t0 = time;
	double t1 = t0 + dt;
	if (t1 > mes.size() - 1)
	{
		t1 = mes.size() - 1;
		t0 = t1 - dt;
	}
	// calculate the position Q
	Matrix Q (2,3);
	get_coord (t.p1, t0, Q(0,0), Q(1,0), is_new);
	get_coord (t.p2, t0, Q(0,1), Q(1,1), is_new);
	get_coord (t.p3, t0, Q(0,2), Q(1,2), is_new);
	// calculate the displacement q
	Matrix q (2,3);
	get_coord (t.p1, t1, q(0,0), q(1,0), is_new);
	get_coord (t.p2, t1, q(0,1), q(1,1), is_new);
	get_coord (t.p3, t1, q(0,2), q(1,2), is_new);
	q.sub (Q);
	// now get displacement per unit time
	q.multiply_by_scalar (1/dt);
	// calculate Nen
	double valsNen[] = {1,0,-1,0,1,-1};
	Matrix Nen (2,3, valsNen);
	// calculate J
	double valsJ[] = {Q(0,0)-Q(0,2),Q(1,0)-Q(1,2),
			  Q(0,1)-Q(0,2),Q(1,1)-Q(1,2)};
	Matrix J (2,2,valsJ);
	// calculate the Nxy
	Matrix Nxy = inverse (J) * Nen;
	// calculate the first component of the strain tensor
	double ex =
		Nxy(0,0) * q(0,0) +
		Nxy(0,1) * q(0,1) +
		Nxy(0,2) * q(0,2);
	double ey = 
		Nxy(1,0) * q(1,0) +
		Nxy(1,1) * q(1,1) +
		Nxy(1,2) * q(1,2);
	double exy =
		Nxy(1,0) * q(0,0) +
		Nxy(1,1) * q(0,1) +
		Nxy(1,2) * q(0,2) +
		Nxy(0,0) * q(1,0) +
		Nxy(0,1) * q(1,1) +
		Nxy(0,2) * q(1,2);
	return GrowthTensor (ex, ey, 0, 0, 0, exy);
}

int Data::load( Tokenizer & tok) 
{
	if( debug)
		fprintf( stderr, "Loading data...\n");
	// read in the number of measurements
	long n_measurements = tok.read_long( "num. of measurements");
	
	// read in all measurements
	for( long i = 0 ; i < n_measurements ; i ++) {
		Measurement m;
		if( debug)
			fprintf( stderr, "\t- loading measurement %ld.\n",
				 i+1);
		m.load( tok);
		mes.push_back( m);
	}
	
	// make sure that the number of points in each measurement is
	// increasing
	for( size_t i = 1 ; i < mes.size() ; i ++) {
		if( mes[i].pts.size() < mes[i-1].pts.size())
			die( "Measurement %ld has 'lost' points!!!\n",
			     long(i+1));
	}
	if( debug)
		fprintf( stderr, "...data loaded.\n");
	return 0;
}

void Data::center( void)
// ======================================================================
// Center all measurements around origin and scale to fit into unit
// bounding box.
// ----------------------------------------------------------------------
{
	double minx = 0, miny = 0, maxx = 0, maxy = 0;
	assert( mes.size() > 0);
	assert( mes[0].pts.size() > 0);

	minx = maxx = mes[0].pts[0].x;
	miny = maxy = mes[0].pts[0].y;
//	for( size_t i = 0 ; i < mes.size() ; i ++) {
	std::vector<Measurement>::iterator m;
	for( m = mes.begin() ; m != mes.end() ; m++)
	{
		PList::iterator p;
		for( p = (*m).pts.begin() ; p != (*m).pts.end() ; p ++)
		{
			minx = std::min<double> ((*p).x, minx);
			miny = std::min<double> ((*p).y, miny);
			maxx = std::max<double> ((*p).x, maxx);
			maxy = std::max<double> ((*p).y, maxy);
		}
	}

	double w = maxx - minx;
	double h = maxy - miny;
	double sc = 2/h;
	if( w > h) sc = 2/w;
	
	for( m = mes.begin() ; m != mes.end() ; m++)
	{
		PList::iterator p;
		minx = maxx = (*m).pts[0].x;
		miny = maxy = (*m).pts[0].y;
		for( p = (*m).pts.begin() ; p != (*m).pts.end() ; p ++)
		{
			minx = std::min<double> ( (*p).x, minx);
			miny = std::min<double> ( (*p).y, miny);
			maxx = std::max<double> ( (*p).x, maxx);
			maxy = std::max<double> ( (*p).y, maxy);
		}
		double cx = (minx + maxx)/2.0;
		double cy = (miny + maxy)/2.0;
		for( p = (*m).pts.begin() ; p != (*m).pts.end() ; p ++)
		{
			(*p).x = ((*p).x - cx)*sc;
			(*p).y = ((*p).y - cy)*sc;
		}
	}
}

static double sqr( double a) { return a*a; }

void Data::extrapolate( void)
// ======================================================================
// extrapolates missing points
// ----------------------------------------------------------------------
{
	fprintf( stderr, "Extrapolation...\n");
	for( size_t i = mes.size() - 1; i >= 1 ; i --)
	{
		fprintf( stderr, "\t- processing measurement %ld\n",
			 long(i+1));
		Measurement & m1 = mes[i];
		Measurement & m2 = mes[i-1];
		fprintf( stderr, "\t\t- extrapolating %ld new points\n",
			 long(m1.pts.size() - m2.pts.size()));
		// if the previous measurement contains the same number
		// of points, we don't have to extrapolate anything
		if( m1.pts.size() == m2.pts.size()) continue;

		// for all extra points in measurement m1, extrapolate their
		// position to measurement m2
		// --------------------------------------------------
		// first, calculate the matrices A1 and A2
		Matrix A1( 3, long(m2.pts.size()));
		Matrix A2( 2, long(m2.pts.size()));
		for( size_t j = 0 ; j < m2.pts.size() ; j ++) {
			A1(0,j) = m1.pts[j].x;
			A1(1,j) = m1.pts[j].y;
			A1(2,j) = 1.0;
			A2(0,j) = m2.pts[j].x;
			A2(1,j) = m2.pts[j].y;
		}
		// calculate the pseudoinverse of A1
		Matrix A1p = pseudo_inverse( A1);
		// for each new point in m1, extrapolate its position in m2
		size_t extra = m2.pts.size();
		for( size_t j = extra ; j < m1.pts.size() ; j ++)
		{
// 			for( size_t k = 0 ; k < extra ; k ++) {
// 				double lsq = sqr(m1.pts[j].x - m1.pts[k].x) +
// 					sqr(m1.pts[j].y - m1.pts[k].y);
// 				A1(2,k) = 1 / sqrt(lsq);
// 			}
// 			Matrix A1p = pseudo_inverse( A1);
			Point & p = m1.pts[j];
			Matrix P(3,1);
			P(0,0) = p.x; P(1,0) = p.y; P(2,0) = 1.0;
			Matrix W = A1p * P;
//			W.print( "w=");
			Matrix Pe = A2*W;
//			Pe.print( "Pe");
			// create a new point
			Point pe;
			pe.x = Pe(0,0); pe.y = Pe(1,0); pe.is_new = 1;
			m2.pts.push_back( pe);
		}
	}
}

static Point
deform(
	const Point & p,
	const Triangle & t,
	const PList & pts1,
	const PList & pts2)
{
	// name the coordinates
	double x1 = pts1[t.p1].x; double y1 = pts1[t.p1].y;
	double x2 = pts1[t.p2].x; double y2 = pts1[t.p2].y;
	double x3 = pts1[t.p3].x; double y3 = pts1[t.p3].y;
/*
	// calculate the angles within this triangle
	double a1; {
		Vector3d v1( x2-x1, y2-y1, 0); v1.normalize();
		Vector3d v2( x3-x1, y3-y1, 0); v2.normalize();
		a1 = fabs( acos( v1*v2)) * 180 / M_PI;
	}
	double a2; {
		Vector3d v1( x3-x2, y3-y2, 0); v1.normalize();
		Vector3d v2( x1-x2, y1-y2, 0); v2.normalize();
		a2 = fabs( acos( v1*v2)) * 180 / M_PI;
	}
	double a3; {
		a3 = 180 - a1 - a2;
//		fprintf( stderr, "a3 = %f\n", a3);
//		assert( a3 >= 0.0);
	}
	// if the triangle is degenerate, don't interpolate anything
	double dangle = 5;
	if( a1 < dangle || a2 < dangle || a3 < dangle) return p;
*/
	// calculate the isoparamteric coordinates (p,q)
	double den = x2*y1-x3*y1-x1*y2+x3*y2+x1*y3-x2*y3;
	if( fabs(den) < 1e-6) return p;
	double n = - (x1*p.y-x3*p.y-p.x*y1+x3*y1+p.x*y3-x1*y3) / den;
	double s = (x1*p.y-x2*p.y-p.x*y1+x2*y1+p.x*y2-x1*y2) / den;
	// name the deformed triangle coordinates
	x1 = pts2[t.p1].x; y1 = pts2[t.p1].y;
	x2 = pts2[t.p2].x; y2 = pts2[t.p2].y;
	x3 = pts2[t.p3].x; y3 = pts2[t.p3].y;
	// calculate the deformed coordinates
	Point res;
	res.x = x1 + n * (x2-x1) + s * (x3-x1);
	res.y = y1 + n * (y2-y1) + s * (y3-y1);
	return res;
}

void Data::extrapolate3( void)
// ======================================================================
// extrapolates missing points
// ----------------------------------------------------------------------
{
	fprintf( stderr, "Extrapolation...\n");
	for( size_t i = mes.size() - 1; i >= 1 ; i --)
	{
		fprintf( stderr, "\t- processing measurement %ld\n",
			 long(i+1));
		Measurement & m1 = mes[i];
		Measurement & m2 = mes[i-1];
		fprintf( stderr, "\t\t- extrapolating %ld new points\n",
			 long(m1.pts.size() - m2.pts.size()));
		// if the previous measurement contains the same number
		// of points, we don't have to extrapolate anything
		if( m1.pts.size() == m2.pts.size()) continue;

		// for all extra points in measurement m1, extrapolate their
		// position to measurement m2
		// --------------------------------------------------
		// first, calculate the matrices A1 and A2
		Matrix A1( 3, long(m2.pts.size()));
		Matrix A2( 2, long(m2.pts.size()));
		for( size_t j = 0 ; j < m2.pts.size() ; j ++) {
			A1(0,j) = m1.pts[j].x;
			A1(1,j) = m1.pts[j].y;
			A1(2,j) = 1.0;
			A2(0,j) = m2.pts[j].x;
			A2(1,j) = m2.pts[j].y;
		}
		// calculate the pseudoinverse of A1
		Matrix A1p = pseudo_inverse( A1);
		// for each new point in m1, extrapolate its position in m2
		size_t extra = m2.pts.size();
		for( size_t j = extra ; j < m1.pts.size() ; j ++)
		{
// 			for( size_t k = 0 ; k < extra ; k ++) {
// 				double lsq = sqr(m1.pts[j].x - m1.pts[k].x) +
// 					sqr(m1.pts[j].y - m1.pts[k].y);
// 				A1(2,k) = 1 / sqrt(lsq);
// 			}
// 			Matrix A1p = pseudo_inverse( A1);
			Point & p = m1.pts[j];
			Matrix P(3,1);
			P(0,0) = p.x; P(1,0) = p.y; P(2,0) = 1.0;
			Matrix W = A1p * P;
//			W.print( "w=");
			Matrix Pe = A2*W;
//			Pe.print( "Pe");
			// create a new point
			Point pe;
			pe.x = Pe(0,0); pe.y = Pe(1,0); pe.is_new = 1;
			m2.pts.push_back( pe);
		}
	}
}

/*
static double
dist_sq( const Point & p, const Triangle & t, const PList & pts)
{
	double cx = (pts[t.p1].x + pts[t.p2].x + pts[t.p3].x) / 3.0;
	double cy = (pts[t.p1].y + pts[t.p2].y + pts[t.p3].y) / 3.0;
	
	return sqr( p.x - cx) + sqr( p.y - cy);
}
*/

static double
dist_sq( const Point & p, const Point & p1, const Point & p2)
// return the distance from point p to edge <p1,p2>
{
	double x = p.x; double y = p.y;
	double x1 = p1.x; double y1 = p1.y;
	double x2 = p1.x; double y2 = p1.y;
	double den = sqr(x2-x1) + sqr(y2-y1);
	if( fabs(den) > 1e-6) {
		double t = -((x - x1)*(x1 - x2) + (y - y1)*(y1 - y2))/den;
		if( t >= 0 && t <= 1) {
			double cx = x1 + t*(x2-x1);
			double cy = y1 + t*(y2-y1);
			return sqr(cx-x) + sqr(cy-y);
		}
	}
	// one of the edge points is the closest one
	double d1 = sqr(x1-x) + sqr(y1-y);
	double d2 = sqr(x2-x) + sqr(y2-y);
	if( d1 < d2) return d1; else return d2;
}

/*
static double
dist_sq( const Point & p, const Triangle & t, const PList & pts)
{
	double x1 = pts[t.p1].x; double y1 = pts[t.p1].y;
	double x2 = pts[t.p2].x; double y2 = pts[t.p2].y;
	double x3 = pts[t.p3].x; double y3 = pts[t.p3].y;
	
	// get iso-coordinates of the point
	double den = x2*y1-x3*y1-x1*y2+x3*y2+x1*y3-x2*y3;
	if( fabs(den) > 1e-6) {
		double n = - (x1*p.y-x3*p.y-p.x*y1+x3*y1+p.x*y3-x1*y3) / den;
		double s = (x1*p.y-x2*p.y-p.x*y1+x2*y1+p.x*y2-x1*y2) / den;
		if( n >= 0 && s >= 0 && n <= 1 && s <= 1 && n+s <= 1)
			return 0;
	}
	// try the 3 edges
	double min = dist_sq( p, pts[t.p1], pts[t.p2]);
	double d = dist_sq( p, pts[t.p2], pts[t.p3]);
	if( d < min) min = d;
	d = dist_sq( p, pts[t.p1], pts[t.p3]);
	if( d < min) min = d;
	return min;
}
*/

static double
isodist_sq( const Point & p, const Triangle & t, const PList & pts)
{
	double x1 = pts[t.p1].x; double y1 = pts[t.p1].y;
	double x2 = pts[t.p2].x; double y2 = pts[t.p2].y;
	double x3 = pts[t.p3].x; double y3 = pts[t.p3].y;
	
	// get iso-coordinates of the point
	double den = x2*y1-x3*y1-x1*y2+x3*y2+x1*y3-x2*y3;
	if( fabs(den) < 1e-6) return 1e3;
	double n = - (x1*p.y-x3*p.y-p.x*y1+x3*y1+p.x*y3-x1*y3) / den;
	double s = (x1*p.y-x2*p.y-p.x*y1+x2*y1+p.x*y2-x1*y2) / den;
	if( n >= 0 && s >= 0 && n <= 1 && s <= 1 && n+s <= 1)
		return 0;
	// try the 3 edges
	double d1 = dist_sq( Point(n,s), Point(0,0), Point(1,0));
	double d2 = dist_sq( Point(n,s), Point(0,0), Point(0,1));
	double d3 = dist_sq( Point(n,s), Point(1,0), Point(0,1));
	double min = d1; if( d2 < min) min = d2; if( d3 < min) min = d3;
	return min;
}

void Data::extrapolate2( void)
// ======================================================================
// extrapolates missing points
// ----------------------------------------------------------------------
{
	fprintf( stderr, "Extrapolation...\n");
	for( size_t i = mes.size() - 1; i >= 1 ; i --)
	{
		fprintf( stderr, "\t- processing measurement %ld\n",
			 long(i+1));
		Measurement & m1 = mes[i];
		Measurement & m2 = mes[i-1];
		fprintf( stderr, "\t\t- extrapolating %ld new points\n",
			 long(m1.pts.size() - m2.pts.size()));
		// if the previous measurement contains the same number
		// of points, we don't have to extrapolate anything
		if( m1.pts.size() == m2.pts.size()) continue;

		// triangulate the previous measurement
		tlist = m2.triangulation();

		// for all extra points in measurement m1, extrapolate their
		// position to measurement m2
		// --------------------------------------------------
		size_t extra = m2.pts.size();
		for( size_t j = extra ; j < m1.pts.size() ; j ++)
		{
			Point & p = m1.pts[j];
			Point np; np.x = np.y = 0; np.is_new = 1;
			double wsum = 0;
			// for each triangle, calculate the deformed
			// position of the point
			TList::iterator t;
			for( t = tlist.begin() ; t != tlist.end() ; t ++) {
				Point dp = deform( p, *t, m1.pts, m2.pts);
				double dsq = isodist_sq( p, *t, m1.pts) + 0.01;
				double w = 1/(pow(dsq,2));
				np.x += w * dp.x;
				np.y += w * dp.y;
				wsum += w;
			}
			np.x /= wsum;
			np.y /= wsum;
			m2.pts.push_back( np);
		}
	}
}

void Data::print( FILE * fp)
{
	fprintf( fp, "# DATA FILE of measurements\n");
	fprintf( fp, "\n");
	fprintf( fp, "%ld # number of all measurements\n", long( mes.size()));
	fprintf( fp, "\n");

	// print all measurements
	for( size_t i = 0 ; i < mes.size() ; i ++)
	{
		fprintf( fp, "# measurement #%ld\n", long(i+1));
		mes[i].print( fp);
		fprintf( fp, "\n");
	}
}

// ======================================================================
// ======================================================================
// ======================================================================
//
//                           USER INTERFACE
//
// ======================================================================
// ======================================================================
// ======================================================================

Data * data = NULL;

class View {

public:
        double eye_x, eye_y, eye_z;
        double center_x, center_y, center_z;
        double up_x, up_y, up_z;
        double alpha;
        double beta;
        double dist;

        void set( double a, double b, double d) {
                alpha = a;
                beta = b;
                dist = d;

                // first deal with alpha
                double vx = sin(alpha * M_PI / 180.0);
                double vy = - cos(alpha * M_PI / 180.0);
                double vz = 0.0;
                // apply beta
                if( beta > 89) beta = 89;
                if( beta < -89) beta = -89;
                vx = vx * cos( beta * M_PI / 180.0);
                vy = vy * cos( beta * M_PI / 180.0);
                vz = sin( beta * M_PI / 180.0);
                // set up vector based on v
                up_x = - vx; up_x = 0;
                up_y = - vy; up_y = 0;
                up_z = 0; up_z = 1;
                // apply distance to v
                vx *= dist;
                vy *= dist;
                vz *= dist;
                // set eye coorinates based on v
                eye_x = center_x + vx;
                eye_y = center_y + vy;
                eye_z = center_z + vz;
        }

	View() {
		center_x = 0;
		center_y = 0;
		center_z = 0;
		set( 45, 40, 2);
	}
} view;

using namespace OGLUI;

class MainWindow
	: public OGLUI::WindowWidget
{
private:
	// some internal variables
	int last_mouse_x, last_mouse_y;
	enum MouseState { ROTATE, ZOOM, TIME_SET, NONE } mouse_state;
	// the label
	LabelWidget * lab;
	
	double tensorScale_;
public:
	// overloaded constructor
	MainWindow (Widget * p_parent);
	// overloaded redraw
	virtual void redraw ();
	// overloaded mouse handler
	virtual bool mouse (const MouseEvent & ev);
	// tensor scale interface
	void setTensorScale (double newVal)
		{
			bool different = tensorScale_ != newVal;
			if (different)
				tensorScale_ = newVal;
			if (draw_tensors && different)
				postRedisplay ();
		}
	double getTensorScale (void) const
		{
			return tensorScale_;
		}
};

// the main window
MainWindow * win = NULL;

static void glDrawText3d(
	double x, double y, double z,
	char * str, ...)
// ---------------------------------------------------------------------------
// draw a string at the specified xyz coordinates
// ---------------------------------------------------------------------------
{
	void * font = GLUT_BITMAP_HELVETICA_12;
	const double fontHeight = 9;

	// prepare the string into 'buff'
	va_list ap;
	va_start( ap, str);
	char buff[ 4096];
	vsprintf( buff, str, ap);
	va_end( ap);

	// set the position of the raster cursor
	glRasterPos3d( x, y, z);
	GLboolean valid;
	glGetBooleanv (GL_CURRENT_RASTER_POSITION_VALID, & valid);
	if (! valid) return;
	// extract the window coordinates
	double pos[4];
	glGetDoublev (GL_CURRENT_RASTER_POSITION, pos);
	// save the state of the opengl
	glPushAttrib (GL_ALL_ATTRIB_BITS);
	glMatrixMode (GL_PROJECTION);
	glPushMatrix ();
	glLoadIdentity ();
	glMatrixMode (GL_MODELVIEW);
	glPushMatrix ();
	glLoadIdentity ();
	double textWidth = glutBitmapLength (font, (unsigned char *) buff);
	double textHeight = fontHeight;
	glViewport (0, 0,
		    GLint(win-> getInGeom().width()),
		    GLint(win-> getInGeom().height()));
	gluOrtho2D (0, win-> getInGeom().width(),
		    0, win-> getInGeom().height());
	glEnable (GL_BLEND);
	glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glDisable (GL_LIGHTING);
	glDisable (GL_DEPTH_TEST);
	glColor4f (0,0,0,0.5);
	OGLUI::glBox (pos[0]-2, pos[1]-2,
		      pos[0]+textWidth+2, pos[1]+textHeight+2);
	glPopAttrib ();
	glMatrixMode (GL_PROJECTION);
	glPopMatrix ();
	glMatrixMode (GL_MODELVIEW);
	glPopMatrix ();
	for( char * c = buff ; * c != '\0' ; c ++)
		glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, *c);
}

static void glDrawTensor2D (
	const Vector3d & c,
	const Vector3d & v1,
	const Vector3d & v2,
	const Color & col1,
	const Color & col2)
{
	glPushAttrib (GL_ALL_ATTRIB_BITS);
	glDisable (GL_DEPTH_TEST);
	glDisable (GL_LIGHTING);
	// draw the filled inside
	long n = 20;
	Color::Green-> glColor ();
	glBegin (GL_TRIANGLE_FAN);
	glColor4f (0,0,0,0);
	glVertex2d (c.x, c.y);
	for (long i = 0 ; i <= n ; i ++)
	{
		double alpha = i * M_PI * 2 / n;
		Color col = fabs (sin (alpha)) * col1
			+ fabs (cos (alpha)) * col2;
		col.glColor ();
		Vector3d p1 = c + sin (alpha) * v1 + cos (alpha) * v2;
		glVertex2d (p1.x, p1.y);
	}
	glEnd ();

	// draw the cross inside
	glColor4f (0,0,0,0.5);
	glBegin (GL_LINES);
	glVertex2d (c.x - v1.x, c.y - v1.y);
	glVertex2d (c.x + v1.x, c.y + v1.y);
	glVertex2d (c.x - v2.x, c.y - v2.y);
	glVertex2d (c.x + v2.x, c.y + v2.y);
	glEnd ();
	glPopAttrib ();
}

class MyToggleListener : public ToggleButtonWidget::Listener
{
private:
	const std::string _str;
public:
	MyToggleListener (const std::string & str) : _str(str) {;}
	void stateChanged ( const ToggleButtonWidget & tb) {
//		std::cerr << "Toggle(" << _str << ") changed state: "
//			  << tb.getState() << "\n";
		if (_str == "draw_cubes")
		{
			draw_cubes = tb.getState();
			win-> postRedisplay();
		}
		else if (_str == "draw_numbers")
		{
			draw_numbers = tb.getState();
			win-> postRedisplay();
		}
		else if (_str == "draw_grid")
		{
			draw_grid = tb.getState();
			win-> postRedisplay();
		}
		else if (_str == "draw_lines")
		{
			draw_lines = tb.getState();
			win-> postRedisplay();
		}
		else if (_str == "autozoom")
		{
			draw_autozoom = tb.getState();
			win-> postRedisplay();
		}
		else if (_str == "antialias")
		{
			draw_antialias = tb.getState();
			win-> postRedisplay();
		}
		else if (_str == "draw_tensors")
		{
			draw_tensors = tb.getState();
			win-> postRedisplay();
		}
		else if (_str == "draw_growth_rates")
		{
			draw_growth_rates = tb.getState();
			win-> postRedisplay();
		}
		else
		{
			std::cerr << "UNKNOWN toggle pressed: '"
				  << _str << "'\n";
		}
	}
};

MainWindow::MainWindow (Widget * p_parent)
	: WindowWidget (p_parent, "mainwindow")
{
	tensorScale_ = 0.01;
	
	setDisplayMode (GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	// set preferred size
	setPrefSize (OGLUI::Size (600,600));
	// set a manager
	OGLUI::SLM::Pointer lm = setLayoutManager (new SLM);
	// add more widgets
	Color::Pointer bg = new Color (0,0,0,0.6);
	Color::Pointer fg = Color::Orange;
	OGLUI::Font::Pointer font = BitmapFont::Helvetica12;
	int width = 130;
	int height = 25;
	int curry = 50;
	int gap = 0;
	int currx = 5;
	double bevel = 0.0;
	double border = 0.0;
	// Label - with current time
	lab = new LabelWidget (this);
	lab-> setBgColor (bg);
	lab-> setFontColor (fg);
	lm-> setGeometry (lab
			  , OGLUI::Geometry (Position (currx,curry)
					     , Size (width,height)));
	lab-> setFont (font);
	lab-> setBevelSize (bevel);
	lab-> setLabel ("Time: 0.00");
	lab-> borderThickness (border);
	// Toggle : draw cubes
	ToggleButtonWidget * tog;
	curry += height + gap;
	tog = new ToggleButtonWidget (this);
	tog-> setBgColor (bg);
	tog-> setFontColor (fg);
	lm-> setGeometry (tog
			  , OGLUI::Geometry (Position (currx,curry)
					     , Size (width,height)));
	tog-> setFont (font);
	tog-> setBevelSize (bevel);
	tog-> setLabel ("Draw cubes");
	tog-> borderThickness (border);
	tog-> setState (draw_cubes);
	tog-> addListener (new MyToggleListener ("draw_cubes"));
	// Toggle : draw numbers
	curry += height + gap;
	tog = new ToggleButtonWidget (this);
	tog-> setBgColor (bg);
	tog-> setFontColor (fg);
	lm-> setGeometry (tog
			  , OGLUI::Geometry (Position (currx,curry)
					     , Size (width,height)));
	tog-> setFont (font);
	tog-> setBevelSize (bevel);
	tog-> setLabel ("Draw numbers");
	tog-> borderThickness (border);
	tog-> setState (draw_numbers);
	tog-> addListener (new MyToggleListener ("draw_numbers"));
	// slider - for scaling tensors
	curry += height + gap;
	SliderWidget <double> * slider = new SliderWidget <double> (this);
	slider-> setBgColor (bg);
	lm-> setGeometry (slider
			  , OGLUI::Geometry (Position (currx,curry)
					     , Size (width,height)));
	slider-> setBevelSize (bevel);
	slider-> minVal (0.0001);
	slider-> maxVal (0.1);
	slider-> currVal ( getTensorScale ());
	slider-> borderThickness (border);
	slider-> valueChangedSignal.connect (
		slot (*this, & MainWindow::setTensorScale));
	// Toggle : draw growth tensors
//	curry += height + gap;
	curry += height;
	tog = new ToggleButtonWidget (this);
	tog-> setBgColor (bg);
	tog-> setFontColor (fg);
	lm-> setGeometry (tog
			  , OGLUI::Geometry (Position (currx,curry)
					     , Size (width,height)));
	tog-> setFont (font);
	tog-> setBevelSize (bevel);
	tog-> setState (draw_tensors);
	tog-> setLabel ("Growth tensors");
	tog-> borderThickness (border);
	tog-> addListener (new MyToggleListener ("draw_tensors"));
	// Toggle : draw grid
	curry += height + gap;
	tog = new ToggleButtonWidget (this);
	tog-> setBgColor (bg);
	tog-> setFontColor (fg);
	lm-> setGeometry (tog
			  , OGLUI::Geometry (Position (currx,curry)
					     , Size (width,height)));
	tog-> setFont (font);
	tog-> setBevelSize (bevel);
	tog-> setLabel ("Grid");
	tog-> borderThickness (border);
	tog-> setState (draw_grid);
	tog-> addListener (new MyToggleListener ("draw_grid"));
	// Toggle : autozoom
	curry += height + gap;
	tog = new ToggleButtonWidget (this);
	tog-> setBgColor (bg);
	tog-> setFontColor (fg);
	lm-> setGeometry (tog
			  , OGLUI::Geometry (Position (currx,curry)
					     , Size (width,height)));
	tog-> setFont (font);
	tog-> setBevelSize (bevel);
	tog-> setLabel ("Autozoom");
	tog-> borderThickness (border);
	tog-> setState (draw_autozoom);
	tog-> addListener (new MyToggleListener ("autozoom"));
	// Toggle : antialias
	curry += height + gap;
	tog = new ToggleButtonWidget (this);
	tog-> setBgColor (bg);
	tog-> setFontColor (fg);
	lm-> setGeometry (tog
			  , OGLUI::Geometry (Position (currx,curry)
					     , Size (width,height)));
	tog-> setFont (font);
	tog-> setBevelSize (bevel);
	tog-> setLabel ("Antialiasing");
	tog-> borderThickness (border);
	tog-> setState (draw_antialias);
	tog-> addListener (new MyToggleListener ("antialias"));
	// Toggle : antialias
	curry += height + gap;
	tog = new ToggleButtonWidget (this);
	tog-> setBgColor (bg);
	tog-> setFontColor (fg);
	lm-> setGeometry (tog
			  , OGLUI::Geometry (Position (currx,curry)
					     , Size (width,height)));
	tog-> setFont (font);
	tog-> setBevelSize (bevel);
	tog-> setLabel ("Growth rates");
	tog-> borderThickness (border);
	tog-> setState (draw_growth_rates);
	tog-> addListener (new MyToggleListener ("draw_growth_rates"));
	// Toggle : draw lines
	curry += height + gap;
	tog = new ToggleButtonWidget (this);
	tog-> setBgColor (bg);
	tog-> setFontColor (fg);
	lm-> setGeometry (tog
			  , OGLUI::Geometry (Position (currx,curry)
					     , Size (width,height)));
	tog-> setFont (font);
	tog-> setBevelSize (bevel);
	tog-> setLabel ("Lines");
	tog-> borderThickness (border);
	tog-> setState (draw_lines);
	tog-> addListener (new MyToggleListener ("draw_lines"));
}

void MainWindow::redraw ()
// ----------------------------------------------------------------------
// called whenever the contents of the OGL window need to be redrawn.
// ----------------------------------------------------------------------
{
	// get inside geometry
	Geometry g = getInGeom();
        glViewport( GLint (g.minx ()), GLint (g.miny ())
		    , GLint(g.width()), GLint(g.height()));
	// figure out min/max of the bounding box
	double minx, maxx, miny, maxy;
	{
		// default values for min/max x/y
		int is_new;
		data-> get_coord (0, curr_time, minx, miny, is_new);
		maxx = minx;
		maxy = miny;
	}
	for (size_t i = 1 ; i < data-> mes[0].pts.size() ; i ++) {
		double x, y; int is_new;
		data-> get_coord( i, curr_time, x, y, is_new);
		minx = std::min<double> (x,minx);
		maxx = std::max<double> (x,maxx);
		miny = std::min<double> (y,miny);
		maxy = std::max<double> (y,maxy);
	}

	if (draw_autozoom)
	{
		// figure out the ortho2d projection
		double cx = (maxx+minx) / 2;
		double cy = (maxy+miny) / 2;
		double r = g.width() / g.height();
		double x1 = minx; double x2 = maxx;
		double h = (maxx - minx) / r;
		double y1 = cy - h/2; double y2 = cy + h/2;
		if (h < maxy - miny)
		{
			y1 = miny; y2 = maxy;
			double w = (maxy - miny) * r;
			x1 = cx - w/2; x2 = cx + w/2;
		}
		// make room around the margins : 10%
		double vx = x2 - cx; double vy = y2 - cy;
		x1 = cx - 1.1 * vx; y1 = cy - 1.1 * vy;
		x2 = cx + 1.1 * vx; y2 = cy + 1.1 * vy;

		// set up projection matrix
		glMatrixMode (GL_PROJECTION );
		glLoadIdentity ();
		gluOrtho2D( x1, x2, y1, y2);
		// set the MODELVIEW matrix
		glMatrixMode( GL_MODELVIEW );
		glLoadIdentity();
	}
	else
	{
		// set up projection matrix
		glMatrixMode( GL_PROJECTION );
		glLoadIdentity();
		gluPerspective( 60.0, g.width()/g.height(),
				0.1, 100.0);
		// set the MODELVIEW matrix
		glMatrixMode( GL_MODELVIEW );
		glLoadIdentity();
		gluLookAt( view.eye_x, view.eye_y, view.eye_z,
			   view.center_x, view.center_y, view.center_z,
			   view.up_x, view.up_y, view.up_z);
	}

	// setup first light
	{
		float light_ambient[] = {0.6, 0.6, 0.6, 0.0 };
		float light_diffuse[] = {0.8, 0.8, 0.8, 0.0 };
		float light_specular[] = {0.8, 0.8, 0.8, 0.0 };
//		float light_position[] = {3.0, 3.0, 3.0, 0.0 };
		float light_position[] = {
			view.eye_x,
			view.eye_y,
			view.eye_z,
			1};
		
		glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
		glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
		glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
		glLightfv(GL_LIGHT0, GL_POSITION, light_position);

		glEnable( GL_LIGHT0 );
	}

        // set the background color to light grey
        glClearColor( 0.8, 0.8, 0.8, 1.0);
        // clear the background
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	// enable transparencies
	glEnable( GL_BLEND);
	glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        // enable anti-aliasing of lines
	if (draw_antialias)
		glEnable( GL_LINE_SMOOTH);
	else
		glDisable( GL_LINE_SMOOTH);
	glLineWidth (1);
        // enable automatic normalization of normals
        glEnable( GL_NORMALIZE);
	// enable Z buffer
	glEnable( GL_DEPTH_TEST);
	// enable lighting
	glEnable( GL_LIGHTING );
	// display a grid
	if (draw_grid)
	{
		glDisable( GL_LIGHTING );
		long ng = 12; double d = 2.0/ng;
		glNormal3d( 0, 0, 1);
		glBegin( GL_QUADS);
		for( long i = 0 ; i < ng ; i ++) {
			for( long j = 0 ; j < ng ; j ++) {
				if( (i+j)%2) glColor4f( 0.5, 0.5, 0.5, 0.3);
				else glColor4f( 0.3, 0.4, 0.3, 0.3);
				double cx = -1 + d * i;
				double cy = -1 + d * j;
				glVertex3d( cx, cy, 0);
				glVertex3d( cx+d, cy, 0);
				glVertex3d( cx+d, cy+d, 0);
				glVertex3d( cx, cy+d, 0);
			}
		}
		glEnd();
		glEnable( GL_LIGHTING );
	}

	// draw the tensors
	Color maxCol (1.0,0.0,0.0,0.8);
	Color neuCol (0.4,0.4,0.4,0.8);
	Color minCol (0.0,0.0,1.0,0.8);
	if (draw_growth_rates)
	{
		// draw the triangles
		glDisable (GL_DEPTH_TEST);
		glDisable( GL_LIGHTING );
		glColor4f( 0, 0, 0, 0.3);
		glBegin (GL_TRIANGLES);
		for( size_t i = 0 ; i < data-> tlist.size() ; i ++) {
			Triangle & t = data-> tlist[i];
			Vector3d pos[3];
			data-> get_triangle_pos (t, curr_time, pos);
			double r = data -> get_triangle_growth_rate
				(t, curr_time, 0.05);
			Color col;
			if (r >= 0)
			{
				double rat = r / data-> max_growth_rate;
				rat = std::min<double> (rat, 1.0);
				col = rat * maxCol + (1-rat) * neuCol;
			}
			else
			{
				double rat = r / data-> min_growth_rate;
				rat = std::min<double> (rat, 1.0);
				col = rat * minCol + (1-rat) * neuCol;
			}
			col.glColor ();
			glVertex3d( pos[0].x, pos[0].y, 0);
			glVertex3d( pos[1].x, pos[1].y, 0);
			glVertex3d( pos[2].x, pos[2].y, 0);
		}
		glEnd();
		glEnable (GL_LIGHTING);
		glEnable (GL_DEPTH_TEST);
	}

	// draw the lines
	if (draw_lines)
	{
		glDisable (GL_DEPTH_TEST);
		glDisable( GL_LIGHTING );
		glColor4f( 0, 0, 0, 1);
		glBegin( GL_LINES);
		for( size_t i = 0 ; i < data-> tlist.size() ; i ++) {
			Triangle & t = data-> tlist[i];
			double x1, y1, x2, y2, x3, y3; int is_new;
			data-> get_coord( t.p1, curr_time, x1, y1, is_new);
			data-> get_coord( t.p2, curr_time, x2, y2, is_new);
			data-> get_coord( t.p3, curr_time, x3, y3, is_new);
			glVertex3d( x1, y1, 0);
			glVertex3d( x2, y2, 0);
			glVertex3d( x2, y2, 0);
			glVertex3d( x3, y3, 0);
			glVertex3d( x3, y3, 0);
			glVertex3d( x1, y1, 0);
		}
		glEnd();
		glEnable( GL_LIGHTING );
		glEnable (GL_DEPTH_TEST);
	}

	// draw growth tensors
	if (draw_tensors)
	{
		glDisable (GL_LIGHTING);
		glDisable (GL_DEPTH_TEST);
		glColor4f (0,0.5,0,0.8);
		for (size_t i = 0 ; i < data-> tlist.size () ; i ++)
		{
			Triangle & t = data-> tlist [i];
			GrowthTensor gt = data-> get_growth_tensor
				(t, curr_time, 0.01);
			// get the eigenvalues
			double a1, a2, a3;
			gt.get_eigenvalues (a1, a2, a3);
			// get the max. abs. value eigenvalues into a1 & a2
			if (fabs (a1) > fabs (a2)) std::swap (a1,a2);
			if (fabs (a2) > fabs (a3)) std::swap (a2,a3);
			if (fabs (a1) > fabs (a2)) std::swap (a1,a2);
			// get the corresponding eigenvectors
			Vector3d v1, v2, v3;
			gt.get_eigenvectors (a3, v1, v2);
			gt.get_eigenvectors (a2, v2, v3);
			v1.normalize ();
			v2.normalize ();
			v1 = tensorScale_ * a3 * v1;
			v2 = tensorScale_ * a2 * v2;
			// get the center of the triangle
			Vector3d c = data-> get_triangle_center
				(t, curr_time);
			Color col1 = a3 >= 0 ? maxCol : minCol;
			Color col2 = a2 >= 0 ? maxCol : minCol;
			glDrawTensor2D (c, v1, v2, col1, col2);
		}
		glEnable (GL_DEPTH_TEST);
		glEnable (GL_LIGHTING);
		
	}
	
	// display all points
	if (draw_cubes)
	{
		// figure out the cube size
		double cube_size = std::min<double> ((maxx-minx)/100,
						     (maxy-miny)/100);
	
		for( size_t i = 0 ; i < data-> mes[0].pts.size() ; i ++) {
			double x, y; int is_new;
			data-> get_coord( i, curr_time, x, y, is_new);
			if( is_new)
				Material_YellowPlastic.set();
			else
				Material_RedPlastic.set();
			glDrawCube( x, y, 0, cube_size);
		}
	}

	// display the names of the points
	if (draw_numbers)
	{
		glDisable (GL_LIGHTING);
		glDisable (GL_DEPTH_TEST);
		for( size_t i = 0 ; i < data-> mes[0].pts.size() ; i ++) {
			double x, y; int is_new;
			data-> get_coord( i, curr_time, x, y, is_new);
			if (is_new)
				glColor4f (1,1,0,1.0);
			else
				glColor4f (1,0.5,0.5,1);
			glDrawText3d( x, y, 0, "%02ld", i+1);
		}
		glEnable (GL_LIGHTING);
		glEnable (GL_DEPTH_TEST);
	}

	// draw the legend for the growth rates
	if (draw_growth_rates)
	{
		glDisable (GL_LIGHTING);
		glDisable (GL_DEPTH_TEST);
		long w = 50;
		glViewport( GLint(g.width()) - w, time_line_height,
			    w, GLint(g.height()-time_line_height));
		glMatrixMode( GL_PROJECTION );
		glLoadIdentity();
		gluOrtho2D( 0, 2,
			    data-> min_growth_rate, data-> max_growth_rate);
		glMatrixMode( GL_MODELVIEW );
		glLoadIdentity();
		// clear the window
		glColor4f (0,0,0,0.3);
		glBegin (GL_QUADS);
		glVertex2d (0,data-> min_growth_rate);
		glVertex2d (2,data-> min_growth_rate);
		glVertex2d (2,data-> max_growth_rate);
		glVertex2d (0,data-> max_growth_rate);
		glEnd ();
		// draw the quad for the positive growth rates
		glBegin (GL_QUADS);
		maxCol.glColor ();
		glVertex2d (0,data-> max_growth_rate);
		glVertex2d (2,data-> max_growth_rate);
		neuCol.glColor ();
		glVertex2d (2,0);
		glVertex2d (0,0);
		glEnd ();
		// draw the quad for the negative growth rates
		glBegin (GL_QUADS);
		minCol.glColor ();
		glVertex2d (0,data-> min_growth_rate);
		glVertex2d (2,data-> min_growth_rate);
		neuCol.glColor ();
		glVertex2d (2,0);
		glVertex2d (0,0);
		glEnd ();
		// draw the legend...
		glMatrixMode( GL_PROJECTION );
		glLoadIdentity();
		gluOrtho2D( 0, 60,
			    0, g.height()-time_line_height-20);
		glMatrixMode( GL_MODELVIEW );
		glLoadIdentity();

		glColor4f (1,1,1,1);
		glDrawText3d (0,2,0,"%5.1f%%",
			      100 * data-> min_growth_rate);
		glDrawText3d (0
			      , g.height()-time_line_height-20-12, 0,
			      "%5.1f%%",
			      100 * data-> max_growth_rate);
		glEnable (GL_DEPTH_TEST);
		glEnable (GL_LIGHTING);
	}

	// draw the time line
	// --------------------------------------------------
	{
		long nmes = long( data-> mes.size());
		double tmax = nmes-1;
		// determine the width of one pixel
		double wp = double(nmes) / g.width();
		glDisable( GL_LIGHTING );
		glDisable( GL_DEPTH_TEST);
		glViewport( 0, 0, GLint(g.width())
			    , time_line_height);
		glMatrixMode( GL_PROJECTION );
		glLoadIdentity();
		double minx = 0; double maxx = tmax;
		double miny = 0; double maxy = 2;
		gluOrtho2D( minx, maxx, miny, maxy);
		glMatrixMode( GL_MODELVIEW );
		glLoadIdentity();
		// draw the background
		glColor4f( 0, 0, 0, 0.5);
		glBegin( GL_QUADS);
		glVertex2d( minx, miny); glVertex2d( maxx, miny);
		glVertex2d( maxx, maxy); glVertex2d( minx, maxy);
		glEnd();
		// draw the background ticks
		glColor4f( 1, 0, 0, 0.8);
		glBegin( GL_LINES);
		for( long i = 0 ; i <= (nmes-1)*10 ; i ++) {
			glVertex2d( i/10.0, 0);
			glVertex2d( i/10.0, i % 10 == 0 ? 1 : 0.4);
		}
		glEnd();
		// draw the tick
		double t = curr_time;
		glColor4f( 1, 1, 0, 1);
		glBegin( GL_TRIANGLES);
		glVertex2d( t, 1);
		double w = 5 * wp;
		glVertex2d( t + w, 2);
		glVertex2d( t - w, 2);
		glEnd();
		glEnable( GL_LIGHTING );
		glEnable( GL_DEPTH_TEST);
	}
}

bool MainWindow::mouse (const MouseEvent & ev)
{
	if (mouse_state == NONE)
		// try propagating to the kids
		if (propagateMouseEvent (ev)) return true;
	// if kids don't want it, it is ours
	if (ev.reason == MouseEvent::Motion)
	{
		// figure out the displacement of the mouse
		int dx = ev.x - last_mouse_x;
		int dy = ev.y - last_mouse_y;
		// remember the mouse coordinates
		last_mouse_x = ev.x;
		last_mouse_y = ev.y;

		// depending on the mode select an action
		if( mouse_state == ROTATE) {
			double alpha = - dx*0.3 + view.alpha;
			double beta = -dy*0.3 + view.beta;
			double dist = view.dist;
			view.set( alpha, beta, dist);
			postRedisplay ();
		} else if( mouse_state == ZOOM) {
			double alpha = view.alpha;
			double beta = view.beta;
			double dist = view.dist * 
				pow(pow( 2.0, dy/180.0), 2.0);
			view.set( alpha, beta, dist);
			postRedisplay ();
		} else if( mouse_state == TIME_SET) {
			double dt = double(dx) / getInGeom().width();
			curr_time += data-> mes.size() * dt;
			if( curr_time < 0)
				curr_time = 0;
			if( curr_time > data-> mes.size()-1)
				curr_time = data-> mes.size()-1;
			char buff [4096];
			sprintf (buff, "Time: %.2f", curr_time);
			lab-> setLabel (buff);
			postRedisplay ();
		} else {
			// mouse is doing nothing
		}
	}
	else if( ev.reason == MouseEvent::Push)
	{
		// capture the mouse
		if (! captureMouse ()) return false;
		if( ev.button == GLUT_LEFT_BUTTON)
		{
			if( ev.y < time_line_height)
			{
				mouse_state = TIME_SET;
			}
			else
			{
				mouse_state = ROTATE;
			}
		}
		else if( ev.button == GLUT_MIDDLE_BUTTON)
		{
			mouse_state = ZOOM;
		}

		last_mouse_x = ev.x;
		last_mouse_y = ev.y;
	}
	else if (ev.reason == MouseEvent::Release)
	{
		releaseMouse ();
		mouse_state = NONE;
	}
	return true;
}

static void usage( void)
{
	fprintf( stderr,
		 "Usage: jo2fem [-ex1|-ex2|-ex3] [-do fname] [fname]\n"
		 "\n"
		 "       - no fname means standard input\n"
		 "       - do fname - specifies output in Dave's leaf format\n"
		);
	exit( -1);
}

// command line parameters
struct CMD_LINE {
	std::string fname;
	std::string progname;
	std::string dave_fname;
	int extrapolation_method;
} cmd_line;

static void parse_command_line( int argc, char ** argv)
// ----------------------------------------------------------------------
// parses command line arguments
// ----------------------------------------------------------------------
{
	// remember how to invoke ourselves
	cmd_line.progname = argv[0];
	cmd_line.fname = "";
	cmd_line.dave_fname = "";
	cmd_line.extrapolation_method = 1;
	long arg = 1;
	while( arg < argc) {
		if( strcasecmp( argv[arg], "-ex1") == 0) {
			cmd_line.extrapolation_method = 0;
			arg ++; continue;
		}
		else if( strcasecmp( argv[arg], "-ex2") == 0) {
			cmd_line.extrapolation_method = 1;
			arg ++; continue;
		} else if( strcasecmp( argv[arg], "-ex3") == 0) {
			cmd_line.extrapolation_method = 2;
			arg ++; continue;
		} else if( strcasecmp( argv[arg], "-do") == 0) {
			arg ++;
			if( arg >= argc) usage();
			cmd_line.dave_fname = argv[arg];
			arg ++; continue;
		} else if( cmd_line.fname == "") {
			cmd_line.fname = argv[arg];
			arg ++; continue;
		}
		usage( );
	}
}

int main( int argc, char ** argv)
{
	parse_command_line( argc, argv);

	// read in the number of measurements
	FILE * fp = stdin;
	if( cmd_line.fname != "") {
		fp = fopen( cmd_line.fname.c_str(), "r");
		if( fp == NULL) die( "Could not open '%s' for reading.\n",
				     cmd_line.fname.c_str());
	}
	Tokenizer tok( fp);
	data = new Data;
	data-> load( tok);
	fclose( fp);
	if( cmd_line.extrapolation_method == 0)
		data-> extrapolate();
	else if( cmd_line.extrapolation_method == 1)
		data-> extrapolate2();
	else if( cmd_line.extrapolation_method == 2)
		data-> extrapolate3();
	else die( "Internal error at %s:%ld", __FILE__, long( __LINE__));
	data-> center();
	data-> triangulate ();
	data-> calculate_growth_rates();

	// if dave output was requested, print it out
	if( cmd_line.dave_fname != "") {
		fprintf( stderr, "Generating data for Dave...\n");
		FILE * fp = fopen( cmd_line.dave_fname.c_str(), "w");
		if( fp == NULL) die( "Cannot open output file %s.",
				     cmd_line.dave_fname.c_str());
		fprintf( fp, "%3ld # number of points\n",
			 long(data-> mes[0].pts.size()));
		fprintf( fp, "%3ld # number of triangles\n",
			 long(data-> tlist.size()));
		fprintf( fp, "%3ld # number of time steps\n",
			 long(data-> mes.size()));
		fprintf( fp,
			 "\n"
			 "# times:\n");
		for( size_t i = 0 ; i < data-> mes.size() ; i ++)
			fprintf( fp, "%.2f ", data-> mes[i].time);
		fprintf( fp, "\n");
		fprintf( fp,
			 "\n"
			 "# point coordinates (x,y) for each time\n");
		size_t n_pts = data-> mes[0].pts.size();
		for( size_t i = 0 ; i < n_pts ; i ++) {
			for( size_t j = 0 ; j < data-> mes.size() ; j ++)
				fprintf( fp, "%9.5f %9.5f ",
					 data-> mes[j].pts[i].x,
					 data-> mes[j].pts[i].y);
			fprintf( fp, "\n");
		}
		fprintf( fp,
			 "\n"
			 "# triangles (indexing points above, "
			 "starting with 0\n");
		for( size_t i = 0 ; i < data-> tlist.size() ; i ++)
			fprintf( fp, "%3ld %3ld %3ld\n",
				 data-> tlist[i].p1,
				 data-> tlist[i].p2,
				 data-> tlist[i].p3);
		fclose( fp);
		fprintf( stderr, "File '%s' was generated.\n",
			 cmd_line.dave_fname.c_str());
		exit(0);
	}
	
	fprintf( stderr, "Number of measurements = %ld\n",
		 long(data-> mes.size()));

        view.set( 0, 89, 2);

	// create an interface
	OGLUI::UI * ui = new OGLUI::UI (argc, argv);
	// create a new window
	win = new MainWindow (ui);
	ui-> run ();
	
	return 0;
}

/*
static void menu_func( int item)
{
	fprintf( stderr, "Menu %d selected.\n", item);
//	curr_measurement = size_t(item);
	glutPostRedisplay();
}
*/
