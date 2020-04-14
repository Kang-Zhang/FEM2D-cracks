#include <cstdlib>
#include <cstdio>
#include <cerrno>
#include <unistd.h>
#include <string.h>
#include <stdarg.h>
#include "TriangulatorDelaunay2D.hh"
#include "die.h"

void
TriangulatorDelaunay2D::add_point( double x, double y)
{
	pts.push_back( Point( x, y));
}

static void execute( char * str, ...)
{
	va_list ap;
	va_start( ap, str);
	char command[ 4096];
	vsprintf( command, str, ap);
	va_end( ap);

	// fprintf( stderr, "execute( '%s')\n", command);
	int res = system( command);
	// fprintf( stderr, "   - system() returned %d\n", res);
	if( res) die( "Could not execute command '%s'", command);
}

std::vector<TriangulatorDelaunay2D::Triangle>
TriangulatorDelaunay2D::triangulate( void)
// TriangulatorDelaunay2D::tlist TriangulatorDelaunay2D::triangulate( void)
{
	fprintf( stderr, "Triangulation starting...\n");

	// special cases:
	if( pts.size() < 3) return std::vector<Triangle>();
	if( pts.size() == 3) {
		Triangle t;
		t.p1 = 0; t.p2 = 1; t.p3 = 2;
		std::vector<Triangle> res;
		res.push_back( t);
		return res;
	}

	// write the set of points to a temporary file
	char fname1[4096];
	sprintf( fname1, "/tmp/triangulator.tmp.XXXXXX");
	int fd = mkstemp( fname1);
	if( fd == -1)
		die( "TriangulatorDelaunay2D::triangulate(): ERROR:\n"
		     "\t- cannot open temporary file '%s'\n"
		     "\t- because unix error = %s\n",
		     fname1, strerror( errno));
	FILE * fp = fdopen( fd, "w");
 	if( fp == NULL) {
 		die( "TriangulatorDelaunay2D::triangulate(): ERROR:\n"
 		     "\t- cannot associate file stream with file descriptor\n"
 		     "\t- unix error = %s\n",
 		     strerror( errno));
	}
	fprintf( fp, "2\n"); // dimension
	fprintf( fp, "%ld\n", long(pts.size())); // number of points
	for( size_t i = 0 ; i < pts.size() ; i ++)
		fprintf( fp, "%.10f %.10f\n", pts[i].x,  pts[i].y);
	fclose( fp);
	// create a temporary file where the result will be stored
	char fname2[4096];
	sprintf( fname2, "/tmp/triangulator-res.tmp.XXXXXX");
	int fd2 = mkstemp( fname2);
	if( fd2 != -1) close( fd2);
	// perform the delaunay triangulation using an external program
	// qhull
//	execute( "qhull d o < %s > %s", fname1, fname2);
	execute( "qhull d QJ o < %s > %s", fname1, fname2);

	// read in the result
	FILE * fp2 = fopen( fname2, "r");
	if( fp2 == NULL)
		die( "TriangulatorDelaunay2d::triangulate(): ERROR\n"
		     "\t- cannot open the resulting triangulation file\n"
		     "\t- '%s'", fname2);
	long dim;
	if( 1 != fscanf( fp2, "%ld", & dim) || dim != 3)
		die( "TriangulatorDelaunay2d::triangulate(): ERROR\n"
		     "\t- the dimension expected to be '3' in file\n"
		     "\t- '%s'", fname2);
	long npts, ntri, nfac;
	if( 3 != fscanf( fp2, "%ld %ld %ld", & npts, & ntri, & nfac))
		die( "TriangulatorDelaunay2d::triangulate(): ERROR\n"
		     "\t- cannot read number of points, triangles, faces in\n"
		     "\t- '%s'", fname2);
	if( npts != long(pts.size()))
		die( "TriangulatorDelaunay2d::triangulate(): ERROR\n"
		     "\t- number of points disagrees (in=%ld out=%ld) in\n"
		     "\t- '%s'", long(pts.size()), npts, fname2);
	// skip the actual points
	for( long i = 0 ; i < npts ; i ++) {
		double x;
		if( 3 != fscanf( fp2, "%lf %lf %lf", &x, &x, &x))
			die( "TriangulatorDelaunay2d::triangulate(): ERROR\n"
			     "\t- coult not read all points in\n"
			     "\t- '%s'", fname2);
	}
	// read in the triangles
	std::vector<Triangle> res;
	for( long i = 0 ; i < ntri ; i ++) {
		Triangle t; long n;
		if( 4 != fscanf( fp2, "%ld %ld %ld %ld",
				 & n, & t.p1, & t.p2, & t.p3))
			die( "TriangulatorDelaunay2d::triangulate(): ERROR\n"
			     "\t- could not read triangle %ld from\n"
			     "\t- '%s'", i+1, fname2);
		if( n != 3)
			die( "TriangulatorDelaunay2d::triangulate(): ERROR\n"
			     "\t- triangle #%ld has '%ld' vertices???\n"
			     "\t- '%s'", i+1, n, fname2);
		if( t.p1 < 0 || t.p1 >= npts ||
		    t.p2 < 0 || t.p2 >= npts ||
		    t.p3 < 0 || t.p3 >= npts)
			die( "TriangulatorDelaunay2d::triangulate(): ERROR\n"
			     "\t- invalid triangle #%ld in\n",
			     "\t- '%s'", i+1, fname2);
		// add the triangle
		res.push_back( t);
	}
	fclose( fp2);

	// remove both temporary files
	unlink (fname1);
	unlink (fname2);

	fprintf( stderr, "\t-end.\n");
	return res;
}
