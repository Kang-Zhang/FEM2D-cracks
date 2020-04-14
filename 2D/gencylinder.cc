#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <malloc.h>
#include <unistd.h>
#include <stdarg.h>
#include <string.h>
#include "Vector3d.hh"

const double E = 100.0;
const double Poisson = 0.45;
static double height = 2;
static double radius_top = 0.4;
static double radius_bottom = 0.4;
static double thick_top = 0.1;
static double thick_bottom = 0.1;
#define MIN(a,b) ( a >= b ? a : b)


double interpolate( double min_x,
		    double max_x,
		    double val1,
		    double val2,
		    double x)
{
	double p = (x - min_x) / (max_x - min_x);
	if( p < -1) p = -1; else if( p > 1) p = 1;
	return (1-p) * val1 + p * val2;
}

double radius_func( double z)
// returns the value of the radius at the specified z
{
	double min_z = - height/2;
	double max_z = height/2;
	return radius_bottom * (z - min_z) / (max_z - min_z) +
		radius_top * (1 - (z - min_z) / (max_z - min_z));
}


class Point {
public:
	double x, y, z;		// coordinates of this point
};

class Triangle {
public:
	long p1, p2, p3;
};

static void usage( void)
{
	fprintf( stderr,
		 "Usage: gencylinder [-h] [-r] n_points\n"
		 "                   [rad_bottom=num] [rad_top=num] \n"
		 "                   [thick_bottom=num] [thick_top=num]\n"
		 "\n"
		 "           -h displays help\n"
		 "           -r applies relaxation so that the density of\n"
		 "              the points is approximately constant\n"
		);
}

static void die( char * str, ...)
{
	va_list ap;
	va_start( ap, str);
	fprintf( stderr, "gencylinder: ");
	vfprintf( stderr, str, ap);
	fprintf( stderr, "\n");
	va_end( ap);
	exit( -1);
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

static int argc;
static char ** argv;
static int relaxation_requested = 0;
static long n_points = -1;
static long n_triangles = -1;
static Point * points = NULL;
static Triangle * triangles = NULL;
static const double EPS = 0.000001;

void parse_arguments( int p_argc, char ** p_argv)
{
	argc = p_argc;
	argv = p_argv;
	int n = 1;

	while( 1)
	{
		if( strcmp( argv[ n], "-h") == 0) {
			usage(); exit(0);
		} else if( strcmp( argv[ n], "-r") == 0) {
			relaxation_requested = 1;
		} else if( n_points == -1) {
			if( sscanf( argv[ n], "%ld", &n_points) != 1)
			{
				usage(); exit(-1);
			}
			if( n_points <= 0)
				die( "number of points must be > 0 !");
		} else {
			char * eq = index( argv[n], '=');
			if( eq == NULL) {
				fprintf( stderr, "%s is not an assignment\n",
					 argv[n]);
				usage(); exit( -1);
			}
			char var[ 4096];
			strncpy( var, argv[n], eq - argv[n]);
			var[ eq - argv[n]] = 0;
			fprintf( stderr, "arg: '%s' = '%s'\n", var, eq+1);
			double val;
			if( sscanf( eq+1, "%lf", & val) != 1)
				die( "assignment %s not a number", argv[n]);
			if( strcmp( "rad_top", var) == 0)
				radius_top = val;
			else if( strcmp( "rad_bottom", var) == 0)
				radius_bottom = val;
			else if( strcmp( "thick_top", var) == 0)
				thick_top = val;
			else if( strcmp( "thick_bottom", var) == 0)
				thick_bottom = val;
			else
				die( "assignment %s not a valid variable",
				     argv[n]);
		}
		n++;
		if( n >= argc) break;
	}
	fprintf( stderr, "%f %f %f %f\n",
		 radius_top, radius_bottom, thick_top, thick_bottom);
}

int main( int argc, char ** argv)
{
	parse_arguments( argc, argv);

	fprintf( stderr, "Number of points requested: %ld\n", n_points);
	fprintf( stderr, "Relaxation requested: %s\n",
		 relaxation_requested ? "yes" : "no");

	// generate random points on a cylinder
	// --------------------------------------------------
	points = new Point[ n_points+2];
	// calculate the minimum distance between points
	double width = 2*M_PI;
	double min_dist = 0.1 * width / n_points ;
	if( 0.1 * height / n_points < min_dist)
		min_dist = 0.1 * height / n_points ;
	double min_dist_sq = min_dist * min_dist;
	fprintf( stderr, "min_dist = %f\n", min_dist);
	// set the index to 0
	long i = 0;
	while( 1) {
		// do we have enough points?
		if( i >= n_points) break;
		// generate a new random point
		points[i].x = drand48()*width;
		points[i].y = drand48()*height - height/2;
		// make sure that this point is at leas some minimum distance
		// away from all other pointr
		int failure = 0;
		for( long j = 0 ; j < i ; j ++) {
			double d_sq =
				(points[i].x - points[j].x) *
				(points[i].x - points[j].x) +
				(points[i].y - points[j].y) *
				(points[i].y - points[j].y);
			if( d_sq < min_dist_sq) {
				// we are not far enough from other points,
				// try again with a new random point
				failure = 1;
				fprintf( stderr, "huh!!!\n");
				break;
			}
		}
		if( ! failure) i ++;
		continue;
	}

	// output the points to a temporary file
	// --------------------------------------------------
	FILE * fp = fopen( "/tmp/gencyl1.dat", "w" ); assert( fp != NULL );
	// bounding box:
	fprintf( fp, "%f %f %f %f\n",
		 0.0,  width,
		 -height/2, height/2);
	// number of points
	fprintf( fp, "%ld\n", n_points);
	// the points themselves
	for( long i = 0 ; i < n_points ; i ++)
		fprintf( fp, "%f %f\n", points[i].x, points[i].y);
	fclose( fp );

	// relax the points if requested
	// --------------------------------------------------
	if( relaxation_requested) {
		execute( "./relax -b < /tmp/gencyl1.dat > /tmp/gencyl2.dat");
	} else {
		execute( "(head -1 /tmp/gencyl1.dat ; "
			 " tail +3 /tmp/gencyl1.dat) > /tmp/gencyl2.dat");
	}

	// read in the points
	// --------------------------------------------------
	fp = fopen( "/tmp/gencyl2.dat", "r"); assert( fp != NULL);
	long tmp; assert( 1 == fscanf( fp, "%ld", & tmp));
	assert( tmp == n_points);
	for( long i = 0 ; i < n_points ; i ++)
		assert( 2 == fscanf( fp, "%lf %lf", & points[i].x,
				     & points[i].y));
	fclose( fp);

	// project the points on a bulged out cylinder
	// and write them out to a file
	// --------------------------------------------------
	fp = fopen( "/tmp/gencyl3.dat", "w"); assert( fp != NULL);
	fprintf( fp, "3\n%ld\n", n_points + 2);
	for( long i = 0 ; i < n_points ; i ++) {
		// calculate the real z
		double z = points[i].y;
		// calculate bulged radius (as a function of z)
		double br = -16*z*z/(5*height*height)+9/5.0;
		// calculate the x and y (bulged)
		double x = br * sin( points[i].x);
		double y = br * cos( points[i].x);
		// print out these in the file
		fprintf( fp, "%f %f %f\n", x, y, z);
		points[i].x = x;
		points[i].y = y;
		points[i].z = z;
	}
	// print extra two points at the top & bottom of the cylinder
	fprintf( fp, "0 0 %f\n", height*3/4);
	fprintf( fp, "0 0 %f\n", -height*3/4);
	fclose( fp);

	// run qhull on the last file
	// --------------------------------------------------
	execute( "qhull -o < /tmp/gencyl3.dat > /tmp/gencyl4.dat");

	// read in the results (the triangles)
	// --------------------------------------------------
	fp = fopen( "/tmp/gencyl4.dat", "r"); assert( fp != NULL);
	// skip dimension (make sure it is 3)
	long dim; assert( 1 == fscanf( fp, "%ld", & dim)); assert( dim == 3);
	// skip number of points
	long np; assert( 1 == fscanf( fp, "%ld", & np));
	assert( np == n_points+2);
	// read in number of triangles
	assert( 1 == fscanf( fp, "%ld", & n_triangles) && n_triangles >= 0);
	// skip number of faces
	long n_faces; assert( 1 == fscanf( fp, "%ld", & n_faces));
	// skip point coordinates
	for( long i = 0 ; i < n_points + 2 ; i ++) {
		double x, y, z;
		assert( 3 == fscanf( fp, "%lf %lf %lf", &x, &y, &z));
	}
	// read in the triangles
	triangles = new Triangle[n_triangles];
	for( long i = 0 ; i < n_triangles ; i ++) {
		// skip vertex count
		long vc; assert( 1 == fscanf( fp, "%ld", &vc) && vc == 3);
		// read in the vertices
		assert( 3 == fscanf( fp, "%ld %ld %ld",
				     & triangles[i].p1,
				     & triangles[i].p2,
				     & triangles[i].p3));
		assert( triangles[i].p1 >= 0 && triangles[i].p1 < n_points+2);
		assert( triangles[i].p2 >= 0 && triangles[i].p2 < n_points+2);
		assert( triangles[i].p3 >= 0 && triangles[i].p3 < n_points+2);
	}

	// project all points back onto the cylinder of unit radius
	// --------------------------------------------------------
	for( long i = 0 ; i < n_points ; i ++) {
		double d = sqrt( points[i].x * points[i].x +
				 points[i].y * points[i].y);
		points[i].x /= d;
		points[i].y /= d;
	}

	// clip out all triangles touching either of the last
	// two points
	// --------------------------------------------------
	long n_removed = 0;
	long ptr = 0;
	for( long i = 0 ; i < n_triangles ; i ++) {
		Triangle * t = triangles + i;
		if( t->p1 < n_points && t->p2 < n_points && t->p3 < n_points) {
			triangles[ ptr] = triangles[ i];
			ptr ++;
		} else {
			n_removed ++;
		}
	}
	fprintf( stderr, "Clipped out %ld triangles\n", n_removed);
	n_triangles -= n_removed;

	// count references to every node
	// --------------------------------------------------
	long * ref = new long[n_points];
	for( long i = 0 ; i < n_points ; i ++) ref[i] = 0;
	for( long i = 0 ; i < n_triangles ; i ++) {
		ref[triangles[i].p1]++;
		ref[triangles[i].p2]++;
		ref[triangles[i].p3]++;
	}
	
	// prepare a table of mappings to help us with
	// deleting unreferenced points
	// --------------------------------------------------
	long * map = new long[n_points];
	long count = 0;
	for( long i = 0 ; i < n_points ; i ++) {
		if( ref[i] != 0) {
			map[i] = count;
			count ++;
		} else {
			map[i] = -1;
		}
	}

	// delete all unreferenced nodes
	// --------------------------------------------------
	long j = 0; n_removed = 0;
	for( long i = 0 ; i < n_points ; i ++) {
		if( ref[i] > 0) {
			points[j] = points[i];
			j ++;
		} else {
			n_removed ++;
		}
	}
	fprintf( stderr, "Removed %ld unreferenced nodes.\n", n_removed);
	n_points -= n_removed;

	// renumber the triangles according to the
	// calculated map
	// --------------------------------------------------
	for( long i = 0 ; i < n_triangles ; i ++) {
		triangles[i].p1 = map[ triangles[i].p1];
		triangles[i].p2 = map[ triangles[i].p2];
		triangles[i].p3 = map[ triangles[i].p3];
	}

	// print out results
	// --------------------------------------------------
	// print the model
	printf( "# Cylindrical model\n");
	// print the model
	printf( "# Number of wedges:  %ld\n", n_triangles);
	printf( "\n");
	printf( "\n");
	printf( "sim_time_total = 100.0\n");
	printf( "min_time_step = 0.1\n");
	printf( "max_time_step = 1.0000\n");
	printf( "sim_time_curr = 0.0\n");
	printf( "growth_x = 0.1\n");
	printf( "growth_y = 0.1\n");
	printf( "growth_z = 0.0\n");
	printf( "shrink_top_t0 = 0.0\n");
	printf( "shrink_top_t1 = 0.0\n");
	printf( "shrink_top_val = 0.0\n");
	printf( "shrink_bot_t0 = 0.0\n");
	printf( "shrink_bot_t1 = 0.0\n");
	printf( "shrink_bot_val = 0.0\n");
	printf( "shrink_height_t0 = 0.0\n");
	printf( "shrink_height_val0 = 0.0\n");
	printf( "shrink_height_t1 = 0.0\n");
	printf( "shrink_height_val1 = 0.0\n");
	printf( "gravity_x = 0.0\n");
	printf( "gravity_y = 0.0\n");
	printf( "gravity_z = 0.0\n");
	printf( "\n");
	printf( "fracture_tip_inertia = 1\n");
	printf( "max_break_size = 0.0001\n");
	printf( "min_refine_size = 0.01\n");
	printf( "crack_tip_element_size = 0.0001\n");
	printf( "min_element_size = -0.01\n");
	printf( "max_element_size = 1\n");
	printf( "precision = 1e-4\n");
	printf( "\n");
	printf( "error_radius = 1\n");
	printf( "min_error_radius = 0.1\n");
	printf( "\n");
	printf( "progressive_save = true\n");
	printf( "progressive_save_fmask = \"/scratch/federl/%%N-%%T.res\"\n");
	printf( "progressive_save_skip = 0\n");
	printf( "\n");
	printf( "# description of how to randomize material properties\n");
	printf( "BEGIN material_properties\n");
	printf( "\tyoung_modulus constant 1000.0\n");
	printf( "\tpoisson_ratio constant 0.3\n");
	printf( "\tyield_stress constant 1000.0\n");
	printf( "\tfracture_toughness constant 0.01\n");
	printf( "END material_properties\n");
	printf( "\n");
	printf( "%ld # number of nodes\n", n_points*2);
	
	for( long i = 0 ; i < n_points ; i ++) {
		double r_in = interpolate( - height/2, height/2,
					   radius_bottom,
					   radius_top,
					   points[i].z);
		double r_out = interpolate( - height/2, height/2,
					    radius_bottom + thick_bottom,
					    radius_top + thick_top,
					    points[i].z);

		Vector3d norm( points[i].x, points[i].y, 0); norm.normalize();
		
		printf( "%.10f %.10f %.10f 0 0 0 %.10f %.10f %.10f 1\n",
			points[i].x * r_in, points[i].y * r_in, points[i].z,
			norm.x, norm.y, norm.z);
 		printf( "%.10f %.10f %.10f 0 0 0 %.10f %.10f %.10f 0\n",
 			points[i].x * r_out, points[i].y * r_out, points[i].z,
			norm.x, norm.y, norm.z);
	}

	printf( "\n");
	printf( "%ld # number of elements\n", n_triangles);
	printf( "\n");
	for( long i = 0 ; i < n_triangles ; i ++) {
		printf( "fivewall %ld %ld %ld %ld %ld %ld %ld\n",
			i,
			triangles[i].p2*2+1,
			triangles[i].p1*2+1,
			triangles[i].p3*2+1,
			triangles[i].p2*2,
			triangles[i].p1*2,
			triangles[i].p3*2
			);
	}

	printf( "\n");
	printf( "0 # number of fracture tips\n");

	return 0;
}
