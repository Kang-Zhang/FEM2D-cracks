#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <malloc.h>
#include <unistd.h>
#include <stdarg.h>
#include "Vector3d.hh"

class Point {
public:
	double x, y, z;		// coordinates of this point
	double fx, fy, fz;	// force vectors
};

class Triangle {
public:
	long p1, p2, p3;
};

static void usage( void)
{
	fprintf( stderr,
		 "Usage: gensphere [-h] [-r] n_points\n"
		 "           -h displays help\n"
		 "           -r applies relaxation so that the density of\n"
		 "              the points is approximately constant\n"
		);
}

static void die( char * str, ...)
{
	va_list ap;
	va_start( ap, str);
	fprintf( stderr, "gensphere: ");
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
static double radius = 1.0;
static long n_points = -1;
static long n_triangles = -1;
static Point * points = NULL;
static Triangle * triangles = NULL;
static const double EPS = 1e-6;
static const double threshold_disp = 1e-3;

void parse_arguments( int p_argc, char ** p_argv)
{
	argc = p_argc;
	argv = p_argv;

	while( 1)
	{
		int c = getopt( argc, argv, "hr");
		if( c == -1) {
			if( optind != argc-1) { usage(); exit(-1); }
			if( sscanf( argv[optind], "%ld", &n_points) != 1)
			{
				usage(); exit(-1);
			}
			if( n_points <= 0)
				die( "number of points must be > 0 !");
			break;
		}
		else if( c == 'h') { usage(); exit( 0); }
		else if( c == 'r') { relaxation_requested = 1; }
		else { usage(); exit(-1); }
	}

	radius = 1.0;
}

void get_force( Point & p2, Point & p1, double & fx, double & fy, double & fz)
{
	double vx = p2.x - p1.x;
	double vy = p2.y - p1.y;
	double vz = p2.z - p1.z;

	// calculate the length of the vector
	double l = sqrt( vx*vx + vy*vy + vz*vz);
	if( l < EPS) {
		fx = fy = fz = 0;
		return;
	}

	// normalize the vector
	double nx = vx / l;
	double ny = vy / l;
	double nz = vz / l;

	// calculate the amount of force, 1 being maximum and then tapering
	// off to 0 as the distance increases
	double f = (l+2)/(exp(l)+1);

	// calculate the force
	fx = f * nx;
	fy = f * ny;
	fz = f * nz;

	return;
}

void relax_points( void)
{
	// keep relaxing until we are happy
	double scale = 1.0;
	double scale_factor = 1.1;
	while( 1)
	{
		// zero out the force for every point
		for( long i = 0 ; i < n_points ; i ++)
			points[i].fx = points[i].fy = points[i].fz = 0.0;

		// calculate the force for every pair of points
		for( long i = 0 ; i < n_points ; i ++)
		{
			for( long j = i+1 ; j < n_points ; j ++)
			{
				double fx, fy, fz;
				get_force( points[i], points[j], fx, fy, fz);
				fx /= scale;
				fy /= scale;
				fz /= scale;
				points[i].fx += fx;
				points[i].fy += fy;
				points[i].fz += fz;
				points[j].fx -= fx;
				points[j].fy -= fy;
				points[j].fz -= fz;
			}
		}

		// calculate and apply the displacement for every point
		// based on the force acting on it. Also calculate for each
		// point by how much it has moved from its original
		// location (after it has been mapped back onto the sphere)
		double max_disp = 0.0;
		for( long i = 0 ; i < n_points ; i++)
		{
			// calculate the displacement
			double dx = points[i].fx;
			double dy = points[i].fy;
			double dz = points[i].fz;
			
			// apply the displacement
			points[i].x += dx;
			points[i].y += dy;
			points[i].z += dz;
		
			// map all points back onto the sphere
			double d = sqrt( points[i].x*points[i].x
					 + points[i].y*points[i].y
					 + points[i].z*points[i].z);
			if( d > EPS) {
				points[i].x *= radius/d;
				points[i].y *= radius/d;
				points[i].z *= radius/d;
			} else {
				double a = (drand48()*2-1)*M_PI;
				double b = drand48()*2*M_PI;
				points[i].x = sin(b)*cos(a)*radius;
				points[i].y = cos(b)*cos(a)*radius;
				points[i].z = sin(a)*radius;
				fprintf( stderr, "RRRRRRRRR........\n");
			}

			// calculate the displacement
			double disp = sqrt( dx*dx + dy*dy + dz*dz);
				
			// remember the maximum displacement
			if( disp > max_disp) max_disp = disp;
		}
		fprintf( stderr, "Max. disp = %.20f\r", max_disp);

		if( max_disp < threshold_disp) break;
		
		scale *= scale_factor;
	}
	fprintf( stderr, "\n");
}

int main( int argc, char ** argv)
{
	parse_arguments( argc, argv);

	fprintf( stderr, "Number of points requested: %ld\n", n_points);
	fprintf( stderr, "Relaxation requested: %s\n",
		 relaxation_requested ? "yes" : "no");

	// generate random points on a sphere
	// --------------------------------------------------
	points = new Point[ n_points];
	// calculate the minimum distance between points
	double min_dist = sqrt( M_PI * radius * radius) * 0.01;
	double min_dist_sq = min_dist * min_dist;
	fprintf( stderr, "min_dist = %f\n", min_dist);
	// set the index to 0
	long i = 0;
	while( 1) {
		// do we have enough points?
		if( i >= n_points) break;
		// generate a new random point
		double alpha = drand48() * M_PI * 2;
		double beta = drand48() * M_PI - M_PI/2;
		points[i].x = sin(alpha) * radius * cos(beta);
		points[i].y = cos(alpha) * radius * cos(beta);
		points[i].z = sin(beta) * radius;
		// make sure that this point is at leas some minimum distance
		// away from all other pointr
		int failure = 0;
		for( long j = 0 ; j < i ; j ++) {
			double d_sq =
				(points[i].x - points[j].x) *
				(points[i].x - points[j].x) +
				(points[i].y - points[j].y) *
				(points[i].y - points[j].y) +
				(points[i].z - points[j].z) *
				(points[i].z - points[j].z);
// 			if( d_sq < min_dist_sq) {
// 				// we are not far enough from other points,
// 				// try again with a new random point
// 				failure = 1;
// 				fprintf( stderr, "huh!!!\n");
// 				break;
// 			}
		}
		if( ! failure) i ++;
		continue;
	}

	// relax the points
	if( relaxation_requested) 
		relax_points();

	// print the points to a file
	// --------------------------------------------------
	FILE * fp = fopen( "/tmp/gencyl3.dat", "w"); assert( fp != NULL);
	fprintf( fp, "3\n%ld\n", n_points);
	for( long i = 0 ; i < n_points ; i ++)
		fprintf( fp, "%f %f %f\n", points[i].x,
			 points[i].y,
			 points[i].z);
	fclose(fp);

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
	assert( np == n_points);
	// read in number of triangles
	assert( 1 == fscanf( fp, "%ld", & n_triangles) && n_triangles >= 0);
	// skip number of faces
	long n_faces; assert( 1 == fscanf( fp, "%ld", & n_faces));
	// skip point coordinates
	for( long i = 0 ; i < n_points ; i ++) {
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
		assert( triangles[i].p1 >= 0 && triangles[i].p1 < n_points);
		assert( triangles[i].p2 >= 0 && triangles[i].p2 < n_points);
		assert( triangles[i].p3 >= 0 && triangles[i].p3 < n_points);
	}

#ifdef DONT_COMILE
	// project all points back onto the sphere (unbulge)
	// --------------------------------------------------
	for( long i = 0 ; i < n_points ; i ++) {
		double d = sqrt( points[i].x * points[i].x +
				 points[i].y * points[i].y);
		points[i].x *= radius/d;
		points[i].y *= radius/d;
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
#endif

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
	long j = 0; long n_removed = 0;
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
	// print the model
	printf( "# Requested number of points: %ld\n", n_points);
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
	printf( "crack_tip_element_size = 0.1\n");
	printf( "crack_tip_sub_element_size = 0.001\n");
	printf( "min_element_size = -0.01\n");
	printf( "max_element_size = 10.1\n");
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
		Vector3d norm( points[i].x, points[i].y, points[i].z); norm.normalize();

		printf( "%.6f %.6f %.6f %.3f %.3f %.3f %.4f %.4f %.4f %d\n",
			points[i].x, points[i].y, points[i].z,
			0.0, 0.0, 0.0,
			norm.x, norm.y, norm.z,
			1);
		printf( "%.6f %.6f %.6f %.3f %.3f %.3f %.4f %.4f %.4f %d\n",
			points[i].x*1.1, points[i].y*1.1, points[i].z*1.1,
			0.0, 0.0, 0.0,
			norm.x, norm.y, norm.z,
			0);
	}

	printf( "\n");

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
	printf( "\n");

	return 0;
}
