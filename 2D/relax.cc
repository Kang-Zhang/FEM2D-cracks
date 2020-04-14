#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>
#include <unistd.h>
#include <GL/glut.h>
#include <GL/gl.h>
#include <assert.h>
#include "BBox.hh"

// --------------------------------------------------------------------------
//
//                         global variables
//
// --------------------------------------------------------------------------
//
// number of points
static long n_points;
// number of iterations for relaxation
static long n_iterations;
// datastructure for a point
struct Point {
	double x, y;		// coordinates
	double tx, ty;		// force vector
};
// the array of points
Point * points = NULL;
double * radii = NULL;
// bounding box
int read_in_bounding_box = 0;
BBox bbox;
// desired distance between points
double desired_d;
// how fast the points move (should be <= 0.5)
double weight;
// command line options
static char display_graphics;
static char * fname;
static int show_grid = 1;
int verbose = 0;
int torus_coord = 0;
// jitter angle (force is rotated by +/- this angle * drand48()
// double jitter_angle = 30.0 * (M_PI/180.0);
double jitter_angle = 1.0 * (M_PI/180.0);
double radius_pert = 0.2;

// animations
long frame_count = 0;

static void die( char * str, ...)
{
	va_list ap;
	va_start( ap, str);
	vfprintf( stderr, str, ap);
	fprintf( stderr, "\n");
	va_end( ap);
	exit( -1);
}

static void read_points( void)
{
	// open input
	FILE * fp = stdin;
	if( fname != NULL) fp = fopen( fname, "r");
	if( fp == NULL) die( "Could not open '%s' for reading.", fname);
	
	// if bounding box should be read in from the file, do so
	if( read_in_bounding_box) {
		fprintf( stderr, "here!\n");
		double minx, miny, maxx, maxy;
		if( 4 != fscanf( fp, "%lf %lf %lf %lf",
				 & minx, & maxx, & miny, & maxy)) {
			// maybe it is the polygonal bounding box
			long n;
			if( 1 != fscanf( fp, "poly: %ld", & n))
				die( "Bounding box description missing "
				     "in file.");
			fprintf( stderr, "nbound = %ld\n", n);
			bbox.reset();
			for( long i = 0 ; i < n ; i ++) {
				double x, y;
				if( 2 != fscanf( fp, "%lf %lf", & x, & y))
					die( "Bounding box invalid.");
				bbox.add( x, y);
			}
		} else {
			bbox.reset();
			bbox.add( minx, miny);
			bbox.add( maxx, miny);
			bbox.add( maxx, maxy);
			bbox.add( minx, maxy);
		}
	} else {
		bbox.reset();
		bbox.add( -0.5, -0.5);
		bbox.add(  0.5, -0.5);
		bbox.add(  0.5,  0.5);
		bbox.add( -0.5,  0.5);
	}

	fprintf( stderr, "The bounding box: ");
	bbox.print( stderr);
	fprintf( stderr, "\n");

	// obtain the number of points
	if( 1 != fscanf( fp, "%ld", &n_points))
		die( "Could not read in number of points.");
	if( n_points < 1)
		die( "There must be more than 0 points.");
	
	// allocate room for points
	points = new Point[ n_points];
	if( points == NULL) die( "Could not allocate enough memory.");
	
	// read in all points
	for( long i = 0 ; i < n_points ; i ++)
		if( 2 != fscanf( fp, "%lf %lf", & points[i].x, & points[i].y))
			die( "Point %ld: could not read in coordinates.", i);
	
	// done
	return;
}

class Square {
public:
	long n;
	long * data;
	Square() {
		n = 0;
		data = NULL;
	}
	void add( long k) {
		if( data == NULL) data = (long *) malloc( sizeof( long));
		else data = (long *) realloc( data, sizeof( long) * (n+1));
		data[n] = k;
		n ++;
	}
	void reset( void) {
		if( data != NULL) {
			free( data);
			data = NULL;
			n = 0;
		}
	}
};

class Grid {
public:
	double minx, maxx, miny, maxy;
	long n_horiz, n_vert;
	Square ** squares;

	Grid( double pminx, double pmaxx, double pminy, double pmaxy,
	      long phoriz, long pvert)
	{
		// copy parameters
		minx = pminx; maxx = pmaxx; miny = pminy; maxy = pmaxy;
		n_horiz = phoriz; n_vert = pvert;

		// allocate squares
		squares = (Square **) malloc( sizeof( Square *) * n_vert);
		for( long i = 0 ; i < n_vert ; i ++)
			squares[i] = new Square[ n_horiz];
	}

	void add( long k) {
		double x = points[k].x;
		double y = points[k].y;
		double square_width = (maxx-minx) / n_horiz;
		double square_height = (maxy-miny) / n_vert;

		long gx = long( (x-minx) / square_width);
		long gy = long( (y-miny) / square_height);

		// make sure the indices are correct
		if( gx < 0) gx = 0;
		if( gx > n_horiz-1) gx = n_horiz - 1;
		if( gy < 0) gy = 0;
		if( gy > n_vert-1) gy = n_vert -1;

		squares[gy][gx].add( k);
	}

	void reset( void) {
		for( long gx = 0 ; gx < n_horiz ; gx ++)
			for( long gy = 0 ; gy < n_vert ; gy ++)
				squares[gy][gx].reset();
	}

	void report_stats( void) {
		long max = squares[0][0].n;
		long min = squares[0][0].n;
		for( long gx = 0 ; gx < n_horiz ; gx ++)
			for( long gy = 0 ; gy < n_vert ; gy ++) {
				if( squares[gy][gx].n > max)
					max = squares[gy][gx].n;
				if( squares[gy][gx].n < min)
					min = squares[gy][gx].n;
			}
		fprintf( stderr, "grid: min / max = %ld / %ld\n", min, max);
	}
};

static Grid * grid = NULL;

static void calculate_force( long i1, long i2, double & d)
// calculates the force between points i1 and i2
{
	// distance between identical points is invalid
	if( i1 == i2) {
		d = -1;
		return;
	}

	// calculate the force between points i1 and i2 and add it to
	// points[i1]
	Point * p = & points[i1];
	Point * q = & points[i2];
	double desired_d = radii[i1] + radii[i2];

	// calculate the vector PQ
	double w = bbox.maxx - bbox.minx;
	double h = bbox.maxy - bbox.miny;

	// calculate vector (dx,dy) frpm p1 to p2, taking into account
	// the torus properties of the plane
	double dx = q->x - p->x;
	if( torus_coord) {
		if( fabs( dx) > w/2) {
			if( dx < 0) dx += w; else dx -= w;
		}
	}
	double dy = q->y - p->y;
	if( torus_coord) {
		if( fabs( dy) > h/2) {
			if( dy < 0) dy += h; else dy -= h;
		}
	}

	// calculate distance PQ
	double d_sq = dx * dx + dy * dy;
	d = sqrt( d_sq);
	// if p and q are more than desired_d apart,
	// they don't have any effect on each other
	// and therefore can be ignored
	if( d > desired_d) return;
	// calculate the force between the two points
	// (in percentage, 1.0 =  max. force, occuring
	// when the points coincide, 0.0 = min force, when
	// the points are desired_d apart.
	double force = 1 - d/desired_d;
	force *= force;
	// compute a unit vector from dx,dy
	double ux, uy;
	if( d / desired_d < 0.001) {
		// points p and q almost coincide, we
		// will therefore create a random vector
		double alpha = drand48() * M_PI * 2;
		ux = sin( alpha) + q->x;
		uy = cos( alpha) + q->y;
	} else {
		ux = dx/d;
		uy = dy/d;
	}
	// jitter the unit vector
	double alpha = (drand48() * 2 - 1) * jitter_angle;
	double vx = ux * cos( alpha) - uy * sin( alpha);
	double vy = ux * sin( alpha) + uy * cos( alpha);
	
	// add force * unit_vector to tx,ty of p
	p-> tx -= force * vx;
	p-> ty -= force * vy;
}

static void calculate_force_total( long k, long gx, long gy, double & min_d)
// calculates the total force acting on point[k] from all neighbors
// - neighbors are points in the current grid[gy][gx] and all 8 adjacent
//   grid squares
{
	min_d = -1;
	// calculate the force from all surrounding points
	for( long dx = -1 ; dx <= 1 ; dx ++) {
		for( long dy = -1 ; dy <= 1 ; dy ++) {
			long hx = (gx + dx + grid->n_horiz) % grid->n_horiz;
			long hy = (gy + dy + grid->n_vert) % grid->n_vert;
			Square * sq = & grid-> squares[hy][hx];
			for( long m = 0 ; m < sq-> n ; m++) {
				double d;
				calculate_force( k, sq-> data[m], d);
				if( d >= 0 && (d < min_d || min_d < 0))
					min_d = d;
			}
		}
	}
}

static double relax_points( void)
{
	// populate the grid
	grid-> reset();
	for( long i = 0 ; i < n_points ; i ++)
		grid-> add( i);
//	grid-> report_stats();
	
	// zero out all forces on all points
	for( long i = 0 ; i < n_points ; i ++)
		points[i].tx = points[i].ty = 0.0;

	// for every point find all neighbors (from its grid square and
	// its surrounding points), calculate the force these generate
	// and add it to the total force
	double min_d = -1;
	for( long gx = 0 ; gx < grid-> n_horiz ; gx ++) {
		for( long gy = 0 ; gy < grid-> n_vert ; gy ++) {
			for( long k = 0 ; k < grid->squares[gy][gx].n;k++) {
				double d;
				calculate_force_total(
					grid-> squares[gy][gx].data[k],
					gx, gy, d);
				if( (d < min_d || min_d < 0) && d >= 0)
					min_d = d;
			}
		}
	}

	// adjust the coordinates by adding the force to x & y, but make
	// sure that the resulting points are in the bounding box
	// the bounding box
	for( long i = 0 ; i < n_points ; i ++)
	{
		points[i].x += weight * points[i].tx;
		points[i].y += weight * points[i].ty;

		if( torus_coord) {
			while( points[i].x < bbox.minx)
				points[i].x += (bbox.maxx-bbox.minx);
			while( points[i].x >= bbox.maxx)
				points[i].x -= (bbox.maxx-bbox.minx);
			while( points[i].y < bbox.miny)
				points[i].y += (bbox.maxy-bbox.miny);
			while( points[i].y >= bbox.maxy)
				points[i].y -= (bbox.maxy-bbox.miny);
		} else {
			bbox.snap_inside( points[i].x, points[i].y);
//			if( points[i].x < bbox.minx) points[i].x = bbox.minx;
//			if( points[i].x > bbox.maxx) points[i].x = bbox.maxx;
//			if( points[i].y < bbox.miny) points[i].y = bbox.miny;
//			if( points[i].y > bbox.maxy) points[i].y = bbox.maxy;
		}
	}
	
	return min_d;
}

static void output_points( void)
{
	printf( "%ld\n\n", n_points);
	for( long i = 0 ; i < n_points ; i ++)
		printf( "%.10f %.10f\n", points[i].x, points[i].y);
	
	delete [] points;
}

static void usage( void)
{
	printf(
		"Usage: relax [-h] [-g] [-f filename] [-i n_iterations] \n"
		"             [-v] [-b] [-t] [-j jitter] [-r radius perturbation]\n"
		"\n"
		"       -h displays this help\n"
		"       -g produces graphical result\n"
		"       -f reads the input from <filename>\n"
		"       -i sets the number of iterations\n"
		"       -b read in bounding box from the file\n"
		"       -v verbose\n"
		"       -j angle jitter in degrees (0..180, default 1)\n"
		"       -r repulsion radius perturbation (0..1, default 0.2)\n"
		"\n"
		);
}


static void parse_arguments( int & argc, char ** argv)
{
	// set defaults
	display_graphics = 0;
	fname = NULL;
	n_iterations = 200;
	
	while( 1)
	{
		int c = getopt( argc, argv, "ghf:i:btj:r:");
		if( c == -1) break;
		if( c == 'g') { display_graphics = 1; continue;}
		if( c == 'h') { usage(); exit(0); }
		if( c == 'f') { fname = optarg; continue; }
		if( c == 'b') { read_in_bounding_box = 1; continue; }
		if( c == 'v') { verbose = 1; continue; }
		if( c == 't') { torus_coord = 1; continue; }
		if( c == 'i') {
			if( 1 != sscanf( optarg, "%ld", & n_iterations)) {
				usage();
				die( "Invalid number of iterations.");
			}
			if( n_iterations < 1) {
				usage();
				die( "Number of iterations < 1 !?!");
			}
			continue;
		}
		if( c == 'j') {
			if( 1 != sscanf( optarg, "%lf", & jitter_angle)) {
				usage();
				die( "Invalid jitter angle.");
			}
			jitter_angle *= M_PI / 180;
			continue;
		}
		if( c == 'r') {
			if( 1 != sscanf( optarg, "%lf", & radius_pert)) {
				usage();
				die( "Invalid radius perturbation.");
			}
			continue;
		}
		usage(); die( "Bad syntax.");
	}
	
	if( display_graphics)
		fprintf( stderr, "Displaying graphics.\n");
	fprintf( stderr, "Number of iterations: %ld\n", n_iterations);
}

static void redraw( void)
{
	// clear the screen
	glDrawBuffer( GL_BACK);
	glClear( GL_COLOR_BUFFER_BIT);

	if( show_grid) {
		double w = (bbox.maxx-bbox.minx) / grid-> n_horiz;
		double h = (bbox.maxy-bbox.miny) / grid-> n_vert;
		for( long gx = 0 ; gx < grid-> n_horiz ; gx ++) {
			for( long gy = 0 ; gy < grid-> n_vert ; gy ++) {
				double g1 = 0.15;
				double g2 = 0.30;
				if( (gx + gy) % 2 == 0)
					glColor3f( g1, g1, g1);
				else
					glColor3f( g2, g2, g2);
				
				double x = bbox.minx + gx * w;
				double y = bbox.miny + gy * h;
				glBegin( GL_POLYGON);
				glVertex2d( x, y);
				glVertex2d( x+w, y);
				glVertex2d( x+w, y+h);
				glVertex2d( x, y+h);
				glEnd();
			}
		}
	}

	glPointSize( 6);
	glColor3f( 1, 1, 1);
	glBegin( GL_POINTS);
	for( long i = 0 ; i < n_points ; i ++)
		glVertex2d( points[i].x,
			    points[i].y);
	glEnd();

	// swap buffers
        glutSwapBuffers();

	// glutPostRedisplay();
}

static void reshape_func( int width, int height)
{	
	// set the viewport to use the entire drawing area
	glViewport( 0, 0, width, height);

	// set up the view matrix to be the identity matrix
	glMatrixMode( GL_MODELVIEW );
	glLoadIdentity( );
	
	// set the projection matrix
	glMatrixMode( GL_PROJECTION );
	glLoadIdentity();
	glOrtho( bbox.minx-0.1*(bbox.maxx-bbox.minx),
		 bbox.maxx+0.1*(bbox.maxx-bbox.minx),
		 bbox.miny-0.1*(bbox.maxy-bbox.miny),
		 bbox.maxy+0.1*(bbox.maxy-bbox.miny),
		 -1, 1);

	redraw();

	if( verbose)
		fprintf( stderr, "reshape %dx%d\n", width, height);
}

static void disp_func( void)
{	
	glDisable( GL_DITHER);
	glDisable( GL_LIGHTING);
	glDisable( GL_DEPTH_TEST);

	// antialised points
	glShadeModel( GL_SMOOTH);
	glPointSize( 2.0);
	glEnable(GL_POINT_SMOOTH);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	
	// pick the color for background
	glClearColor( 0.1, 0, 0, 0);

	redraw();
}

static void idle_func( void)
{
	double d = relax_points();
	printf( "Frame %06ld: min dist = %f\n", frame_count, d);
	frame_count ++;
	
	redraw();
}

enum Menu { MENU_STEP, MENU_RUN, MENU_STOP, MENU_QUIT, MENU_GRID };

void menu_func( int menu)
// ----------------------------------------------------------------------
// called when some of the menus are selected
// ----------------------------------------------------------------------
{
	double d;

	switch( menu) {
	case MENU_STEP:
		d = relax_points();
		printf( "Relaxed points. min dist = %f\n", d);
		redraw();
		break;
	case MENU_RUN:
		glutIdleFunc( idle_func);
		break;
	case MENU_STOP:
		glutIdleFunc( NULL);
		break;
	case MENU_QUIT:
		exit( 0);
		break;
	case MENU_GRID:
		show_grid = ! show_grid;
		glutPostRedisplay();
		break;
	}
}

static void create_ui( int & argc, char ** argv)
{
	glutInitWindowSize( 500, 500);
	glutInit( & argc, argv);
        glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowSize( 768, 768);
        char title[ 4096];
        sprintf( title, "%s %s", argv[0], fname);
        int main_window = glutCreateWindow( title);
        assert( main_window);
        glutSetWindow( main_window);

        glutDisplayFunc( disp_func);
        glutReshapeFunc( reshape_func);
        
        // create a popup menu
        int m = glutCreateMenu( menu_func);
        glutSetMenu( m);
	glutAddMenuEntry( "Step", MENU_STEP);
	glutAddMenuEntry( "Run", MENU_RUN);
	glutAddMenuEntry( "Stop", MENU_STOP);
	glutAddMenuEntry( "Quit", MENU_QUIT);
	glutAddMenuEntry( "Grid", MENU_GRID);
        glutAttachMenu( GLUT_RIGHT_BUTTON);

        glutMainLoop();
}


int main( int argc, char ** argv)
{
	// parse command line arguments
	parse_arguments( argc, argv);
	
	// read in the points
	read_points();
	
	// calculate the desired radius
	double area = bbox.get_area();
	desired_d = 2*sqrt(area/(n_points*M_PI));
	desired_d = 1.6 * desired_d;
	weight = 0.3333333 * desired_d;
	weight = 1.0 * weight;
	if( verbose)
		fprintf( stderr, "# desired distance: %f\n", desired_d);

	// calculate the radii of each point:
	//    - optimum is desired_d/2
	//    - perturbed with error +/-err %
	radii = new double[n_points];
	for( long i = 0 ; i < n_points ; i ++) {
		radii[i] = desired_d / 2;
		radii[i] *= 1 + (drand48() * 2 - 1) * radius_pert;
	}

	// create a grid
	long grid_h = long( (bbox.maxx-bbox.minx) / desired_d + 1);
	long grid_v = long( (bbox.maxy-bbox.miny) / desired_d + 1);
	grid = new Grid( bbox.minx, bbox.maxx, bbox.miny, bbox.maxy,
			 grid_h, grid_v);
	
	// if graphics has to be displayed - create the widgets
	if( display_graphics) create_ui( argc, argv);
	
	// relax points
	for( long i = 0 ; i < n_iterations ; i ++) {
		double d = relax_points();
		if( verbose)
			fprintf( stderr, "%4ld) %f\n", i, d);
	}
	
	// output the points
	output_points();
	
	// return success
	return 0;
}
