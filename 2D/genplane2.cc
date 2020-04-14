#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdarg.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include "Distribution.hh"
#include "MapTexture2D.hh"
#include "die.hh"

const double E = 100.0;
const double poisson = 0.3;
const double HEIGHT = 0.085;
const double thick_bl = HEIGHT/1;
const double thick_br = HEIGHT/1;
const double thick_tr = HEIGHT/1;
const double thick_tl = HEIGHT/1;
// bounding box width and height
double bb_width=-1, bb_height=-1;
// growth directions
double growth_r1 = 1;
double growth_r2 = 0;
double growth_alpha = (M_PI / 8) * 1;
MapTexture2D * hmap = NULL;
MapTexture2D * bmap = NULL;
const double thick_rnd = 0.0;

double get_random_young_mod( void) {
	return Distribution::get_normal_rnd( 100, 5);
}

double get_random_poisson_mod( void) {
	return Distribution::get_normal_rnd( 0.3, 0.03);
}

double get_random_yield_stress( void) {
	return Distribution::get_normal_rnd( 10, 0.5);
}

class Node {
public:
	double x, y, z;
	double dx, dy, dz;
	double gx, gy, gz;
	int type;
	long num;

	Node( long pnum, double px, double py, double pz, int ptype) {
		num = pnum;
		x = px; y = py; z = pz;
		type = ptype;
		dx = 0; dy = 0; dz = 0;
		gx = 0; gy = 0; gz = 0;
	}
};

class Wedge {
public:
	Node * n1, * n2, * n3, * n4, * n5, * n6;
	double young_mod, poisson_mod, yield_stress;
	long num;

	Wedge( long pnum, double young, double poisson, double yield,
	       Node * pn1, Node * pn2, Node * pn3,
	       Node * pn4, Node * pn5, Node * pn6) {
		num = pnum;
		young_mod = young; poisson_mod = poisson; yield_stress = yield;
		n1 = pn1; n2 = pn2; n3 = pn3;
		n4 = pn4; n5 = pn5; n6 = pn6;
	}
};

static long req_n_wedges = -1;
static long n_layers = -1;
static long n_wedges = 0;
static Wedge ** wedges = NULL;
static long n_nodes = 0;
static Node ** nodes = NULL;

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

static void usage( void)
{
	fprintf( stderr,
		 "Usage: genplane [-n n_ponts] [-m n_layers] [-t type]\n"
		 "                [-d width height] [-gvf alpha r1 r2]\n"
		 "                [-hmap fname min max]\n"
		 "                [-bmap fname min max]\n"
		 "\n"
		 "      -n    specifies number of points to be used,\n"
		 "            default is 100\n"
		 "      -m    specifies number of layers, defaults is 1\n"
		 "      -t    requests type of point distrd. subdivision:\n"
		 "               random - random distribution\n"
		 "               repelled - random & repelled (default)\n"
		 "               ut2 - two triangles per square\n"
		 "               ut2 - two triangles per square\n"
		 "               ut4 - four triangles per square\n"
		 "               hex - hexagonal subdivision\n"
		 "      -d    width & height of the patch (default 1x1)\n"
		 "     -gv    growth vector parameters (default 0 1 1)\n"
		 "               alpha = angle in degrees\n"
		 "                  r1 = rate of growth in direction alpha\n"
		 "                  r2 = rate of growth in perpend. dir.\n"
		 "   -hmap    defines the height of the surface through a\n"
		 "            texture map (default is constant 0.05)\n"
		 "             fname = name of the file with the texture\n"
		 "                     (must be gray-scale SGI RGB image\n"
		 "               min = value for pixel value 0\n"
		 "               max = value for pixel value 255\n"
		 "   -bmap    defines the texture map for the profile of\n"
		 "            of the bakcground (default = constant 0)\n"
		 "             fname = name of the file with the texture\n"
		 "                     (must be gray-scale SGI RGB image\n"
		 "               min = value for pixel value 0\n"
		 "               max = value for pixel value 255\n"
		 "      -h    displays this help\n"
		);
	exit( -1);
}

enum SubdivisionType { Unspecified, UnifTri2, UnifTri4, Hexagonal, Rand };
static SubdivisionType subdivision_type = Unspecified;
int relax = 1;

void parse_arguments( int argc, char ** argv)
{
	// default parameters
	req_n_wedges = 100;
	n_layers = 1;
	subdivision_type = Rand; relax = 1;
	growth_alpha = 0; growth_r1 = 1; growth_r2 = 1;
	char * hmap_fname=NULL; double hmap_min = 0.05; double hmap_max = 0.05;
	char * bmap_fname=NULL; double bmap_min = 0.0; double bmap_max = 0.0;
	bb_width = 1;
	bb_height = 1;
	
	int ag = 1;

	while( 1) {
		if( ag >= argc) break;
//		fprintf( stderr, "Processing '%s'\n", argv[ag]);

		if( strcmp( argv[ag], "-n") == 0) {
			// number of points
			if( ag + 1 >= argc ||
			    1 != sscanf( argv[ag+1], "%ld", & req_n_wedges) ||
			    req_n_wedges < 1)
				usage();
			ag += 2; continue;
		}
		else if( strcmp( argv[ag], "-m") == 0) {
			// number of material layers
			if( ag + 1 >= argc ||
			    1 != sscanf( argv[ag+1], "%ld", & n_layers) ||
			    n_layers < 1)
				usage();
			ag += 2; continue;
		}
		else if( strcmp( argv[ag], "-d") == 0) {
			// width of the patch
			if( ag + 2 >= argc ||
			    1 != sscanf( argv[ag+1], "%lf", & bb_width) ||
			    1 != sscanf( argv[ag+2], "%lf", & bb_height) ||
			    bb_width <= 0 ||
			    bb_height <= 0)
				usage();
			ag += 3; continue;
		}
		else if( strcmp( argv[ag], "-t") == 0) {
			// subdivision type
			if( ag + 1 >= argc) usage();
			ag ++;
			if( strcmp( "random", argv[ag]) == 0) {
				subdivision_type = Rand; relax = 0;
			} else if( strcmp( "repelled", argv[ag]) == 0) {
				subdivision_type = Rand; relax = 1;
			} else if( strcmp( "ut2", argv[ag]) == 0) {
				subdivision_type = UnifTri2; relax = 0;
			} else if( strcmp( "ut4", argv[ag]) == 0) {
				subdivision_type = UnifTri4; relax = 0;
			} else if( strcmp( "hex", argv[ag]) == 0) {
				subdivision_type = Hexagonal; relax = 0;
			} else {
				usage();
			}
			ag ++; continue;
		}
		else if( strcmp( argv[ag], "-gv") == 0) {
			// growth vector
			if( ag + 3 >= argc ||
			    1 != sscanf( argv[ag+1], "%lf", & growth_alpha) ||
			    1 != sscanf( argv[ag+2], "%lf", & growth_r1) ||
			    1 != sscanf( argv[ag+3], "%lf", & growth_r2))
				usage();
			ag += 4; continue;
		}
		else if( strcmp( argv[ag], "-hmap") == 0) {
			// growth vector
			if( ag + 3 >= argc) usage();
			hmap_fname = argv[ag+1];
			if( 1 != sscanf( argv[ag+2], "%lf", & hmap_min) ||
			    1 != sscanf( argv[ag+3], "%lf", & hmap_max))
				usage();
			ag += 4; continue;
		}
		else if( strcmp( argv[ag], "-bmap") == 0) {
			// growth vector
			if( ag + 3 >= argc) usage();
			bmap_fname = argv[ag+1];
			if( 1 != sscanf( argv[ag+2], "%lf", & bmap_min) ||
			    1 != sscanf( argv[ag+3], "%lf", & bmap_max))
				usage();
			ag += 4; continue;
		}
		else if( strcmp( argv[ag], "-h") == 0) {
			usage();
		}


		// unknown argument
		usage();
	}

	fprintf( stderr, "Req. # = %ld\n", req_n_wedges);
	fprintf( stderr, "n. layers = %ld\n", n_layers);
	fprintf( stderr, "Subdivision type = %d\n", int( subdivision_type));
	fprintf( stderr, "gv= %f %f %f\n", growth_alpha, growth_r1, growth_r2);
	fprintf( stderr, "relax = %d\n", relax);
	fprintf( stderr, "dimensions = %f x %f\n", bb_width, bb_height);
	fprintf( stderr, "hmap = %s %f %f\n", hmap_fname, hmap_min, hmap_max);
	fprintf( stderr, "bmap = %s %f %f\n", bmap_fname, bmap_min, bmap_max);

	hmap = new MapTexture2D( hmap_fname, hmap_min, hmap_max,
				 - bb_width/2, bb_width/2,
				 - bb_height/2, bb_height/2);
	bmap = new MapTexture2D( bmap_fname, bmap_min, bmap_max,
				 - bb_width/2, bb_width/2,
				 - bb_height/2, bb_height/2);
}

static double surface_height( double x, double y, long layer)
{
	double iso_z = double(layer) / n_layers;
	double z = iso_z * hmap-> get_value( x, y) +
		(1-iso_z) * bmap-> get_value( x, y);
	return z;

	double minx = - bb_width / 2; double maxx = bb_width /2;
	double miny = - bb_height/ 2; double maxy = bb_height/2;
	double p = (x - minx) / (maxx-minx);
	double q = (y - miny) / (maxy-miny);
	if( p < 0) p = 0; if( p > 1) p = 1;
	if( q < 0) q = 0; if( q > 1) q = 1;
	double height =
		thick_bl * (1-p)*(1-q) +
		thick_br * (  p)*(1-q) +
		thick_tr * (  p)*(  q) +
		thick_tl * (1-p)*(  q);
	double ret = height * (layer / double(n_layers));
	// paraboloid adjustment
	double h1 = 1.0; double h2 = 1.0;
	ret = ret * ((x*x+y*y)*(0.5*(h2-h1))+h1);
	// randomization adjustment
	ret = ret * ( 1 + (drand48()-0.5)*thick_rnd);
	return ret;
}


static void add_node( double x, double y, long layer)
{
	nodes = (Node **) realloc( nodes, sizeof( Node *) * (n_nodes+1));
	assert( nodes != NULL);
	double z = surface_height( x, y, layer);
	nodes[ n_nodes] = new Node( n_nodes, x, y, z, layer == 0 ? 1 : 0);
	n_nodes ++;
}

static void add_wedge( long num, long n1, long n2, long n3,
		       long n4, long n5, long n6,
		       double E, double poisson, double yield)
{
	assert( n1 < n_nodes);
	assert( n2 < n_nodes);
	assert( n3 < n_nodes);
	assert( n4 < n_nodes);
	assert( n5 < n_nodes);
	assert( n6 < n_nodes);

	wedges = (Wedge **) realloc( wedges, sizeof( Wedge *) * (n_wedges+1));
	assert( wedges != NULL);
	wedges[ n_wedges] = new Wedge( num, E, poisson, yield,
				       nodes[n1], nodes[n2], nodes[n3],
				       nodes[n4], nodes[n5], nodes[n6]);
	n_wedges ++;
}

static int is_degenerate_triangle( double x1, double y1, double z1,
				   double x2, double y2, double z2,
				   double x3, double y3, double z3)
{
	// calculate the three vectors of the triangle
	double v1x = x2 - x1;
	double v1y = y2 - y1;
	double v1z = z2 - z1;
	double v2x = x3 - x1;
	double v2y = y3 - y1;
	double v2z = z3 - z1;
	double v3x = x3 - x2;
	double v3y = y3 - y2;
	double v3z = z3 - z2;
	// calculate the sides
	double a = sqrt( v1x*v1x + v1y*v1y + v1z*v1z);
	double b = sqrt( v2x*v2x + v2y*v2y + v2z*v2z);
	double c = sqrt( v3x*v3x + v3y*v3y + v3z*v3z);
	// calculate the cosine of alpha
	double cosa = (b*b + c*c - a*a) / (2*b*c);
	double cosb = (a*a + c*c - b*b) / (2*a*c);
	double cosc = (a*a + b*b - c*c) / (2*a*b);
	// calculate angles
	if( cosa < -1) cosa = -1; if( cosa > 1) cosa = 1;
	if( cosb < -1) cosb = -1; if( cosb > 1) cosb = 1;
	if( cosc < -1) cosc = -1; if( cosc > 1) cosc = 1;
	double alpha = acos( cosa) * 180 / M_PI;
	double beta = acos( cosb) * 180 / M_PI;
	double gamma = acos( cosc) * 180 / M_PI;

	double min_angle = 3;
	int is_degenerate = 
		isnan( alpha) || isnan( beta) || isnan( gamma) ||
		alpha < min_angle || beta < min_angle || gamma < min_angle;

//	fprintf( stderr, "angles: %f %f %f (%f) %s\n", 
//		 alpha, beta, gamma, alpha + beta + gamma,
//		 is_degenerate ? "degenerate" : "OK"
//		 );
	
	return is_degenerate;
}

int main( int argc, char ** argv)
{
	if( 0) {
		while( 1) {
			fprintf( stderr, "%f %f %f\n",
				 get_random_young_mod(),
				 get_random_poisson_mod(),
				 get_random_yield_stress());
		}
	}

	parse_arguments( argc, argv);

	if( subdivision_type == UnifTri2) {
		// determine the size of a square
		// ==================================================
		// each square will contain 2 triangles. So:
		// n_wedges = n_squares * 2
		// But n_squares has to be a power of 2
		long sub = 2;
		while( sub * sub * 2 <= req_n_wedges) sub += 1;
		sub -= 1;
		
		double square_size = 2.0 / sub;

		for( long i = 0 ; i <= sub ; i ++)
			for( long j = 0 ; j <= sub ; j ++) {
				add_node( (i*square_size-1) / 2,
					  (j*square_size-1) / 2,
					  0);
				add_node( (i*square_size-1) / 2,
					  (j*square_size-1) / 2,
					  1);
		}

		for( long i = 0 ; i < sub ; i ++) {
			for( long j = 0 ; j < sub ; j ++) {
				long w = (sub+1)*2;
				long p = i*w + j*2;
				add_wedge( 2*(i*sub+j),
					   p, p+w, p+2,
					   p+1, p+w+1, p+2+1,
					   get_random_young_mod(),
					   get_random_poisson_mod(),
					   get_random_yield_stress()
					);
				add_wedge( 2*(i*sub+j)+1,
					   p+2,p+w,p+w+2,
					   p+2+1,p+w+1,p+w+2+1,
					   get_random_young_mod(),
					   get_random_poisson_mod(),
					   get_random_yield_stress()
					);
			}
		}
	} else if( subdivision_type == UnifTri4) {
		long sub = 2;
		while( sub * sub * 4 <= req_n_wedges) sub ++;
		sub --;
		// printf( "# sub = %ld\n", sub);
		double size = 2.0 / sub;

		for( long row = 0 ; row <= sub ; row ++) {
			for( long col = 0 ; col <= sub ; col ++) {
				add_node( (col*size-1.0)/2,
					  (row*size-1.0)/2,
					  0);
				add_node( (col*size-1.0)/2,
					  (row*size-1.0)/2,
					  1);
			}
			if( row < sub) {
				for( long col = 0 ; col < sub ; col ++) {
					add_node( (col*size-1.0+size/2)/2,
						  (row*size-1.0+size/2)/2,
						  0);
					add_node( (col*size-1.0+size/2)/2,
						  (row*size-1.0+size/2)/2,
						  1);
				}
			}
		}

		for( long r = 0 ; r < sub ; r ++) {
			for( long c = 0 ; c < sub ; c ++) {
				long pr = 2*(sub + 1 + sub);
				long ref = r * pr + c*2;
				add_wedge( 4*sub*r+4*c+0,
					   ref, ref+2, ref+sub*2+2,
					   ref+1, ref+2+1, ref+sub*2+2+1,
					   get_random_young_mod(),
					   get_random_poisson_mod(),
					   get_random_yield_stress());
				add_wedge( 4*sub*r+4*c+1,
					   ref+2,ref+2+pr,ref+sub*2+2,
					   ref+2+1,ref+2+pr+1,ref+sub*2+2+1,
					   get_random_young_mod(),
					   get_random_poisson_mod(),
					   get_random_yield_stress());
				add_wedge( 4*sub*r+4*c+2,
					   ref+2+pr,ref+pr,ref+sub*2+2,
					   ref+2+pr+1,ref+pr+1,ref+sub*2+2+1,
					   get_random_young_mod(),
					   get_random_poisson_mod(),
					   get_random_yield_stress());
				add_wedge( 4*sub*r+4*c+3,
					   ref+pr,ref,ref+sub*2+2,
					   ref+pr+1,ref+1,ref+sub*2+2+1,
					   get_random_young_mod(),
					   get_random_poisson_mod(),
					   get_random_yield_stress());
			}
		}

	} else if( subdivision_type == Hexagonal) {
		// determine the number of rows and columns
		long n_rows = 2;
		long n_cols = 2;
		while( 1) {
			long n_t = (n_rows-1) * (n_cols - 1) * 2;
			if( n_t > req_n_wedges) break;
			n_rows ++;
			n_cols ++;
		}
		n_rows --; n_cols --;
		fprintf( stderr, "HEX: %ldx%ld\n", n_rows, n_cols);

		// calculate the length of a triangle side
		double s = 2.0 / (n_cols - 0.5);
		// calculate the height
		double h = sqrt( 3*s*s/4);
		// calculate total height
		double th = (n_rows-1) * h;
		// calculate initial y offset
		double y_off = - th / 2;
		
		// generate the nodes
		for( long row = 0 ; row < n_rows ; row ++) {
			double x_off = -1;
			if( row % 2 == 1) x_off += s/2;
			for( long col = 0 ; col < n_cols ; col ++) {
				add_node( x_off + col * s,
					  y_off + row * h,
					  0);
				add_node( x_off + col * s,
					  y_off + row * h,
					  1);
			}
		}

		// generate wedges
		long num = 0;
		for( long row = 1 ; row < n_rows ; row ++) {
			for( long col = 1 ; col < n_cols ; col ++) {
				long ind = row * n_cols + col;
				if( row % 2 == 1) {
					add_wedge(
						num ++,
						0 + 2*(ind - n_rows - 1),
						0 + 2*(ind - n_rows),
						0 + 2*(ind -1),
						1 + 2*(ind - n_rows - 1),
						1 + 2*(ind - n_rows),
						1 + 2*(ind -1),
						get_random_young_mod(),
						get_random_poisson_mod(),
						get_random_yield_stress()
						);
					add_wedge(
						num ++,
						0 + 2*(ind - n_rows),
						0 + 2*(ind),
						0 + 2*(ind -1),
						1 + 2*(ind - n_rows),
						1 + 2*(ind),
						1 + 2*(ind -1),
						get_random_young_mod(),
						get_random_poisson_mod(),
						get_random_yield_stress()
						);
				} else {
					add_wedge(
						num ++,
						0 + 2*(ind - n_rows - 1),
						0 + 2*(ind - n_rows),
						0 + 2*(ind),
						1 + 2*(ind - n_rows - 1),
						1 + 2*(ind - n_rows),
						1 + 2*(ind),
						get_random_young_mod(),
						get_random_poisson_mod(),
						get_random_yield_stress()
						);
					add_wedge(
						num ++,
						0 + 2*(ind - n_rows - 1),
						0 + 2*(ind),
						0 + 2*(ind -1),
						1 + 2*(ind - n_rows - 1),
						1 + 2*(ind),
						1 + 2*(ind -1),
						get_random_young_mod(),
						get_random_poisson_mod(),
						get_random_yield_stress()
						);
				}
			}
		}
	} else if( subdivision_type == Rand) {
		{
			FILE * fp = fopen( "/tmp/rndpts2", "w");
			assert( fp != NULL);
			long n_pts = req_n_wedges + 6;
			 // bounding box
			fprintf( fp, "%.15f %.15f %.15f %.15f\n",
				 - bb_width / 2, bb_width / 2,
				 - bb_height / 2, bb_height / 2);
			// number of points
			fprintf( fp, "%ld\n", n_pts);
			// the points themselves
			for( long i = 0 ; i < n_pts ; i ++)
				fprintf( fp, "%.15f %.15f\n",
					 drand48() * bb_width - bb_width/2,
					 drand48() * bb_height - bb_height/2);
			fclose( fp);
		}
//		execute( "rbox %ld > /tmp/rndpts", req_n_wedges+6);
//		execute( "tail +2 < /tmp/rndpts > /tmp/rndpts2");
		if( relax)
			execute( "./relax -b < /tmp/rndpts2 > /tmp/rndpts3");
		else
			execute( "tail +2 < /tmp/rndpts2 > /tmp/rndpts3");
		execute( "(echo 2;cat /tmp/rndpts3) > /tmp/rndpts4");
		execute( "qhull d o < /tmp/rndpts4 > /tmp/rndpts5");

		// open the last file and get all the info about nodes and
		// wedges
		FILE * fp = fopen( "/tmp/rndpts5", "r");
		if( fp == NULL) die( "could not fopen(/tmp/rndpts5):%s",
				     strerror( errno));
		// skip the dimension
		long dim; assert( 1 == fscanf( fp, "%ld", & dim) && dim == 3);
		// get number of points & triangles & skip n_facets
		long n_points, n_triangles, n_facets;
		assert( 3 == fscanf( fp, "%ld %ld %ld", & n_points,
				     & n_triangles, & n_facets));
		assert( n_points == req_n_wedges + 6);
		// read in all points
		for( long i = 0 ; i < n_points ; i ++) {
			double x, y, z;
			assert( 3 == fscanf( fp, "%lf %lf %lf", &x,&y,&z));
			// add a node for every layer
			for( long j = 0 ; j <= n_layers ; j ++)
				add_node( x, y, j);
		}
		// read in all triangles
		for( long i = 0 ; i < n_triangles ; i ++) {
			long n, n1, n2, n3;
			assert( 4 == fscanf( fp, "%ld %ld %ld %ld",
					     &n, &n1, &n2, &n3));
			assert( n == 3);
			assert( n1 >= 0 && n1 < n_points);
			assert( n2 >= 0 && n2 < n_points);
			assert( n3 >= 0 && n3 < n_points);
			// skip degenerate triangles
			if( ! relax || ! is_degenerate_triangle( 
				    nodes[n1*(n_layers+1)]->x,
				    nodes[n1*(n_layers+1)]->y,
				    0,
				    nodes[n2*(n_layers+1)]->x,
				    nodes[n2*(n_layers+1)]->y,
				    0,
				    nodes[n3*(n_layers+1)]->x,
				    nodes[n3*(n_layers+1)]->y,
				    0))
			{
				// add a wedge for every layer
				for( long j = 0 ; j < n_layers ; j ++)
					add_wedge(
						i,
						n1*(n_layers+1)+j+1,
						n2*(n_layers+1)+j+1,
						n3*(n_layers+1)+j+1,
						n1*(n_layers+1)+j,
						n2*(n_layers+1)+j,
						n3*(n_layers+1)+j,
						get_random_young_mod(),
						get_random_poisson_mod(),
						get_random_yield_stress());
			}
		}
		fclose( fp);
	}

	// generate growth vectors
	for( long i = 0 ; i < n_nodes ; i ++) {
		double x = nodes[i]-> x;
		double y = nodes[i]-> y;
		double z = nodes[i]-> z;
		double r1 = growth_r1;
		double r2 = growth_r2;
		double a = growth_alpha * (M_PI/180);
		double gx = 0.5*x*(r1+r2+(r2-r1)*cos(2*a))+y*(r1-r2)*sin(a)*cos(a);
		double gy=  0.5*y*(r1+r2+(r1-r2)*cos(2*a))+x*(r1-r2)*sin(a)*cos(a);
		double gz = 0;
		
		nodes[i]-> gx = gx;
		nodes[i]-> gy = gy;
		nodes[i]-> gz = gz;

		// tapered growth (none on top)
//		 nodes[i]-> gx = nodes[i]-> x * (bb_height/2-nodes[i]-> y);
		// tapered growth (half on top)
//		nodes[i]-> gx = nodes[i]-> x *
//			(1.0+(bb_height/2-nodes[i]-> y))/2.0;
//		 nodes[i]-> gy = gy*2;
//		 nodes[i]-> gz = 0 ;
	}

	// print the model
	printf( "# Requested number of wedges: %ld\n", req_n_wedges);
	printf( "# Number of wedges:  %ld\n", n_wedges);
	printf( "\n");
	printf( "\n");
	printf( "sim_time_total = 1.0\n");
	printf( "min_time_step = 0.001\n");
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
	printf( "min_refine_size = 10.01\n");
	printf( "crack_tip_element_size = 0.1\n");
	printf( "crack_tip_sub_element_size = 0.001\n");
	printf( "min_element_size = -0.01\n");
	printf( "max_element_size = 0.1\n");
	printf( "precision = 1e-6\n");
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
	printf( "\tyield_stress constant 70.0\n");
	printf( "\tfracture_toughness constant 0.01\n");
	printf( "END material_properties\n");

	printf( "\n");
	printf( "%ld # number of nodes\n", n_nodes);
	for( long i = 0 ; i < n_nodes ; i ++)
		printf( "%.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %.3f %d\n",
			nodes[i]-> x, nodes[i]-> y, nodes[i]-> z,
			nodes[i]-> dx, nodes[i]-> dy, nodes[i]-> dz,
			nodes[i]-> gx, nodes[i]-> gy, nodes[i]-> gz,
			nodes[i]-> type);

	printf( "\n");
	printf( "%ld # number of elements\n", n_wedges);
	printf( "\n");
	for( long i = 0 ; i < n_wedges ; i ++) {
		printf( "fivewall %ld %ld %ld %ld %ld %ld %ld\n",
			wedges[i]-> num,
			wedges[i]-> n1-> num,
			wedges[i]-> n2-> num,
			wedges[i]-> n3-> num,
			wedges[i]-> n4-> num,
			wedges[i]-> n5-> num,
			wedges[i]-> n6-> num);
	}

	printf( "\n0 # number of fracture tips\n");

	return 0;
}
