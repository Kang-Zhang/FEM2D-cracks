#include <cstdio>
#include <cstdlib>
#include <stdarg.h>
#include <cmath>
#include <sys/types.h>
#include <unistd.h>
#include <GL/glut.h>
#include <GL/gl.h>
#include <cassert>
#include <iostream>
#include <vector>
#include <fstream>
#include <cerrno>
#include <cstring>
#include <sstream>
#include "BBox.hh"

#include "Vector.hpp"


using VectLib::V2D;
using VectLib::V3D;


#ifdef DONT_COMPILE


int main ()
{
    int n = 100;
    double minx = -0.5;
    double maxx = 0.5;
    double miny = -0.5;
    double maxy = 0.5;
    std::cout << n + 4 << " 2 1 0\n";
    std::cout << "1 " << minx << " " << miny << " 1\n";
    std::cout << "2 " << maxx << " " << miny << " 1\n";
    std::cout << "3 " << maxx << " " << maxy << " 1\n";
    std::cout << "4 " << minx << " " << maxy << " 1\n";
    for (int i = 0 ; i < n ; i ++)
	std::cout << 5 + i << " "
		  << drand48() * (maxx - minx) + minx << " "
		  << drand48() * (maxy - miny) + miny << " 1\n";
    std::cout << "4 0\n";
    std::cout << "1 1 2\n"
	      << "2 2 3\n"
	      << "3 3 4\n"
	      << "4 4 1\n";
    std::cout << "0\n";
    return 0;
}

#endif

// --------------------------------------------------------------------------
//
//                         global variables
//
// --------------------------------------------------------------------------
//
// number of desired points
static size_t n_points_wanted;
// number of iterations for relaxation
static long n_iterations;
// datastructure for a point
struct Point {
    Point (const V2D & pos, double desired_d, double radius_pert)
	: cur (pos)
	, force (0,0)
	, sum (0)
	, die_count (0)
    {
	double r = desired_d / 2;
	rad = r * (1 + (drand48() * 2 - 1) * radius_pert);
    }
    Point ()
	: cur (0,0)
	, force (0,0)
	, sum (0)
	, rad (0) 
	, die_count (0)
	{;}

    VectLib::V2D cur; // position of the point
    VectLib::V2D force; // force acting on the point
    double sum; // sum of all force magnitudes (used for splitting points)
    double rad; // radius of point's force field
    long die_count;
};
// the array of points
std::vector <Point> points;
//double * radii = NULL;
// bounding box
BBox bbox;
// desired distance between points
double desired_d;
// how fast the points move (should be <= 0.5)
double weight;
double weight_mod = 0.9;
// what is the total maximum force allowed on a point
double max_force;
// if the sum of force magnitudes on a particle is below the following threshold, the particle
// will split
double force_split_threshold = 1.1;
// threshold above which a particle will die
double force_die_threshold = 30.0;
// how many times in a row a particle has to be above die threshold before it actually dies
double max_die_count = 5;
// command line options
static bool display_graphics = false;
static bool show_grid = true;
int verbose = 0;
// jitter angle (force is rotated by +/- this angle * drand48()
//double jitter_angle = 1.0 * (M_PI/180.0);
// jitter radius (points are jittered within this radius * desired_d)
double jitter_radius = 0.0001;
// force radius perurbation (percentage of desired_d)
double radius_pert = 0.25;

// animations
long step_count = 0;

static void die(const char * str, ...)
{
    va_list ap;
    va_start(ap, str);
    vfprintf(stderr, str, ap);
    fprintf(stderr, "\n");
    va_end(ap);
    exit(-1);
}

std::istream & operator >> (std::istream & is, V2D & p)
{
    if (! (is >> p.x()))
	die ("Expected x coordinate.");
    if (! (is >> p.y()))
	die ("Expected y coordinate.");
    return is;
}

std::istream & operator >> (std::istream & is, BBox & b)
{
    size_t n;
    if (! (is >> n))
	die ("Number of points on a polygon expected.");
    for (size_t i = 0 ; i < n ; i ++) {
	V2D p;
	if (! (is >> p))
	    die ("Expected x,y coordinates for a point %ld on the polygon.", i+1);
	bbox.add (p);
    }
    std::cerr << "Read in " << n << "-polygon\n";
    std::cerr << "Bounding box of the polygon: (" << bbox.minx << "," << bbox.miny << ") x ("
	      << bbox.maxx << "," << bbox.maxy << ")\n";
    return is;
}

/*
class Square {
public:
    long n;
    long * data;
    Square() {
	n = 0;
	data = NULL;
    }
    void add(long k) {
	if(data == NULL) data = (long *) malloc(sizeof(long));
	else data = (long *) realloc(data, sizeof(long) * (n+1));
	data[n] = k;
	n ++;
    }
    void reset(void) {
	if(data != NULL) {
	    free(data);
	    data = NULL;
	    n = 0;
	}
    }
};
*/
class Square {
private:
    std::vector <long> _data;
public:
    Square() {
    }
    void add (long k) {
	_data.push_back (k);
    }
    void reset(void) {
	_data.clear();
    }
    size_t n () const {
	return _data.size();
    }
    size_t data (size_t ind) const {
	return _data[ind];
    }
};

class Grid {
public:
    double minx, maxx, miny, maxy;
    long n_horiz, n_vert;
    Square ** squares;
    
    Grid(double pminx, double pmaxx, double pminy, double pmaxy,
	  long phoriz, long pvert)
	{
	    // copy parameters
	    minx = pminx; maxx = pmaxx; miny = pminy; maxy = pmaxy;
	    n_horiz = phoriz; n_vert = pvert;
	    
	    // allocate squares
	    squares = (Square **) malloc(sizeof(Square *) * n_vert);
	    for(long i = 0 ; i < n_vert ; i ++)
		squares[i] = new Square[ n_horiz];
	}
    
    void add(long k) {
	double x = points[k].cur.x();
	double y = points[k].cur.y();
	double square_width = (maxx-minx) / n_horiz;
	double square_height = (maxy-miny) / n_vert;
	
	long gx = long((x-minx) / square_width);
	long gy = long((y-miny) / square_height);
	
	// make sure the indices are correct
	if(gx < 0) gx = 0;
	if(gx > n_horiz-1) gx = n_horiz - 1;
	if(gy < 0) gy = 0;
	if(gy > n_vert-1) gy = n_vert -1;
	
	squares[gy][gx].add(k);
    }
    
    void reset(void) {
	for(long gx = 0 ; gx < n_horiz ; gx ++)
	    for(long gy = 0 ; gy < n_vert ; gy ++)
		squares[gy][gx].reset();
    }
    
    void report_stats(void) {
	size_t max = squares[0][0].n();
	size_t min = squares[0][0].n();
	for(long gx = 0 ; gx < n_horiz ; gx ++)
	    for(long gy = 0 ; gy < n_vert ; gy ++) {
		if(squares[gy][gx].n() > max)
		    max = squares[gy][gx].n();
		if(squares[gy][gx].n() < min)
		    min = squares[gy][gx].n();
	    }
	fprintf(stderr, "grid: min / max = %ld / %ld\n", min, max);
    }
};

static Grid * grid = NULL;

static void calculate_force(size_t i1, size_t i2, double & d)
// calculates the force between points i1 and i2
{
    // distance between identical points is invalid
    if(i1 == i2) {
	d = 0;
	return;
    }
    
    // calculate the force between points i1 and i2 and add it to
    // points[i1]
    Point & p = points[i1];
    Point & q = points[i2];
    double desired_d = p.rad + q.rad;
    
    // calculate vector (dx,dy) frpm p1 to p2
    double dx = q.cur.x() - p.cur.x();
    double dy = q.cur.y() - p.cur.y();
    
    // calculate distance PQ
    d = sqrt (dx * dx + dy * dy);
    // if p and q are more than desired_d apart,
    // they don't have any effect on each other
    // and therefore can be ignored
    if(d > desired_d) return;
    // calculate the force between the two points
    // (in percentage, 1.0 =  max. force, occuring
    // when the points coincide, 0.0 = min force, when
    // the points are desired_d apart.
/*
    double force = 1 - d/desired_d;
    force *= force;
*/
    double force = - d * d / (desired_d * desired_d) + 1;
/*
    double dr = d / desired_d;
    double force = (cos (dr * M_PI) + 1) / 2;
*/
    // compute a unit vector from dx,dy
//    double ux, uy;
    V2D u;
    if(d / desired_d < 0.001) {
	// points p and q almost coincide, we
	// will therefore create a random vector
	double alpha = drand48() * M_PI * 2;
	u = V2D (sin(alpha) + q.cur.x(), cos(alpha) + q.cur.y());
    } else {
	u = V2D (dx/d, dy/d);
    }
//     // jitter the unit vector
//     double alpha = (drand48() * 2 - 1) * jitter_angle;
//     V2D v = V2D (u.x() * cos(alpha) - u.y() * sin(alpha)
// 		 ,u.x() * sin(alpha) + u.y() * cos(alpha));
/*
    double vx = ux * cos(alpha) - uy * sin(alpha);
    double vy = ux * sin(alpha) + uy * cos(alpha);
*/
    
    // add force * unit_vector to tx,ty of p
//    p.force -= force * v;
    p.force -= force * u;
    p.sum += force;
/*
    p.force.x() -= force * vx;
    p.force.y() -= force * vy;
*/
}

static void calculate_force_total(long k, long gx, long gy, double & min_d)
// calculates the total force acting on point[k] from all neighbors
// - neighbors are points in the current grid[gy][gx] and all 8 adjacent
//   grid squares
{
    points[k].sum = 0;
    min_d = -1;
    // calculate the force from all surrounding points
    for(long dx = -1 ; dx <= 1 ; dx ++) {
	long hx = gx + dx;
	if (hx < 0 || hx >= grid-> n_horiz) continue;
	for(long dy = -1 ; dy <= 1 ; dy ++) {
	    long hy = gy + dy;
	    if (hy < 0 || hy >= grid-> n_vert) continue;
	    Square * sq = & grid-> squares[hy][hx];
	    for(size_t m = 0 ; m < sq-> n() ; m++) {
		double d;
		calculate_force(k, sq-> data(m), d);
		if(d >= 0 && (d < min_d || min_d < 0))
		    min_d = d;
	    }
	}
    }
}

static double relax_points(void)
{
    // populate the grid
    grid-> reset();
    for(size_t i = 0 ; i < points.size() ; i ++)
	grid-> add(i);
    
    // zero out all forces on all points
    for(size_t i = 0 ; i < points.size() ; i ++) {
	points[i].force = V2D (0,0);
	points[i].sum = 0;
    }
    
    // for every point find all neighbors (from its grid square and
    // its surrounding points), calculate the force these generate
    // and add it to the total force
    double min_d = -1;
    for(long gx = 0 ; gx < grid-> n_horiz ; gx ++) {
	for(long gy = 0 ; gy < grid-> n_vert ; gy ++) {
	    for(size_t k = 0 ; k < grid->squares[gy][gx].n() ; k++) {
		double d;
		calculate_force_total(
		    grid-> squares[gy][gx].data(k),
		    gx, gy, d);
		if((d < min_d || min_d < 0) && d >= 0)
		    min_d = d;
	    }
	}
    }

    // clamp the forces to max_force
    for (size_t i = 0 ; i < points.size () ; i ++) {
	double nsq = V2D::normSq (points[i].force);
	double len = sqrt (nsq);
	if (len > max_force)
	    points[i].force = max_force * V2D::normalize (points[i].force, 1e-6);
    }
    
    // adjust the coordinates by adding the force to x & y, but make
    // sure that the resulting points are in the bounding box
    // the bounding box
    for(size_t i = 0 ; i < points.size() ; i ++)
    {
	V2D x;
	x = 2.0 * x;
	points[i].cur += weight * points[i].force;
//	points[i].cur = points[i].cur + points[i].cur;
	points[i].cur = bbox.snap_inside (points[i].cur);
    }

    // pick a particle with the highest sum and locate it somewhere near a particle with the
    // lowest sum
    if (0) {
	size_t ind_min = 0;
	size_t ind_max = 0;
	for (size_t i = 0 ; i < points.size() ; i ++) {
	    if (points[ind_min].sum > points[i].sum) ind_min = i;
	    if (points[ind_max].sum < points[i].sum) ind_max = i;
	}
	double a = drand48() * 2 * M_PI;
	double r = drand48() * desired_d / 4;
	points[ind_max].cur = bbox.snap_inside
	    (points[ind_min].cur + V2D (sin(a) * r, cos(a) * r));
    }
    
    // split all particles whose sum of force magnitues is below a threshold
    {
	std::vector <Point> newPoints;
	for (size_t i = 0 ; i < points.size() ; i ++) {
//	    if (newPoints.size() + points.size() >= n_points_wanted) break;
	    if (points[i].sum > force_split_threshold) continue;
	    Point p (points[i].cur, desired_d, radius_pert);
	    V2D fwd = drand48() * 5 * desired_d * V2D::normalize (points[i].force, 1e-6); 
	    p.cur = bbox.snap_inside (p.cur + fwd);
	    newPoints.push_back (p);
	}
	// append all new points to the end of points
	for (size_t i = 0 ; i < newPoints.size() ; i ++) {
	    points.push_back (newPoints[i]);
	}
    }

    // remove all particles whose sum of force magnitudes is above a threshold
    {
	std::vector <Point> newPoints;
	for (size_t i = 0 ; i < points.size() ; i ++) {
	    if (points[i].sum > force_die_threshold)
		points[i].die_count ++;
	    else
		points[i].die_count = 0;
	    if (points[i].die_count < max_die_count)
		newPoints.push_back (points[i]);
	}
	// if the points to be kept is 0, then take at least one point
	if (newPoints.empty()) newPoints.push_back (points[0]);
	if (newPoints.size() < points.size())
	    std::cerr << "\tdying " << points.size() - newPoints.size() << " particles.\n";
	points = newPoints;
    }

    // report some statistics
    {
	double maxf = points[0].sum;
	double minf = points[0].sum;
	for (size_t i = 0 ; i < points.size() ; i ++) {
	    if (maxf < points[i].sum) maxf = points[i].sum;
	    if (minf > points[i].sum) minf = points[i].sum;
	}
	printf("Step %06ld: n:%ld min.dist = %9.6f force [%10.8f,%10.8f]\n"
	       , step_count, points.size(), min_d, minf, maxf);
	step_count ++;
    }

    // jitter the particles
    if (1) {
//	double p = 0.05;
	for (size_t i = 0 ; i < points.size() ; i ++) {
//	    if (drand48() > p) continue;
	    double a = drand48() * 2 * M_PI;
	    double r = drand48() * jitter_radius;
	    V2D d (sin(a) * r, cos(a) * r);
	    points[i].cur = bbox.snap_inside (points[i].cur + d);
	}
    }

    return min_d;
}

static bool output_poly_file (const std::string & fname)
{
    // determine which points are not on the boundary
    std::vector <V2D> pts;
    for (size_t i = 0 ; i < points.size() ; i ++)
	if (bbox.dist_to_boundary (points[i].cur) > desired_d * 0.1)
	    pts.push_back (points[i].cur);
    std::cerr << "Discarding " << points.size() - pts.size() << " points on the boundary.\n";
    std::ofstream out (fname.c_str());
    out << "# generated by relax.cc\n";
    out << bbox.pts().size() + pts.size() << " 2 0 0 1\n";
    out << "# the bounding polygon:\n";
    for (size_t i = 0 ; i < bbox.pts().size() ; i ++)
	out << i << " "
	    << bbox.pts()[i].x() << " "
	    << bbox.pts()[i].y() << " 1\n";
    out << "\n# the points inside the polygon\n";
    for(size_t i = 0 ; i < pts.size() ; i ++)
	out << bbox.pts().size() + i << " "
	    << pts[i].x() << " " << pts[i].y() << " 0\n";
    out << "\n# segments:\n";
    out << bbox.pts().size() << " 1\n";
    for (size_t i = 0 ; i < bbox.pts().size() ; i ++)
	out << i << " "
	    << i << " " << (i + 1) % bbox.pts().size() << " 1\n";
    out << "\n# no holes\n";
    out << "0\n";
    out << "\n# no regional attributes\n";
    out << "0\n";
    return true;
}

static bool triangulate ()
{
    // generate the prefix for temporary files
    std::string prefix;
    int fd;
    {
	char fname [4096];
	sprintf (fname, "/tmp/relax-data-XXXXXX");
	fd = mkstemp (fname);
	if (fd == -1) {
	    std::cerr << "Cannot create a temporary files with prefix " << fname << "\n"
		      << "\t- system error: " << strerror (errno) << "\n";
	    return false;
	}
	prefix = fname;
    }
    // output a poly file
    output_poly_file (prefix + ".poly");
    // call triangle to generate the output files
    double area = bbox.get_area ();
    std::ostringstream os;
    os << "triangle -v -a" << area / n_points_wanted << " -J -B " << prefix << ".poly";
    std::string cmd = os.str();
    std::cerr << "Command to be executed is: " << cmd << "\n";
    system (cmd.c_str());
    // call showme
    cmd = std::string ("showme") + " " + prefix + ".poly";
    // now show the results using showme
    int res = fork ();
    if (res == 0) {
	std::cerr << "Command to be executed is: " << cmd << "\n";
	system (cmd.c_str());
    }
    if (res == 0 || res == -1) {
	// delete temporary files
	close (fd);
	unlink (prefix.c_str());
	unlink ((prefix + ".poly").c_str());
	unlink ((prefix + ".1.v.node").c_str());
	unlink ((prefix + ".1.v.edge").c_str());
	unlink ((prefix + ".1.poly").c_str());
	unlink ((prefix + ".1.node").c_str());
	unlink ((prefix + ".1.ele").c_str());
    }
    if (res == 0) {
	// exit child
	exit (0);
    }
    return true;
}

static void output_points(void)
{
    output_poly_file ("out.poly");
}

static void usage(void)
{
    printf(
	"Usage: relax [-h] [-g] [-i n_iterations] \n"
	"             [-v] [-b] [-t] [-j jitter] [-r radius perturbation]\n"
	"             [-s thr] [-d thr] [-w weight]\n"
	"\n"
	"       -h displays this help\n"
	"       -g produces graphical result\n"
	"       -i sets the number of iterations\n"
	"       -v verbose\n"
	"       -j radius (relative to desired_d)\n"
	"       -r repulsion radius perturbation (0..1, default 0.2)\n"
	"       -s threshold for splitting\n"
	"       -d threshold for dying\n"
	"       -w sensitivity to the force (default 1)\n"
	"\n"
	);
}


static void parse_arguments(int & argc, char ** argv)
{
    // set defaults
    display_graphics = false;
    n_iterations = 200;
    
    while(1)
    {
	int c = getopt(argc, argv, "ghf:i:btj:r:s:d:w:");
	if(c == -1) break;
	if(c == 'g') { display_graphics = true; continue;}
	if(c == 'h') { usage(); exit(0); }
	if(c == 'v') { verbose = 1; continue; }
	if(c == 'i') {
	    if(1 != sscanf(optarg, "%ld", & n_iterations)) {
		usage();
		die("Invalid number of iterations.");
	    }
	    if(n_iterations < 1) {
		usage();
		die("Number of iterations < 1 !?!");
	    }
	    continue;
	}
	if(c == 'j') {
	    if(1 != sscanf(optarg, "%lf", & jitter_radius)) {
		usage();
		die("Invalid jitter radius.");
	    }
	    continue;
	}
	if(c == 'w') {
	    if(1 != sscanf(optarg, "%lf", & weight_mod)) {
		usage();
		die("Invalid weight modifier.");
	    }
	    continue;
	}
	if(c == 'r') {
	    if(1 != sscanf(optarg, "%lf", & radius_pert)) {
		usage();
		die("Invalid radius perturbation.");
	    }
	    continue;
	}
	if(c == 's') {
	    if(1 != sscanf(optarg, "%lf", & force_split_threshold)) {
		usage();
		die("Invalid split threshold.");
	    }
	    continue;
	}
	if(c == 'd') {
	    if(1 != sscanf(optarg, "%lf", & force_die_threshold)) {
		usage();
		die("Invalid die threshold.");
	    }
	    continue;
	}
	usage(); die("Bad syntax.");
    }
    
    if(display_graphics)
	fprintf(stderr, "Displaying graphics.\n");
    fprintf(stderr, "Number of iterations: %ld\n", n_iterations);
}

static void redraw(void)
{
    // clear the screen
    glDrawBuffer(GL_BACK);
    glClear(GL_COLOR_BUFFER_BIT);
    
    if(show_grid) {
	double w = (bbox.maxx-bbox.minx) / grid-> n_horiz;
	double h = (bbox.maxy-bbox.miny) / grid-> n_vert;
	for(long gx = 0 ; gx < grid-> n_horiz ; gx ++) {
	    for(long gy = 0 ; gy < grid-> n_vert ; gy ++) {
		double g1 = 0.15;
		double g2 = 0.30;
		if((gx + gy) % 2 == 0)
		    glColor3f(g1, g1, g1);
		else
		    glColor3f(g2, g2, g2);
		
		double x = bbox.minx + gx * w;
		double y = bbox.miny + gy * h;
		glBegin(GL_POLYGON);
		glVertex2d(x, y);
		glVertex2d(x+w, y);
		glVertex2d(x+w, y+h);
		glVertex2d(x, y+h);
		glEnd();
	    }
	}
    }
    
    glPointSize(6);
    glBegin(GL_POINTS);
//    std::cerr << "desired_d = " << desired_d << "\n";
    for(size_t i = 0 ; i < points.size() ; i ++) {
	double r = points[i].sum / force_die_threshold;
	glColor3f(1, 1-r, 1-r);
	if (1) {
	    double d = bbox.dist_to_boundary (points[i].cur);
//	    std::cerr << i << " - " << d << "\n";
	    if (d < desired_d * 0.1)
		glColor3f (0, 1, 0);
	}
	glVertex2d(points[i].cur.x(),
		   points[i].cur.y());
    }
    glEnd();
    
    // swap buffers
    glutSwapBuffers();
    
    // glutPostRedisplay();
}

static void reshape_func(int width, int height)
{	
    // set the viewport to use the entire drawing area
    glViewport(0, 0, width, height);
    
    // set up the view matrix to be the identity matrix
    glMatrixMode(GL_MODELVIEW );
    glLoadIdentity();
    
    // set the projection matrix
    glMatrixMode(GL_PROJECTION );
    glLoadIdentity();
    glOrtho(bbox.minx-0.1*(bbox.maxx-bbox.minx),
	     bbox.maxx+0.1*(bbox.maxx-bbox.minx),
	     bbox.miny-0.1*(bbox.maxy-bbox.miny),
	     bbox.maxy+0.1*(bbox.maxy-bbox.miny),
	     -1, 1);
    
    redraw();
    
    if(verbose)
	fprintf(stderr, "reshape %dx%d\n", width, height);
}

static void disp_func(void)
{	
    glDisable(GL_DITHER);
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);
    
    // antialised points
    glShadeModel(GL_SMOOTH);
    glPointSize(2.0);
    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    // pick the color for background
    glClearColor(0.1, 0, 0, 0);
    
    redraw();
}

static void idle_func(void)
{
    relax_points();
    
    redraw();
}

enum Menu { MENU_STEP, MENU_RUN, MENU_STOP, MENU_QUIT, MENU_GRID, MENU_SAVE, MENU_TRIANGULATE };

void menu_func(int menu)
// ----------------------------------------------------------------------
// called when some of the menus are selected
// ----------------------------------------------------------------------
{
    switch(menu) {
    case MENU_STEP:
	relax_points();
	redraw();
	break;
    case MENU_RUN:
	glutIdleFunc(idle_func);
	break;
    case MENU_STOP:
	glutIdleFunc(NULL);
	break;
    case MENU_QUIT:
	exit(0);
	break;
    case MENU_GRID:
	show_grid = ! show_grid;
	glutPostRedisplay();
	break;
    case MENU_SAVE:
	output_points ();
	break;
    case MENU_TRIANGULATE:
	triangulate ();
	break;
    }
}

static void create_ui(int & argc, char ** argv)
{
    glutInitWindowSize(500, 500);
    glutInit(& argc, argv);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutInitWindowSize(768, 768);
    std::string title;
    for (int i = 0 ; i < argc ; i ++)
	title = title + argv[i] + " ";
    int main_window = glutCreateWindow(title.c_str ());
    assert(main_window);
    glutSetWindow(main_window);
    
    glutDisplayFunc(disp_func);
    glutReshapeFunc(reshape_func);
    
    // create a popup menu
    int m = glutCreateMenu(menu_func);
    glutSetMenu(m);
    glutAddMenuEntry("Step", MENU_STEP);
    glutAddMenuEntry("Run", MENU_RUN);
    glutAddMenuEntry("Stop", MENU_STOP);
    glutAddMenuEntry("Quit", MENU_QUIT);
    glutAddMenuEntry("Grid", MENU_GRID);
    glutAddMenuEntry("Save", MENU_SAVE);
    glutAddMenuEntry("Triangulate", MENU_TRIANGULATE);
    glutAttachMenu(GLUT_RIGHT_BUTTON);
    
    glutMainLoop();
}


int main(int argc, char ** argv)
{
    if (0) {
	V2D p (3, 9);
	V2D p1 (3, 8);
	V2D p2 (4, 9);
	std::cerr << "d = " << VectLib::dist_to_line_seg (p, p1, p2) << "\n";
	return 0;
    }
    // parse command line arguments
    parse_arguments(argc, argv);
    
    // read in the bounding box
    if (! (std::cin >> bbox))
	die ("Bounding box could not be read in from the file.");
    // consistency check - to make sure all vertices of the polygon are 'inside'
    for (size_t i = 0 ; i < bbox.pts().size() ; i ++) {
	if (! bbox.is_inside (bbox.pts()[i]))
	    std::cerr << "Error! Vertex " << i << " not inside the polygon.\n";
    }
    std::cerr << "Area of the bounding box is " << bbox.get_area() << "\n";

    if (0) {
	for (size_t i = 0 ; i < bbox.pts().size(); i ++)
	    std::cerr << i << " - " << bbox.dist_to_boundary (bbox.pts()[i]) << "\n";
    }

    // read in number of points
    if (! (std::cin >> n_points_wanted))
	die ("Number of points could not be read in from the file.");

    // calculate the desired radius
    double area = bbox.get_area();
    desired_d = 2*sqrt(area/(n_points_wanted*M_PI));
    desired_d = 1.6 * desired_d;
    weight = 0.133333333 * weight_mod * desired_d;
    max_force = 2.0;
    if(verbose)
	fprintf(stderr, "# desired distance: %f\n", desired_d);

    // create a grid
    long grid_h = long((bbox.maxx-bbox.minx) / desired_d + 1);
    long grid_v = long((bbox.maxy-bbox.miny) / desired_d + 1);
    grid = new Grid(bbox.minx, bbox.maxx, bbox.miny, bbox.maxy,
		     grid_h, grid_v);
    
    // generate the points
//    for(size_t i = 0 ; i < n_points_wanted ; i ++) {
    for(size_t i = 0 ; i < 1 ; i ++) {
	V2D loc = bbox.random_point ();
	Point p (loc, desired_d, radius_pert);
	points.push_back (p);
	std::cerr << "Generated point " << i << " with rad = " << p.rad << "\n";
    }

    // if graphics has to be displayed - create the widgets
    if(display_graphics) create_ui(argc, argv);
    
    // relax points
    for(long i = 0 ; i < n_iterations ; i ++) {
	double d = relax_points();
	if(verbose)
	    fprintf(stderr, "%4ld) %f\n", i, d);
    }
    
    // output the points
    output_points();
    
    // return success
    return 0;
}
