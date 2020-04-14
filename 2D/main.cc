#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <stdarg.h>
#include <ctype.h>
#include <GL/glut.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include <malloc.h>
#include <time.h>
#include <algorithm>
#include <string>
#include "Matrix.hh"
#include "SparseMatrix.hh"
#include "conjugate_gradient.hh"
//#include "cubic.hh"
//#include "Timer.hh"
#include "Tokenizer.hh"
#include "die.hh"
#include "Distribution.hh"
#include "Map.hh"
#include "Vector3d.hh"
#include "StressTensor.hh"
#include "plane_split.hh"
#include "projected_area.hh"
#include "Stack.hh"
// #include "LSQF.hh"
#include "glDrawCube.hh"
#include "glDrawTetra.hh"
#include "glDrawWedge.hh"
#include "Materials.hh"
#include "Vector.hh"
#include "Fifo.hh"
#include "linfunc.hh"
#include "xmath.hh"

using namespace std;

// int ns_debug_flag = 0;

// global variables affecting rendering
// -------------------------------------------------------
int nodal_stresses_from_elements = 0;
int antialiasing = 1;
int draw_axes = 0;
int draw_nodal_stresses = 0;
int draw_elemental_stresses = 0;
int draw_wireframe = 1;
int draw_solid = 1;
int draw_nodes = 0;
double draw_node_cube_size = 0.01;
int draw_font_num = 0;
int draw_fog = 0;
int draw_growth_vectors = 0;
int draw_original_shapes = 0;
int draw_deformed_model = 1;
int optimize_memory_flag = 0;
int fast_ke_removal_flag = 0;
bool precise_nodal_stress_flag = true;
int use_local_relaxation_flag = 1;

// keep track of time spent on relaxing and stress calculation
long long relax_time = 0;
long long stress_time = 0;

long ke_recalc_n = 0;

// global simulation parameters
// --------------------------------------------------
// when an element is split into two triangles by the stress plane, it
// won't be allowed to do so if the two new angles would yield a new angle
// smaller than this constant (given in degrees)
double min_split_angle = 0.1;
// minimum angle between fracture plane and a boundary node
double min_ftbn_angle = 30;
//double min_ftbe_angle = 10;
double min_ftbe_angle = 5;
// when an element is split into two triangles by the stress plane, it
// won't be allowed to do so if either one of the two new edges would be
// smaller than this constant
double min_split_length = 0.01;
// choose one of the above methods
int split_limit_by_angle = 1;
// what is a degenerate triangle? we can define it through the smallest
// possible angle inside a good triangle. Degenerate triangles will
// be eliminated (hopefuly) from the mesh
//double degenerate_triangle_min_angle = 10.0;
double degenerate_triangle_min_angle = 3;
// what is an acceptable triangle? non-acceptable triangle will try to
// be prettified if it does not affect the mesh topology
double pretty_triangle_min_angle = 20.0;

void * draw_fonts[] = {
    GLUT_BITMAP_TIMES_ROMAN_10,
    GLUT_BITMAP_HELVETICA_10,
    GLUT_BITMAP_HELVETICA_12,
    GLUT_BITMAP_8_BY_13,
    GLUT_BITMAP_9_BY_15,
    GLUT_BITMAP_HELVETICA_18,
    GLUT_BITMAP_TIMES_ROMAN_24
};
long draw_n_fonts = 7;
long color_display = 0;

// autosave global variables
// ----------------------------------------------------------------------
long autosave_interval = -1;	// default - autosave off
long autosave_n_keep = 5;	// number of autosaves to keep
char * autosave_fname_mask = "/tmp/fem-autosave%ld";

// animation global variables
// ----------------------------------------------------------------------
char * animation_fname_mask = NULL;
long animation_count = 0;

int win_width, win_height;	// window dimensions

static int execute( char * str, ...)
{
    va_list ap;
    va_start( ap, str);
    char command[ 4096];
    vsprintf( command, str, ap);
    va_end( ap);

    // fprintf( stderr, "execute( '%s')\n", command);
    int res = system( command);
    // fprintf( stderr, "   - system() returned %d\n", res);
    if( res) {
	fprintf( stderr, "Could not execute command '%s'", command);
	return -1;
    } else {
	return 0;
    }
}

static void usage( const char * pname)
{
    die (
	"Usage: %s [-p]\n"
	"          [-autosave fname_mask interval history]\n"
	"          [-anim fname_mask count]\n"
	"          [-break Element | Node | Edge ] \n"
	"          [-memopt]\n"
	"          [-no_local_relax]\n"
	"          [-fast_ke_removal]\n"
	"          [-calculate_nodal_stress]"
	"          [-no-precise-nodal-stress]"
	"          filename.\n", pname);
}

static double linterpolate( double val1, double val2, double r)
{
    return val1 + r * (val2-val1);
}

int where_is_triangle_wrt_plane(
    double nx, double ny, double nz,
    double x0, double y0, double z0,
    double x1, double y1, double z1,
    double x2, double y2, double z2,
    double & s
    )
{
    if( split_limit_by_angle)
	return where_is_triangle_wrt_plane_angle(
	    nx, ny, nz,
	    x0, y0, z0,
	    x1, y1, z1,
	    x2, y2, z2,
	    min_split_angle,
	    s);
    else
	return where_is_triangle_wrt_plane_length(
	    nx, ny, nz,
	    x0, y0, z0,
	    x1, y1, z1,
	    x2, y2, z2,
	    min_split_length,
	    s);
}

static double triangle_quality( Vector3d & p1,
    Vector3d & p2,
    Vector3d & p3)
{
    Vector3d v1, v2;
    // calculate angle at point 1:
    v1 = p2-p1; v1.normalize();
    v2 = p3-p1; v2.normalize();
    double a1 = fabs( acos( v1*v2)) * 180 / M_PI;

    // calculate angle at point 2:
    v1 = p1-p2; v1.normalize();
    v2 = p3-p2; v2.normalize();
    double a2 = fabs( acos( v1*v2)) * 180 / M_PI;

    // calculate angle at point 3:
    v1 = p1-p3; v1.normalize();
    v2 = p2-p3; v2.normalize();
    double a3 = fabs( acos( v1*v2)) * 180 / M_PI;

    // get the max:
    double max = a1;
    if( a2 > max) max = a2;
    if( a3 > max) max = a3;

    // get the min
    double min = a1;
    if( a2 < min) min = a2;
    if( a3 < min) min = a3;

    // calculate difference
    double diff = max - min;

    // divide by 180
    double r = diff / 180;

    return 1 - r;
}
				

static void draw_text( double x, double y, double z,
    char * str, ...)
// ---------------------------------------------------------------------------
// draw a string at the specified xyz coordinates
// ---------------------------------------------------------------------------
{
    // prepare the string into 'buff'
    va_list ap;
    va_start( ap, str);
    char buff[ 4096];
    vsprintf( buff, str, ap);
    va_end( ap);

    // render the string
    glRasterPos3d( x, y, z);
    for( char * c = buff ; * c != '\0' ; c ++)
	glutBitmapCharacter( draw_fonts[ draw_font_num], * c);
}

// main window
int main_window = 0;
enum Menu { MENU_SIMULATION_STEP,
	    MENU_START_SIMULATION,
	    MENU_STOP_SIMULATION,
	    MENU_SAVE_MODEL,
	    MENU_ORIENT_ELEMENTS,
	    MENU_VIEW_MODEL
};

class View {

 public:

    Vector3d eye;
    Vector3d center;
    Vector3d upv;

    void set( Vector3d p_eye,
	Vector3d p_center,
	Vector3d p_upv) {
	eye = p_eye;
	center = p_center;
	upv = p_upv;
    }

    View() {
	eye.set( 1, 1, 1);
	center.set( 0, 0, 0);
	upv.set( 0, 0, 1);
    }

    void zoom( double r) {
	Vector3d v = center - eye;
	eye = eye + r * v;
	center = center + r * v;
	return ;
	if( r < -0.999) r = -0.999;
	Vector3d new_eye = center + (eye-center) * (1 + r);
	if( (new_eye-center).length() > 1e-3)
	    eye = new_eye;
    }

    void rotate_center( double alpha, double beta) {
	double radius = (center-eye).length();
	Vector3d e1 = center - eye; e1.normalize(); e1 = e1 * radius;
	Vector3d e2 = cross_product( e1, upv); e2.normalize(); e2 = e2 * radius;
	Vector3d e3 = cross_product( e1, e2); e3.normalize(); e3 = e3 * radius;
	double c1 = cos( alpha / 360) * cos( beta / 360);
	double c2 = sin( alpha / 360) * cos( beta / 360);
	double c3 = sin( beta / 360);
	center = eye + c1 * e1 + c2 * e2 + c3 * e3;
	upv = -1 * e3;
    }

    void rotate_eye( double alpha, double beta) {
	double radius = (center-eye).length();
	Vector3d e1 = eye - center; e1.normalize(); e1 = e1 * radius;
	Vector3d e2 = cross_product( e1, upv); e2.normalize(); e2 = e2 * radius;
	Vector3d e3 = cross_product( e1, e2); e3.normalize(); e3 = e3 * radius;
	double c1 = cos( alpha / 360) * cos( beta / 360);
	double c2 = sin( alpha / 360) * cos( beta / 360);
	double c3 = sin( beta / 360);
	eye = center + c1 * e1 + c2 * e2 + c3 * e3;
	upv = -1 * e3;
    }

    void translate( double x, double y) {
	double radius = (center-eye).length();
	Vector3d e1 = eye - center; e1.normalize(); e1 = e1 * radius;
	Vector3d e2 = cross_product( e1, upv); e2.normalize(); e2 = e2 * radius;
	Vector3d e3 = cross_product( e1, e2); e3.normalize(); e3 = e3 * radius;

	center = center + x * e2 + y * e3;
	eye = eye + x * e2 + y * e3;
    }

} view;

class WedgeElement;
class Node {
 private:
    void _init( void) {
	ps_val = 0;
	ps_nx = ps_ny = ps_nz = 0;
	n_refs = 0;
	refs = NULL;
	is_fixed = 0;
	is_in_pzone = 0;
	selected = 0;
	has_large_error = 0;
	young_modulus = 0;
	poisson_ratio = 0;
	yield_stress = 0;
	ftoughness = 0;
    }
 public:
    double ox, oy, oz;	// original location of the node
    double px, py, pz;	// current location of the node
    char is_fixed;          // is the node attached to the background?
    // 0 - no, it is free to move
    // 1 - yes, and it grows with background
    // 2 - no, but its displacement is fixed
    char selected;		// is the node selected
    char has_large_error;	// is there a large error at this node
    char is_in_pzone;	// is this node in a plastic zone?
    long n_refs;		// number of elements using this node
    WedgeElement ** refs;	// list of elements using this node
    Vector3d gv;		// background growth vector (only
				// used for fixed nodes, oh well)
	
    // principal stress at the node: the value and the corresponding
    // normal vector of the associated stress plane
    double ps_val, ps_nx, ps_ny, ps_nz;
    // material properties at the node, calculated as a (possibly
    // weighted) average from the elements surrounding the node
    double young_modulus;
    double poisson_ratio;
    double yield_stress;
    double ftoughness;
    StressTensor stensor;	// stress tensor at the node

    void read_data( Tokenizer & tok) {
	ox = tok.read_double( "node orig. X");
	oy = tok.read_double( "node orig. Y");
	oz = tok.read_double( "node orig. Z");
	double dx = tok.read_double( "node disp. X");
	double dy = tok.read_double( "node disp. Y");
	double dz = tok.read_double( "node disp. Z");
	px = ox + dx;
	py = oy + dy;
	pz = oz + dz;
	// read growth vector
	gv.x = tok.read_double( "normal X");
	gv.y = tok.read_double( "normal Y");
	gv.z = tok.read_double( "normal Z");

	// read fixed status / pzone status
	is_fixed = tok.read_long( "is_fixed");
	if( is_fixed >= 4) {
	    is_fixed -= 4;
	    is_in_pzone = 1;
	}
	assert( is_fixed == 0 || is_fixed == 1 || is_fixed == 2);
    }

    void add_ref( WedgeElement * e) {
	refs = (WedgeElement **) realloc(
	    refs, sizeof( WedgeElement *) * (n_refs + 1));
	refs[n_refs] = e;
	n_refs ++;
    }

    void del_ref( void * e) {
	long ind = -1;
	for( long i = 0 ; i < n_refs ; i ++)
	    if( refs[i] == e) {
		ind = i;
		break;
	    }
	assert( ind != -1);

	if( n_refs > 0)
	    refs[ind] = refs[n_refs-1];

	// 		if( ind < n_refs - 1)
	// 			memmove( & refs[ ind], & refs[ ind+1],
	// 				 sizeof( void *) * (n_refs-1-ind));
	n_refs --;
    }

    static double dist( Node * n1, Node * n2) {
	return sqrt( (n1-> ox - n2-> ox)*(n1-> ox - n2-> ox) +
	    (n1-> oy - n2-> oy)*(n1-> oy - n2-> oy) +
	    (n1-> oz - n2-> oz)*(n1-> oz - n2-> oz));
    }

    static double dist_p( Node * n1, Node * n2) {
	return sqrt( (n1-> px - n2-> px)*(n1-> px - n2-> px) +
	    (n1-> py - n2-> py)*(n1-> py - n2-> py) +
	    (n1-> pz - n2-> pz)*(n1-> pz - n2-> pz));
    }

    static double angle( Node * n1, Node * n2, Node * n3) {
	Vector3d v1( n1->ox-n2->ox, n1->oy-n2->oy, n1->oz-n2->oz);
	Vector3d v2( n3->ox-n2->ox, n3->oy-n2->oy, n3->oz-n2->oz);
	v1.normalize();
	v2.normalize();
	return fabs( acos( v1*v2)) * 180 / M_PI;
    }

    Node( ) {
	_init();
    }

    Node( Node & n) {
	_init();
	ox = n.ox; oy = n.oy; oz =n.oz;
	px = n.px; py = n.py; pz =n.pz;
	gv = n.gv;
	is_fixed = n.is_fixed;
	is_in_pzone = n.is_in_pzone;
    }
    Node( double pox, double poy, double poz,
	double ppx, double ppy, double ppz,
	long pis_fixed)
    {
	assert( ! isnan(pox));
	assert( ! isnan(poy));
	assert( ! isnan(poz));
	assert( ! isnan(ppx));
	assert( ! isnan(ppy));
	assert( ! isnan(ppz));
	_init();
	ox = pox; oy = poy; oz = poz;
	px = ppx; py = ppy; pz = ppz;
	gv.set( 0, 0, 0);
	is_fixed = pis_fixed;
    }

    Node( Node * n1, Node * n2, double r)
    // create an interpolated node of n1 and n2 at isoparametric
    // coordinate r
    {
	_init();
	// interpolate original coordinates
	ox = n1-> ox * (1-r) + n2-> ox * r;
	oy = n1-> oy * (1-r) + n2-> oy * r;
	oz = n1-> oz * (1-r) + n2-> oz * r;
	// interpolate current coordinates
	px = n1-> px * (1-r) + n2-> px * r;
	py = n1-> py * (1-r) + n2-> py * r;
	pz = n1-> pz * (1-r) + n2-> pz * r;
	// the new type is:
	//    fixed - if both n1 and n2 are fixed
	//     free - otherwise
	if( n1-> is_fixed == n2-> is_fixed)
	    is_fixed = n1-> is_fixed;
	else if( n1-> is_fixed == 0 || n2-> is_fixed == 0)
	    is_fixed = 0;
	else
	    is_fixed = 1;
	// plastic zone of the node is inherited from the neighbors
	if( n1-> is_in_pzone || n2-> is_in_pzone)
	    is_in_pzone = 1;
	else
	    is_in_pzone = 0;
	// interpolate the growth vector
	gv = n1-> gv * (1-r) + n2-> gv * r;
    }

    ~Node() {
	if( refs != NULL) free( refs);
    }

    void grow( double time,
	double growth_x, double growth_y, double growth_z)
    {
	if( is_fixed == 1) {
	    px = ox + gv.x * growth_x * time;
	    py = oy + gv.y * growth_y * time;
	    pz = oz + gv.z * growth_z * time;
	}
    }
};

static long n_nodes = 0;
static long n_nodes_alloc = 0;
static Node ** nodes = NULL;

class WedgeElement {
 private:
    Matrix * Ke;		// the element's stiffness matrix
    Vector_double * RHSe;	// the right hand side correponding to no
				// deformation
    Vector3d zyz;
 protected:

    double size;		// volume of the element (actually, this is
				// really only the are of the triangle

    double ox[6];		// the reference shape of the element
    double oy[6];
    double oz[6];
	
#ifdef CACHE_B_S
    Matrix ** B_cache;	// cached values of B
    char * B_cache_req;	// B requests
#endif

    StressTensor st;	// stress tensor inside the element
    // (average of the stress tesors near
    //  the surface)
 public:
    long p[6];		// references to the nodes

    double young_mod,	// material properties
	poisson_mod,
	yield_stress,
	fracture_toughness;


    double max_pstress;	// (the max. eigenvalue of st)

    char selected;		// whether this element is selected (used
				// for interactive debugging)
    char mark;		// needed to keep track of which elements
    // are chosen for refinement

    char is_degenerate_flag; // inidicates whether a trianlge is
    // degenerate or not, set in the
    // constructor
    char is_pretty_flag;	// indicates whether a triangle is pretty
				// or not
    char is_below_fp;	// is this element below fracture plane?

    // gaussian points needed for numerical integration
    static struct GP { double w, e, n, s; } * gps;
    static long n_gauss;
    static long n_gauss_surface;
    static char cache_Ke;
    static char fast_Ke_removal;

    enum DrawType { WIREFRAME_TOP, WIREFRAME_ALL,
		    WALLS_TOP, WALLS_ALL};

    WedgeElement() {
	_init();
    }

    WedgeElement( long n1, long n2, long n3,
	long n4, long n5, long n6,
	Map * ymm, Map * prm, Map * ysm, Map* ftm)
    {
	assert( ymm != NULL);
	assert( prm != NULL);
	assert( ysm != NULL);
	assert( ftm != NULL);

	_init();
	p[0] = n1; p[1] = n2; p[2] = n3;
	p[3] = n4; p[4] = n5; p[5] = n6;
	young_mod = ymm-> get_value(
	    nodes[p[0]]-> ox,
	    nodes[p[0]]-> oy,
	    nodes[p[0]]-> oz,
	    nodes[p[1]]-> ox,
	    nodes[p[1]]-> oy,
	    nodes[p[1]]-> oz,
	    nodes[p[2]]-> ox,
	    nodes[p[2]]-> oy,
	    nodes[p[2]]-> oz);
	poisson_mod = prm-> get_value(
	    nodes[p[0]]-> ox,
	    nodes[p[0]]-> oy,
	    nodes[p[0]]-> oz,
	    nodes[p[1]]-> ox,
	    nodes[p[1]]-> oy,
	    nodes[p[1]]-> oz,
	    nodes[p[2]]-> ox,
	    nodes[p[2]]-> oy,
	    nodes[p[2]]-> oz);
	yield_stress = ysm-> get_value(
	    nodes[p[0]]-> ox,
	    nodes[p[0]]-> oy,
	    nodes[p[0]]-> oz,
	    nodes[p[1]]-> ox,
	    nodes[p[1]]-> oy,
	    nodes[p[1]]-> oz,
	    nodes[p[2]]-> ox,
	    nodes[p[2]]-> oy,
	    nodes[p[2]]-> oz);
	fracture_toughness = ftm-> get_value(
	    nodes[p[0]]-> ox,
	    nodes[p[0]]-> oy,
	    nodes[p[0]]-> oz,
	    nodes[p[1]]-> ox,
	    nodes[p[1]]-> oy,
	    nodes[p[1]]-> oz,
	    nodes[p[2]]-> ox,
	    nodes[p[2]]-> oy,
	    nodes[p[2]]-> oz);

	double min_angle = get_top_face_min_angle();
	is_degenerate_flag = min_angle < degenerate_triangle_min_angle;
	is_pretty_flag = min_angle >= pretty_triangle_min_angle;

	// if a triangle is degenerate, set the young modulus to
	// zero - so that the element automatically does not have
	// any effect on the mesh
	if( is_degenerate_flag)
	    young_mod = 0;
		
	// some error checking
	// --------------------------------------------------
	for( long i = 0 ; i < 6 ; i ++)
	    assert( p[i] >= 0 && p[i] < n_nodes);
	for( long i = 0 ; i < 6 ; i ++)
	    for( long j = i + 1 ; j < 6 ; j ++)
		assert( p[i] != p[j]);
    }

    void grow( double time,
	double tf_t0, double tf_t1, double tf_val,
	double bf_t0, double bf_t1, double bf_val,
	double h_t0, double h_val0, double h_t1, double h_val1)
    {
	// name the bottom face points of the reference triangle
	Vector3d p3( nodes[p[3]]-> ox,
	    nodes[p[3]]-> oy,
	    nodes[p[3]]-> oz);
	Vector3d p4( nodes[p[4]]-> ox,
	    nodes[p[4]]-> oy,
	    nodes[p[4]]-> oz);
	Vector3d p5( nodes[p[5]]-> ox,
	    nodes[p[5]]-> oy,
	    nodes[p[5]]-> oz);

	// calculate the current shrinkage amounts
	double top_shrinkage = linfunc( tf_t0, 0, tf_t1, tf_val, time);
	double bot_shrinkage = linfunc( bf_t0, 0, bf_t1, bf_val, time);
	double height_shrinkage =
	    linfunc( h_t0, h_val0, h_t1, h_val1, time);
		
	// figure out the original height
	double oheight = (p3 - Vector3d( nodes[p[0]]-> ox, 
		nodes[p[0]]-> oy, 
		nodes[p[0]]-> oz)).length();

	// figure out the current height
	double height = oheight * ( 1 - height_shrinkage);
		
	// figure out the normal of the reference trinagle
	Vector3d norm = cross_product( p3-p4, p5-p4); norm.normalize();

	// flip the normal if necessary (no guarantee the
	// background triangles will be in counterclockwise order)
	if( norm * (Vector3d( nodes[p[0]]->ox,
		    nodes[p[0]]->oy,
		    nodes[p[0]]->oz)-p3) < 0)
	    norm = -1 * norm;
		
	// extrude the top face
	Vector3d p0 = p3 + norm * height;
	Vector3d p1 = p4 + norm * height;
	Vector3d p2 = p5 + norm * height;

	// figure out the center of the top & bottom face
	Vector3d tc = (p0 + p1 + p2) / 3;
	Vector3d bc = (p3 + p4 + p5) / 3;

	// apply shrinkages to the top face
	p0 = tc + (p0-tc) * (1 - top_shrinkage);
	p1 = tc + (p1-tc) * (1 - top_shrinkage);
	p2 = tc + (p2-tc) * (1 - top_shrinkage);
	// apply shrinkages to the bottom face
	p3 = bc + (p3-bc) * (1 - bot_shrinkage);
	p4 = bc + (p4-bc) * (1 - bot_shrinkage);
	p5 = bc + (p5-bc) * (1 - bot_shrinkage);

	// set the results
	p0.get( ox[0], oy[0], oz[0]);
	p1.get( ox[1], oy[1], oz[1]);
	p2.get( ox[2], oy[2], oz[2]);
	p3.get( ox[3], oy[3], oz[3]);
	p4.get( ox[4], oy[4], oz[4]);
	p5.get( ox[5], oy[5], oz[5]);

	// invalidate cached Ke and RHSe
	if( cache_Ke) {
	    if( Ke != NULL) { delete Ke; Ke = NULL;}
	    if( RHSe != NULL) {delete RHSe; RHSe = NULL; }
	}
		
    }

    virtual ~WedgeElement() {
	// get rid of the stifness matrix
	if( Ke != NULL) delete Ke;
	if( RHSe != NULL) delete RHSe;
		
#ifdef CACHE_B_S
	if( ! optimize_memory_flag) {
	    // delete the cached B's
	    for( long i = 0 ; i < n_gauss ; i ++)
		if( B_cache[i] != NULL) delete B_cache[i];
	    delete [] B_cache;
	    delete [] B_cache_req;
	}
#endif
    }

 private:

    void _init( void) {
	// prepare K
	Ke = NULL;
	RHSe = NULL;
	size = -1;
	p[0] = p[1] = p[2] = p[3] = p[4] = p[5] = -1;
	mark = 0;
	selected = 0;
	max_pstress = -1;

	// prepare gaussian points
	prepare_gaussian_points();

#ifdef CACHE_B_S
	if( ! optimize_memory_flag) {
	    // prepare a cache for B's - an empty one
	    B_cache = new (Matrix *)[n_gauss];
	    B_cache_req = new char[n_gauss];
	    for( long i = 0 ; i < n_gauss ; i ++) {
		B_cache[i] = NULL;
		B_cache_req[i] = 0;
	    }
	} else {
	    B_cache = NULL;
	    B_cache_req = NULL;
	}
#endif
    }

 public:

    double get_size( void) {
	if( size >= 0) return size;
	// get the lenghts of the sides
	double l1 = Node::dist(nodes[p[0]],nodes[p[1]]);
	double l2 = Node::dist(nodes[p[1]],nodes[p[2]]);
	double l3 = Node::dist(nodes[p[2]],nodes[p[0]]);
	if( l1 >= l2 && l1 >= l3) size = l1;
	else if( l2 >= l1 && l2 >= l3) size = l2;
	else size = l3;
	return size;

	// get me the side lengths
	double x1 = nodes[p[0]]-> ox;
	double y1 = nodes[p[0]]-> oy;
	double z1 = nodes[p[0]]-> oz;
	double x2 = nodes[p[1]]-> ox;
	double y2 = nodes[p[1]]-> oy;
	double z2 = nodes[p[1]]-> oz;
	double x3 = nodes[p[2]]-> ox;
	double y3 = nodes[p[2]]-> oy;
	double z3 = nodes[p[2]]-> oz;
	double a = (x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)+(z2-z1)*(z2-z1);
	double b = (x2-x3)*(x2-x3)+(y2-y3)*(y2-y3)+(z2-z3)*(z2-z3);
	double c = (x3-x1)*(x3-x1)+(y3-y1)*(y3-y1)+(z3-z1)*(z3-z1);
	a = sqrt(a); b = sqrt(b); c = sqrt(c);
	double s = (a+b+c)/2.0;
	size = sqrt(s*(s-a)*(s-b)*(s-c));
	return size;
    }

    double get_height_at_node( long n0) const {
	long in0 = get_internal_index( n0);
	assert( in0 >= 0);
	long in0_op = in0 + 3;
	if( in0 > 2) in0_op = in0 - 3;
	double res = sqrt( sqr(ox[in0]-ox[in0_op]) +
	    sqr(oy[in0]-oy[in0_op]) +
	    sqr(oz[in0]-oz[in0_op]));
	return res;
    }

    int has_marked_node( void)
    // returns: 1 if one of the element's nodes is marked
    //	    0 otherwise
    {
	for( long i = 0 ; i < 6 ; i ++)
	    if( nodes[p[i]]-> selected) return 1;
	return 0;
    }

    double get_top_face_min_angle( void)
    // ------------------------------------------------------------------
    // returns the min. angle of the top face
    // ------------------------------------------------------------------
    {
	// name the points of the top face
	Vector3d P[3];
	for( long i = 0 ; i < 3 ; i ++)
	    P[i].set( nodes[p[i]]-> ox,
		nodes[p[i]]-> oy,
		nodes[p[i]]-> oz);
	// angle 1 - at node 0
	Vector3d v1 = P[1]-P[0]; v1.normalize();
	Vector3d v2 = P[2]-P[0]; v2.normalize();
	double a1 = fabs( acos( v1*v2)) * 180 / M_PI;
	// angle 2 - at node 1
	v1 = P[2]-P[1]; v1.normalize();
	v2 = P[0]-P[1]; v2.normalize();
	double a2 = fabs( acos( v1*v2)) * 180 / M_PI;
	// angle 3 - at node 2
	double a3 = 180 - a1 - a2;
		
	// return the minimum
	if( a1 <= a2 && a1 <= a3) return a1;
	if( a2 <= a3 && a2 <= a1) return a2;
	return a3;
    }

    long get_internal_index( long ind) const
    // ------------------------------------------------------------
    // returns the internal index (0..5) of the global node ind
    // ------------------------------------------------------------
    {
	for( long i = 0 ; i < 6 ; i ++)
	    if( p[i] == ind) return i;
	return -1;
    }

    double get_top_face_max_angle( void)
    // ------------------------------------------------------------------
    // returns the min. angle of the top face
    // ------------------------------------------------------------------
    {
	// name the points of the top face
	Vector3d P[3];
	for( long i = 0 ; i < 3 ; i ++)
	    P[i].set( nodes[p[i]]-> ox,
		nodes[p[i]]-> oy,
		nodes[p[i]]-> oz);
	// angle 1 - at node 0
	Vector3d v1 = P[1]-P[0]; v1.normalize();
	Vector3d v2 = P[2]-P[0]; v2.normalize();
	double a1 = fabs( acos( v1*v2)) * 180 / M_PI;
	// angle 2 - at node 1
	v1 = P[2]-P[1]; v1.normalize();
	v2 = P[0]-P[1]; v2.normalize();
	double a2 = fabs( acos( v1*v2)) * 180 / M_PI;
	// angle 3 - at node 2
	double a3 = 180 - a1 - a2;
		
	// return the minimum
	if( a1 >= a2 && a1 >= a3) return a1;
	if( a2 >= a3 && a2 >= a1) return a2;
	return a3;
    }

    void get_shortest_top_edge( long & n1, long & n2)
    // ------------------------------------------------------------------
    // returns global node indecies of the shortest top edge
    // ------------------------------------------------------------------
    {
	double l1 = Node::dist(nodes[p[0]],nodes[p[1]]);
	double l2 = Node::dist(nodes[p[1]],nodes[p[2]]);
	double l3 = Node::dist(nodes[p[2]],nodes[p[0]]);
	if( l1 <= l2 && l1 <= l3) {
	    n1 = p[0]; n2 = p[1];
	} else if( l2 <= l1 && l2 <= l3) {
	    n1 = p[1]; n2 = p[2];
	} else {
	    n1 = p[2]; n2 = p[0];
	}
    }

    StressTensor get_nodal_stress_global_ind( long ind)
    // ------------------------------------------------------------------
    // calculates the stress at node ind (ind is the global number)
    // ------------------------------------------------------------------
    {
	// find the index of this node
	long i = -1;
	for( long k = 0 ; k < 6 ; k ++)
	    if( p[k] == ind) {
		i = k;
		break;
	    }
	assert( i != -1);

	if( 1) {
	    long transl[6] = {0,1,2,6,7,8};
	    i = transl[i];
	    StressTensor res;
	    calculate_stress_tensor_at_gp( res, i);
	    return res;
	} else {
	    // calculate the stress at this node
	    double iso[6][3] = {
		// e,   n,    s
		{  1,   0,   -1},
		{  0,   1,   -1},
		{  0,   0,   -1},
		{  1,   0,    1},
		{  0,   1,    1},
		{  0,   0,    1}};
	    StressTensor res;
	    calculate_stress_tensor_at_iso( res,
		iso[i][0],
		iso[i][1],
		iso[i][2]);
	    return res;
	}
    }

    double get_quality( void)
    // ------------------------------------------------------------------
    // returns quality of the element
    //  - quality = average of top and bottom face triangle quality
    //  - triangle quality:
    //      - all 60 degree angles --> quality = 1
    //      - one angle 180, the other 0 --> quality = 0
    //      - algorithm:
    //          - get the difference between the min. and max. angles
    //            in a triangle
    //          - divide the difference by 180, obtaining a number
    //            between 0 and 1 (0 indicating all angles equal,
    //            1 indicating max. difference in angles
    // ------------------------------------------------------------------
    {
	Vector3d P[6];
	for( long i = 0 ; i < 6 ; i ++)
	    P[i].set( ox[i], oy[i], oz[i]);

	double q1 = triangle_quality( P[0], P[1], P[2]);
	double q2 = triangle_quality( P[3], P[4], P[5]);
	double q = (q1 + q2) / 2;

	// 		q *= q;

	return q;
    }

    double get_angle_at_global_ind( long ind)
    // ------------------------------------------------------------------
    // calculates the angle of the triangle at node 'ind'
    // ------------------------------------------------------------------
    {
	// find the index of this node
	long i = -1;
	for( long k = 0 ; k < 6 ; k ++)
	    if( p[k] == ind) {
		i = k;
		break;
	    }
	assert( i != -1);

	// assume top face
	long n1 = (i+0)%3;
	long n2 = (i+1)%3;
	long n3 = (i+2)%3;
	// if bottom face, add 3
	if( i > 2) {
	    n1 += 3;
	    n2 += 3;
	    n3 += 3;
	}
		
	Vector3d v1( ox[n2]-ox[n1], oy[n2]-oy[n1], oz[n2]-oz[n1]);
	Vector3d v2( ox[n3]-ox[n1], oy[n3]-oy[n1], oz[n3]-oz[n1]);
	v1.normalize();
	v2.normalize();
	double alpha = fabs( acos( v1*v2)) * 180 / M_PI;

	return alpha;
    }

    double get_projected_area( Vector3d & pn)
    // ------------------------------------------------------------------
    // calculates the area of an element projected onto a plane
    // whose normal is 'pn'
    // ------------------------------------------------------------------
    {
	// if the normal is degenerate (i.e. too close to zero)
	// just return zero
	if( pn.length() < 1e-10) return 0;

	// name the points of the elements
	Vector3d pt[6];
	for( long i = 0 ; i < 6 ; i ++)
	    pt[i].set( nodes[p[i]]-> ox,
		nodes[p[i]]-> oy,
		nodes[p[i]]-> oz);
	// create a triangle list - a closed manifold (hopefuly)
	// with proper orientation of elements
	TriangleList list( 8);
	list.add( pt[0], pt[1], pt[2]); // top face
	list.add( pt[4], pt[3], pt[5]); // bottom face
	list.add( pt[0], pt[4], pt[1]);
	list.add( pt[0], pt[3], pt[4]);
	list.add( pt[1], pt[5], pt[2]);
	list.add( pt[1], pt[4], pt[5]);
	list.add( pt[2], pt[3], pt[0]);
	list.add( pt[2], pt[5], pt[3]);
	// get the projection
	return projected_area( pn, list);
    }
	
    void prepare_gaussian_points( void)
    {
	// if gaussian points were already prepared, skip this
	if( gps != NULL) return;
		
	// __________initialze gauss points__________
	if( 0) {
	    n_gauss = 6;
	    n_gauss_surface = 3;
	    gps = new GP[ n_gauss];
	    // 1st GAUSS point
	    gps[0].w = 1/6.0 * 1.0;
	    gps[0].e = 2/3.0; gps[0].n = 1/6.0; gps[0].s = -0.5773502692;
	    // 2nd GAUSS point
	    gps[1].w = 1/6.0 * 1.0;
	    gps[1].e = 1/6.0; gps[1].n = 2/3.0; gps[1].s = -0.5773502692;
	    // 3rd GAUSS point
	    gps[2].w = 1/6.0 * 1.0;
	    gps[2].e = 1/6.0; gps[2].n = 1/3.0; gps[2].s = -0.5773502692;
	    // 4th GAUSS point
	    gps[3].w = 1/6.0 * 1.0;
	    gps[3].e = 2/3.0; gps[3].n = 1/6.0; gps[3].s = 0.5773502692;
	    // 5th GAUSS point
	    gps[4].w = 1/6.0 * 1.0;
	    gps[4].e = 1/6.0; gps[4].n = 2/3.0; gps[4].s = 0.5773502692;
	    // 6th GAUSS point
	    gps[5].w = 1/6.0 * 1.0;
	    gps[5].e = 1/6.0; gps[5].n = 1/3.0; gps[5].s = 0.5773502692;

	} else {

	    n_gauss = 9;
	    n_gauss_surface = 3;
	    gps = new GP[n_gauss];
	    long i = 0;
	    // 1st GAUSS point
	    gps[i].w = 1/6.0 * 0.5555555556;
	    gps[i].e = 2/3.0; gps[i].n = 1/6.0; gps[i].s = -0.7745966692;
	    i ++;
	    // 2nd GAUSS point
	    gps[i].w = 1/6.0 * 0.5555555556;
	    gps[i].e = 1/6.0; gps[i].n = 2/3.0; gps[i].s = -0.7745966692;
	    i ++;
	    // 3rd GAUSS point
	    gps[i].w = 1/6.0 * 0.5555555556;
	    //		gps[i].e = 1/6.0; gps[i].n = 1/3.0; gps[i].s = -0.7745966692;
	    gps[i].e = 1/6.0; gps[i].n = 1/6.0; gps[i].s = -0.7745966692;
	    i ++;
	    // 4th GAUSS point
	    gps[i].w = 1/6.0 * 0.8888888889;
	    gps[i].e = 2/3.0; gps[i].n = 1/6.0; gps[i].s = 0;
	    i ++;
	    // 5th GAUSS point
	    gps[i].w = 1/6.0 * 0.8888888889;
	    gps[i].e = 1/6.0; gps[i].n = 2/3.0; gps[i].s = 0;
	    i ++;
	    // 6th GAUSS point
	    gps[i].w = 1/6.0 * 0.8888888889;
	    //		gps[i].e = 1/6.0; gps[i].n = 1/3.0; gps[i].s = 0;
	    gps[i].e = 1/6.0; gps[i].n = 1/6.0; gps[i].s = 0;
	    i ++;
	    // 7th GAUSS point
	    gps[i].w = 1/6.0 * 0.5555555556;
	    gps[i].e = 2/3.0; gps[i].n = 1/6.0; gps[i].s = 0.7745966692;
	    i ++;
	    // 8th GAUSS point
	    gps[i].w = 1/6.0 * 0.5555555556;
	    gps[i].e = 1/6.0; gps[i].n = 2/3.0; gps[i].s = 0.7745966692;
	    i ++;
	    // 9th GAUSS point
	    gps[i].w = 1/6.0 * 0.5555555556;
	    //		gps[i].e = 1/6.0; gps[i].n = 1/3.0; gps[i].s = 0.7745966692;
	    gps[i].e = 1/6.0; gps[i].n = 1/6.0; gps[i].s = 0.7745966692;
	    i ++;
	}

    }

    void read_data( Tokenizer & tok) {
	(void) tok.read_long( "ID");
	p[0] = tok.read_long( "element-> p[0]");
	p[1] = tok.read_long( "element-> p[1]");
	p[2] = tok.read_long( "element-> p[2]");
	p[3] = tok.read_long( "element-> p[3]");
	p[4] = tok.read_long( "element-> p[4]");
	p[5] = tok.read_long( "element-> p[5]");
		
	// some error checking
	// --------------------------------------------------
	for( long i = 0 ; i < 6 ; i ++)
	    if( p[i] < 0 || p[i] >= n_nodes)
		die( "Element's p[%ld] invalid", i);
	for( long i = 0 ; i < 6 ; i ++)
	    for( long j = i + 1 ; j < 6 ; j ++)
		if( p[i] == p[j])
		    die( "Element: duplicate nodes.");

	// 		// adjust references to the nodes
	// 		// --------------------------------------------------
	// 		for( long i = 0 ; i < 6 ; i ++)
	// 			nodes[p[i]]-> add_ref( this);
    }

    void add_to_global( SparseMatrix & Kg, Vector_double & RHSg) {
	assert( Kg.n_rows == n_nodes*3);
	assert( Kg.n_cols == n_nodes*3);
	assert( RHSg.size == n_nodes*3);

	// obtain stiffness matrix and RHS
	Matrix K(18,18);
	Vector_double RHS( 18);
	calculate_K( K, RHS);

	for( long rn = 0 ; rn < 6 ; rn ++)
	    for( long cn = 0 ; cn < 6 ; cn ++)
		for( long i = 0 ; i < 3 ; i ++)
		    for( long j = 0 ; j < 3 ; j ++)
			Kg.add(p[rn]*3+i,
			    p[cn]*3+j,
			    K.data[rn*3+i][cn*3+j]);
	// now add the RHS of the element to the appropriate places
	// of the global RHS
	for( long rn = 0 ; rn < 6 ; rn ++) {
	    RHSg.array[p[rn]*3+0] += RHS.array[rn*3+0];
	    RHSg.array[p[rn]*3+1] += RHS.array[rn*3+1];
	    RHSg.array[p[rn]*3+2] += RHS.array[rn*3+2];
	}
    }

    void remove_from_global(
	SparseMatrix & Kg, Vector_double & RHSg)

    {
	if( fast_Ke_removal) {
	    _remove_from_global_fast( Kg, RHSg);
	} else {
	    _remove_from_global_precise( Kg, RHSg);
	}
    }

    void _remove_from_global_fast(
	SparseMatrix & Kg, Vector_double & RHSg)
    {
	assert( Kg.n_rows == n_nodes*3);
	assert( Kg.n_cols == n_nodes*3);

	// get element's stiffness matrix and RHS
	Matrix K(18,18);
	Vector_double RHS( 18);
	calculate_K( K, RHS);

	for( long rn = 0 ; rn < 6 ; rn ++)
	    for( long cn = 0 ; cn < 6 ; cn ++)
		for( long i = 0 ; i < 3 ; i ++)
		    for( long j = 0 ; j < 3 ; j ++)
			Kg.add(p[rn]*3+i,
			    p[cn]*3+j,
			    -K.data[rn*3+i][cn*3+j]);
	// now remove the RHS of the element from the appropriate
	// places of the global RHS
	for( long rn = 0 ; rn < 6 ; rn ++) {
	    RHSg.array[p[rn]*3+0] -= RHS.array[rn*3+0];
	    RHSg.array[p[rn]*3+1] -= RHS.array[rn*3+1];
	    RHSg.array[p[rn]*3+2] -= RHS.array[rn*3+2];
	}
    }

    void _remove_from_global_precise(
	SparseMatrix & Kg, Vector_double & RHSg)
    // removes the element's Ke from Kg
    //
    // Remark:
    //    a naive way to do this would be to to just subtract
    //    Ke from Kg, but that would eventually produce numerical errors
    //    (consider summing a list of doubles of various magnitudes,
    //     and then subtracting them in a random order - you don't end
    //     up with a ZERO!!!)
    //
    // A better (correct) way is to rebuild the entries of Kg that
    // the element affected.
    //
    // Method:
    //    - reset all elements of Kg affected by this element to ZERO
    //    - construct a list of all elements sharing two nodes with
    //      this element (an edge or face neighbor)
    //    - recalculate all entries of Kg affected by both this element
    //      and the elements on the list - by adding their appropriate
    //      entries in elemental Ke's
    {
	assert( Kg.n_rows == n_nodes*3);
	assert( Kg.n_cols == n_nodes*3);
	assert( RHSg.size == n_nodes*3);

	// zero out entries affected by this element
	for( long rn = 0 ; rn < 6 ; rn ++) {
	    long grn = p[rn];
	    for( long cn = 0 ; cn < 6 ; cn ++) {
		long gcn = p[cn];
		for( long i = 0 ; i < 3 ; i ++)
		    for( long j = 0 ; j < 3 ; j ++)
			Kg.set( grn*3+i, gcn*3+j, 0);
	    }
	    RHSg.array[grn*3+0] = 0;
	    RHSg.array[grn*3+1] = 0;
	    RHSg.array[grn*3+2] = 0;
	}

	// create a list of elements which share a node with
	// this element
	Stack<WedgeElement *> elist;
	for( long i = 0 ; i < 6 ; i ++) {
	    for( long j = 0 ; j < nodes[p[i]]-> n_refs ; j ++) {
		WedgeElement * e = nodes[p[i]]-> refs[j];
		if( e == this) continue;
		elist.push_unique( e);
	    }
	}

	// for each element on the list, add its appropriate
	// contributions to the erased parts of Kg and RHSg
	Matrix K(18,18); Vector_double RHS(18);
	while( 1) {
	    if( elist.is_empty()) break;
	    WedgeElement * e = elist.pop();
	    // get element's K and RHS
	    e-> calculate_K( K, RHS);
	    // process each entry of K
	    for( long rn = 0 ; rn < 6 ; rn ++) {
		long grn = e-> p[rn];
		// grn has to be in this element, if not
		// this row can be skipped
		if( grn != p[0] &&
		    grn != p[1] &&
		    grn != p[2] &&
		    grn != p[3] &&
		    grn != p[4] &&
		    grn != p[5]) continue;
		for( long cn = 0 ; cn < 6 ; cn ++) {
		    long gcn = e-> p[cn];
		    // gcn also has to be in this
		    // element
		    if( gcn != p[0] &&
			gcn != p[1] &&
			gcn != p[2] &&
			gcn != p[3] &&
			gcn != p[4] &&
			gcn != p[5]) continue;
		    // add the element's K[rn,cn] to
		    // Kg[grn,gcn]
		    for( long i = 0 ; i < 3 ; i ++)
			for( long j = 0 ; j < 3 ; j ++)
			    Kg.add( grn*3+i, gcn*3+j,
				K.data[rn*3+i][cn*3+j]);
		}
				
		// also process RHS
		RHSg.array[grn*3+0] += RHS.array[rn*3+0];
		RHSg.array[grn*3+1] += RHS.array[rn*3+1];
		RHSg.array[grn*3+2] += RHS.array[rn*3+2];
	    }
	}
    }

    void calculate_J( Matrix & J, double e, double n, double s) {
	assert( J.n_cols == 3);
	assert( J.n_rows == 3);
		
	Matrix Nens(3, 6); calculate_Nens( Nens, e, n, s);
	// Nens.print( "Nens= ", "%7.2f");
	Matrix Q( 6, 3); calculate_Q( Q);
	// J = Nens x Q
	Nens.multiply( Q, J);
	// J.print( "J= ", "%7.2f");
    }

	
    void calculate_D( Matrix & D) {
	assert( D.n_cols == 6);
	assert( D.n_rows == 6);
	// calculate D
	D.reset();
	D.data[  0][ 0] =  1.0 - poisson_mod;
	D.data[  1][ 1] =  1.0 - poisson_mod;
	D.data[  2][ 2] =  1.0 - poisson_mod;
	D.data[  3][ 3] =  0.5 - poisson_mod;
	D.data[  4][ 4] =  0.5 - poisson_mod;
	D.data[  5][ 5] =  0.5 - poisson_mod;
	D.data[  0][ 1] =  poisson_mod;
	D.data[  0][ 2] =  poisson_mod;
	D.data[  1][ 0] =  poisson_mod;
	D.data[  1][ 2] =  poisson_mod;
	D.data[  2][ 0] =  poisson_mod;
	D.data[  2][ 1] =  poisson_mod;
	D.multiply_by_scalar( young_mod /
	    ((1+poisson_mod)*(1-2*poisson_mod)));
    }

    // shape function derivatives w.r.t. isoparametric coordinates
    double Nae( double e, double n, double s) { return 0.5*(1-s);  }
    double Nbe( double e, double n, double s) { return 0;          }
    double Nce( double e, double n, double s) { return -0.5*(1-s); }
    double Nde( double e, double n, double s) { return 0.5*(1+s);  }
    double Nee( double e, double n, double s) { return 0;          }
    double Nfe( double e, double n, double s) { return -0.5*(1+s); }

    double Nan( double e, double n, double s) { return 0;          }
    double Nbn( double e, double n, double s) { return 0.5*(1-s);  }
    double Ncn( double e, double n, double s) { return -0.5*(1-s); }
    double Ndn( double e, double n, double s) { return 0;          }
    double Nen( double e, double n, double s) { return 0.5*(1+s);  }
    double Nfn( double e, double n, double s) { return -0.5*(1+s); }

    double Nas( double e, double n, double s) { return -0.5*e;     }
    double Nbs( double e, double n, double s) { return -0.5*n;     }
    double Ncs( double e, double n, double s) { return -0.5*(1-e-n);}
    double Nds( double e, double n, double s) { return 0.5*e;      }
    double Nes( double e, double n, double s) { return 0.5*n;      }
    double Nfs( double e, double n, double s) { return 0.5*(1-e-n);}
	
    // shape functions
    double Na( double e, double n, double s) { return 0.5*e*(1-s); }
    double Nb( double e, double n, double s) { return 0.5*n*(1-s); }
    double Nc( double e, double n, double s) { return 0.5*(1-e-n)*(1-s); }
    double Nd( double e, double n, double s) { return 0.5*e*(1+s); }
    double Ne( double e, double n, double s) { return 0.5*n*(1+s); }
    double Nf( double e, double n, double s) { return 0.5*(1-e-n)*(1+s); }
    // shortcuts
    double Na( GP & gp) { return Na( gp.e, gp.n, gp.s); }
    double Nb( GP & gp) { return Nb( gp.e, gp.n, gp.s); }
    double Nc( GP & gp) { return Nc( gp.e, gp.n, gp.s); }
    double Nd( GP & gp) { return Nd( gp.e, gp.n, gp.s); }
    double Ne( GP & gp) { return Ne( gp.e, gp.n, gp.s); }
    double Nf( GP & gp) { return Nf( gp.e, gp.n, gp.s); }


    void calculate_Nens( Matrix & Nens, double e, double n, double s) {
	assert( Nens.n_cols == 6);
	assert( Nens.n_rows == 3);

	Nens.data[0][0] = Nae(e,n,s);
	Nens.data[0][1] = Nbe(e,n,s);
	Nens.data[0][2] = Nce(e,n,s);
	Nens.data[0][3] = Nde(e,n,s);
	Nens.data[0][4] = Nee(e,n,s);
	Nens.data[0][5] = Nfe(e,n,s);

	Nens.data[1][0] = Nan(e,n,s);
	Nens.data[1][1] = Nbn(e,n,s);
	Nens.data[1][2] = Ncn(e,n,s);
	Nens.data[1][3] = Ndn(e,n,s);
	Nens.data[1][4] = Nen(e,n,s);
	Nens.data[1][5] = Nfn(e,n,s);

	Nens.data[2][0] = Nas(e,n,s);
	Nens.data[2][1] = Nbs(e,n,s);
	Nens.data[2][2] = Ncs(e,n,s);
	Nens.data[2][3] = Nds(e,n,s);
	Nens.data[2][4] = Nes(e,n,s);
	Nens.data[2][5] = Nfs(e,n,s);
    }

    void calculate_Q( Matrix & Q) {
	assert( Q.n_cols == 3);
	assert( Q.n_rows == 6);
	Q.data[0][0] = ox[0];
	Q.data[0][1] = oy[0];
	Q.data[0][2] = oz[0];
	Q.data[1][0] = ox[1];
	Q.data[1][1] = oy[1];
	Q.data[1][2] = oz[1];
	Q.data[2][0] = ox[2];
	Q.data[2][1] = oy[2];
	Q.data[2][2] = oz[2];
	Q.data[3][0] = ox[3];
	Q.data[3][1] = oy[3];
	Q.data[3][2] = oz[3];
	Q.data[4][0] = ox[4];
	Q.data[4][1] = oy[4];
	Q.data[4][2] = oz[4];
	Q.data[5][0] = ox[5];
	Q.data[5][1] = oy[5];
	Q.data[5][2] = oz[5];
    }

    void get_orig_shape( Vector3d (&p)[6]) {
	for( long i = 0 ; i < 6 ; i ++)
	    p[i].set( ox[i], oy[i], oz[i]);
    }

    Vector3d iso_to_xyz( GP & gp)
    // -----------------------------------------------------------------
    // calculate the cartesian coordinates of the gaussian point
    // in the undeformed element (xyz = [N(gp)].[Q]
    // -----------------------------------------------------------------
    {
	// calculate shape function values at gp
	double na = Na( gp);
	double nb = Nb( gp);
	double nc = Nc( gp);
	double nd = Nd( gp);
	double ne = Ne( gp);
	double nf = Nf( gp);
	// name the points
	Vector3d P[6];
	for( long i = 0 ; i < 6 ; i ++)
	    P[i].set( nodes[p[i]]-> ox,
		nodes[p[i]]-> oy,
		nodes[p[i]]-> oz);
	Vector3d res =
	    na*P[0] +
	    nb*P[1] +
	    nc*P[2] +
	    nd*P[3] +
	    ne*P[4] +
	    nf*P[5];
	return res;
    }


    void calculate_B_at_gp( Matrix & B, long ind) {
	assert( ind >= 0 && ind < n_gauss);

	calculate_B_at_iso( B,
	    gps[ind].e,
	    gps[ind].n,
	    gps[ind].s);
#ifdef CACHE_B_S
	if( optimize_memory_flag) {
	} else {
	    if( B_cache[ind] == NULL) {
		calculate_B_at_iso( B,
		    gps[ind].e,
		    gps[ind].n,
		    gps[ind].s);
		// if B at this point was requested more
		// than once,
		// cache it
		if( B_cache_req[ind] < 10)
		    B_cache_req[ind] ++;
		else
		    B_cache[ind] = new Matrix( B);
	    } else {
		B.copy( * B_cache[ind]);
	    }
	}
#endif
    }

    void calculate_B_at_iso( Matrix & B, double e, double n, double s) {
	// calculate J
	Matrix J( 3, 3); calculate_J( J, e, n, s);
	// calculate Jm = inverse( J)
	Matrix Jm( 3, 3); if( ! is_degenerate()) J.inverse( Jm);
	// calculate Nens
	Matrix Nens( 3, 6); calculate_Nens( Nens, e, n, s);
	// calculate Nxyz = J^-1 x Nens
	Matrix Nxyz( 3, 6);
	Jm.multiply( Nens, Nxyz);

	// now substitute all values from Nxyz to B
	double nax = Nxyz.data[0][0];
	double nbx = Nxyz.data[0][1];
	double ncx = Nxyz.data[0][2];
	double ndx = Nxyz.data[0][3];
	double nex = Nxyz.data[0][4];
	double nfx = Nxyz.data[0][5];
	double nay = Nxyz.data[1][0];
	double nby = Nxyz.data[1][1];
	double ncy = Nxyz.data[1][2];
	double ndy = Nxyz.data[1][3];
	double ney = Nxyz.data[1][4];
	double nfy = Nxyz.data[1][5];
	double naz = Nxyz.data[2][0];
	double nbz = Nxyz.data[2][1];
	double ncz = Nxyz.data[2][2];
	double ndz = Nxyz.data[2][3];
	double nez = Nxyz.data[2][4];
	double nfz = Nxyz.data[2][5];
		
	double a[6][18] = {
	    { nax,0,0,nbx,0,0,ncx,0,0,ndx,0,0,nex,0,0,nfx,0,0},
	    { 0,nay,0,0,nby,0,0,ncy,0,0,ndy,0,0,ney,0,0,nfy,0},
	    { 0,0,naz,0,0,nbz,0,0,ncz,0,0,ndz,0,0,nez,0,0,nfz},
	    { 0,naz,nay,0,nbz,nby,0,ncz,ncy,0,ndz,ndy,0,nez,ney,0,nfz,nfy},
	    { naz,0,nax,nbz,0,nbx,ncz,0,ncx,ndz,0,ndx,nez,0,nex,nfz,0,nfx},
	    { nay,nax,0,nby,nbx,0,ncy,ncx,0,ndy,ndx,0,ney,nex,0,nfy,nfx,0}
	};

	B.reset();
	for( long r = 0 ; r < 6 ; r ++)
	    for( long c = 0 ; c < 18 ; c ++)
		B.data[r][c] = a[r][c];
    }

    void calculate_K( Matrix & K, Vector_double & RHS) {

	// if Ke and RHSe are being cached, and they are cached
	// right now, return the cached values
	if( cache_Ke) {
	    if( Ke != NULL) {
		assert( RHSe != NULL);
		K.copy( * Ke);
		RHS.copy( * RHSe);
		return;
	    }
	}

	ke_recalc_n ++;

	// zero out the K matrix
	K.reset();

	// calculate D for this material - it will be the same
	// anywhere in the loop
	Matrix D(6,6); calculate_D( D);

	// make temporary variables, so that we don't have
	// to create them every time
	Matrix J(3,3);
	Matrix B(6,18);
	Matrix Bt(18,6);
	Matrix BtD(18,6);
	Matrix BtDB(18,18);

	// for every gauss point add the corresponding B D B detJ
	for( long i = 0 ; i < n_gauss ; i ++) {
	    // calculate B at a gauss point[i]
	    //			Matrix B(6,18);
	    calculate_B_at_gp( B, i);
	    // calculate Bt = transpose B
	    //			Matrix Bt(18,6);
	    B.transpose( Bt);
	    // calculate D for this material
	    //			Matrix D(6,6); calculate_D( D);
	    // calculate Bt x D
	    //			Matrix BtD(18,6);
	    Bt.multiply( D,BtD);
	    // calculate Bt x D x B
	    //			Matrix BtDB(18,18);
	    BtD.multiply(B,BtDB);
	    // calculate det J at a gauss point
	    //			Matrix J(3,3);
	    calculate_J( J, gps[i].e, gps[i].n, gps[i].s);
	    double detJ = fabs( J.determinant());
	    // 			fprintf( stderr, "detJ(%d) = %.30f\n",
	    // 				 int( is_degenerate_flag),
	    // 				 detJ);
	    if( isnan( detJ)) {
		fprintf( stderr, "detJ = NaN !!!");
		fprintf( stderr,
		    "p[] = %ld %ld %ld %ld %ld %ld\n",
		    p[0],p[1],p[2],p[3],p[4],p[5]);
		fprintf( stderr,
		    "coords:\n");
		for( long i = 0 ; i < 6 ; i ++)
		    fprintf( stderr,
			"%ld) o: %.10f %.10f %.10f "
			"p: %.10f %.10f %.10f\n",
			i,
			nodes[p[i]]-> ox,
			nodes[p[i]]-> oy,
			nodes[p[i]]-> oz,
			nodes[p[i]]-> px,
			nodes[p[i]]-> py,
			nodes[p[i]]-> pz);
		detJ = 0;
	    }
	    // multiply BtDB by detJ
	    BtDB.multiply_by_scalar( -detJ);
	    // multiply by a gauss weight
	    BtDB.multiply_by_scalar( gps[i].w);
	    // add the result to K
	    K.add(BtDB);
	}

	// now construct the RHS = K.Q
	Vector_double Q(18);
	for( long i = 0 ; i < 6 ; i ++) {
	    Q.array[i*3+0] = ox[i];
	    Q.array[i*3+1] = oy[i];
	    Q.array[i*3+2] = oz[i];
	}
	K.mult( Q, RHS);

	// cache the results
	if( cache_Ke) {
	    if( Ke == NULL) {
		Ke = new Matrix( K);
		assert( RHSe == NULL);
		RHSe = new Vector_double( RHS);
	    } else {
		Ke-> copy( K);
		assert( RHSe != NULL);
		RHSe-> copy( RHS);
	    }
	}
    }

    void calculate_Disp( Matrix & Disp) {
	assert( Disp.n_rows == 18);
	assert( Disp.n_cols = 1);
	for( long i = 0 ; i < 6 ; i ++)	{
	    Disp.data[ i*3+0][ 0] = nodes[p[i]]->px - ox[i];
	    Disp.data[ i*3+1][ 0] = nodes[p[i]]->py - oy[i];
	    Disp.data[ i*3+2][ 0] = nodes[p[i]]->pz - oz[i];
	}
    }

    void calculate_strain_tensor( StressTensor & s,
	double iso_e,
	double iso_n,
	double iso_s)
    // --------------------------------------------------
    // calculates stress tensor into Strain
    // --------------------------------------------------
    {
	// Get the strain matrix
	Matrix B( 6, 18);
	// calculate B at the requested isoparametric coordinate
	calculate_B_at_iso( B, iso_e, iso_n, iso_s);
	// get the displacement matrix Disp
	Matrix Disp( 18, 1);
	calculate_Disp( Disp);
	// calculate Strain = B * Disp
	Matrix tmpstrain(6,1);
	B.multiply( Disp, tmpstrain);
	// extract the stress from the Stress matrix
	s.set( tmpstrain.data[0][0], // sxx
	    tmpstrain.data[1][0], // syy
	    tmpstrain.data[2][0], // szz
	    tmpstrain.data[3][0], // syz
	    tmpstrain.data[4][0], // sxz
	    tmpstrain.data[5][0]); // sxy
    }

    void calculate_stress_tensor_at_gp( StressTensor & s, long ind) {
	// --------------------------------------------------
	// calculates stress tensor into Stress
	// --------------------------------------------------
	{
	    // preconditions
	    assert( ind >= 0 && ind < n_gauss);

	    // get [B]
	    Matrix B( 6, 18);
	    calculate_B_at_gp( B, ind);
	    // get the displacement matrix Disp
	    Matrix Disp( 18, 1);
	    calculate_Disp( Disp);
	    // calculate Strain = B * Disp
	    Matrix Strain( 6, 1);
	    B.multiply( Disp, Strain);
	    // calculate the material matrix D for this element
	    Matrix D(6,6);
	    calculate_D( D);
	    // calculate Stress = D * Strain
	    Matrix Stress(6,1);
	    D.multiply( Strain, Stress);
		
	    // extract the stress from the Stress matrix
	    s.set( Stress.data[0][0], // sxx
		Stress.data[1][0], // syy
		Stress.data[2][0], // szz
		Stress.data[3][0], // syz
		Stress.data[4][0], // sxz
		Stress.data[5][0]); // sxy
	}
		
    }

    void calculate_stress_tensor_at_iso( StressTensor & s,
	double iso_e,
	double iso_n,
	double iso_s)
    // --------------------------------------------------
    // calculates stress tensor into Stress
    // --------------------------------------------------
    {
	// Get the strain matrix
	Matrix B( 6, 18);
	// calculate B at the requested isoparametric coordinate
	calculate_B_at_iso( B, iso_e, iso_n, iso_s);
	// get the displacement matrix Disp
	Matrix Disp( 18, 1);
	calculate_Disp( Disp);
	// calculate Strain = B * Disp
	Matrix Strain( 6, 1);
	B.multiply( Disp, Strain);
	// calculate the material matrix D for this element
	Matrix D(6,6);
	calculate_D( D);
	// calculate Stress = D * Strain
	Matrix Stress(6,1);
	D.multiply( Strain, Stress);
		
	// extract the stress from the Stress matrix
	s.set( Stress.data[0][0], // sxx
	    Stress.data[1][0], // syy
	    Stress.data[2][0], // szz
	    Stress.data[3][0], // syz
	    Stress.data[4][0], // sxz
	    Stress.data[5][0]); // sxy
    }

    StressTensor & get_stress_tensor( void) 
    {
	return st;
    }
	
    void calculate_stress_tensor( void)
    // --------------------------------------------------
    // calculate the current stress tensor of the element into 'st'
    // calculate max. princ. stress (max. eigenvalue of 'st')
    {
	// calculate the stress tensor 'st'
	st.set( 0, 0, 0, 0, 0, 0);
	StressTensor s;
	for( long i = 0 ; i < n_gauss_surface ; i ++) {
	    calculate_stress_tensor_at_gp( s, i);
	    // add s to Stress
	    st.add( s);
	}
	st.multiply_by_scalar( 1.0 / n_gauss_surface);

	// calculate the max. principal stress
	double s1, s2, s3;
	s.get_eigenvalues( s1, s2, s3);
	max_pstress = s3;
    }
	
    void draw_triangle_normal( double x0, double y0, double z0,
	double x1, double y1, double z1,
	double x2, double y2, double z2)
    {
	// figure out the normal
	double ux = x1 - x0; double uy = y1 - y0; double uz = z1 - z0;
	double vx = x2 - x0; double vy = y2 - y0; double vz = z2 - z0;
	double nx = uy * vz - uz * vy;
	double ny = uz * vx - ux * vz;
	double nz = ux * vy - uy * vx;

	// set the normal
	glNormal3d( -nx, -ny, -nz);
    }

    Vector3d get_top_face_normal( void) {
	double x0 = nodes[p[0]]-> ox;
	double y0 = nodes[p[0]]-> oy;
	double z0 = nodes[p[0]]-> oz;
	double x1 = nodes[p[1]]-> ox;
	double y1 = nodes[p[1]]-> oy;
	double z1 = nodes[p[1]]-> oz;
	double x2 = nodes[p[2]]-> ox;
	double y2 = nodes[p[2]]-> oy;
	double z2 = nodes[p[2]]-> oz;
	// figure out the normal
	double ux = x1 - x0; double uy = y1 - y0; double uz = z1 - z0;
	double vx = x2 - x0; double vy = y2 - y0; double vz = z2 - z0;
	Vector3d res( 
	    uy * vz - uz * vy,
	    uz * vx - ux * vz,
	    ux * vy - uy * vx);
	res.normalize();
	return res;
		
    }

    Vector3d get_bottom_face_normal( void)
    // -------------------------------------------------------------------
    // calculate the normal of the bottom surface (the triangle), using
    // deformed coordinates
    // -------------------------------------------------------------------
    {
	Vector3d P[3] = {
	    Vector3d( nodes[p[3]]-> px,
		nodes[p[3]]-> py,
		nodes[p[3]]-> pz),
	    Vector3d( nodes[p[4]]-> px,
		nodes[p[4]]-> py,
		nodes[p[4]]-> pz),
	    Vector3d( nodes[p[5]]-> px,
		nodes[p[5]]-> py,
		nodes[p[5]]-> pz)
	};

	Vector3d res = cross_product(P[1]-P[0],P[2]-P[0]);
	res.normalize();
	return res;
    }

    void calculate_nodal_forces( Matrix & F)
    // --------------------------------------------------------------
    // Calculate nodal forces exerted by this element on all of
    // its nodes ( [F] = [K][disp] )
    // --------------------------------------------------------------
    {
	assert( F.n_rows == 18 && F.n_cols == 1);
	Matrix K(18,18); Vector_double RHS(18);
	calculate_K( K, RHS);
	//		Matrix disp(18,1);
	//		calculate_Disp( disp);
	//		K.multiply( disp, F);
	Matrix pos(18,1);
	for( long i = 0 ; i < 6 ; i ++) {
	    pos.data[i*3+0][0] = nodes[p[i]]-> px;
	    pos.data[i*3+1][0] = nodes[p[i]]-> py;
	    pos.data[i*3+2][0] = nodes[p[i]]-> pz;
	}
	K.multiply( pos, F);
	//		F.sub( RHS);
	for( long i = 0 ; i < 18 ; i ++)
	    F.data[i][0] -= RHS.array[i];
    }

    void calculate_nodal_unit_area_forces(
	Vector3d & f1, Vector3d & f2, Vector3d & f3)
    // -------------------------------------------------------------
    // Calculates the forces acting on p[0..2] in per/unit area
    //   - first calculate the total forces acting on each node 0-3
    //   - calculate the forces acting on each of the three walls
    //   - divide the wall forces by the wall area
    //   - compute the result (each node force is half of one wall
    //     force + half of the other wall force)	
    // -------------------------------------------------------------
    {
	// calculate the total nodal forces (tnf1..3)
	Matrix F(18,1);
	calculate_nodal_forces( F);
	Vector3d tnf1( F.get(0,0), F.get(1,0), F.get(2,0));
	Vector3d tnf2( F.get(3,0), F.get(4,0), F.get(5,0));
	Vector3d tnf3( F.get(6,0), F.get(7,0), F.get(8,0));
	// name the points
	Vector3d P[6];
	for( long i = 0 ; i < 6 ; i ++)
	    P[i].set( nodes[p[i]]-> px,
		nodes[p[i]]-> py,
		nodes[p[i]]-> pz);
	// calculate normals of 3 walls
	//  Wall 1 defined by p[0],p[1],p[4],p[3]
	//  Wall 2 defined by p[1],p[2],p[5],p[4]
	//  Wall 3 defined by p[2],p[0],p[3],p[5]
	Vector3d w1n = cross_product(P[1]-P[0],P[3]-P[0]);
	Vector3d w2n = cross_product(P[2]-P[1],P[4]-P[1]);
	Vector3d w3n = cross_product(P[0]-P[2],P[5]-P[2]);
	w1n.normalize();
	w2n.normalize();
	w3n.normalize();
	fprintf( stderr, "w1n = %s\n", w1n.to_str());
	fprintf( stderr, "w2n = %s\n", w2n.to_str());
	fprintf( stderr, "w3n = %s\n", w3n.to_str());
	// calculate the total wall forces (wall1 = node1-2,
	// wall2 = node2-3, wall3 = node3,1)
	// twf1 = project( tnf1, wall_1_normal) +
	//        project( tnf2, wall_1_normal)
	// twf2 = project( tnf2, wall_2_normal) +
	//        project( tnf3, wall_2_normal)
	// twf3 = project( tnf3, wall_3_normal) +
	//        project( tnf1, wall_3_normal)
	Vector3d twf1 = project( tnf1, w1n) + project( tnf2, w1n);
	Vector3d twf2 = project( tnf2, w2n) + project( tnf3, w2n);
	Vector3d twf3 = project( tnf3, w3n) + project( tnf1, w3n);
	// get areas of walls
	double wa1 =
	    cross_product(P[1]-P[0],P[3]-P[0]).length()/2.0 +
	    cross_product(P[1]-P[3],P[4]-P[3]).length()/2.0;
	double wa2 =
	    cross_product(P[2]-P[1],P[4]-P[1]).length()/2.0 +
	    cross_product(P[2]-P[5],P[4]-P[5]).length()/2.0;
	double wa3 =
	    cross_product(P[0]-P[2],P[3]-P[2]).length()/2.0 +
	    cross_product(P[2]-P[5],P[3]-P[5]).length()/2.0;
	fprintf( stderr, "wa1 = %f\n", wa1);
	fprintf( stderr, "wa2 = %f\n", wa2);
	fprintf( stderr, "wa3 = %f\n", wa3);
	// divide total wall forces by the corresponding areas
	// of the wall, to get the force per unit areas
	twf1.scale( 1/wa1);
	twf2.scale( 1/wa2);
	twf3.scale( 1/wa3);
	// compute the results
	f1 = (twf1 + twf3); f1.scale( 0.5);
	f2 = (twf1 + twf2); f2.scale( 0.5);
	f3 = (twf2 + twf3); f3.scale( 0.5);
	fprintf( stderr, "f1=%s\n", f1.to_str());
	fprintf( stderr, "f2=%s\n", f2.to_str());
	fprintf( stderr, "f3=%s\n", f3.to_str());
    }

    void
	draw( DrawType draw_type, double ys_min, double ys_max,
	    int draw_deformed)
    {
	if( draw_type == WALLS_ALL || draw_type == WALLS_TOP) {
	    // calculate the color based on stresses
	    double ts = 0.0;
	    if( max_pstress > 0)
		ts = max_pstress / yield_stress;
	    if( ts > 1) ts = 1;
	    double cs = 0.0;
	    if( max_pstress < 0)
		cs = - max_pstress / yield_stress;
	    if( cs > 1) cs = 1;
	    ts = pow( ts, 2);
	    cs = pow( cs, 2);
			
	    // set the material
	    if( color_display == 0 && ! is_degenerate() &&
		! selected) {
		float amb[] = {0.3, 0.13, 0.0, 1.0 };
		float dif[] = {0.7, 0.3, 0.0, 1.0 };
		float spec[] = {0.02, 0.02, 0.02, 1.0 };
		float emis[] = {0.0, 0.0, 0.0, 1.0 };
		float shin[] = {0.0};
				
		glMaterialfv(GL_FRONT, GL_AMBIENT, amb);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, dif);
		glMaterialfv(GL_FRONT, GL_SPECULAR, spec);
		glMaterialfv(GL_FRONT, GL_EMISSION, emis);
		glMaterialfv(GL_FRONT, GL_SHININESS, shin);

	    } else if( color_display == 0 && selected) {
		float amb[] = {0.3, 0.13, 0.0, 1.0 };
		float dif[] = {0.7, 0.7, 0.0, 1.0 };
		float spec[] = {0.02, 0.02, 0.02, 1.0 };
		float emis[] = {0.0, 0.0, 0.0, 1.0 };
		float shin[] = {0.0};
				
		glMaterialfv(GL_FRONT, GL_AMBIENT, amb);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, dif);
		glMaterialfv(GL_FRONT, GL_SPECULAR, spec);
		glMaterialfv(GL_FRONT, GL_EMISSION, emis);
		glMaterialfv(GL_FRONT, GL_SHININESS, shin);

	    } else if( color_display == 0 && is_degenerate()) {
		float amb[] = {0.3, 0.13, 0.0, 1.0 };
		float dif[] = {0.3, 0.7, 0.0, 1.0 };
		float spec[] = {0.02, 0.02, 0.02, 1.0 };
		float emis[] = {0.0, 0.0, 0.0, 1.0 };
		float shin[] = {0.0};
				
		glMaterialfv(GL_FRONT, GL_AMBIENT, amb);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, dif);
		glMaterialfv(GL_FRONT, GL_SPECULAR, spec);
		glMaterialfv(GL_FRONT, GL_EMISSION, emis);
		glMaterialfv(GL_FRONT, GL_SHININESS, shin);

	    } else if( color_display == 1) {
		double no_stress_red = 0.4;
		double no_stress_green = 0.4;
		double no_stress_blue = 0.4;
		double full_tensile_red = 1.0;
		double full_tensile_green = 0.0;
		double full_tensile_blue = 0.0;
		double full_compressive_red = 0.0;
		double full_compressive_green = 0.0;
		double full_compressive_blue = 1.0;
		double full_both_red = 1.0;
		double full_both_green = 0;
		double full_both_blue = 0.0;
				
		double red =
		    (1-cs)*(1-ts)*no_stress_red
		    + (1-cs)*(  ts)*full_tensile_red
		    + (  cs)*(1-ts)*full_compressive_red
		    + (  cs)*(  ts)*full_both_red;
		double green =
		    (1-cs)*(1-ts)*no_stress_green
		    + (1-cs)*(  ts)*full_tensile_green
		    + (  cs)*(1-ts)*full_compressive_green
		    + (  cs)*(  ts)*full_both_green;
		double blue =
		    (1-cs)*(1-ts)*no_stress_blue
		    + (1-cs)*(  ts)*full_tensile_blue
		    + (  cs)*(1-ts)*full_compressive_blue
		    + (  cs)*(  ts)*full_both_blue;

		if( selected) {
		    red = 1.0; green = 1.0; blue = 0.0;
		}

		glDisable( GL_LIGHTING);
		glColor4f( red, green, blue, 1.0);
	    } else if( color_display == 2) {
		double red = (yield_stress - ys_min)
		    /(ys_max - ys_min);
		double green = red;
		double blue = red;
				
		float amb[] = { red * 0.3, green * 0.3,
				blue * 0.3, 1.0 };
		float dif[] = { red * 0.7, green * 0.7,
				blue * 0.7, 1.0 };
		float spec[] = {0.02, 0.02, 0.02, 1.0 };
		float emis[] = {0.0, 0.0, 0.0, 1.0 };
		float shin[] = {0.0};

		glMaterialfv(GL_FRONT, GL_AMBIENT, amb);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, dif);
		glMaterialfv(GL_FRONT, GL_SPECULAR, spec);
		glMaterialfv(GL_FRONT, GL_EMISSION, emis);
		glMaterialfv(GL_FRONT, GL_SHININESS, shin);

		glDisable( GL_LIGHTING);
		glColor4f( red, green, blue, 1.0);
	    }
				
	    long triangles[24] = {
		0, 2, 1,
		3, 5, 4,
		3, 4, 1,
		3, 1, 0,
		4, 5, 2,
		4, 2, 1,
		5, 3, 0,
		5, 0, 2};
	    glBegin( GL_TRIANGLES);
	    for( long i = 0 ; i < 8 ; i ++) {
		if( ! selected &&
		    (draw_type == WALLS_TOP && i > 0))
		    continue;
		Node * n0 = nodes[p[triangles[i*3+0]]];
		Node * n1 = nodes[p[triangles[i*3+1]]];
		Node * n2 = nodes[p[triangles[i*3+2]]];
		if( draw_deformed) {
		    draw_triangle_normal(
			n0-> px, n0-> py, n0-> pz,
			n1-> px, n1-> py, n1-> pz,
			n2-> px, n2-> py, n2-> pz);

		    // draw the triangle
		    glVertex3d( n0-> px, n0-> py, n0-> pz);
		    glVertex3d( n1-> px, n1-> py, n1-> pz);
		    glVertex3d( n2-> px, n2-> py, n2-> pz);
		} else {
		    draw_triangle_normal(
			n0-> ox, n0-> oy, n0-> oz,
			n1-> ox, n1-> oy, n1-> oz,
			n2-> ox, n2-> oy, n2-> oz);
		    // draw the triangle
		    glVertex3d( n0-> ox, n0-> oy, n0-> oz);
		    glVertex3d( n1-> ox, n1-> oy, n1-> oz);
		    glVertex3d( n2-> ox, n2-> oy, n2-> oz);
		}
	    }
	    glEnd();
	}
	else if( draw_type == WIREFRAME_ALL ||
	    draw_type == WIREFRAME_TOP) {

	    glBegin( GL_LINES);
	    long pairs[18] = {
		0, 1, 1, 2, 2, 0, 3, 4, 4, 5, 5, 3,
		0, 3, 1, 4, 2, 5 };
	    for( long i = 0 ; i < 9 ; i ++) {
		if( draw_type == WIREFRAME_TOP && i > 2)
		    continue;
		Node * n1 = nodes[p[pairs[i*2]]];
		Node * n2 = nodes[p[pairs[i*2+1]]];
		if( draw_deformed) {
		    glVertex3d( n1-> px, n1-> py, n1-> pz);
		    glVertex3d( n2-> px, n2-> py, n2-> pz);
		} else {
		    glVertex3d( n1-> ox, n1-> oy, n1-> oz);
		    glVertex3d( n2-> ox, n2-> oy, n2-> oz);
		}
	    }
	    glEnd();
	}
    }

    // 	int shares_edge_with( WedgeElement & e) {
    // 		long in_common_count = 0;
    // 		for( long i = 0 ; i < 3 ; i ++)
    // 			for( long j = 0 ; j < 3 ; j ++)
    // 				if( p[i] == e.p[j]) in_common_count ++;
    // 		if( in_common_count > 1) return 1; else return 0;
    // 	}

    long get_num_shared_nodes( const WedgeElement * e) const {
	long count = 0;
	for( long i = 0 ; i < 6 ; i ++) {
	    for( long j = 0 ; j < 6 ; j ++) {
		if( p[i] != e-> p[j]) continue;
		count ++;
		break;
	    }
	}
	return count;
    }

    void iso_to_cartesian_orig( double e, double n, double s,
	double & x, double & y, double & z) {
	// convert isoparametric coordinate e, n, s to the cartesian
	// coordinate (using the original xyz)
	double w[6];
	w[0] = Na(e,n,s);
	w[1] = Nb(e,n,s);
	w[2] = Nc(e,n,s);
	w[3] = Nd(e,n,s);
	w[4] = Ne(e,n,s);
	w[5] = Nf(e,n,s);
	x = y = z = 0;
	for( long i = 0 ; i < 6 ; i ++) {
	    x += w[i] * nodes[p[i]]-> ox;
	    y += w[i] * nodes[p[i]]-> oy;
	    z += w[i] * nodes[p[i]]-> oz;
	}
    }

    int is_degenerate( void) {
	// returns whether this element is degenerate
	return is_degenerate_flag;
    }

    int is_pretty( void) {
	// returns whether this element is pretty or not
	return is_pretty_flag;
    }

    double plane_intersection( long ind, Vector3d pn)
    // -----------------------------------------------------------------
    // Calculates the area of the intersection of the plane given by
    // plane normal [pn] and node[ind], with the element.
    // The intersection area is calculated for the undeformed element,
    // but the coordinates of the actual intersection are calculated for
    // the deformed element. Also, the true plane intersection is only
    // calculated for the top face, then it is assumed the plane is
    // perpendicular to the top face (i.e. the intersection with the
    // bottom plane is calculated using an interpolation. It is assumed
    // 'ind' is one of p[0], p[1] or p[3].
    // -----------------------------------------------------------------
    {
	// calculate n1=ind, n2 and n3 are the other two top face
	// nodes
	Node * n1 = nodes[ind];
	Node * n2 = NULL;
	Node * n3 = NULL;
	Node * b1 = NULL;
	Node * b2 = NULL;
	Node * b3 = NULL;
	if( p[0] == ind) {
	    n2 = nodes[p[1]];
	    n3 = nodes[p[2]];
	    b1 = nodes[p[0+3]];
	    b2 = nodes[p[1+3]];
	    b3 = nodes[p[2+3]];
	} else if( p[1] == ind) {
	    n2 = nodes[p[2]];
	    n3 = nodes[p[0]];
	    b1 = nodes[p[1+3]];
	    b2 = nodes[p[2+3]];
	    b3 = nodes[p[0+3]];
	} else if( p[2] == ind) {
	    n2 = nodes[p[0]];
	    n3 = nodes[p[1]];
	    b1 = nodes[p[2+3]];
	    b2 = nodes[p[0+3]];
	    b3 = nodes[p[1+3]];
	} else {
	    // ind is not a top node - return no intersection
	    return 0;
	}
	// name the deformed coordinates
	double n1x = n1-> px;
	double n1y = n1-> py;
	double n1z = n1-> pz;
	double n2x = n2-> px;
	double n2y = n2-> py;
	double n2z = n2-> pz;
	double n3x = n3-> px;
	double n3y = n3-> py;
	double n3z = n3-> pz;
	// calculate vector from n2 to n3
	double vx = n3x - n2x;
	double vy = n3y - n2y;
	double vz = n3z - n2z;
	// calculate the equation of the plane, so that it is given by
	// a*x + b*y + c*z + d = 0
	double a = pn.x;
	double b = pn.y;
	double c = pn.z;
	double d = - (a*n1x + b*n1y + c*n1z);
	// the line from n2 to n3 is a parametric line
	// x[t] = vx*t + nx2
	// y[t] = vy*t + ny2
	// z[t] = vz*t + nz2
	// calculate t for which the line intersects the plane
	double den = (a*vx + b*vy + c*vz);
	// if vector n2-n3 is parallel to the plane, report 0
	if( fabs( den) < 1e-10) return 0.0;
	double t = -((d + a*n2x + b*n2y + c*n2z)/den);
	// if t is outside of interval 0,1, report 0 intersection
	if( t < 0 || t > 1) return 0.0;
	//		fprintf( stderr, "t=%f\n", t);
	// find the intersection points with the wedge element,
	// but in the undeformed coordinates
	//
	// p1 = n1
	Vector3d p1( n1-> ox, n1-> oy, n1-> oz);
	// p2 = intersection of n2 and n3 (top face)
	Vector3d p2( linterpolate( n2-> ox, n3-> ox, t),
	    linterpolate( n2-> oy, n3-> oy, t),
	    linterpolate( n2-> oz, n3-> oz, t));
	// p3 = intersection of b2 and b3 (bottom face)
	Vector3d p3( linterpolate( b2-> ox, b3-> ox, t),
	    linterpolate( b2-> oy, b3-> oy, t),
	    linterpolate( b2-> oz, b3-> oz, t));
	// p4 = b1
	Vector3d p4( b1-> ox, b1-> oy, b1-> oz);
	//		fprintf( stderr, "P1=%s\n", p1.to_str());
	//		fprintf( stderr, "P2=%s\n", p2.to_str());
	//		fprintf( stderr, "P3=%s\n", p3.to_str());
	//		fprintf( stderr, "P4=%s\n", p4.to_str());
	// calculate the area of intersection = the sum of triangle
	// areas (p1,p2,p3) and (p1,p3,p4)
	double a1 = cross_product(p2-p1,p3-p1).length()/2;
	double a2 = cross_product(p3-p1,p4-p1).length()/2;
	//		fprintf( stderr, "a1 = %f , a2 = %f\n", a1, a2);
	double area = a1 + a2;
	return area;
    }
};

/*
  template <class T> void swap( T & a, T & b) {
  T tmp = a;
  a = b;
  b = tmp;
  }
*/

class Edge {
 public:	long n1, n2;
    Edge( long pn1, long pn2) {
	assert( pn1 != pn2);
	n1 = pn1;
	n2 = pn2;
    }
    Edge(void) {
	n1 = -1; n2 = -1;
    }
    Edge( const Edge & se) {
	n1 = se.n1;
	n2 = se.n2;
    }
    int operator == ( const Edge & e) const {
	if( n1 == e.n1 && n2 == e.n2) return 1;
	if( n2 == e.n1 && n1 == e.n2) return 1;
	return 0;
    }
    double length( void) const {
	return sqrt( sqr(nodes[n1]->ox-nodes[n2]->ox) +
	    sqr(nodes[n1]->oy-nodes[n2]->oy) +
	    sqr(nodes[n1]->oz-nodes[n2]->oz));
    }
};

// static variables of WedgeElement
WedgeElement::GP * WedgeElement::gps = NULL;
long WedgeElement::n_gauss = -1;
long WedgeElement::n_gauss_surface = -1;
char WedgeElement::cache_Ke = 0;
char WedgeElement::fast_Ke_removal = 0;

class Model {

 public:
    // simulation constants:
    // -------------------------------------------------------

    double time_total;	// total simulation time
    double time_curr;	// current simulation time
    double min_time_step;   // minimum time step allowed
    double max_time_step;	// maximum time step allowed
    //	double t0;		// last time known when there are no fracture
    //				// candidates
    //	double t1;		// known time when there are more than one
    //				// fracture candidates
    //	double dt;		// current time step
    // variables used to determine next <time> when applying growth
    long min_dtc;
    long max_dtc;
    long dtc;
    double t1;
    long fc_count;

    double growth_x;	// X-growth per 1 time unit
    double growth_y;	// Y-growth per 1 time unit
    double growth_z;	// Z-growth per 1 time unit

    double shrink_top_t0;	// Top face shrink start time
    double shrink_top_t1;	// Top face shrink end time
    double shrink_top_val;	// Top face total end shrinkage
    double shrink_bot_t0;	// Bottom face shrink start time
    double shrink_bot_t1;	// Bottom face shrink end time
    double shrink_bot_val;	// Bottom face total end shrinkage
    double shrink_height_t0; // Time when height starts to change
    double shrink_height_t1; // time when height stops to change
    double shrink_height_val0; // height before t0
    double shrink_height_val1; // height after t1

    double fracture_tip_inertia; // inertia of fracture tip propagation
    // max. pstress at a fracture tip will be
    // multiplied by this coefficient

    double gravity_x;	// GRAVITY CONSTANTS
    double gravity_y;	//      || 
    double gravity_z;	//      ||

    double max_break_size;	// the maximum allowed size for an 
				// element
				// before it is allowed to be broken
    double min_refine_size;	// when automatic (adaptive) mesh refinement
				// is used, elements below this size are
				// not refined

    // what should be the element size around crack tips?
    double crack_tip_element_size;

    // what should be the element size around crack tip in the sub-mode?
    double crack_tip_sub_element_size;

    // the minimum element size (elements are not allowed to be split
    // any further than this value)
    double min_element_size;
	
    // max. element size - this is used at the beggining, when a model
    // is read in - all elements below this value are subdivided, until
    // none is larger than this value
    double max_element_size;

    // size of the plastic zone
    double pzone_size;

    double precision;	// precision to be used in numerical solutions
    double error_radius;	// the radius of error for local discretization
    double min_error_radius; // minimum error radius

    // material property - randomization mappings
    Map * ym_map , * pr_map, * ys_map, * ft_map;

    // refinement map
    Map * refinement_map;

    // name of the input file
    char * input_fname;

    SparseMatrix * Kg;	// model's global Kg
    Vector_double * RHSg;	// model's global right hand side
    SparseMatrix * KgC;	// model's global Kg + global constraints
    Vector_double * RHSgC;	// model's RightHandSide + global constraints
				// KgC and RHSgC are used to pass the values
				// between relax_local() & relax_global()
    double relax_local_err;	// what error did relax_local() achieve?
				// (relax_global() will decide based on this
				//  error whether to do anything)
    long n_fractures;	// count the number of fractures here
    Fifo<long> ftips;	// here we store fracture tips (when nodal
    // breaking is used)

    long get_number_of_fracture_candidates( void);
    double get_mat_thickness_at_node( const long n0) const;
    int edge_exists( long n1, long n2);
    int face_edge_exists( long n1, long n2);
    void grow( double time);
    int refine_neighborhood( WedgeElement * e, double max_size);
    void recalculate_Kg( void);
    void insert_constraints( SparseMatrix & k, Vector_double & Sol);
    void insert_constraints_local( SparseMatrix & k, Vector_double & Sol);
    long add_node( Node * n);
	
    long add_element( WedgeElement * e);
    void remove_element( WedgeElement * e);
    int delaunay_edge_flip( long tn1, long tn2);
    void insert_sink( long ind);
    void subdivide_element_bisection( long ind);
    long find_background_node( long ind);
    Vector3d get_surface_normal( long ind);
    Vector3d get_face_normal( long n);
	
    typedef Stack<long> NodeList;
    NodeList * get_face_neighbors_circle( long n0);
    NodeList * get_face_neighbors( long n0);
    NodeList get_node_column( long ind);

    typedef Stack<WedgeElement *> ElementList;
    ElementList get_elements_sharing_edge( Edge & e);

    // progressive save information
    char progressive_save_on;
    char * progressive_save_fmask;
    double progressive_save_skip;
    double progressive_save_last;
    void progressive_save( void);
    void auto_save( void);
    int split_single_point_node( long ind);
    long split_edge( long n1, long n2, double r);
    int refine_mesh( void);
    int _refine_mesh_von_mises( void);
    int _refine_mesh_max_stress( void);
    void interpolate_growth_vector( long n1, long n2, long n3, double t);

 public:
    double get_Kg_error( void);
    int check_ref_consistency( void);

    enum BreakMethod {
	BreakElement,
	BreakNode,
	BreakEdge
    };

    Model() {
	Kg = new SparseMatrix( 1, 1);
	RHSg = new Vector_double( 1);
	KgC = new SparseMatrix( 1, 1);
	RHSgC = new Vector_double( 1);
	relax_local_err = 1e100;
	ym_map = pr_map = ys_map = ft_map = NULL;
	refinement_map = NULL;
	input_fname = NULL;
	progressive_save_on = 0;
	progressive_save_fmask = NULL;
	progressive_save_skip = 0.0;
	progressive_save_last = time_curr;
	break_method = BreakElement;
	n_fractures = 0;
    }
	
    ~Model() {
	if( Kg != NULL) delete Kg;
	if( RHSg != NULL) delete RHSg;
	if( KgC != NULL) delete KgC;
	if( RHSgC != NULL) delete RHSgC;
	if( ym_map != NULL) delete ym_map;
	if( pr_map != NULL) delete pr_map;
	if( ys_map != NULL) delete ys_map;
	if( ft_map != NULL) delete ft_map;
	if( refinement_map != NULL) delete refinement_map;
	if( input_fname != NULL) free( input_fname);
	if( progressive_save_fmask != NULL)
	    free( progressive_save_fmask);
    }

    void draw( void);
    void draw_orig_shapes( void);
    void load( const char * fname);
    int simulation_step( void);
    void save( const char * fname);
    void _save_header( FILE * fp);
    void fracture_extend( void);
    long fracture_new( int count_mode = 0);
    void orient_elements( void);
    void subdivide_element( long ind);
    void subdivide_marked_elements( void);
    void select_nodes( long n_it);
    long select_elements( long n, long n_it);
    long mark_local_change_nodes( double init_value,
	int skip_fixed_nodes = 1);
    void mark_plastic_zone( long ind);
    long flag_nodes_with_error( void);
    void relax( int forced = 0);
    void relax_local( void);
    void _relax_local( void);
    void relax_global( int forced = 0);
    void _relax_global( int forced = 0);
    void calculate_stresses( void);
    void _calculate_elemental_stresses( void);
    void _calculate_nodal_stresses( void);
    void _calculate_nodal_stresses_tensor_weighing( void);
    void _calculate_nodal_stresses_tensor_LSQF_interpolation( void);
    void _calculate_nodal_stresses_JH( void);
    void calculate_precise_stress_at_node( long n0);
    double calculate_potential_energy( void);
    //	void calculate_nodal_forces( void);
    void set_break_method( BreakMethod m) {
	break_method = m;
    }
    void print_stresses( FILE * fp);
    int is_node_split_valid( long tn0);
    long get_element_index( WedgeElement * e);
    int is_node_surrounded( long ind);
    int is_node_moveable( long ind);
    int is_node_moveable_raw( long ind);
    void collapse_edge( long tn1, long tn2);
    long get_neighbor_count( WedgeElement * e);
    int prettify_mesh( void);
    int angle_smooth_mesh( void);
    int repell_smooth_mesh( void);
    int smooth_mesh_around_crack_tips( void);
    int smooth_mesh_around_nodes( const NodeList & nlist, long depth);
    int refine_mesh_around_node( long ind, long depth, double req_size);
    //	int is_background_node( long ind);
    int peel( double peel_stress);

 private:
 public:
    BreakMethod break_method;

    int fracture_element( int count_mode = 0);
    int _fracture_at_node( long n0,
	NodeList & new_ftips,
	double & farea);
    int fracture_edge( void);

};

class GlobalVariables {
 public:
    char * fname;
};
GlobalVariables glob;

Model * model;
static char * fname = NULL;
static string progname;
static string viewer_progname;
static int profile_on = 0;
static int calculate_nodal_stress_flag = 0;
static long n_elements = 0;
static WedgeElement ** elements = NULL;
static enum { BreakElementMethod,
	      BreakNodeMethod,
	      BreakEdgeMethod
} b_method = BreakNodeMethod;


static void parse_command_line( int argc, char ** argv)
// ----------------------------------------------------------------------
// parses command line arguments
// ----------------------------------------------------------------------
{
    // remember how to invoke ourselves
    progname = argv[0];
    // figure out how to invoke viewer
    {
	viewer_progname = progname;
	string::size_type pos = viewer_progname.rfind( '/');
	if( pos == string::npos) pos = 0; else pos++;
	viewer_progname.erase( pos);
	viewer_progname.append( "viewer");
    }

    fname = NULL;
    long arg = 1;
    while( arg < argc) {
	if( strcmp( argv[arg], "-p") == 0 && profile_on == 0) {
	    profile_on = 1;
	    arg ++;	continue;
	} else if( strcmp( argv[arg], "-calculate_nodal_stress") == 0){
	    arg ++;
	    calculate_nodal_stress_flag = 1;
	    continue;
	} else if( strcmp( argv[arg], "-break") == 0) {
	    arg ++;
	    if( arg >= argc) usage( argv[0]);
	    if( strcmp( argv[ arg], "Element") == 0)
		b_method = BreakElementMethod;
	    else if( strcmp( argv[ arg], "Node") == 0)
		b_method = BreakNodeMethod;
	    else if( strcmp( argv[ arg], "Edge") == 0)
		b_method = BreakEdgeMethod;
	    else
		usage( argv[0]);
	    arg ++; continue;
	} else if( strcmp( argv[arg], "-memopt") == 0) {
	    optimize_memory_flag = 1;
	    arg ++; continue;
	} else if( strcmp( argv[arg], "-fast_ke_removal") == 0) {
	    fast_ke_removal_flag = 1;
	    arg ++; continue;
	} else if( strcmp( argv[arg], "-no_local_relax") == 0) {
	    use_local_relaxation_flag = 0;
	    arg ++; continue;
	} else if( strcmp( argv[arg], "-anim") == 0) {
	    arg ++;
	    if( arg >= argc) usage( argv[0]);
	    animation_fname_mask = strdup( argv[arg]);
	    arg ++;
	    if( arg >= argc) usage( argv[0]);
	    if( 1 != sscanf( argv[ arg], "%ld",
		    & animation_count))
		usage( argv[0]);
	    arg ++;
	    continue;
	} else if( strcmp( argv[arg], "-autosave") == 0) {
	    arg ++;
	    if( arg >= argc) usage( argv[0]);
	    autosave_fname_mask = strdup( argv[arg]);
	    arg ++;
	    if( arg >= argc) usage( argv[0]);
	    if( 1 != sscanf( argv[ arg], "%ld",
		    & autosave_interval)) usage( argv[0]);
	    arg ++;
	    if( arg >= argc) usage( argv[0]);
	    if( 1 != sscanf( argv[ arg], "%ld",
		    & autosave_n_keep)) usage( argv[0]);
	    arg ++; continue;
	} else if (strcmp (argv[arg], "-no-precise-nodal-stress") == 0) {
	    precise_nodal_stress_flag = false;
	    arg ++; continue;
	} else if( fname == NULL) {
	    fname = argv[arg];
	    arg ++; continue;
	}
	usage( argv[0]);
    }
    if( fname == NULL) usage( argv[0]);

    fprintf( stderr, "Parameters:\n");
    fprintf( stderr, "   fname = %s\n", fname);
    fprintf( stderr, "   profile = %s\n", profile_on ? "YES" : "NO");
    fprintf( stderr, "   break method = ");
    if( b_method == BreakElementMethod)
	fprintf( stderr, "Element");
    else if( b_method == BreakNodeMethod)
	fprintf( stderr, "Node");
    else if( b_method == BreakEdgeMethod)
	fprintf( stderr, "Edge");
    fprintf( stderr, "\n");
	

    fprintf( stderr, "\n");
    fprintf( stderr, "\n");
}

double Model::get_mat_thickness_at_node( const long n0) const
// determines the height of the material at this node
{
    assert( n0 >= 0 || n0 < n_nodes);

    if( nodes[n0]-> n_refs <= 0) return -1;

    double h_sum = 0;
    for( long i = 0 ; i < nodes[n0]-> n_refs ; i ++)
	h_sum += nodes[n0]-> refs[i]-> get_height_at_node( n0);
    return h_sum / nodes[n0]-> n_refs;
}

long Model::get_number_of_fracture_candidates( void)
// ----------------------------------------------------------------------
// return the number of places a fracture could occur
// (depends on the method of breaking)
// ----------------------------------------------------------------------
{
    if( break_method == BreakElement)
	return fracture_element( 1);
    else if( break_method == BreakNode)
	return fracture_new( 1);

    fprintf( stderr, "TBD: candidate count not implemented!!!\n");
    assert( 0);
    return 0;
}

int Model::edge_exists( long n1, long n2)
// ----------------------------------------------------------------------
// returns: 1 - if edge n1,n2 can be found in the model
//          0 - otherwise
// ----------------------------------------------------------------------
{
    assert( n1 >= 0 && n1 < n_nodes);
    assert( n2 >= 0 && n2 < n_nodes);
    for( long i = 0 ; i < nodes[n1]-> n_refs ; i ++) {
	WedgeElement & e = * nodes[n1]-> refs[i];
	if( e.get_internal_index(n2) >= 0)
	    return 1;
    }
    return 0;
}

int Model::face_edge_exists( long n1, long n2)
// ----------------------------------------------------------------------
// returns: 1 - if edge n1,n2 can be found in the model on some element
//              either on its top or bottom face
//          0 - otherwise
// ----------------------------------------------------------------------
{
    assert( n1 >= 0 && n1 < n_nodes);
    assert( n2 >= 0 && n2 < n_nodes);
    for( long i = 0 ; i < nodes[n1]-> n_refs ; i ++) {
	WedgeElement & e = * nodes[n1]-> refs[i];
	long in1 = e.get_internal_index(n1);
	long in2 = e.get_internal_index(n2);
	if( in1 == in2) continue;
	if( in2 == -1) continue;
	if( in1 < 3 && in2 < 3) return 1;
	if( in1 >= 3 && in2 >= 3) return 1;
    }
    return 0;
}

void Model::grow( double time)
// ----------------------------------------------------------------------
// Applies changes to the model to bring its shape to time
//    - adjusts positions of nodes
//    - adjusts shapes of elements
//    - recalculates the global stiffness matrix Kg and RHSg
// ----------------------------------------------------------------------
{
    fprintf( stderr, "\t- growing for time %f\n", time);
    // =============================================================
    // adjust node positions
    // =============================================================
    // use surface normals to determine growth
    for( long i = 0 ; i < n_nodes ; i ++)
	nodes[i]-> grow( time, growth_x, growth_y, growth_z);
    /*
      for( long i = 0 ; i < n_nodes ; i ++) {
      Node * n = nodes[i];
      // only apply growth to the fixed nodes
      if( ! n-> is_fixed) continue;
      n-> px = n-> ox + n-> gv.x * growth_x * time;
      n-> py = n-> oy + n-> gv.y * growth_y * time;
      n-> pz = n-> oz + n-> gv.z * growth_z * time;
      assert( ! isnan( n-> px));
      assert( ! isnan( n-> py));
      assert( ! isnan( n-> pz));
      }
    */

    // apply growth to individual elements
    for( long i = 0 ; i < n_elements ; i ++)
	elements[i]-> grow( time,
	    shrink_top_t0,
	    shrink_top_t1,
	    shrink_top_val,
	    shrink_bot_t0,
	    shrink_bot_t1,
	    shrink_bot_val,
	    shrink_height_t0,
	    shrink_height_val0,
	    shrink_height_t1,
	    shrink_height_val1
			    );

    // recalculate Kg and RHSg
    recalculate_Kg();
}

long Model::get_element_index( WedgeElement * e)
// ---------------------------------------------------------------------------
// given an element pointer, find its index in the elements[] array
// ---------------------------------------------------------------------------
{
    for( long i = 0 ; i < n_elements ; i ++)
	if( elements[i] == e) return i;
    return -1;
}

void Model::subdivide_marked_elements( void)
{
    while( 1) {
	int go_again = 0;
	for( long i = 0 ; i < n_elements ; i ++) {
	    WedgeElement * e = elements[i];
	    if( e-> mark) {
		subdivide_element( i);
		go_again = 1;
		continue;
	    }
	}
	if( ! go_again) break;
    }
}

void Model::select_nodes( long n_it)
// ----------------------------------------------------------------------
// for a number of iterations requested:
//    - select all nodes in all elements which have at least one node
//      selected
// ----------------------------------------------------------------------
{
    for( long it = 0 ; it < n_it ; it ++) {
	// make a list of all selected nodes
	NodeList lst;
	for( long i = 0 ; i < n_nodes ; i ++)
	    if( nodes[i]-> selected)
		lst.push(i);
	// now select all nodes in elements that touch nodes
	// on the list
	while( ! lst.is_empty()) {
	    Node & n = * nodes[ lst.pop()];
	    for( long i = 0 ; i < n.n_refs ; i ++)
		for( long j = 0 ; j < 6 ; j ++) {
		    long ind = n.refs[i]->p[j];
		    if( nodes[ind]-> is_fixed) continue;
		    nodes[ ind]->selected = 1;
		}
	}
    }
}

long Model::select_elements( long ind, long n_it)
// ----------------------------------------------------------------------
// for n_it = 0 - select no elements
// for n_it = 1 - select all elements touching node <ind>
// for n_it = 2 - select all elements touching <n0> and also select
// selects elements that touch
// for a number of iterations requested:
//    - 
// ----------------------------------------------------------------------
{
    assert( 0); // not implemented
    return 0;
}

long
    Model::mark_local_change_nodes( double init_rad, int skip_fixed_nodes)
// ----------------------------------------------------------------------
// Input: some nodes are selected
// Output: all nodes withing 'init_rad' of the originally selected
//         nodes will be selected
// ----------------------------------------------------------------------
{
    // set an error value for each node
    double err[ n_nodes];
    for( long i = 0 ; i < n_nodes ; i ++) 
	if( nodes[i]-> selected)
	    err[i] = init_rad;
	else
	    err[i] = 0;
	
    // initialize stack with nodes that have error > 0
    Stack<long> stack;
    // also, prepare an array of nodes, so that we can tell very
    // quickly whether a node is in the stack or not, to avoid
    // duplicates
    char flag[n_nodes];
    for( long i = 0 ; i < n_nodes ; i ++) {
	if( err[i] > 0) {
	    stack.push(i);
	    flag[i] = 1;
	} else {
	    flag[i] = 0;
	}
    }
    // while stack is not empty, pick a node from the stack that has
    // the max. error associated with it, and try to diffuse its error
    // to the around space, putting into stack all nodes affected
    while( ! stack.is_empty()) {
	long ind = -1;
	if( 0) {
	    ind = stack.pop_max();
	} else {
	    ind = stack.array[0];
	    for( long i = 1 ; i < stack.n ; i ++)
		if( err[stack.array[i]] > err[ind])
		    ind = stack.array[i];
	    stack.remove_val(ind);
	}
	assert( ind >= 0 && ind < n_nodes);
	flag[ind] = 0;
	Node * n1 = nodes[ind];
	// skip fixed nodes
	if( skip_fixed_nodes && n1->is_fixed) continue;
	for( long j = 0 ; j < n1-> n_refs ; j ++) {
	    WedgeElement & e = *(WedgeElement*)n1->refs[j];
	    for( long k = 0 ; k < 6 ; k ++) {
		Node * n2 = nodes[e.p[k]];
		// do not propagate into fixed nodes
		// and onto yourself
		if( skip_fixed_nodes && n2->is_fixed) continue;
		if( n1 == n2) continue;
		// calculate length between n1 and n2
		double d = Node::dist(n1,n2);
		// calculate the new value at n2
		double val2 = err[ind] - d;
		if( val2 < 0) val2 = 0;
		// is the node n2 affected?
		if( val2 <= err[e.p[k]]) continue;
		// yes it is:
		err[e.p[k]] = val2;
		if( ! flag[e.p[k]]) {
		    stack.push(e.p[k]);
		    flag[e.p[k]] = 1;
		}
	    }
	}
    }

    // mark all nodes whose error is > 0
    long n_marked = 0;
    for( long i = 0 ; i < n_nodes ; i ++){
	nodes[i]-> selected = 0;
	if( nodes[i]-> n_refs == 0) continue;
	if( skip_fixed_nodes && nodes[i]-> is_fixed) continue;
	if( err[i] > 0)
	{
	    nodes[i]-> selected = 1;
	    n_marked ++;
	}
    }

    return n_marked;
}

void Model::mark_plastic_zone( long ind)
// ----------------------------------------------------------------------
// all nodes connected to <ind> through top faces which have yield stress
// exceeded are marked as plastic zone nodes
// ----------------------------------------------------------------------
{
    assert( ind >= 0 && ind < n_nodes);

    // mark all nodes in the plastic zone
    for( long i = 0 ; i < n_nodes ; i ++)
	nodes[i]-> selected = 0;
    nodes[ ind]-> selected = 1;
    mark_local_change_nodes( pzone_size);

    // set all selected nodes as plastic nodes
    for( long i = 0 ; i < n_nodes ; i ++)
	nodes[i]-> is_in_pzone = nodes[i]-> selected ||
	    nodes[i]-> is_in_pzone;

    return;

    Node * n = nodes[ind];
    // if this node does not have yield stress exceeded, we are done
    if( n-> ps_val < n-> yield_stress) return;
    // if this node has been processed already, we are done
    if( n-> is_in_pzone) return;
    // this node is in plastic zone
    n-> is_in_pzone = 1;
    // get a list of all nodes connected to node <ind>
    NodeList * lst = get_face_neighbors( ind);
    // process each neighbor
    while( ! lst-> is_empty())
	mark_plastic_zone( lst-> pop());
    delete lst;
}

#ifdef DONT_COMPILE
int Model::split_single_point_node( long ind)
// ---------------------------------------------------------------------------
// checks to see if at this there is a triangle that touches other material
// only at the node, if so, the triangle is given its own copy of this node
//
// Returns: 1 = if something was split
//          0 = if nothing was split
// ---------------------------------------------------------------------------
{
    return 0;

    assert( 0 <= ind && ind < n_nodes);

    Node & n = * nodes[ ind];

    assert( n.is_fixed == 0);

    // check every element touching node n
    for( long i = 0 ; i < n.n_refs ; i ++) {
	WedgeElement & e = * n.refs[i];
	// if element e neighbors with any other element connected
	// to the node we are checking, it is left alone
	int has_edge_neighbor = 0;
	for( long j = 0 ; j < n.n_refs ; j ++) {
	    if( i == j) continue;
	    WedgeElement & e2 = * n.refs[j];
	    if( e.shares_edge_with( e2)) {
		has_edge_neighbor = 1;
		break;
	    }
	}
	if( has_edge_neighbor) continue;

	fprintf( stderr, "\t- splitting node %ld\n", ind);
	// this element does not have an edge neighbor, it
	// is therefore separated from node n

	// duplicate node n
	long cn0 = add_node( new Node( n.ox, n.oy, n.oz,
		n.px, n.py, n.pz,
		0));
		
	// find out n1 and n2
	long tn0, bn0, tn1, bn1, tn2, bn2;
	if( e.p[0] == ind) {
	    tn0 = e.p[0]; tn1 = e.p[1]; tn2 = e.p[2];
	    bn0 = e.p[3]; bn1 = e.p[4]; bn2 = e.p[5];
	} else if( e.p[1] == ind) {
	    tn0 = e.p[1]; tn1 = e.p[2]; tn2 = e.p[0];
	    bn0 = e.p[4]; bn1 = e.p[5]; bn2 = e.p[3];
	} else {
	    tn0 = e.p[2]; tn1 = e.p[0]; tn2 = e.p[1];
	    bn0 = e.p[5]; bn1 = e.p[3]; bn2 = e.p[4];
	}

	remove_element( & e);
	delete & e;
	add_element( 
	    new WedgeElement(  cn0, tn1, tn2,
		bn0, bn1, bn2,
		ym_map, pr_map, ys_map, ft_map));
		
	// recursively finish off this node
	split_single_point_node( ind);
	return 1;
    }

    // nothing was split:
    return 0;
}
#endif

int Model::split_single_point_node( long ind)
// ---------------------------------------------------------------------------
// splits single point contact nodes
//
// Method:
//   - create an initial list of all elements connected to node <ind>
//   - while list not empty do
//       - extract top element from the list into list2
//       - repeat
//           - find an element in list that is connected to anything in
//             list 2
//           - if such element exists, put it into list2 and continue
//           - otherwise break
//       - if list is empty (we have extracted the last group of elements
//         so we are done), break
//       - otherwise, create a copy of ind and assign it to all elements
//         in list 2
//   - return
//
// Returns: 1 = if something was split
//          0 = if nothing was split
// ---------------------------------------------------------------------------
{
    assert( 0 <= ind && ind < n_nodes);

    Node & n = * nodes[ ind];

    // create a list of all elements connected to node <ind>
    Stack<WedgeElement *> list1;
    for( long i = 0 ; i < n.n_refs ; i ++)
	list1.push( n.refs[i]);

    int res = 0;

    // process all groups in list1:
    while( ! list1.is_empty()) {
	Stack<WedgeElement *> list2;
	// extract the first element
	list2.push( list1.pop());
	// now extract the group
	while( 1) {
	    int something_extracted = 0;
	    // check all elements in list1
	    for( long i = 0 ; i < list1.n ; i ++) {
		WedgeElement * e = list1(i);
		for( long j = 0 ; j < list2.n ; j ++) {
		    WedgeElement * e2 = list2(j);
		    long k = e-> get_num_shared_nodes( e2);
		    if( k < 4) continue;
		    list1.remove_val( e);
		    list2.push( e);
		    something_extracted = 1;
		    e = NULL;
		    break;
		}
		if( e == NULL) break;
	    }
	    if( ! something_extracted) break;
	}

	// list2 now contains another group of elements. If list1 is
	// empty, we are done
	if( list1.is_empty()) break;
		
	// otherwise make a copy of ind, and assign this copy to
	// all elements in list2
	res = 1;
	long cind = add_node( new Node( n));
	while ( ! list2.is_empty()) {
	    WedgeElement * e = list2.pop();
	    long in0 = e-> get_internal_index( ind);
	    assert( in0 >= 0 && in0 < 6);
	    long p[6];
	    for( long k = 0 ; k < 6 ; k ++) p[k] = e-> p[k];
	    p[in0]=cind;
	    remove_element( e); delete e;
	    add_element( new WedgeElement(
			     p[0], p[1], p[2],
			     p[3], p[4], p[5],
			     ym_map, pr_map, ys_map, ft_map));
	}
    }
	
    return res;
}

class SplitEdge {
 public:	long n1, n2, n3;
    SplitEdge( long pn1, long pn2, long pn3 = -1) {
	assert( pn1 != pn2);
	n1 = pn1;
	n2 = pn2;
	n3 = pn3;
    }
    SplitEdge(void) {
	n1 = -1; n2 = -1; n3 = -1;
    }
    SplitEdge( const SplitEdge & se) {
	n1 = se.n1;
	n2 = se.n2;
	n3 = se.n3;
    }
    int operator == ( const SplitEdge & e) const {
	if( n1 == e.n1 && n2 == e.n2) return 1;
	if( n2 == e.n1 && n1 == e.n2) return 1;
	return 0;
    }
};

long
    Model::split_edge( long n1, long n2, double r)
// ----------------------------------------------------------------------
// splits an edge by inserting a node, which means splitting all wedges
// that share this edge. Also, the multi-layer structure of the wedges
// is observed, i.e. if other edges have to be split, all of their wedges
// are also split
//   n1 = first node of the edge
//   n2 = second node of the edge
//    r = where to split the edge (0 = at node n1, 1 = at node n2)
// returns:
//    -1 = nothing was split (the edge does not exist)
//     n = the index of the new node inserted between <n1> and <n2>
// algorithm:
//     - initialize elist1 with n1, n2
//     - in elist1 we keep a list of edges that have to be split
//     - initialize elist2 to be empty
//     - in elist2 we keep a list of edges that have been split, including
//       the inserted node
//     - while edge list is not empty:
//         - take an edge out of elist1 and put it into n1,n2
//         - for each element e in which n1,n2 forms an edge either on top
//           or bottom face:
//             - find the oposite edge in the element to n1,n2 and put
//               it into n3, n4
//             - if edge n3,n4 is not in elist1, add it there
//             - subdivide n1,n2, creating n12 (but check elist2 first)
//             - subdivide n3,n4, creating n34 (but check elist2 first)
//             - remove element e
//             - add two new elements
// ----------------------------------------------------------------------
{
    //	fprintf( stderr, "split_edge(%ld,%ld,%f)\n", n1, n2, r);
    assert( n1 >= 0);
    assert( n1 < n_nodes);
    assert( n2 >= 0);
    assert( n2 < n_nodes);
    assert( r > 0);
    assert( r < 1);
    // initialize a list of edges to be processed, intially containing
    // only n1,n2
    Stack<SplitEdge> elist1;
    elist1.push( SplitEdge( n1, n2));
	
    // creat a list of split edges - initialy empty
    Stack<SplitEdge> elist2;

    // while the unprocessed edge list is not empty:
    while( ! elist1.is_empty()) {
	// extract an unprocessed edge
	SplitEdge tmp = elist1.pop();
	long n1 = tmp.n1;
	long n2 = tmp.n2;
	//		fprintf( stderr, "\t popped %ld %ld\n", n1, n2);
	for( long i = nodes[n1]-> n_refs-1 ; i >= 0 ; i --) {
	    WedgeElement * e = nodes[n1]->refs[i];
	    // get internal indices of n1 and n2
	    long in1 = e-> get_internal_index( n1);
	    assert( in1 >= 0);
	    long in2 = e-> get_internal_index( n2);
	    // if n2 is not part of element 'e', skip this element
	    if( in2 < 0) continue;
	    // figure out whether n1,n2 are on top face or
	    // bottom face, or two different faces. If on two
	    // different faces, this element is skipped,
	    // although it should never happen!!!
	    int n12_on_top;
	    if( in1 < 3 && in2 < 3) n12_on_top = 1;
	    else if( in1 >= 3 && in2 >= 3) n12_on_top = 0;
	    else {
		fprintf( stderr,
		    "--- Weird Al Yankowic ---\n");
		continue;
	    }
	    // if n1,n2 are on the bottom face, but both are
	    // fixed nodes, just don't do anything
	    if( nodes[n1]-> is_fixed &&
		nodes[n2]-> is_fixed &&
		! n12_on_top) continue;
	    // figure out whether n1,n2 are clockwise or
	    // counter-clockwise in the element
	    int n12_clockwise = 0;
	    if( n12_on_top) {
		if( (in1+1)%3 == in2)
		    n12_clockwise = 1;
	    } else {
		if( (in1-3+1)%3 == in2-3)
		    n12_clockwise = 1;
	    }
	    // figure out n3 and n4 - the mirrors of n1,n2 on
	    // the opposite face. Also determine on12 = the
	    // opposite on n1,n2 and on34 = the opposite of n3,
	    // n4
	    long n3 = -1;
	    long n4 = -1;
	    long on12 = -1;
	    long on34 = -1;
	    if( n12_on_top) {
		n3 = e-> p[in1+3];
		n4 = e-> p[in2+3];
		if( (in1 == 0 && in2 == 1) ||
		    (in1 == 1 && in2 == 0)) {
		    on12 = e->p[2];	on34 = e->p[2+3];
		} else if( (in1 == 1 && in2 == 2) ||
		    (in1 == 2 && in2 == 1)) {
		    on12 = e->p[0];	on34 = e->p[0+3];
		} else {
		    on12 = e->p[1];	on34 = e->p[1+3];
		}
	    } else {
		n3 = e-> p[in1-3];
		n4 = e-> p[in2-3];
		if( (in1 == 3 && in2 == 4) ||
		    (in1 == 4 && in2 == 3)) {
		    on12 = e->p[5];	on34 = e->p[5-3];
		} else if( (in1 == 4 && in2 == 5) ||
		    (in1 == 5 && in2 == 4)) {
		    on12 = e->p[3];	on34 = e->p[3-3];
		} else {
		    on12 = e->p[4];	on34 = e->p[4-3];
		}
	    }
	    // now add n3,n4 to the list of unprocessed edges, but
	    // only if that edge is not there yet
	    assert(n3 != n4);
	    elist1.push_unique( SplitEdge( n3, n4));
	    // 			fprintf( stderr,
	    // 				 "\tAdding unprocessed edge [%ld,%ld]\n",
	    // 				 n3, n4);
	    // check from elist2 whether n12 has been already
	    // created or not
	    long n12 = -1;
	    for( long j = 0 ; j < elist2.n ; j ++) {
		if( (elist2.array[j].n1 == n1 &&
			elist2.array[j].n2 == n2) ||
		    (elist2.array[j].n2 == n1 &&
			elist2.array[j].n1 == n2)) {
		    n12 = elist2.array[j].n3;
		    break;
		}
	    }
	    // if n12 has not yet been created, create it
	    // and put it into elist2
	    if( n12 == -1) {
		assert( n1 != n2);
		n12 = add_node(
		    new Node( nodes[n1], nodes[n2], r));
		elist2.push( SplitEdge( n1, n2, n12));
		//  				fprintf( stderr,
		//  					 "\tAdding split edge12 "
		//  					 "[%ld,%ld,%ld]\n",
		//  					 n1, n2, n12);
	    }
	    // check from elist2 whether n34 has been already
	    // created or not
	    long n34 = -1;
	    for( long j = 0 ; j < elist2.n ; j ++) {
		if( (elist2.array[j].n1 == n3 &&
			elist2.array[j].n2 == n4) ||
		    (elist2.array[j].n2 == n3 &&
			elist2.array[j].n1 == n4)) {
		    n34 = elist2.array[j].n3;
		    break;
		}
	    }
	    // if n34 has not yet been created, create it
	    // and put it into elist2
	    if( n34 == -1) {
		assert( n3 != n4);
		n34 = add_node(
		    new Node( nodes[n3], nodes[n4], r));
		elist2.push( SplitEdge( n3, n4, n34));
		//  				fprintf( stderr,
		//  					 "\tAdding split edge34 "
		//  					 "[%ld,%ld,%ld]\n",
		//  					 n3, n4, n34);
	    }
	    // remove element e from the model
	    remove_element( e);
	    delete e; e = NULL;
	    // replace with two new elements, preserving
	    // orientation of faces
	    if( n12_on_top && n12_clockwise) {
		add_element( new WedgeElement(
				 n1, n12, on12,
				 n3, n34, on34,
				 ym_map, pr_map, ys_map,
				 ft_map));
		add_element( new WedgeElement(
				 n12, n2, on12,
				 n34, n4, on34,
				 ym_map, pr_map, ys_map,
				 ft_map));
	    } else if( n12_on_top && ! n12_clockwise) {
		add_element( new WedgeElement(
				 n12, n1, on12,
				 n34, n3, on34,
				 ym_map, pr_map, ys_map,
				 ft_map));
		add_element( new WedgeElement(
				 n2, n12, on12,
				 n4, n34, on34,
				 ym_map, pr_map, ys_map,
				 ft_map));
	    } else if( ! n12_on_top && n12_clockwise) {
		add_element( new WedgeElement(
				 n3, n34, on34,
				 n1, n12, on12,
				 ym_map, pr_map, ys_map,
				 ft_map));
		add_element( new WedgeElement(
				 n34, n4, on34,
				 n12, n2, on12,
				 ym_map, pr_map, ys_map,
				 ft_map));
	    } else {
		add_element( new WedgeElement(
				 n34, n3, on34,
				 n12, n1, on12,
				 ym_map, pr_map, ys_map,
				 ft_map));
		add_element( new WedgeElement(
				 n4, n34, on34,
				 n2, n12, on12,
				 ym_map, pr_map, ys_map,
				 ft_map));
	    }
	}
    }

    // now we have to return the result
    if( elist2.is_empty()) return -1;
    // find n1,n2 in elist2 (it should actually be the first one
    // on the list
    while( ! elist2.is_empty()) {
	SplitEdge edge = elist2.pop();
	if( edge.n1 == n1 && edge.n2 == n2)
	    return edge.n3;
    }
    // hmmm - report an error, edge <n1,n2> not on the list of split
    // edges
    assert( 0);
    return -1;
}

void Model::recalculate_Kg( void)
// ----------------------------------------------------------------------
// assemble a global stiffness matrix from all elements, as well as
// the right hand side
// ----------------------------------------------------------------------
{
    if( Kg != NULL) delete Kg;
    if( RHSg != NULL) delete RHSg;

    Kg = new SparseMatrix( n_nodes*3, n_nodes*3);
    RHSg = new Vector_double( n_nodes * 3);
    long last_shown_value = -1;
    for( long i = 0 ; i < n_elements ; i ++) {
	elements[i]-> add_to_global( * Kg, * RHSg);
	long new_value = (100 * (i+1)) / n_elements;
	if( last_shown_value != new_value) {
	    fprintf( stderr, "Calculating global Kg... %ld%%\r",
		new_value);
	    last_shown_value = new_value;
	}
	// 		if( n_elements <= 100 || (i % (n_elements/100) == 0))
	// 			fprintf( stderr, "Calculating global Kg... %6.2f%%\r",
	// 				 double(i+1)/n_elements*100);
    }
    fprintf( stderr, "\n");
    //	printf( "global K is %ld x %ld\n",
    //		Kg-> n_rows, Kg-> n_cols);

    //	if( KgC != NULL) delete KgC; KgC = NULL;
    //	if( RHSgC != NULL) delete RHSgC; RHSgC = NULL;
}

void Model::insert_constraints_local( SparseMatrix & k, Vector_double & Sol)
// ----------------------------------------------------------------------
// eliminates entries from K where the solution is known
//    - solution is known at the nodes that are fixed & and not-selected
// ----------------------------------------------------------------------
{
    //#define is_knownl(n) ((n)->is_fixed) || ((n)->n_refs==0) || (!(n)->selected))
#define is_knownl(n) (!(n)->selected)
    // for every fixed node zero out the appropriate rows, then set
    // the proper diagonal entries and Sol[] entries
    for( long i = 0 ; i < n_nodes ; i ++) {
	if( ! is_knownl(nodes[i])) continue;
	// zero out the corresponding rows
	k.reset_row( i*3+0); k.reset_row( i*3+1); k.reset_row( i*3+2);
	// set corresponding diagonal entries to 1.0
	k.set( i*3+0, i*3+0, 1.0);
	k.set( i*3+1, i*3+1, 1.0);
	k.set( i*3+2, i*3+2, 1.0);
	// and setup the right hand side
	Sol.array[i*3+0] = nodes[i]->px;
	Sol.array[i*3+1] = nodes[i]->py;
	Sol.array[i*3+2] = nodes[i]->pz;
    }

    // from each row eliminate entries for which we know the value
    for( long row = 0 ; row < n_nodes*3 ; row ++) {
	// skip rows corresponding to known variables, since
	// these were processed above
	if( is_knownl( nodes[row/3])) {
	    assert( k.rows[row]-> n_entries == 1);
	    continue;
	}
	// name the row for convenience
	SparseMatrix::Row * r = k.rows[row];
	// go through entries in this row, and if you find one
	// that corresponds to a known variable, eliminate it
	// by subtracting its product with the known value from
	// the solution side
	for( long i = 0 ; i < r-> n_entries ; i ++) {
	    // name the node this column corresponds to
	    long col = r-> entries[i].col;
	    // name the node this column corresponds to
	    Node * n = nodes[col/3];
	    // if the column corresponds to a node that is not
	    // fixed, then we do not know the value of the
	    // variable, so we don't have to adjust the entry
	    if( ! is_knownl( n)) continue;
	    // otherwise, we extract the value from the matrix,
	    // and the value of the known variable, adjust
	    // the sol. vector, and mark the entry for
	    // deletion
	    double val_mat = r-> entries[i].val;
	    double val_known = Sol.array[col];
	    Sol.array[row] -= val_mat * val_known;
	    r-> entries[i].col = -1;
	}
	// delete the marked entries from the row
	r-> delete_marked_entries();
    }
}

void Model::insert_constraints( SparseMatrix & k, Vector_double & Sol)
// ----------------------------------------------------------------------
// eliminates entries from K where the solution is known
//    - solution is known at the nodes that are fixed
// ----------------------------------------------------------------------
{
#define is_known(n) (((n)->is_fixed) || ((n)->n_refs==0))
    // for every fixed node zero out the appropriate rows, then set
    // the proper diagonal entries and Sol[] entries
    for( long i = 0 ; i < n_nodes ; i ++) {
	//		if( ! nodes[i]-> is_fixed) continue;
	// 		if( (! nodes[i]-> is_fixed) &&
	// 		    ( nodes[i]-> n_refs > 0)) continue;
	if( ! is_known(nodes[i])) continue;
	// zero out the corresponding rows
	k.reset_row( i*3+0); k.reset_row( i*3+1); k.reset_row( i*3+2);
	// set corresponding diagonal entries to 1.0
	k.set( i*3+0, i*3+0, 1.0);
	k.set( i*3+1, i*3+1, 1.0);
	k.set( i*3+2, i*3+2, 1.0);
	// and setup the right hand side
	Sol.array[i*3+0] = nodes[i]->px;
	Sol.array[i*3+1] = nodes[i]->py;
	Sol.array[i*3+2] = nodes[i]->pz;
    }

    // from each row eliminate entries for which we know the value
    for( long row = 0 ; row < n_nodes*3 ; row ++) {
	// skip rows corresponding to known variables, since
	// these were processed above
	if( is_known( nodes[row/3])) {
	    continue;
	}
	// name the row for convenience
	SparseMatrix::Row * r = k.rows[row];
	// go through entries in this row, and if you find one
	// that corresponds to a known variable, eliminate it
	// by subtracting its product with the known value from
	// the solution side
	for( long i = 0 ; i < r-> n_entries ; i ++) {
	    // name the node this column corresponds to
	    long col = r-> entries[i].col;
	    // name the node this column corresponds to
	    Node * n = nodes[col/3];
	    // if the column corresponds to a node that is not
	    // fixed, then we do not know the value of the
	    // variable, so we don't have to adjust the entry
	    //			if( ! n-> is_fixed)) continue;
	    //			if( (! n-> is_fixed) && (n-> n_refs > 0)) continue;
	    if( ! is_known( n)) continue;
	    // otherwise, we extract the value from the matrix,
	    // and the value of the known variable, adjust
	    // the sol. vector, and mark the entry for
	    // deletion
	    double val_mat = r-> entries[i].val;
	    double val_known = Sol.array[col];
	    Sol.array[row] -= val_mat * val_known;
	    r-> entries[i].col = -1;
	}
	// delete the marked entries from the row
	r-> delete_marked_entries();
    }
}

/*
  int Model::is_background_node( long ind)
  // ----------------------------------------------------------------------
  // finds out whether a node is on the surface
  //   - node is a surface node if it is located
  //     in any element on its bottom face
  //   - special consideration: if n_refs of the node is 0, it is not
  //     a surface node
  // ----------------------------------------------------------------------
  {
  if( nodes[ ind]-> n_refs == 0) return 0;
  // pick any element sharing this node, i.e. the first one
  WedgeElement & e = * nodes[ind]-> refs[0];
  if( e.p[3] == ind || e.p[4] == ind || e.p[5] == ind)
  return 1;
  else
  return 0;
  }
*/

Vector3d Model::get_face_normal( long n)
// ----------------------------------------------------------------------
// calculate the surface normal at node <n>, which is the average
// of normals of all faces connected to this element
// ----------------------------------------------------------------------
{
    Vector3d res;
    for( long i = 0 ; i < nodes[n]-> n_refs ; i ++) {
	WedgeElement & e = * nodes[n]-> refs[i];
	long ni = e.get_internal_index( n);
	assert( ni > -1);
	if( ni < 3)
	    res = res + e.get_top_face_normal();
	else
	    res = res + e.get_bottom_face_normal();
    }
    res.normalize();
    if( res.length() < 0.5) {
	fprintf( stderr, "TBD: degenerate face normal !!!\n");
	res.x = res.y = res.z = 0;
    }
    return res;
}

Vector3d Model::get_surface_normal( long ind)
// ----------------------------------------------------------------------
// calculate the surface normal at this node, which is the average
// of normals of all elements connected to this node
// ----------------------------------------------------------------------
{
    Vector3d res;
    for( long i = 0 ; i < nodes[ind]-> n_refs ; i ++) {
	WedgeElement & e = * nodes[ind]-> refs[i];
	res = res + e.get_bottom_face_normal();
    }
    res.normalize();
    return res;
}

int Model::peel( double peel_stress)
// ----------------------------------------------------------------------
// finds a surface(background) node with the most stress acting in the
// direction of the surface normal
//
// if the stress is greater than some threshold, the node is made
// 'free', instead of fixed, which will have the effect of peeling
// the material from the surface
//
// Returns: 1 if something was peeled
//          0 if nothing was peeled
// ----------------------------------------------------------------------
{
    fprintf( stderr, "Peel():\n");
    long round = 1;
    while( 1) {
	fprintf( stderr, "\t-round %ld\n", round);
	relax();
	calculate_stresses();
	// consider every node
	double max_stress = -1;
	char something_peeled = 0;
	for( long i = 0 ; i < n_nodes ; i ++) {
	    Node & node = * nodes[i];
	    // if the node is not fixed, skip it
	    if( node.is_fixed != 1) continue;
	    //			// if the node is not a surface node, skip it
	    //			if( ! is_background_node(i)) continue;
	    // find the normal of the surface at the node
	    Vector3d sn = get_surface_normal( i);
	    // find the stress in the direction of the normal
	    Vector3d st = node.stensor * sn;
	    // project st on the surface normal (TBD: is this
	    // right to project the stress onto the normal?)
	    st = project( st, sn);
	    double val = st.length();
	    if( max_stress < val) max_stress = val;
	    if( val > peel_stress) {
		nodes[i]-> is_fixed = 0;
		something_peeled = 1;
	    }
	}
	fprintf( stderr, "\tmax_stress = %.10f\n", max_stress);
	if( ! something_peeled) break;
	round ++;
    }
    return 0;
}

void Model::relax( int forced)
{
    relax_local();
    (void) flag_nodes_with_error();
    relax_global( forced);
    (void) flag_nodes_with_error();
}

long Model::flag_nodes_with_error( void)
// ----------------------------------------------------------------------
// all nodes that can be moved: (! is_fixed and n_refs > 0)
// - figure out the error, if it is large, set large_err = 1
// ----------------------------------------------------------------------
{
    clock_t start_time = clock();

    // current solution vector:
    Vector_double X(3*n_nodes);
    for( long i = 0 ; i < n_nodes ; i ++) {
	X.array[i*3+0] = nodes[i]-> px;
	X.array[i*3+1] = nodes[i]-> py;
	X.array[i*3+2] = nodes[i]-> pz;
    }
    // error vector = KgC * X - RHSgC
    Vector_double err(3*n_nodes);
    KgC->mult(X,err);
    err.sub(*RHSgC);
    // mark all nodes above with large eror
    long n_flagged = 0;
    double max_err = -1;
    for( long i = 0 ; i < n_nodes ; i ++) {
	nodes[i]-> has_large_error = 0;
	if( nodes[i]-> is_fixed || nodes[i]-> n_refs == 0) continue;
	if( fabs(err(i*3+0)) >= precision ||
	    fabs(err(i*3+1)) >= precision ||
	    fabs(err(i*3+2)) >= precision ) {
	    nodes[i]-> has_large_error = 1;
	    n_flagged ++;
	}
	if( fabs(err(i*3+0)) > max_err) max_err = fabs(err(i*3+0));
	if( fabs(err(i*3+1)) > max_err) max_err = fabs(err(i*3+1));
	if( fabs(err(i*3+2)) > max_err) max_err = fabs(err(i*3+2));
    }

    clock_t end_time = clock();

    fprintf( stderr, "\t- flagged %ld nodes in %.4fs max_err = %.15f\n",
	n_flagged,
	double( end_time - start_time) / CLOCKS_PER_SEC,
	max_err);

    return n_flagged;
		 
}

void Model::relax_global (int forced)
{
    clock_t start = clock ();
    _relax_global (forced);
    clock_t end = clock ();
    relax_time += end - start;
}

void Model::relax_local ()
{
    clock_t start = clock ();
    _relax_local ();
    clock_t end = clock ();
    relax_time += end - start;
}

void Model::_relax_global( int forced)
// ----------------------------------------------------------------------
// WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
//
// - this function must be preceded immediately by relax_global()
//
// WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING WARNING
//
// - if the error obtained by relax_local() is already acceptable
//   then just calculate elemental stresses
// - otherwise perform relaxation on all free nodes
// - free nodes are: non_fixed, and with n_refs > 0
// - calculate principal stresses for every element
// - mark all free nodes as marked
// ----------------------------------------------------------------------
{
    // make sure we have KgC and RHSgC
    // ---------------------------------
    assert( KgC != NULL);
    assert( RHSgC != NULL);

    // only perform global relaxation if relax_local() left some nodes
    // with large error, or if it is forced
    if( relax_local_err >= precision || forced)
    {
	// create an initial guess vector
	// ------------------------------
	Vector_double X( 3*n_nodes);
	for( long i = 0 ; i < n_nodes ; i ++) {
	    X.array[i*3+0] = nodes[i]->px;
	    X.array[i*3+1] = nodes[i]->py;
	    X.array[i*3+2] = nodes[i]->pz;
	}
		
	// solve the system of equations using conjugate gradient
	// ------------------------------------------------------
	int res = conjugate_gradient(
	    *KgC, *RHSgC, X, precision);
			
	if( res) {
	    fprintf( stderr, "conjugate-gradient.res = %d\n", res);
	    assert( res == 0);
	}
		
	// retrieve the solution - displacements for the free nodes
	for( long i = 0 ; i < n_nodes ; i ++) {
	    Node * n = nodes[i];
	    n-> selected = 0;
	    if( n-> is_fixed || n-> n_refs == 0) continue;
	    n-> px = X.array[i*3+0];
	    n-> py = X.array[i*3+1];
	    n-> pz = X.array[i*3+2];
	    n-> selected = 1;
	    assert( ! isnan( n-> px));
	    assert( ! isnan( n-> py));
	    assert( ! isnan( n-> pz));
	}
    }

    // calculate stresses
    //	calculate_stresses();
}

void Model::_relax_local( void)
// ----------------------------------------------------------------------
// perform local relaxation - in the neighborhood of the nodes where the
// error of the current solution is large
//
// Steps:
//
//    - find the error of the current solution for all nodes
//    - find out where the error is large
//    - mark the nodes in a certain radius around large errors
//    - perform local relaxation - only relaxing nodes in the neighborhood
//      of large error
//    - recalculate error of the current solution
//    - if the error is below threshold, decrease radius
//    - otherwise increase radius
//    - relax_global() will finish the relaxation if the error was too
//                      large
//
// Output:
//    - KgC and RHSgC will contain the global system with constraints
//      inserted, so that relax_global() does not have to recompute them
//    - modified nodes will be marked
//    - relax_local_err will be set with the error of the final solution
// ----------------------------------------------------------------------
{
    assert( Kg != NULL);
    assert( RHSg != NULL);

    // create a right hand side vector
    // -------------------------------
    // resize Sol if necessary
    if( RHSgC-> size != 3 * n_nodes) {
	delete RHSgC; RHSgC = new Vector_double( 3*n_nodes);
    }
    // populate RHSgC
    for( long i = 0 ; i < n_nodes ; i ++) {
	RHSgC->array[i*3+0] = gravity_x + RHSg-> array[i*3+0];
	RHSgC->array[i*3+1] = gravity_y + RHSg-> array[i*3+1];
	RHSgC->array[i*3+2] = gravity_z + RHSg-> array[i*3+2];
    }
	
    // calculate KgC = Kg + constraints
    // --------------------------------
    KgC-> copy( * Kg);
    insert_constraints( * KgC, * RHSgC);

    // if local relaxation was not requested, just return, but
    // tell relax_global() to do its relaxation (by setting
    // the error to an 'unacceptable value'
    if( ! use_local_relaxation_flag) {
	relax_local_err = precision + 1;
	return;
    }

    // create the current solution vector
    // ----------------------------------
    Vector_double X( 3*n_nodes);
    for( long i = 0 ; i < n_nodes ; i ++) {
	X.array[i*3+0] = nodes[i]->px;
	X.array[i*3+1] = nodes[i]->py;
	X.array[i*3+2] = nodes[i]->pz;
    }

    // calculate how much different it is from the solution
    //   err = (KgC.X - RHSgC)
    Vector_double err( X.size);
    KgC->mult( X, err);
    err.sub( * RHSgC);
    // report the error of the solution (same calculation that CG uses!!!)
    double solerr = fabs(err.max_abs());
    fprintf( stderr, "\t- before local relax sol. error = %.20f\n",
	solerr);

    // set error for all nodes with a large error
    long n_marked = 0;
    for( long i = 0 ; i < n_nodes ; i ++) {
	Node * n = nodes[i];
	n-> selected = 0;
	// nodes that are fixed & nodes with no references
	// are marked as not-selected
	if( n-> is_fixed || n-> n_refs == 0) {
	    continue;
	}
	if( fabs( err.array[i*3+0]) >= precision ||
	    fabs( err.array[i*3+1]) >= precision ||
	    fabs( err.array[i*3+2]) >= precision)
	{
	    nodes[i]-> selected = 1;
	    n_marked ++;
	}
    }
    fprintf( stderr,
	"\t- local relax: marked nodes: %ld (%.2f%%)\n",
	n_marked,
	100 * double(n_marked) / n_nodes);
    clock_t start_time = clock();
    n_marked = mark_local_change_nodes( error_radius);
    clock_t end_time = clock();
    fprintf( stderr,
	"\t- after error diffusion, marked: %ld (%.2f%%) %.3fs\n",
	n_marked,
	100 * double(n_marked) / n_nodes,
	double( end_time - start_time) / CLOCKS_PER_SEC
	     );

    // debug - print out a big warning if nothing was marked, but the
    //         solution error is larger than precision (because that
    //         would be a conceptual error - if relax_local() does not
    //         do anything, then relax_global() should do nothing either
    if( solerr >= precision && n_marked == 0) {
	fprintf( stderr,
	    "  #     #   #   ####  #   # # #   #  ###     \n"
	    "  #     #  # #  #   # ##  # # ##  # #   #    \n"
	    "   # # #  #   # #   # ##  # # ##  # #        \n"
	    "   # # #  ##### ####  # # # # # # # # ###    \n"
	    "    # #   #   # #   # #  ## # #  ## #   #    \n"
	    "    # #   #   # #   # #   # # #   #  ###     \n"
	    "                                             \n"
	    " solerr(%.20f) >= precision(%.20f)\n",
	    solerr, precision);
    }

    // if nothing was marked, just tell global relax what the error
    // was
    if( n_marked == 0) {
	relax_local_err = solerr;
	return;
    }

    // perform local relaxation - treating only selected nodes as
    // the free nodes
    // ----------------------------------------------------------

    // create the right-hand side vector
    Vector_double SolL( 3*n_nodes);
    for( long i = 0 ; i < n_nodes ; i ++) {
	SolL.array[i*3+0] = gravity_x + RHSg-> array[i*3+0];
	SolL.array[i*3+1] = gravity_y + RHSg-> array[i*3+1];
	SolL.array[i*3+2] = gravity_z + RHSg-> array[i*3+2];
    }
    // create the coefficient matrix with removed known values
    SparseMatrix KgCL( * Kg);
    insert_constraints_local( KgCL, SolL);

    // create a compressed equivalents of SolL, X and KgCL

    // prepare translation vectors (back and forth)
    long transl[3*n_nodes];
    long transl_inv[n_marked];
    long ind = 0;
    for( long i = 0; i < n_nodes ; i ++) {
	if( nodes[i]-> selected) {
	    assert(i*3+2 < n_nodes * 3);
	    transl[i*3+0] = ind * 3 + 0;
	    transl[i*3+1] = ind * 3 + 1;
	    transl[i*3+2] = ind * 3 + 2;
	    assert( ind*3+2 < n_marked * 3);
	    transl_inv[ind] = i;
	    ind += 1;
	} else {
	    transl[i*3+0] = -1;
	    transl[i*3+1] = -1;
	    transl[i*3+2] = -1;
	}
    }
    assert( ind == n_marked);

    // using the translation vectors, create a compressed
    // version of:
    //      KgCL -> K_comp
    //      SolL -> Sol_comp
    //      X    -> X_comp
    SparseMatrix K_comp( n_marked * 3, n_marked * 3);
    Vector_double Sol_comp( n_marked * 3);
    Vector_double X_comp( n_marked * 3);
    for( long row = 0 ; row < 3 * n_nodes ; row ++) {
	long trow = transl[row];
	if( trow == -1) continue;
	SparseMatrix::Row * r = KgCL.rows[row];
	for( long j = 0 ; j < r-> n_entries ; j ++) {
	    long tcol = transl[r-> entries[j].col];;
	    if( tcol == -1) continue;
	    K_comp.set( trow, tcol, r-> entries[j].val);
	}
	Sol_comp.array[trow] = SolL.array[row];
	X_comp.array[trow] = X.array[row];
    }

    // solve the compressed system using conjugate gradient
    int res = conjugate_gradient(
	K_comp,Sol_comp,X_comp,precision);
    if( res) die( "conjugate-gradient-l.res = %d\n", res);
	
    // retrieve results from compressed X
    for( long row = 0 ; row < n_marked ; row ++) {
	long tind = transl_inv[row];
	assert( tind >= 0 && tind < n_nodes);
	Node * n = nodes[ tind];
	assert( n-> selected);
	n-> px = X_comp.array[row*3+0];
	n-> py = X_comp.array[row*3+1];
	n-> pz = X_comp.array[row*3+2];
	assert( ! isnan( n-> px));
	assert( ! isnan( n-> py));
	assert( ! isnan( n-> pz));
    }

    // recalculate error after local relaxation
    // -----------------------------------------------------------
    // current solution vector:
    for( long i = 0 ; i < n_nodes ; i ++) {
	X.array[i*3+0] = nodes[i]->px;
	X.array[i*3+1] = nodes[i]->py;
	X.array[i*3+2] = nodes[i]->pz;
    }
    // calculate the error of the solution
    //   err = (KgC.X - RHSgC)
    KgC-> mult( X, err);
    err.sub( * RHSgC);
    // report the error of the solution (same calculation that CG uses)
    solerr = fabs(err.max_abs());
    fprintf( stderr, "\t- after local relax sol. error = %.20f\n", solerr);

    // set the error so that relax_global() knows what we did here
    relax_local_err = solerr;

    // dynamically adjust the error radius
    //   - if global error is smaller than requested precision,
    //     try to make the radius smaller, otherwiser larger
    if( solerr < precision) {
	error_radius *= 0.99;
	if( error_radius < min_error_radius)
	    error_radius = min_error_radius;
	fprintf( stderr,
	    "\t- trying smaller local change radius "
	    "(r=%.20f\n",
	    error_radius
		 );
    } else {
	// try a larger increment if our solution was bad
	error_radius *= 1.5;
	fprintf( stderr,
	    "\t- trying LARGER local change radius "
	    "(r=%.20f\n",
	    error_radius
		 );
    }
}

void Model::draw( void)
{
    // if walls are not requested, still draw all selected elements
    if( draw_solid == 0) {
	for( long i = 0 ; i < n_elements ; i ++) {
	    glLoadName( i + 1);
	    WedgeElement & e = * elements[i];
	    if( ! e.selected) continue;
	    e.draw( WedgeElement::WALLS_ALL,
		ys_map-> get_min(),
		ys_map-> get_max(),
		draw_deformed_model);
	}
    }

    glDisable( GL_BLEND);
    for( long i = 0 ; i < n_elements ; i ++) {
	glLoadName( i + 1);
	WedgeElement & e = * elements[i];
	if( draw_solid == 1)
	    e.draw( WedgeElement::WALLS_TOP,
		ys_map-> get_min(),
		ys_map-> get_max(),
		draw_deformed_model);
	else if( draw_solid == 2)
	    e.draw( WedgeElement::WALLS_ALL,
		ys_map-> get_min(),
		ys_map-> get_max(),
		draw_deformed_model);
    }

    glLoadName( 0); // no more picking

    if( draw_wireframe) {
	glDisable( GL_LIGHTING);
	if( draw_solid) {
	    //			glDisable( GL_BLEND);
	    //			glLineWidth( 0.1);
	    //			glColor4f( 0, 0, 0, 0.5);
	    glEnable( GL_BLEND);
	    glColor4f( 0, 0, 0, 0.5);
	    glLineWidth( 1);
	    glDisable( GL_DEPTH_TEST);
	} else {
	    glEnable( GL_BLEND);
	    glColor4f( 0, 0, 0, 0.2);
	    glLineWidth( 0.1);
	    glDisable( GL_DEPTH_TEST);
	}
	for( long i = 0 ; i < n_elements ; i ++) {
	    WedgeElement & e = * elements[i];
	    if( draw_wireframe == 1)
		e.draw( WedgeElement::WIREFRAME_TOP,
		    ys_map-> get_min(),
		    ys_map-> get_max(),
		    draw_deformed_model);
	    else if( draw_wireframe == 2)
		e.draw( WedgeElement::WIREFRAME_ALL,
		    ys_map-> get_min(),
		    ys_map-> get_max(),
		    draw_deformed_model);
	}
	glEnable( GL_DEPTH_TEST);
	glEnable( GL_LIGHTING);
    }

    // draw nodes - how they are drawn depends on the value
    // of draw_nodes
    //
    // 0 = nothing gets drawn
    // 1 = only top nodes are drawn
    // 2 = top & bottom nodes are drawn
    // 3 = top nodes are drawn, toghether with numbers
    // 4 = top & bottom nodes are drawn, together with numbers
    if( draw_nodes > 0) {
	// draw the nodes
	glEnable( GL_LIGHTING);
	glEnable( GL_BLEND);
	for( long i = 0 ; i < n_nodes ; i ++) {
	    double x, y, z;
	    if( draw_deformed_model) {
		x = nodes[i]-> px;
		y = nodes[i]-> py;
		z = nodes[i]-> pz;
	    } else {
		x = nodes[i]-> ox;
		y = nodes[i]-> oy;
		z = nodes[i]-> oz;
	    }
	    glLoadName( n_elements + i + 1);
	    if( nodes[i]-> n_refs == 0) continue;
	    if( nodes[i]-> selected)
		if( nodes[i]-> has_large_error)
		    Material_SelectedErrorNode.set();
		else
		    Material_SelectedNode.set();
	    else
		if( nodes[i]-> has_large_error)
		    Material_UnselectedErrorNode.set();
		else
		    Material_UnselectedNode.set();
	    if( nodes[i]-> is_fixed)
		glDrawTetra( x, y, z, draw_node_cube_size*1.5);
	    else
		glDrawCube( x, y, z, draw_node_cube_size);
	    glLoadName( 0);
	}
	// now draw fracture tips (bigger cubes)
	for( long j = 0 ; j < ftips.n ; j ++) {
	    long i = ftips.array[j];
	    double x, y, z;
	    if( draw_deformed_model) {
		x = nodes[i]-> px;
		y = nodes[i]-> py;
		z = nodes[i]-> pz;
	    } else {
		x = nodes[i]-> ox;
		y = nodes[i]-> oy;
		z = nodes[i]-> oz;
	    }
	    glLoadName( n_elements + i + 1);
	    if( nodes[i]-> n_refs == 0) continue;
	    if( nodes[i]-> selected)
		if( nodes[i]-> has_large_error)
		    Material_SelectedErrorNode.set();
		else
		    Material_SelectedNode.set();
	    else
		if( nodes[i]-> has_large_error)
		    Material_UnselectedErrorNode.set();
		else
		    Material_UnselectedNode.set();
	    if( nodes[i]-> is_fixed)
		glDrawTetra( x, y, z, draw_node_cube_size*3);
	    else
		glDrawCube( x, y, z, draw_node_cube_size*2);
	    glLoadName( 0);
	}
	if( draw_nodes > 1) {
	    glDisable( GL_LIGHTING);
	    glDisable( GL_DEPTH_TEST);
	    glColor4f( 0, 0, 0, 1);
	    for( long i = 0 ; i < n_nodes ; i ++) {
		double x, y, z;
		if( draw_deformed_model) {
		    x = nodes[i]-> px;
		    y = nodes[i]-> py;
		    z = nodes[i]-> pz;
		} else {
		    x = nodes[i]-> ox;
		    y = nodes[i]-> oy;
		    z = nodes[i]-> oz;
		}
		glLoadName( n_elements + i + 1);
		// skip nodes with no references
		if( nodes[i]-> n_refs == 0) continue;
		draw_text( x, y, z, "%ld", i);
		glLoadName(0);
	    }
	    glEnable( GL_LIGHTING);
	    glEnable( GL_DEPTH_TEST);
	}
    }
}

void Model::draw_orig_shapes( void)
// ----------------------------------------------------------------------
// draw shapes of original elements
// ----------------------------------------------------------------------
{
    for( long i = 0 ; i < n_elements ; i ++) {
	WedgeElement * e = elements[i];
	Vector3d p[6];
	e-> get_orig_shape( p);
	// calculate the center of the deformed shape (the
	// center of the deformed bottom shape)
	Vector3d cd =
	    (
		Vector3d( nodes[e->p[3]]-> px,
		    nodes[e->p[3]]-> py,
		    nodes[e->p[3]]-> pz) +
		Vector3d( nodes[e->p[4]]-> px,
		    nodes[e->p[4]]-> py,
		    nodes[e->p[4]]-> pz) +
		Vector3d( nodes[e->p[5]]-> px,
		    nodes[e->p[5]]-> py,
		    nodes[e->p[5]]-> pz)) * (1.0/3);
	// calculate the center of the undeformed shape (the
	// center of the undeformed bottom shape)
	Vector3d cu = (p[3]+p[4]+p[5]) * (1.0/3);
	// shift all 6 points to be centered around the deformed
	// center
	for( long i = 0 ; i < 6 ; i ++)
	    p[i] = p[i] - cu + cd;

	struct edge {
	    long p1, p2;
	} edges[9] = { {0,1}, {1,2}, {2,0},
		       {3,4}, {4,5}, {5,3},
		       {0,3}, {1,4}, {2,5}};

	if( draw_original_shapes == 1) {
	    // draw the wireframe
	    glDisable( GL_LIGHTING);
	    glEnable( GL_BLEND);
	    glLineWidth( 0.1);
	    glColor4f( 0, 0, 0.7, 0.2);
	    glBegin( GL_LINES);
	    for( long i = 0 ; i < 9 ; i ++) {
		glVertex3d( p[edges[i].p1].x,
		    p[edges[i].p1].y,
		    p[edges[i].p1].z);
		glVertex3d( p[edges[i].p2].x,
		    p[edges[i].p2].y,
		    p[edges[i].p2].z);
	    }
	    glEnd();
	}

	if( draw_original_shapes == 2) {
	    // draw the faces
	    glEnable( GL_LIGHTING);
	    glDisable( GL_BLEND);
	    if( e-> selected)
		Material_YellowPlastic.set();
	    else
		Material_BluePlastic.set();
	    glDrawWedge( p);
	}
    }
	
}

static void disp_func( void)
// ----------------------------------------------------------------------
// called whenever the contents of the OGL window need to be redrawn.
// ----------------------------------------------------------------------
{
    // set the MODELVIEW matrix
    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();

    gluLookAt( view.eye.x, view.eye.y, view.eye.z,
	view.center.x, view.center.y, view.center.z,
	view.upv.x, view.upv.y, view.upv.z);

    // setup first light
    {
	float light_ambient[] = {0.6, 0.6, 0.6, 0.0 };
	float light_diffuse[] = {0.8, 0.8, 0.8, 0.0 };
	float light_specular[] = {0.8, 0.8, 0.8, 0.0 };
	//		float light_position[] = {3.0, 3.0, 3.0, 0.0 };
	float light_position[] = {0.5, 1.5, 3.0, 0.0 };
		
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);

	glEnable( GL_LIGHT0 );
    }

    // enable second light
    {
	float light_ambient[] = {0.0, 0.0, 0.0, 0.0 };
	float light_diffuse[] = {0.3, 0.3, 0.3, 0.0 };
	float light_specular[] = {0.0, 0.0, 0.0, 0.0 };
	float light_position[] = {-3.0, -3.0, -3.0, 0.0 };
		
	glLightfv(GL_LIGHT1, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT1, GL_SPECULAR, light_specular);
	glLightfv(GL_LIGHT1, GL_POSITION, light_position);
		
	//		glEnable( GL_LIGHT1 );
    }

    // set the background color to light grey
    glClearColor( 0.8, 0.8, 0.8, 1.0);
    //        glClearColor( 0.0, 0.0, 0.0, 1.0);
        
    // clear the background
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // enable automatic normalization of normals
    glEnable( GL_NORMALIZE);

    // enable Z buffer
    glEnable( GL_DEPTH_TEST);

    // enable / disable fog
    if( draw_fog > 0 ) {
	GLfloat fogColor[4] = { 0.8f, 0.8f, 0.8f, 1.0f};
	//		glClearColor( fogColor[0],fogColor[1],fogColor[2],fogColor[3]);
	if( draw_fog == 1)
	    glFogi( GL_FOG_MODE, GL_EXP);
	else if( draw_fog == 2)
	    glFogi( GL_FOG_MODE, GL_EXP2);
	else
	    glFogi( GL_FOG_MODE, GL_EXP);
	glFogfv( GL_FOG_COLOR, fogColor);
	glFogf( GL_FOG_DENSITY, 0.35f);
	glHint( GL_FOG_HINT, GL_DONT_CARE);
	glFogf( GL_FOG_START, 0.01f);
	glFogf( GL_FOG_END, 0.1f);
	glEnable( GL_FOG);
    } else {
	glDisable( GL_FOG);
    }

    // enable lighting
    glEnable( GL_LIGHTING );

    // draw axes
    if( draw_axes) {
	glDisable( GL_LIGHTING);
	glLineWidth( 2.0);
	glBegin( GL_LINES);
	glColor3f( 0.5, 0, 0);
	glVertex3d( 0, 0, 0); glVertex3d( 1, 0, 0);
	glColor3f( 0, 0.5, 0);
	glVertex3d( 0, 0, 0); glVertex3d( 0, 1, 0);
	glColor3f( 0, 0.4, 0.4);
	glVertex3d( 0, 0, 0); glVertex3d( 0, 0, 1);
	glEnd();
	glColor3f( 0.5, 0, 0);
	draw_text( 1, 0, 0, "X");
	glColor3f( 0, 0.5, 0);
	draw_text( 0, 1, 0, "Y");
	glColor3f( 0, 0.4, 0.4);
	draw_text( 0, 0, 1, "Z");
	glEnable( GL_LIGHTING);
    }

    // enable anti-aliasing of lines
    if( antialiasing) {
	glEnable( GL_BLEND);
	glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable( GL_LINE_SMOOTH);
	glLineWidth( 1);
    } else {
	glDisable( GL_BLEND);
	glLineWidth( 1);
    }

    // draw the model
    glLineWidth( 0.1);
    model-> draw();

    // draw surface normals
    if( draw_growth_vectors) {
	glDisable( GL_LIGHTING);
	glEnable( GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	//		glColor4f( 1, 0, 0, 0.5);
	glColor4f( 0, 0, 0, 1);
	glBegin( GL_LINES);
	for( long i = 0 ; i < n_nodes ; i ++) {
	    if( ! nodes[i]-> is_fixed) continue;
	    glVertex3d(
		nodes[i]-> px,
		nodes[i]-> py,
		nodes[i]-> pz);
	    glVertex3d(
		nodes[i]-> px + 0.3 * nodes[i]-> gv.x,
		nodes[i]-> py + 0.3 * nodes[i]-> gv.y,
		nodes[i]-> pz + 0.3 * nodes[i]-> gv.z);
	}
	glEnd();
	glDisable( GL_BLEND);
    }

    // draw elemental stresses
    if( draw_elemental_stresses) {
	glDisable( GL_LIGHTING);
	glEnable( GL_BLEND);
	glColor4f( 0, 0.5, 0, 1);
	glLineWidth( 1.5);
	glBegin( GL_LINES);
	for( long i = 0 ; i < n_elements ; i ++) {
	    WedgeElement * e = elements[i];
	    // for each surface gaussian point
	    for( long j=0 ; j<WedgeElement::n_gauss_surface ; j++)
	    {
		double se = e-> gps[j].e;
		double sn = e-> gps[j].n;
		double ss = e-> gps[j].s;
		// retrieve stress at the gaussian point
		StressTensor s;
		e-> calculate_stress_tensor_at_gp( s, j);
		// calculate eigen values from s
		double e1, e2, e3;
		s.get_eigenvalues( e1, e2, e3);
		// calculate eigenvectors corrsponding to the
		// extracted eigenvalue
		Vector3d v1, v2;
		int res = s.get_eigenvectors( e3, v1, v2);
		if( res != 1) continue;
		v1.normalize();
		if( e3 < 0) e3 = 0;
		v1.scale( e3*0.0001);

		// render results
		double na = e-> Na( se, sn, ss);
		double nb = e-> Nb( se, sn, ss);
		double nc = e-> Nc( se, sn, ss);
		double nd = e-> Nd( se, sn, ss);
		double ne = e-> Ne( se, sn, ss);
		double nf = e-> Nf( se, sn, ss);
		Node * n0 = nodes[ e-> p[0]];
		Node * n1 = nodes[ e-> p[1]];
		Node * n2 = nodes[ e-> p[2]];
		Node * n3 = nodes[ e-> p[3]];
		Node * n4 = nodes[ e-> p[4]];
		Node * n5 = nodes[ e-> p[5]];
		double x0 = n0-> px;
		double x1 = n1-> px;
		double x2 = n2-> px;
		double x3 = n3-> px;
		double x4 = n4-> px;
		double x5 = n5-> px;
		double y0 = n0-> py;
		double y1 = n1-> py;
		double y2 = n2-> py;
		double y3 = n3-> py;
		double y4 = n4-> py;
		double y5 = n5-> py;
		double z0 = n0-> pz;
		double z1 = n1-> pz;
		double z2 = n2-> pz;
		double z3 = n3-> pz;
		double z4 = n4-> pz;
		double z5 = n5-> pz;
		double tx =
		    na * x0 + nb * x1 + nc * x2 +
		    nd * x3 + ne * x4 + nf * x5;
		double ty =
		    na * y0 + nb * y1 + nc * y2 +
		    nd * y3 + ne * y4 + nf * y5;
		double tz =
		    na * z0 + nb * z1 + nc * z2 +
		    nd * z3 + ne * z4 + nf * z5;
		glVertex3d( tx-v1.x/2,ty-v1.y/2,tz-v1.z/2);
		glVertex3d( tx+v1.x/2,ty+v1.y/2,tz+v1.z/2);
	    }
	}
	glEnd();
	glDisable( GL_BLEND);
    }

    // draw max. principal stress vectors as requested
    if( draw_nodal_stresses > 0) {
	glEnable( GL_BLEND);
	glDisable( GL_LIGHTING);
	for( long i = 0 ; i < n_nodes ; i ++) {
	    if( draw_nodal_stresses < 3)
		if( nodes[i]-> is_fixed != 0) continue;
	    if( nodes[i]-> is_fixed == 0)
		glColor4f( 1, 0, 0, 0.5);
	    else
		glColor4f( 0.8, 0.4, 0, 0.5);
	    Node & n = * nodes[i];
	    glLineWidth( 0.5);
	    if( n.selected || 1) {
		glColor4f( 1, 0, 0, 1);
		glLineWidth( 3.0);
	    }
	    glBegin( GL_LINES);
	    double nx = n.ps_nx * n.ps_val * 0.001;
	    double ny = n.ps_ny * n.ps_val * 0.001;
	    double nz = n.ps_nz * n.ps_val * 0.001;
	    glVertex3d( n.px - nx/2,
		n.py - ny/2,
		n.pz - nz/2);
	    glVertex3d( n.px + nx/2,
		n.py + ny/2,
		n.pz + nz/2);
	    glEnd();
	}
	glDisable( GL_BLEND);
    }

    // draw stress planes as requested
    if( draw_nodal_stresses > 1) {
	glEnable (GL_BLEND);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glDisable( GL_LIGHTING);
	for( long i = 0 ; i < n_nodes ; i ++) {
	    if( draw_nodal_stresses < 3)
		if( nodes[i]-> is_fixed != 0) continue;
	    if( nodes[i]-> is_fixed == 0)
		glColor4f( 0, 0, 1, 0.2);
	    else
		glColor4f( 0, 0.3, 1, 0.2);
	    Node & n = * nodes[i];
	    Vector3d nor( n.ps_nx, n.ps_ny, n.ps_nz);
	    // calculate u = nor x (0,0,1)
	    Vector3d u( nor);
	    Vector3d tmp( 0, 0, 1);
	    u.cross_product( tmp);
	    u.normalize();
	    // calculate v = nor x u
	    Vector3d v( nor);
	    v.cross_product( u);
	    v.normalize();
	    // draw circle
	    glBegin( GL_TRIANGLE_FAN);
	    long n_fans = 10;
	    glVertex3d( n.px, n.py, n.pz);
	    for( long k = 0 ; k <= n_fans ; k ++) {
		double alpha = k * M_PI * 2 / n_fans;
		double s = cos( alpha);
		double t = sin( alpha);
		double x = s * u.x + t * v.x;
		double y = s * u.y + t * v.y;
		double z = s * u.z + t * v.z;
		double scale =
		    model-> crack_tip_element_size / 10;
		x *= scale;
		y *= scale;
		z *= scale;
		glVertex3d( n.px + x,
		    n.py + y,
		    n.pz + z);
	    }
	    glEnd();
	}
		
    }

    if( draw_original_shapes)
	model-> draw_orig_shapes();

    // set the title of the window
    char buff[ 4096];
    double etr = (clock()/model-> time_curr) *
	(model-> time_total- model-> time_curr);
    long etr_s = long(etr / CLOCKS_PER_SEC);
    long etr_d = etr_s / (24*60*60);
    etr_s = etr_s % (24*60*60);
    long etr_h = etr_s / (60*60);
    etr_s = etr_s % (60*60);
    long etr_m = etr_s / (60);
    etr_s = etr_s % (60);
    sprintf( buff,
	"%s: t=%.5f dt=%.6f ne=%ld nn=%ld "
	"nf=%ld(%.4fs/f) "
	"ETR=%ldd%ldh%ldm%lds "
	"rt=%.1f st=%.1f",
	model-> input_fname, model-> time_curr,
	model-> dtc * model-> min_time_step,
	n_elements, n_nodes,
	model-> n_fractures, 
	(double(clock())/ CLOCKS_PER_SEC)/model-> n_fractures,
	etr_d, etr_h, etr_m, etr_s,
	double (relax_time) / CLOCKS_PER_SEC,
	double (stress_time) / CLOCKS_PER_SEC
    );
    glutSetWindowTitle( buff);

    // swap buffers
    glutSwapBuffers();
}

static void reshape_func( int width, int height)
// ----------------------------------------------------------------------
// called whenever the OGL window is reshaped
// ----------------------------------------------------------------------
{
    win_width = width;
    win_height = height;
    glViewport( 0, 0, width, height );
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    gluPerspective( 60.0, width/double(height), 0.25, 4.0);
}

int mouse_x, mouse_y;
enum MouseState { ROTATE, ZOOM, TRANSLATE, NONE } mouse_state = NONE;

void mouse_func( int button, int state, int x, int y)
// ----------------------------------------------------------------------
// called when mouse buttons are pressed
// ----------------------------------------------------------------------
{
    if( state == GLUT_UP) {
	// mouse was released
	mouse_state = NONE;
	return;
    }

    fprintf( stderr, "modifiers = %d\n", int( glutGetModifiers()));

    // get the keyboard modifiers
    int shift_down = glutGetModifiers() & 1;
    int ctrl_down = glutGetModifiers() & 2;
    int alt_down = glutGetModifiers() & 4;
    if( button == GLUT_LEFT_BUTTON && ! shift_down && ! alt_down &&
	! ctrl_down)
	mouse_state = ROTATE;
    else if( button == GLUT_LEFT_BUTTON && shift_down && ! alt_down)
	mouse_state = ZOOM;
    else if( button == GLUT_LEFT_BUTTON && ! shift_down && alt_down)
	mouse_state = TRANSLATE;
    else if( button == GLUT_LEFT_BUTTON && ! shift_down && ctrl_down)
	mouse_state = TRANSLATE;
    else if( button == GLUT_MIDDLE_BUTTON)
	mouse_state = ZOOM;
    else
	mouse_state = NONE;

    mouse_x = x;
    mouse_y = y;
}

void mouse_motion_func( int x, int y)
// ----------------------------------------------------------------------
// called when mouse moves
// ----------------------------------------------------------------------
{
    // figure out the displacement of the mouse
    int dx = x - mouse_x;
    int dy = y - mouse_y;
    // remember the mouse coordinates
    mouse_x = x;
    mouse_y = y;


    // depending on the mode select an action
    if( mouse_state == ROTATE) {
	view.rotate_eye( dx, -dy);
	glutPostRedisplay();
    } else if( mouse_state == ZOOM) {
	view.zoom( dy / 100.0);
	glutPostRedisplay();
    } else if( mouse_state == TRANSLATE) {
	view.translate( dx / 500.0, - dy / 500.0);
	glutPostRedisplay();
    }
}

long Model::add_node( Node * n)
// ----------------------------------------------------------------------
// adds a node to the list, returns the index of the new node
// ----------------------------------------------------------------------
{
    // preconditions:
    assert( n_nodes <= n_nodes_alloc);

    // reallocate the array of nodes if needed
    if( n_nodes == n_nodes_alloc) {
	n_nodes_alloc = long( n_nodes_alloc * 1.5 + 10);
	nodes = (Node **) realloc( nodes,
	    sizeof( Node *) * (n_nodes_alloc));
    }

    // put the new node at the end
    assert( nodes != NULL);
    nodes[ n_nodes] = n;
    n_nodes ++;

    // expand the global stiffness matrix accordingly
    Kg-> expand( n_nodes*3, n_nodes*3);
    RHSg-> resize( n_nodes * 3);

    return n_nodes - 1;
}

long Model::add_element( WedgeElement * e)
// ----------------------------------------------------------------------
// adds an element to the list, returns the index of the new element
// ----------------------------------------------------------------------
{
    elements = (WedgeElement **)
	realloc( elements,
	    sizeof( WedgeElement *) *(n_elements + 1));
    assert( elements != NULL);
    elements[ n_elements] = e;
    n_elements ++;

    // calculate the shape of the element
    e-> grow( time_curr,
	shrink_top_t0,
	shrink_top_t1,
	shrink_top_val,
	shrink_bot_t0,
	shrink_bot_t1,
	shrink_bot_val,
	shrink_height_t0,
	shrink_height_val0,
	shrink_height_t1,
	shrink_height_val1);

    // add the stiffness matrix of this element to the global K
    e-> add_to_global( * Kg, * RHSg);
	
    // adjust references to the nodes
    // --------------------------------------------------
    for( long i = 0 ; i < 6 ; i ++)
	nodes[e->p[i]]-> add_ref( e);

    return n_elements-1;
}

int Model::check_ref_consistency( void)
// ----------------------------------------------------------------------
// makes sure that the references in nodes are correct
//
// - returns: 0 = everything is OK
//            1 = there is a reference in a node to an element which
//                does not contain the node
//            2 = an element references a node but the node does not
//                have a reference to it
{
    // make sure all nodes have only references to elements which
    // contain them
    for( long i = 0 ; i < n_nodes ; i ++) {
	Node & n = * nodes[i];
	for( long j = 0 ; j < n.n_refs ; j ++) {
	    WedgeElement & e = * n.refs[j];
	    if( e.p[0] != i &&
		e.p[1] != i &&
		e.p[2] != i &&
		e.p[3] != i &&
		e.p[4] != i &&
		e.p[5] != i) return 1;
	}
    }

    // make sure that all elements have a representation in all their
    // nodes
    for( long i = 0 ; i < n_elements ; i ++) {
	WedgeElement & e = * elements[i];
	for( long j = 0 ; j < 6 ; j ++) {
	    assert( e.p[j] >= 0 && e.p[j] < n_nodes);
	    Node & n = * nodes[ e.p[j]];
	    int found = 0;
	    for( long k = 0 ; k < n.n_refs ; k ++)
		if( n.refs[k] == & e) {
		    found = 1;
		    break;
		}
	    if( ! found) return 2;
	}
    }

    // everything checked out
    return 0;
}

double Model::get_Kg_error( void)
// returns the sum square differences between the current stiffness matrix
// and the correct one
{
    SparseMatrix Ktmp( Kg-> n_rows, Kg-> n_cols);
    Ktmp.copy( * Kg);
    recalculate_Kg();

    double err_sum = 0;
    for( long i = 0 ; i < Kg-> n_rows ; i ++) {
	SparseMatrix::Row * row1 = Ktmp.rows[i];
	for( long j = 0 ; j < row1-> n_entries ; j ++) {
	    double err = row1->entries[j].val -
		Kg->get(i,row1->entries[j].col);
	    err = err * err;
	    if( err > 0.001) {
		fprintf( stderr,
		    "Kg: big error at [%ld,%ld]\n",
		    i/3, j/3);
	    }
	    err_sum += err;
	}
    }
    // put the old matrix back
    Kg-> copy( Ktmp);
    return err_sum;
}

void Model::remove_element( WedgeElement * e)
// ----------------------------------------------------------------------
// remove element ind from the array
// ----------------------------------------------------------------------
{
    // find the index of this element
    long ind = get_element_index( e);
    assert( ind != -1);

    // put the last element into the place of this one
    if( n_elements > 0)
	elements[ ind] = elements[ n_elements-1];
    /*
    // shift the array of element pointers
    if( ind < n_elements - 1)
    memmove( & elements[ind], & elements[ind+1],
    sizeof( elements[ind]) * (n_elements-1-ind));
    */
    n_elements --;

    // remove the stiffness matrix of this element from the global K
    e-> remove_from_global( * Kg, * RHSg);

    // adjust references to the nodes
    for( long i = 0 ; i < 6 ; i ++)
	nodes[e->p[i]]-> del_ref( e);
}

WedgeElement * 
    find_element( WedgeElement * exclude_e, long n1_ind, long n2_ind)
// ----------------------------------------------------------------------
// finds an element that has nodes n1 & n2, but is different from e
//
// returns NULL if element not found
// ----------------------------------------------------------------------
{
    Node & n1 = * nodes[n1_ind];
    for( long i = 0 ; i < n1.n_refs ; i ++) {
	WedgeElement * e = n1.refs[i];
	if( e == exclude_e) continue;
	if( e-> p[0] == n2_ind ||
	    e-> p[1] == n2_ind ||
	    e-> p[2] == n2_ind)
	{
	    return e;
	}
    }
    return NULL;
}

long find_element_ind( WedgeElement * exclude_e, long n1_ind, long n2_ind)
// ----------------------------------------------------------------------
// finds an element that has nodes n1 & n2, but is different from e
//
// returns -1 if element not found
// ----------------------------------------------------------------------
{
    long res = model-> get_element_index(
	find_element( exclude_e, n1_ind, n2_ind));
    if( res < 0) return res;
    assert( res < n_elements);
    assert( exclude_e != elements[ res]);
    return res;
}

long find_oriented_element( long n1_ind, long n2_ind)
// ----------------------------------------------------------------------
// finds an element that has nodes n1 & n2
//
// returns -1 if element not found
// ----------------------------------------------------------------------
{
    Node & n1 = * nodes[ n1_ind];
    for( long i = 0 ; i < n1.n_refs ; i ++) {
	WedgeElement * e = n1.refs[i];
	if( (e-> p[0] == n1_ind && e-> p[1] == n2_ind) ||
	    (e-> p[1] == n1_ind && e-> p[2] == n2_ind) ||
	    (e-> p[2] == n1_ind && e-> p[0] == n2_ind))
	{
	    return model-> get_element_index( e);
	}
    }
    return -1;
}

static int opposite_sides( long i1, long i2, long i3, long i4)
// ----------------------------------------------------------------------
// returns TRUE, if n3 and n4 are distinctly on the opposite sides
// of a line going through n1 and n2, otherwise returns false
// ----------------------------------------------------------------------
{
    assert( i1 >= 0 && i1 < n_nodes);
    assert( i2 >= 0 && i2 < n_nodes);
    assert( i3 >= 0 && i3 < n_nodes);
    assert( i4 >= 0 && i4 < n_nodes);

    Node * n1 = nodes[i1];
    Node * n2 = nodes[i2];
    Node * n3 = nodes[i3];
    Node * n4 = nodes[i4];

    // calculate vector v = n2 - n1
    double vx = n2-> ox - n1-> ox;
    double vy = n2-> oy - n1-> oy;

    // calculate vector u = n3 - n1
    double ux = n3-> ox - n1-> ox;
    double uy = n3-> oy - n1-> oy;

    // calculate vector w = n4 - n1
    double wx = n4-> ox - n1-> ox;
    double wy = n4-> oy - n1-> oy;

    // calculate the Z component of the vector product vu and vw
    double vuz = vx*uy - vy*ux;
    double vwz = vx*wy - vy*wx;

    // TBD: adjust vuz and vwz to ZEROS if really small
    if( fabs( vuz) < 0.0001) vuz = 0;
    if( fabs( vwz) < 0.0001) vwz = 0;

    int res = vuz*vwz < 0;
    return res;
}

Edge find_edge_above( const Edge ed)
// ----------------------------------------------------------------------
// given an edge <ed>, try to find an edge above it (assume <ed> is on
// some bottom face of an element, so return the corresponding edge
// on the top face)
// 
// if such edge cannot be found, the result will have <n1> = -1
// ----------------------------------------------------------------------
{
    Edge res(-1,0);
    for( long i = 0 ; i < nodes[ed.n1]-> n_refs ; i ++) {
	WedgeElement * e = nodes[ed.n1]-> refs[i];
	long in1 = e-> get_internal_index( ed.n1);
	assert( in1 != -1);
	if( in1 < 3) continue;
	long in2 = e-> get_internal_index( ed.n2);
	if( in2 == -1) continue;
	if( in2 < 3) continue;
	res.n1 = e-> p[in1-3];
	res.n2 = e-> p[in2-3];
	break;
    }
    return res;
}

int Model::delaunay_edge_flip( long tn1, long tn2)
// ----------------------------------------------------------------------
// - finds both triangles this edge belongs to (this oriented edge)
// - if the distance between the opposite vertices of the two triangles
//   is smaller than the edge, the edge is deleted and replaced by
//   new, shorter edge, thus replacing the two triangles with two new
//   ones
// - if flip has been performed, recursively check the edges of the
//   new two triangles
// - if the edge does not belong to two triangles, nothing happens
//
// - returns: 0 - if flip was not performed
//            1 - if flip was performed
// ----------------------------------------------------------------------
{
    assert( tn1 >= 0 && tn1 < n_nodes);
    assert( tn2 >= 0 && tn2 < n_nodes);

    // find both triangles
    long ind1 = find_oriented_element( tn1, tn2);
    if( ind1 == -1) return 0;
    WedgeElement * e1 = elements[ ind1];
    long ind2 = find_oriented_element( tn2, tn1);
    if( ind2 == -1) return 0;
    WedgeElement * e2 = elements[ ind2];

    // calculate n3 (the third vertex from e1) and n4 (the third
    // vertex from e2)
    long tn3, tn4;
    if( e1->p[0]!=tn1 && e1->p[0]!=tn2) tn3 = e1->p[0];
    else if( e1->p[1]!=tn1 && e1->p[1]!=tn2) tn3 = e1->p[1];
    else tn3 = e1-> p[2];
    if( e2->p[0]!=tn1 && e2->p[0]!=tn2) tn4 = e2->p[0];
    else if( e2->p[1]!=tn1 && e2->p[1]!=tn2) tn4 = e2->p[1];
    else tn4 = e2-> p[2];

    // calculate bn1, bn2, bn3 and bn4
    long bn1=-1, bn2=-1, bn3=-1, bn4=-1;
    if( e1-> p[0] == tn1) bn1 = e1-> p[3];
    else if( e1-> p[1] == tn1) bn1 = e1-> p[4];
    else bn1 = e1-> p[5];
    if( e1-> p[0] == tn2) bn2 = e1-> p[3];
    else if( e1-> p[1] == tn2) bn2 = e1-> p[4];
    else bn2 = e1-> p[5];
    if( e1-> p[0] == tn3) bn3 = e1-> p[3];
    else if( e1-> p[1] == tn3) bn3 = e1-> p[4];
    else bn3 = e1-> p[5];
    if( e2-> p[0] == tn4) bn4 = e2-> p[3];
    else if( e2-> p[1] == tn4) bn4 = e2-> p[4];
    else bn4 = e2-> p[5];

    // calculate the length of the existing edge (n1-n2)
    double l12 = Node::dist( nodes[tn1], nodes[tn2]);
    // calculate the length of the replacing edge (n3-n4)
    double l34 = Node::dist( nodes[tn3], nodes[tn4]);

    // if the replacing edge would be longer, do nothing
    if( l34 >= l12) return 0;

    // if the new edge does not intersect the old edge, do nothing
    if( ! opposite_sides( tn3, tn4, tn1, tn2)) return 0;

    // delete the elements e1 and e2 from the list
    remove_element( e1);
    remove_element( e2);

    // add two new elements, (n1, n3, n4) and (n2, n3, n4)
    add_element( new WedgeElement(
		     tn1, tn4, tn3, bn1, bn4, bn3,
		     ym_map, pr_map, ys_map, ft_map));
    add_element( new WedgeElement(
		     tn2, tn3, tn4, bn2, bn3, bn4,
		     ym_map, pr_map, ys_map, ft_map));
	
    // delete elements e1 and e2 from memory
    delete e1; e1 = NULL;
    delete e2; e2 = NULL;

    // try to flip also the edge on the bottom
    delaunay_edge_flip( bn1, bn2);

    // try to flip also the edge on the top
    Edge ea = find_edge_above( Edge( tn1, tn2));
    if( ea.n1 != -1)
	delaunay_edge_flip( ea.n1, ea.n2);

    // recursive check
    delaunay_edge_flip( tn1, tn3);
    delaunay_edge_flip( tn1, tn4);
    delaunay_edge_flip( tn2, tn3);
    delaunay_edge_flip( tn2, tn4);

    // we have flipped the edge
    return 1;
}

#ifdef DONT_COMPILE
int Model::delaunay_edge_flip( long tn1, long tn2)
// ----------------------------------------------------------------------
// - finds both wedges this edge belongs to (so that the edge
//   is completely on the top face)
// - if such two wedges do not exist, nothing is done
// - finds tn3 and tn4 - the oposite nodes of the two wedges on top faces
// - finds a spot (r) on tn1,tn2 where the combined distance to both
//   tn3 and tn4 is the smallest
// - calculates the combined distance
// - if the combined distance is larger than the distance tn1,tn2, do nothing
// - split edge tn1,tn2 at (r)
// - recursively call yourself on the outside 4 edges
//
// - returns: 0 - if flip was not performed
//            1 - if flip was performed
// ----------------------------------------------------------------------
{
    return 0;

    assert( tn1 >= 0 && tn1 < n_nodes);
    assert( tn2 >= 0 && tn2 < n_nodes);

    // find the first wedge
    WedgeElement * e1 = find_element( NULL, tn1, tn2);
    if( e1 == NULL) return 0;
    // find the second wedge
    WedgeElement * e2 = find_element( e1, tn1, tn2);
    if( e2 == NULL) return 0;
    // find tn3
    long tn3 = -1;
    if( e1->p[0] != tn1 && e1->p[0] != tn2) tn3 = e1->p[0];
    else if( e1->p[1] != tn1 && e1->p[1] != tn2) tn3 = e1->p[1];
    else tn3 = e1->p[2];
    // find tn4
    long tn4 = -1;
    if( e2->p[0] != tn1 && e2->p[0] != tn2) tn4 = e2->p[0];
    else if( e2->p[1] != tn1 && e2->p[1] != tn2) tn4 = e2->p[1];
    else tn4 = e2->p[2];
    // find r - the isoparametric coordinate on edge tn1,tn2 where
    // the combined distance from r to tn3 and tn4 is the smallest
    double x1 = nodes[tn1]-> ox;
    double y1 = nodes[tn1]-> oy;
    double z1 = nodes[tn1]-> oz;
    double x2 = nodes[tn2]-> ox;
    double y2 = nodes[tn2]-> oy;
    double z2 = nodes[tn2]-> oz;
    double x3 = nodes[tn3]-> ox;
    double y3 = nodes[tn3]-> oy;
    double z3 = nodes[tn3]-> oz;
    double x4 = nodes[tn4]-> ox;
    double y4 = nodes[tn4]-> oy;
    double z4 = nodes[tn4]-> oz;
    double r = ((x1 - x2)*(2*x1 - x3 - x4) + (y1 - y2)*(2*y1 - y3 - y4) + 
	(z1 - z2)*(2*z1 - z3 - z4))
	/ (2*(sqr(x1 - x2) + sqr(y1 - y2) + sqr(z1 - z2)));
    // if 'r' cannot be calculated, something is fishy (i.e. tn1 =
    // tn2???, or too small and numerical precision of the computer is
    // exhausted) so no adjustment will be done
    if( ! finite( r)) return 0;
    // if r is not on the edge tn1-tn2, no adjustment will be done either
    if( r <= 0 || r >= 1) return 0;
    // calculate the cartesian coordinates of the point <r>
    double px = x1 * (1-r) + x2 * r;
    double py = y1 * (1-r) + y2 * r;
    double pz = z1 * (1-r) + z2 * r;
    // calculate the distance between n1 and n2
    double dist12 = Node::dist( nodes[tn1], nodes[tn2]);
    // calculate the distance between tn3 and tn4, which is calculated
    // as the sum of the distances between tn3 and <r> and tn4 and <r>
    double dist34 =
	sqrt( sqr(x3-px) + sqr(y3-py) + sqr(z3-pz)) +
	sqrt( sqr(x4-px) + sqr(y4-py) + sqr(z4-pz));
    // if the distance between tn1-tn2 is smaller than tn3-tn4, no
    // adjustment is needed
    if( dist12 < dist34) return 0;
    // adjustment is needed, so do it:
    split_edge( tn1, tn2, r);
    // recursively adjust the outside edges
    delaunay_edge_flip( tn1, tn3);
    delaunay_edge_flip( tn1, tn4);
    delaunay_edge_flip( tn2, tn3);
    delaunay_edge_flip( tn2, tn4);

    // we have adjusted the edge
    return 1;
}
#endif

void Model::subdivide_element_bisection( long ind)
// ----------------------------------------------------------------------
// subdivides the selected element
//
// - finds the longest edge of this triangle
// - repeat:
//     - find the neighboring triangle (sharing the longest edge)
//     - if the neighboring triangle's longest edge is not the
//       shared edge, recursively subdivide the neighboring triangle
//     - otherwise split this triangle and the neighboring one
// ----------------------------------------------------------------------
{
    //	fprintf( stderr, "bisecting el. %ld\n", ind);
    assert( ind >= 0 && ind < n_elements);

    // get a pointer to the element we are subdividing
    WedgeElement * e = elements[ ind];

    // find out which side is the longest one, so that n1 is the
    // opposite node while n2 and n3 reflect the longest edge
    double l0 = Node::dist(nodes[e->p[1]],nodes[e->p[2]]);
    double l1 = Node::dist(nodes[e->p[0]],nodes[e->p[2]]);
    double l2 = Node::dist(nodes[e->p[0]],nodes[e->p[1]]);
    //  	fprintf( stderr, "\t-l0 = %.20f\n", l0);
    //  	fprintf( stderr, "\t-l1 = %.20f\n", l1);
    //  	fprintf( stderr, "\t-l2 = %.20f\n", l2);
    double maxl=-1;
    long tn1, tn2, tn3, bn1, bn2, bn3;
    if( l0 >= l1 && l0 >= l2) {
	tn1=e->p[0]; tn2=e->p[1]; tn3=e->p[2];
	bn1=e->p[3]; bn2=e->p[4]; bn3=e->p[5];
	maxl = l0;
    } else if( l1 >= l0 && l1 >= l2) {
	tn1=e->p[1]; tn2=e->p[2]; tn3=e->p[0];
	bn1=e->p[4]; bn2=e->p[5]; bn3=e->p[3];
	maxl = l1;
    } else {
	tn1=e->p[2]; tn2=e->p[0]; tn3=e->p[1];
	bn1=e->p[5]; bn2=e->p[3]; bn3=e->p[4];
	maxl = l2;
    }

    // delete the triangle we are subdividing from the list
    remove_element( e);

    long neigh_ind = find_element_ind( e, tn2, tn3);
    int had_neighbor_before = neigh_ind != -1;
    //	fprintf( stderr, "\t- had neighbor = %ld\n", neigh_ind);
    //	fprintf( stderr, "\t- wanting to split edge %ld %ld\n", tn2, tn3);

    // find which triangle is on the opposite side of n2-n3
 loop:
    long ind_opp = find_element_ind( e, tn2, tn3);
    //	fprintf( stderr, "\t- neighbor is %ld\n", ind_opp);
    WedgeElement * ep = NULL;
    if( ind_opp != -1) ep = elements[ ind_opp];

    int has_neighbor_now = ep != NULL;
    assert( had_neighbor_before == has_neighbor_now);

    // find node n4 (the opposite node of n2-n3)
    long tn4 = -1, bn4 = -1;
    if( ep != NULL) {
	if(ep->p[0]!=tn2&&ep->p[0]!=tn3){
	    tn4=ep->p[0];bn4=ep->p[3];
	} else if(ep->p[1]!=tn2&&ep->p[1]!=tn3){
	    tn4=ep->p[1];bn4=ep->p[4];
	} else {
	    tn4=ep->p[2];bn4=ep->p[5];
	}
    }
 
    // if n2-n4 or n3-n4 is bigger than maxl, then recursively
    // subdivide the neighbor and try again
    if( ep != NULL
	&& ( Node::dist( nodes[tn2], nodes[tn4]) * 0.9999 > maxl
	    || Node::dist( nodes[tn3], nodes[tn4]) * 0.9999 > maxl))
    {
	//		fprintf( stderr, "\t- recursive subdivision %ld\n", ind_opp);
	subdivide_element_bisection( ind_opp);
	goto loop;
    } else {
	//		fprintf( stderr, "\t- no recursion\n");
    }

    // put the original triangle back in the list
    add_element( e);
    // subdivide the edge tn2, tn3
    //	fprintf( stderr, "\t- finally splitting edge %ld %ld\n", tn2, tn3);
    split_edge( tn2, tn3, 0.5);

    // perform delaunay adjustments
    delaunay_edge_flip( tn1, tn2);
    delaunay_edge_flip( tn1, tn3);
    if( ep != NULL) {
	delaunay_edge_flip( tn4, tn2);
	delaunay_edge_flip( tn4, tn3);
    }
}

static int circumcenter(
    const Vector3d & p1,
    const Vector3d & p2,
    const Vector3d & p3,
    double & t,
    double & s)
{
    double x2 = p2.x - p1.x;
    double y2 = p2.y - p1.y;
    double z2 = p2.z - p1.z;

    double x3 = p3.x - p1.x;
    double y3 = p3.y - p1.y;
    double z3 = p3.z - p1.z;

    double d = 2*(sqr(x3)*(sqr(y2) + sqr(z2)) + 
	sqr(y3*z2 - y2*z3) - 2*x2*x3*(y2*y3 + z2*z3) + 
	sqr(x2)*(sqr(y3) + sqr(z3)));
    if( fabs(d) < 1e-10) return -1;

    s = -(sqr(x2) + sqr(y2) + sqr(z2))*
	(x2*x3 - sqr(x3) + y2*y3 - sqr(y3) + z2*z3 - sqr(z3));
    s /= d;

    t = (sqr(x2) - x2*x3 + sqr(y2) - y2*y3 + z2*(z2 - z3))*
	(sqr(x3) + sqr(y3) + sqr(z3));
    t /= d;

    fprintf( stderr, "t=%f s=%f\n", t, s);

    // self check
    Vector3d c = (p2-p1)*t + (p3-p1)*s+p1;

    fprintf( stderr, "d(c,1) = %.10f\n", (c-p1).length());
    fprintf( stderr, "d(c,2) = %.10f\n", (c-p2).length());
    fprintf( stderr, "d(c,3) = %.10f\n", (c-p3).length());

    return 0;
}

void Model::insert_sink( long ind)
// inserts a sink into the element
{
    assert( ind >= 0 && ind < n_elements);
    WedgeElement * e = elements[ind];
    // name the nodes
    long n0 = e-> p[0];
    long n1 = e-> p[1];
    long n2 = e-> p[2];
    long n3 = e-> p[3];
    long n4 = e-> p[4];
    long n5 = e-> p[5];
    Node * n[6];
    for( long i = 0 ; i < 6 ; i ++) n[i] = nodes[e->p[i]];
    // name the points
    Vector3d p0(n[0]->ox,n[0]->oy,n[0]->oz);
    Vector3d p1(n[1]->ox,n[1]->oy,n[1]->oz);
    Vector3d p2(n[2]->ox,n[2]->oy,n[2]->oz);
    //	Vector3d p3(n[3]->ox,n[3]->oy,n[3]->oz);
    //	Vector3d p4(n[4]->ox,n[4]->oy,n[4]->oz);
    //	Vector3d p5(n[5]->ox,n[5]->oy,n[5]->oz);
    // calculate the circumcenter of the top face in local coordinates
    // (t,s), so that the center is defined as (p1-p0)*t + (p2-p0)*s
    double t, s;
    if( circumcenter( p0, p1, p2, t, s)) {
	fprintf( stderr,
	    "ERROR: Model::insert_sink(%ld) failed.\n"
	    "       -top face is degenerate\n", ind);
	return;
    }
	
    fprintf( stderr, "s=%f t=%f\n", s, t);
    // is the circumcenter inside the triangle?
    if( t < 0 || s < 0 || t+s > 1) {
	fprintf( stderr,
	    "WARNING: Model::insert_sink(%ld) sink is outsize\n"
	    "         of the triangle.\n", ind);
	return;
    }

    // remove the element
    remove_element( e); delete e; e = NULL;

    // insert a sink (one on the top face, one on the bottom face
    long tops = add_node( new Node(
			      n[1]->ox*t + n[0]->ox*(1-t-s) + n[2]->ox*s,
			      n[1]->oy*t + n[0]->oy*(1-t-s) + n[2]->oy*s,
			      n[1]->oz*t + n[0]->oz*(1-t-s) + n[2]->oz*s,
			      n[1]->px*t + n[0]->px*(1-t-s) + n[2]->px*s,
			      n[1]->py*t + n[0]->py*(1-t-s) + n[2]->py*s,
			      n[1]->pz*t + n[0]->pz*(1-t-s) + n[2]->pz*s,
			      0));
    long bots = add_node( new Node(
			      n[4]->ox*t + n[3]->ox*(1-t-s) + n[5]->ox*s,
			      n[4]->oy*t + n[3]->oy*(1-t-s) + n[5]->oy*s,
			      n[4]->oz*t + n[3]->oz*(1-t-s) + n[5]->oz*s,
			      n[4]->px*t + n[3]->px*(1-t-s) + n[5]->px*s,
			      n[4]->py*t + n[3]->py*(1-t-s) + n[5]->py*s,
			      n[4]->pz*t + n[3]->pz*(1-t-s) + n[5]->pz*s,
			      1));

    // insert 3 new elements
    add_element( new WedgeElement(
		     n0, n1, tops,
		     n3, n4, bots,
		     ym_map, pr_map, ys_map, ft_map));
    add_element( new WedgeElement(
		     n1, n2, tops,
		     n4, n5, bots,
		     ym_map, pr_map, ys_map, ft_map));
    add_element( new WedgeElement(
		     n2, n0, tops,
		     n5, n3, bots,
		     ym_map, pr_map, ys_map, ft_map));
				      
    delaunay_edge_flip( n0, n1);
    delaunay_edge_flip( n1, n2);
    delaunay_edge_flip( n2, n0);
    // delete the old element
    delete e;
}

void Model::subdivide_element( long ind)
{
    subdivide_element_bisection( ind);
}

int Model::refine_neighborhood( WedgeElement * e, double max_size)
// - all elements sharing nodes with 'e' will be checked to make sure that
//   their size is smaller than 'max_size'
// - if any element is bigger than 'max_size', the element is subdivided
// - returns: 0 - if all neighboring elements have size smaller than max_size
//            1 - if some elements had to be subdivided (more calls to this
//                function might be necessary to completely refine the
//                neighborhood)
{
    WedgeElement * we = e;
    long count = 0;
    // mark all elements that are touching 'we' and are bigger than
    // max_size
    for( long i = 0 ; i < n_elements ; i ++) {
	WedgeElement * te = elements[i];
	te-> mark = 0;
	// check if this element shares any nodes with 'e'
	if(te->p[0]==we->p[0]||te->p[0]==we->p[1]||te->p[0]==we->p[2]||
	    te->p[1]==we->p[0]||te->p[1]==we->p[1]||te->p[1]==we->p[2]||
	    te->p[2]==we->p[0]||te->p[2]==we->p[1]||te->p[2]==we->p[2])
	{
	    if( te-> get_size() > max_size)
	    {
		te-> mark = 1;
		count ++;
	    }
	}
    }
    // if no elements were marked, return 0
    if( count == 0) return 0;

    // subdivide all marked elements
    subdivide_marked_elements();
    return 1;
}

int Model::fracture_element( int count_mode)
// ----------------------------------------------------------------------
// - find the element with the highest ratio of principal stress over
//   yield stress
// - check if it exceeds yield stress
// - if yield stress is not exceeded, nothing happens
// - if the area of the element is bigger than some minimal area,
//   the element is subdivided and nothing else happens
// - if yield stress is exceeded, the element's stiffness matrix is
//   removed form the global stiffness matrix
// - the broken element is removed from memory
// - returns:
//     - if count_mode == 0
//         1 if something was broken
//         2 if nothing was broken, but something was subdivided
//         0 if nothing was broken (no element exceeded a yield stress)
//     - if count_mode == 1
//         returns number of elements exceeding yield stress
// 
// ----------------------------------------------------------------------
{
    // if the count mode is selected, just count the number of elements
    // that exceed yield stress
    if( count_mode) {
	long count = 0;
	for( long i = 0 ; i < n_elements ; i++)
	    if( elements[i]-> max_pstress >
		elements[i]-> yield_stress)
		count ++;
	return count;
    }

    // find the first non-broken element
    if( n_elements <= 0) return 0;
    long max_id = 0;

    // find an element with the highest principal stress / yield stress
    // ratio
    for( long i = max_id + 1 ; i < n_elements ; i ++) {
	//		if( elements[i]-> pstress3 > elements[max_id]-> pstress3)
	if( elements[i]-> max_pstress / elements[i]-> yield_stress >
	    elements[max_id]-> max_pstress /
	    elements[max_id]-> yield_stress)
	{
	    max_id = i;
	}
    }

    // check if this max. stress element should be broken (an element
    // should be broken if its maximum principal stress exceeds the
    // yield stress of the material from which it is built)
    if( elements[max_id]-> max_pstress <
	elements[max_id]-> yield_stress) return 0;

    // The element max_id has exceeded yield stress.  Check to make
    // sure that everything touching this triangle has area smaller
    // than max_break_size, and if not, subdivide the offender and
    // return '2'
    if( refine_neighborhood( elements[ max_id], max_break_size)) {
	return 2;
    }

    // replace element with 3 new elements
    WedgeElement * e = elements[ max_id];
    // remove the element from the list
    remove_element( e);
    // delete the element from memory
    delete e; e = NULL;
    n_fractures ++;
    return 1;
}

int Model::is_node_split_valid( long n0)
// ---------------------------------------------------------------------------
// - this function determines whether a split by the stress plane at this
//   node would result in split of the material
// - if the stress plane would result in no splitting, the function
//   returns 0, otherwise it returns 1
// ---------------------------------------------------------------------------
{
    // make sure this split would generate something above and something
    // below the plane - only then is the split valid
    int something_above = 0;
    int something_below = 0;

    Node & n = * nodes[ n0];
    for( long i = 0 ; i < n.n_refs ; i ++) {
	WedgeElement & e = * n.refs[i];
	assert( get_element_index( & e) != -1);
		
	// into n1 and n2 calculate the indices of the
	// other 2 nodes that belong to the element on the
	// same face
	long in0 = e.get_internal_index( n0);
	assert( in0 > -1);
	long n1 = -1; long n2 = -1;
	if( in0 < 3) {
	    // <n0> is on the top face of <e>
	    n1 = e.p[(in0+1)%3];
	    n2 = e.p[(in0+2)%3];
	} else {
	    n1 = e.p[(in0+1)%3+3];
	    n2 = e.p[(in0+2)%3+3];
	}
		
	Vector3d nv( nodes[n0]-> ps_nx,
	    nodes[n0]-> ps_ny,
	    nodes[n0]-> ps_nz);
	Vector3d p0( nodes[n0]-> ox,
	    nodes[n0]-> oy,
	    nodes[n0]-> oz);
	Vector3d p1( nodes[n1]-> ox,
	    nodes[n1]-> oy,
	    nodes[n1]-> oz);
	Vector3d p2( nodes[n2]-> ox,
	    nodes[n2]-> oy,
	    nodes[n2]-> oz);
		
	// calculate 'r' the isoparametric coordinate of the
	// intersection
	double r = split_edge_by_plane_raw( nv, p0, p1, p2);
	// and the point p3 at the split
	Vector3d p3 = p1 * (1-r) + p2 * r;
		
	// find out the angles
	double angle1 = Vector3d::angle( p1-p0,p3-p0);
	double angle2 = Vector3d::angle( p2-p0,p3-p0);
		
	// Determine whether the triangle is more above or
	// below.
	double a1 = 90 - Vector3d::angle( p1-p0, nv);
	double a2 = 90 - Vector3d::angle( p2-p0, nv);
	double above_angle = 0;
	if( a1 > above_angle) above_angle = a1;
	if( a2 > above_angle) above_angle = a2;
	double below_angle = 0;
	if( fabs(a1) > below_angle) below_angle = fabs(a1);
	if( fabs(a2) > below_angle) below_angle = fabs(a2);
	int more_below = above_angle < below_angle;
	double p01l = (p0-p1).length();
	double p02l = (p0-p2).length();
	double p03l = (p0-p3).length();
	double p12l = (p1-p2).length();
	double p13l = (p1-p3).length();
	double p23l = (p2-p3).length();
	int too_small =
	    (p01l < min_element_size) ||
	    (p02l < min_element_size) ||
	    (p03l < min_element_size) ||
	    (p12l < min_element_size) ||
	    (p13l < min_element_size) ||
	    (p23l < min_element_size);
		
		
	// should the whole triangle be considered above
	// or below?
	if( r <= 0 ||
	    r >= 1 ||
	    ! finite( r) ||
	    angle1 < min_split_angle ||
	    angle2 < min_split_angle ||
	    too_small)
	{
	    if( more_below)
		something_below = 1;
	    else
		something_above = 1;
	    continue;
	}

	// if the fracture plane goes too close to a
	// boundary node, we are avoiding back-cracking
	int n1_is_bn = ! is_node_moveable( n1);
	int n2_is_bn = ! is_node_moveable( n2);
	int n12_is_be = n1_is_bn && n2_is_bn;
	int is_too_close_to_bn =
	    (n1_is_bn && angle1<min_ftbn_angle) ||
	    (n2_is_bn && angle2<min_ftbn_angle);
	int is_too_close_to_be =
	    (n1_is_bn && angle1<min_ftbe_angle) ||
	    (n2_is_bn && angle2<min_ftbe_angle);
	if( (is_too_close_to_bn && ! n12_is_be)||
	    (is_too_close_to_be && n12_is_be))
	{
	    // if the element is more above the plane
	    if( more_below) {
		something_below = 1;
		continue;
	    } else {
		something_above = 1;
		continue;
	    }
	}

	// the plane splits the element
	something_above = 1;
	something_below = 1;
	break;
    }

    if( something_above && something_below)
	return 1;
    else
	return 0;
}

long Model::find_background_node( long ind)
// ---------------------------------------------------------------------------
// given a node, find (one of) its corresponding background node
// ---------------------------------------------------------------------------
{
    assert( ind >= 0 && ind < n_nodes && nodes[ind]-> is_fixed == 0);
    // go through all elements connected to the node
    for( long i = 0 ; i < nodes[ind]-> n_refs ; i ++) {
	WedgeElement & e = * nodes[ind]-> refs[i];
	if( e.p[0] == ind) return e.p[3];
	if( e.p[1] == ind) return e.p[4];
	if( e.p[2] == ind) return e.p[5];
    }
    assert( 0); // this should never happen
    return 0;
}

int Model::is_node_moveable( long n0)
// ----------------------------------------------------------------------
// determines whether a node can be freely repositioned without changing
// the the total size of the surrounded elements, taking into consideration
// the multilayer composition of the elements
//
// Returns: 1 - if it can be moved
//          0 - otherwise
//
// Method:
//    - make a list <nlist> of all nodes 'above' and 'below' <n0>,
//      including <n0> itself
//    - if every node in <nlist> is moveable with respect to its
//      own layer, return 1, otherwiser return 0
{
    assert( n0 >= 0 && n0 < n_nodes);

    Stack<long> nlist;
    nlist.push_unique( n0);
    // add all nodes 'above' <n0>
    long last_added = n0;
    while( 1) {
	int found_some = 0;
	Node & n = * nodes[last_added];
	for( long i = 0 ; i < n.n_refs ; i ++) {
	    WedgeElement & e = * n.refs[i];
	    long in = e.get_internal_index( last_added);
	    if( in < 3) continue;
	    found_some = 1;
	    last_added = e.p[in-3];
	    nlist.push( e.p[in-3]);
	    break;
	}
	if( ! found_some) break;
    }
    // add all nodes 'below' <n0>
    last_added = n0;
    while( 1) {
	int found_some = 0;
	Node & n = * nodes[last_added];
	for( long i = 0 ; i < n.n_refs ; i ++) {
	    WedgeElement & e = * n.refs[i];
	    long in = e.get_internal_index( last_added);
	    if( in >= 3) continue;
	    found_some = 1;
	    last_added = e.p[in+3];
	    nlist.push( e.p[in+3]);
	    break;
	}
	if( ! found_some) break;
    }

    // make sure every node in the nlist is surrounded
    while( ! nlist.is_empty())
	if( ! is_node_moveable_raw( nlist.pop()))
	    return 0;
    return 1;
	
}

int Model::is_node_moveable_raw( long n0)
// ----------------------------------------------------------------------
// Determines whether a node can be freely repositioned without changing
// the the total size of the surrounded elements, taking into consideration
// the multilayer composition of the elements.
//
// Returns: 1 - if it can be moved
//          0 - otherwise
//
// Method:
//    - initialize an empty list <top_list>
//    - initialize an empty list <bot_list>
//    - for each element <e> adjacent to the node <n0>:
//        - if node <n0> is on top face of <e>:
//             - let <n1> and <n2> be the other two nodes on the top face
//             - if <n1> is in <top_list>, remove it, otherwise add it
//             - if <n2> is in <top_list>, remove it, otherwise add it
//        - if node <n0> is on bottom face of <e>:
//             - let <n1> and <n2> be the other two nodes on the bottom face
//             - if <n1> is in <bot_list>, remove it, otherwise add it
//             - if <n2> is in <bot_list>, remove it, otherwise add it
//   - if both <top_list> and <bot_list> are empty, return 1
//     otherwise return 0
// ----------------------------------------------------------------------
{
    assert( n0 >= 0 && n0 < n_nodes);

    Stack<long> top_list;
    Stack<long> bot_list;
    Node & n = * nodes[n0];
    for( long i = 0 ; i < n.n_refs ; i ++) {
	WedgeElement & e = * n.refs[ i];
	long in0 = e.get_internal_index( n0);
	if( in0 < 3) {
	    // n0 is on top face of e
	    long n1 = e.p[(in0+1)%3];
	    long n2 = e.p[(in0+2)%3];
	    if( top_list.contains(n1))
		top_list.remove_val(n1);
	    else
		top_list.push(n1);
	    if( top_list.contains(n2))
		top_list.remove_val(n2);
	    else
		top_list.push(n2);
	} else {
	    // n0 is on bottom face of e
	    long n1 = e.p[(in0+1)%3+3];
	    long n2 = e.p[(in0+2)%3+3];
	    if( bot_list.contains(n1))
		bot_list.remove_val(n1);
	    else
		bot_list.push(n1);
	    if( bot_list.contains(n2))
		bot_list.remove_val(n2);
	    else
		bot_list.push(n2);
	}
    }

    // return result
    if( top_list.is_empty() && bot_list.is_empty())
	return 1;
    else
	return 0;
}

int Model::is_node_surrounded( long ind)
// ----------------------------------------------------------------------
// Determines whether the given node is surrounded by material on the
// surface. Returns 1 if yes, returns 0 otherwise.
// ----------------------------------------------------------------------
{
    assert( ind >= 0);
    assert( ind < n_nodes);

    // name the node for convenience
    Node & n = * nodes[ ind];
    // if this node is connected to less than 2 elements, it cannot be
    // surrounded by material
    if( n.n_refs < 2) return 0;
    // create a list of all nodes this node is connected to
    long n_list = 0;
    long * list = new long [ n.n_refs * 2];
    for( long i = 0 ; i < n.n_refs ; i ++) {
	WedgeElement & e = * n.refs[ i];
	// insert 2 nodes from 'e' other than 'ind'
	for( long j = 0 ; j < 3 ; j ++) {
	    if( e.p[j] == ind) continue;
	    assert( n_list < n.n_refs * 2);
	    list[ n_list] = e.p[j];
	    n_list ++;
	}
    }
    // internal check - n_list should be n.n_refs * 2
    assert( n_list == n.n_refs * 2);
    // make sure that each node is represented in the list exactly twice
    for( long i = 0 ; i < n_list ; i ++) {
	if( list[i] == -1) continue;
	long count = 0;
	long x = list[i];
	for( long j = i ; j < n_list ; j ++) {
	    if( x != list[j]) continue;
	    list[j] = -1;
	    count ++;
	}
	if( count != 2) {
	    delete [] list;
	    return 0;
	}
    }

    delete [] list;

    // test passed = return 1
    return 1;
}

Model::ElementList Model::get_elements_sharing_edge( Edge & edge)
// ----------------------------------------------------------------------
// Returns a list of elements sharing edge <edge>
// ----------------------------------------------------------------------
{
    // preconditions
    assert( edge.n1 >= 0 && edge.n1 < n_nodes);
    assert( edge.n2 >= 0 && edge.n2 < n_nodes);
    assert( edge.n1 != edge.n2);

    // here we'll store the result
    ElementList list;

    // process all elements containing node <n1>
    for( long i = 0 ; i < nodes[edge.n1]-> n_refs ; i ++) {
	WedgeElement * e = nodes[edge.n1]-> refs[i];
	long in1 = e-> get_internal_index( edge.n1);
	assert(in1 >= 0);
	long in2 = e-> get_internal_index( edge.n2);
	if( in2 < 0) continue;
	if( (in1 <  3 && in2 <  3) ||
	    (in1 >= 3 && in2 >= 3))
	    list.push_unique( e);
    }

    // 	fprintf( stderr, "elements sharing edge %ld,%ld:\n", edge.n1, edge.n2);
    // 	for( long i = 0 ; i < list.n ; i ++)
    // 		fprintf( stderr, "%2ld) %p\n", i, list(i));
    return list;
	
}

void Model::collapse_edge( long n1, long n2)
// ----------------------------------------------------------------------
// collapses an edge <n1,n2> (merging tn1 into tn2), while considering
// the multi-layer structure of the model.
//
// Algorithm:
//     All elements sharing the edge <n1,n2> will be removed, and all other
//     elements referencing <n1> will have that reference changed to <n2>.
// ----------------------------------------------------------------------
{
    //	fprintf( stderr, "Collapse edge %ld %ld\n", n1, n2);
    // preconditions:
    assert( n1 >= 0);
    assert( n1 < n_nodes);
    assert( n2 >= 0);
    assert( n2 < n_nodes);

    // prepare a list of edges 'above' and 'below' tn1,tn2
    Stack<Edge> e_list; // the final list of edges
    Stack<Edge> ue_list; // the unprocessed edge list
    ue_list.push( Edge(n1, n2));
    while( ! ue_list.is_empty()) {
	// get an uprocessed edge
	Edge e = ue_list.pop();
	// push it onto the final list stack
	e_list.push_unique( e);
	// get a list of elements sharing the edge <e>
	ElementList el_list = get_elements_sharing_edge( e);
	// 		fprintf( stderr, "el_list:\n");
	// 		for( long i = 0 ; i < el_list.n ; i ++)
	// 			fprintf( stderr, "    %p\n", el_list(i));
	// from each of these elements, extract all edges above or
	// below <e>, and if they are not in <e_list>, add them
	// into <ue_list>
	while( ! el_list.is_empty()) {
	    WedgeElement * el = el_list.pop();
	    // find the internal indices in element <el>
	    long in1 = el-> get_internal_index( e.n1);
	    long in2 = el-> get_internal_index( e.n2);
	    assert( in1 >= 0 && in2 >= 0);
	    if( in1 < 3 && in2 >= 3) continue;
	    if( in2 < 3 && in1 >= 3) continue;
	    // extract the edge 'below' or 'above'
	    Edge ne;
	    if( in1 < 3) {
		ne.n1 = el-> p[in1+3];
		ne.n2 = el-> p[in2+3];
	    } else {
		ne.n1 = el-> p[in1-3];
		ne.n2 = el-> p[in2-3];
	    }
	    // put <ne> into unprocessed list, provided it is
	    // not already in the final list
	    if( ! e_list.contains( ne))
		ue_list.push_unique( ne);
	}
    }

    // unmark all elements
    for( long i = 0 ; i < n_elements ; i ++)
	elements[i]-> selected = 0;

    // We now have in <e_lsit> a list of edges to be collapsed. Since we
    // are removing <n1>, process all elements touching <n1>.
    //
    // Start by creating a list of all elements to be processed: <el>.
    // Mark all elements which share a whole edge.
    ElementList el;
    for( long i = 0 ; i < e_list.n ; i ++) {
	Node & n = * nodes[e_list(i).n1];
	for( long j = 0 ; j < n.n_refs ; j ++) {
	    WedgeElement * e = n.refs[j];
	    if( e-> get_internal_index( e_list(i).n2) > -1)
		e-> selected = 1;
	    el.push_unique( e);
	}
    }
    // Prepare a translation table.
    long tt[n_nodes];
    for( long i = 0 ; i < n_nodes ; i ++)
	tt[i] = i;
    for( long i = 0 ; i < e_list.n ; i ++)
	tt[e_list(i).n1] = e_list(i).n2;

    // Process all elements on <el>
    while( ! el.is_empty()) {
	WedgeElement * e = el.pop();
	// remove the element
	remove_element( e);
	// add a new element
	if( ! e-> selected)
	    add_element( new WedgeElement(
			     tt[e-> p[0]], tt[e-> p[1]], tt[e-> p[2]],
			     tt[e-> p[3]], tt[e-> p[4]], tt[e-> p[5]],
			     ym_map, pr_map, ys_map, ft_map));
	delete e;
    }
}

int
    Model::refine_mesh( void)
{
    if( break_method == BreakElement)
	return _refine_mesh_max_stress();
    else
	return _refine_mesh_von_mises();
}

int
    Model::_refine_mesh_max_stress( void)
// ----------------------------------------------------------------------
// refine mesh so that elements whose stress is near yield stress are
// small
//
// - consider each element
// - skip element if it is 'too' small already
// - if element's max. pstress is close to its yield stress,
//   mark the element for refinement
// - return: 0 = no elements needed to be refined
//           1 = an element needed to be refined
//
// ----------------------------------------------------------------------
{
    long mark_count = 0;
    for( long i = 0 ; i < n_elements ; i ++) {
	WedgeElement * e = elements[i];
	e-> mark = 0;
	// skip all elements that are too small to be subdivided
	if( e-> get_size() < min_refine_size) continue;
	// calculate the ratio of this element's max. pstress
	// and its yield stress
	double ratio = e-> max_pstress / e-> yield_stress;
	if( ratio < 0.999) continue;
	e-> mark = 1;
	mark_count ++;
    }
    fprintf( stderr, "\t- %ld marked\n",
	mark_count);
	
    // if no elements were marked for subdivision, return 0
    if( mark_count == 0) return 0;
	
    // otherwise subdivide all elements marked
    subdivide_marked_elements();
    return 1;
}

int
    Model::_refine_mesh_von_mises( void)
// ----------------------------------------------------------------------
// adaptively refines elements to minimize the difference in calculated
// von mises stresses at the nodes
//
// Returns: 0 if nothing was subdivided
//          1 if something was subdivided (a relaxation should follow,
//            followed by another call to this function)
//
// Algorithm:
//   - unmark all elements
//   - for each node:
//       - calculate extrapolated stress tensor at the node from each
//         element the node belongs to
//       - from each tensor, calculate von-mises-stress
//       - calculate the difference between min and max obtained
//       - calculate the ratio of the difference to the max.
//         yield stress
//       - if the ratio is larger than some predefined threshold, mark
//         all elements this node belongs to
//   - subdivide all marked elements
//
// Optimization:
//   - use the knowledge from local relaxation - the modified nodes
//     are all marked
//   - therefore, there is no need to check nodes which are not connected
//     to at least one marked node
//   - also, there is no need to check fixed nodes, and nodes with zero
//     references
// ----------------------------------------------------------------------
{
    // unmark all elements
    for( long i = 0 ; i < n_elements ; i ++)
	elements[i]-> mark = 0;

    double max_ratio = -1;
    long n_marked = 0;
    for( long i = 0 ; i < n_nodes ; i ++) {
	// for each node, calculate the von misses error
	Node & node = * nodes[i];
	// skip fixed nodes
	if( node.is_fixed || node.n_refs == 0) continue;
	// see if this node is connected to an unmarked node
	int is_connected_to_marked_node = 0;
	for( long j = 0 ; j < node.n_refs ; j ++) {
	    WedgeElement * e = node.refs[j];
	    if( e-> has_marked_node()) {
		is_connected_to_marked_node = 1;
		break;
	    }
	}
	if( ! is_connected_to_marked_node) continue;
	// OK, the node is connected to a marked node, we have to
	// do the calculations
	int has_min = 0; double min_vm = 0;
	int has_max = 0; double max_vm = 0;
	// for all elements belonging to this node
	for( long j = 0 ; j < node.n_refs ; j ++) {
	    WedgeElement & e = * node.refs[j];
	    // ask the element to calculate the stress at this
	    // node
	    StressTensor s = e.get_nodal_stress_global_ind( i);
	    // get eigenvalues
	    double s1, s2, s3;
	    s.get_eigenvalues( s1, s2, s3);
	    // calculate von mises stress
	    double vm = sqrt( ((s1-s2) * (s1-s2) +
		    (s1-s3) * (s1-s3) +
		    (s3-s2) * (s3-s2)) / 2.0);
	    // record min & max values
	    if( ! has_min || vm < min_vm) {
		has_min = 1;
		min_vm = vm;
	    }
	    if( ! has_max || vm > max_vm) {
		has_max = 1;
		max_vm = vm;
	    }
	}
	// calculate the difference
	double diff_vm = max_vm - min_vm;
	// calculate the ratio
	double ratio = diff_vm / ys_map-> get_max();
	// if the ratio is bigger than some threshold, mark
	// all elements connected to this node for subdivision
	if( ratio > max_ratio) max_ratio = ratio;
	if( ratio > 0.1)
	    for( long j = 0 ; j < node.n_refs ; j ++) {
		WedgeElement & e =
		    * (node.refs[j]);
		if( e.get_size() > min_refine_size) {
		    e.mark = 1;
		    n_marked ++;
		}
	    }
    }
    fprintf( stderr, "\t- vm: max_ratio = %f marked = %ld\n",
	max_ratio, n_marked);

    // subdivide all marked elements
    if( n_marked > 0)
	subdivide_marked_elements();
    return n_marked > 0;
}

long Model::get_neighbor_count( WedgeElement * e)
// ----------------------------------------------------------------------
// returns number of neighbors this element has on the top face
// ----------------------------------------------------------------------
{
    long count = 0;
    for( long i = 0 ; i < 3 ; i ++) {
	long n1 = e-> p[i];
	long n2 = e-> p[(i+1)%3];
	if( NULL != find_element( e, n1, n2))
	    count ++;
    }
    return count;
}


Model::NodeList *
    Model::get_face_neighbors( long n0)
// ----------------------------------------------------------------------
// finds and returns a list of all nodes connected to this node
// through top or bottom faces of elements.
// ----------------------------------------------------------------------
{
    assert( n0 >= 0 && n0 < n_nodes);

    // name the node for convenience
    Node & n = * nodes[n0];

    NodeList * res = new NodeList();
    for( long i = 0 ; i < n.n_refs ; i ++) {
	WedgeElement & e = * n.refs[i];
	long in0 = e.get_internal_index( n0);
	assert( in0 != -1);
	if( in0 < 3) {
	    res-> push_unique( e.p[(in0+1)%3]);
	    res-> push_unique( e.p[(in0+2)%3]);
	} else {
	    res-> push_unique( e.p[(in0+1)%3+3]);
	    res-> push_unique( e.p[(in0+2)%3+3]);
	}
    }
    return res;
}

Model::NodeList *
    Model::get_face_neighbors_circle( long n0)
// ----------------------------------------------------------------------
// finds and returns a circular list of all nodes connected to this node
// through top or bottom faces of elements. The nodes on the list are all
// circularly connected through face edges.
// ----------------------------------------------------------------------
{
    assert( n0 >= 0 && n0 < n_nodes);

    // name the node for convenience
    //	Node & n = * nodes[n0];
    // make an empty list
    NodeList * res = new NodeList();
    // get the list of all nodes connected to node <n0> through faces
    NodeList * nlist = get_face_neighbors( n0);
    // if this node has no neighbors, return an empty list
    if( nlist-> is_empty()) {
	delete nlist;
	return res;
    }
    // put the first node into the result
    long last = nlist-> pop();
    res-> push( last);

    while( 1) {
	// find a node <n1> in nlist, such that a face edge
	// <last>,<n1> exists
	//		long n1 = -1;
		
    }
    while( ! nlist-> is_empty()) {
	res-> push( nlist-> pop());
    }

    return res;
}

class DoubleEdge {
 public:
    long n0, n1, n2; // edge n0,n1 and n1,n2

    DoubleEdge( void) {
	n0 = -1;
	n1 = -2;
	n2 = -3;
    }

    // calculate the best position for node <n>
    Vector3d get_opt_pos( long n,
	const Vector3d coord[],
	const Plane & plane) const
    {
	Vector3d  nv( coord[ n].x,coord[ n].y,coord[ n].z);
	Vector3d n0v( coord[n0].x,coord[n0].y,coord[n0].z);
	Vector3d n1v( coord[n1].x,coord[n1].y,coord[n1].z);
	Vector3d n2v( coord[n2].x,coord[n2].y,coord[n2].z);

	// project all 4 points onto the plane of <n>
	plane.project_onto( nv);
	plane.project_onto( n0v);
	plane.project_onto( n1v);
	plane.project_onto( n2v);

	// calculate vectors
	Vector3d v1 = n0v - n1v; v1.normalize();
	Vector3d v2 = n2v - n1v; v2.normalize();
	Vector3d solv = v1 + v2;

	if( solv.length() <= 1e-10) {
	    // n0, n1 and n2 on the same line - special case
	    // calculate normal of triangle n0, n2 and n
	    Vector3d v12 = n0v - n1v;
	    Vector3d norm = cross_product( nv-n0v, v12);
	    if( norm.length() < 1e-10) {
		// special case: n is also on the line,
		// so don't do anything at all
		return nv;
	    } else {
		// get cross product of norm and v12
		solv = cross_product( norm, v12);
	    }
	}
	if( solv.length() <= 1e-10) {
	    // if the sol. vector is still too small, do nothing
	    return nv;
	}
	solv.normalize();

	Vector3d oldv = nv - n1v;
	solv = solv * oldv.length();
	Vector3d p1 = n1v + solv;
	Vector3d p2 = n1v - solv;
	// as the desired position, pick the one closer to <n>
	if( (p1-nv).length() < (p2-nv).length())
	    return p1;
	else
	    return p2;
    }
};

class SmoothNode {
 public:
    long n;
    Stack<DoubleEdge> de;
    Vector3d opos;
    Plane plane;

    void calc_opt_pos( Vector3d coord[]) {
	// get optimum position for each double edge, add it into
	// the sum, then average these
	Vector3d sum;
	for( long i = 0 ; i < de.n ; i ++) {
	    Vector3d pos = de(i).get_opt_pos( n, coord, plane);
	    sum = sum + pos;
	}
	if( de.n > 0)
	    opos = sum / de.n;
	else
	    opos.set( coord[n].x,
		coord[n].y,
		coord[n].z);
    }

    double calc_max_angle_diff( Vector3d coord[])
    // --------------------------------------------------------------
    // get the max. difference
    // --------------------------------------------------------------
    {
	// for each double edge, calculate the angle between
	// the current position of <n> and the desired position,
	// and return the maximum
	double max_angle = 0;
	for( long i = 0 ; i < de.n ; i ++) {
	    Vector3d curr_pos = coord[ n];
	    Vector3d desired_pos = de(i).get_opt_pos(
		n, coord, plane);
	    Vector3d curr_v = curr_pos - coord[de(i).n1];
	    Vector3d desired_v = desired_pos - coord[de(i).n1];
	    double angle = Vector3d::get_angle( curr_v, desired_v);
	    if( angle > max_angle) max_angle = angle;
	}
	return max_angle;
    }
};

class Triangle {
 public:
    long n0, n1, n2;
    Triangle( void) {
	n0 = n1 = n2 = -1;
    }
    Triangle( long p_n0, long p_n1, long p_n2) {
	n0 = p_n0;
	n1 = p_n1;
	n2 = p_n2;
    }
    int operator == ( const Triangle & t) {
	if( (n0 == t.n0) && 
	    ( (n1 == t.n1 && n2 == t.n2) ||
		(n1 == t.n2 && n2 == t.n1)))
	    return 1;
	else
	    return 0;
    }
    void get_opt_pos( Vector3d p[3], Vector3d coord[]) {
	Vector3d v01 = coord[n1] - coord[n0];
	Vector3d v12 = coord[n2] - coord[n1];
	Vector3d v20 = coord[n0] - coord[n2];
	double l01 = v01.length();
	double l12 = v12.length();
	double l20 = v20.length();
	double lavg = (l01+l12+l20) / 3;
	double t01 = -0.5 * (lavg / l01 - 1);
	double t12 = -0.5 * (lavg / l12 - 1);
	double t20 = -0.5 * (lavg / l20 - 1);
	// new points
	p[0] = coord[n0] + (v01 * t01 - v20 * t20) / 2;
	p[1] = coord[n1] + (v12 * t12 - v01 * t01) / 2;
	p[2] = coord[n2] + (v20 * t20 - v12 * t12) / 2;
    }
};

int
    Model::repell_smooth_mesh( void)
// ----------------------------------------------------------------------
// smooth mesh (the selected nodes) by particle repelling
//   - each node will attempt to move into its perfect position with
//     respect to every triangle it belongs to (every triangle tries
//     to achieve 60 degree angles)
// ----------------------------------------------------------------------
{
    // make a list of the nodes to be adjusted
    NodeList nlist;
    for( long i = 0 ; i < n_nodes ; i ++)
	if( nodes[i]-> selected && nodes[i]-> n_refs > 0)
	    nlist.push(i);
    if( nlist.n == 0) return 0;
	
    fprintf( stderr, "Repelling %ld nodes.\n", nlist.n);
	
    // copy all nodal coordinates into a temporary array
    Vector3d coord[ n_nodes];
    for( long i = 0 ; i < n_nodes ; i ++)
	coord[i].set( nodes[i]-> ox, nodes[i]-> oy, nodes[i]-> oz);
	
    // construct a list of triangles connected to the selected nodes
    Stack <Triangle> tlist;
    for( long i = 0 ; i < n_elements ; i ++) {
	WedgeElement * e = elements[i];
	if( nodes[ e-> p[0]]-> selected ||
	    nodes[ e-> p[1]]-> selected ||
	    nodes[ e-> p[2]]-> selected)
	    tlist.push_unique(
		Triangle( e-> p[0], e-> p[1], e-> p[2]));
	if( nodes[ e-> p[3]]-> selected ||
	    nodes[ e-> p[4]]-> selected ||
	    nodes[ e-> p[5]]-> selected)
	    tlist.push_unique(
		Triangle( e-> p[3], e-> p[4], e-> p[5]));
    }

    // relaxation loop:
    long n_it = 10;
    for( long it = 0 ; it < n_it ; it ++) {
	fprintf( stderr, "repell iteration %ld\n", it);
	// calculate the plane of each node on the list <nlist>
	Plane planes[ n_nodes];
	for( long i = 0 ; i < nlist.n ; i ++) {
	    long n = nlist(i);
	    planes[n].norm = get_face_normal( n);
	    planes[n].pt = Vector3d(
		nodes[ n]-> ox,
		nodes[ n]-> oy,
		nodes[ n]-> oz);
	}
	// calculate optimum positions of all nodes of all triangles
	// on tlist
	Vector3d pos_sum[ n_nodes];
	long pos_n[ n_nodes];
	for( long i = 0 ; i < n_nodes ; i ++) {
	    pos_sum[i].set(0,0,0);
	    pos_n[i] = 0;
	}
	for( long i = 0 ; i < tlist.n ; i ++) {
	    Vector3d p[3];
	    tlist(i).get_opt_pos( p, coord);
	    // add these positions to the accumulation arrays
	    pos_sum[ tlist(i).n0] = pos_sum[ tlist(i).n0] + p[0];
	    pos_n[ tlist(i).n0] ++;
	    pos_sum[ tlist(i).n1] = pos_sum[ tlist(i).n1] + p[1];
	    pos_n[ tlist(i).n1] ++;
	    pos_sum[ tlist(i).n2] = pos_sum[ tlist(i).n2] + p[2];
	    pos_n[ tlist(i).n2] ++;
	}
	// calculate the averages (only for the selecte nodes), project
	// these onto the planes of the node and put them into coord[]
	for( long i = 0 ; i < nlist.n ; i ++) {
	    long n = nlist(i);
	    if( pos_n[n] > 0)
		coord[n] = pos_sum[n] / pos_n[n];
	    planes[n].project_onto( coord[n]);
	    nodes[n]-> ox = coord[n].x;
	    nodes[n]-> oy = coord[n].y;
	    nodes[n]-> oz = coord[n].z;
	}
    }

    // apply growth to every node modified
    for( long i = 0 ; i < nlist.n ; i ++)
	nodes[ nlist(i)]-> grow( time_curr,
	    growth_x,
	    growth_y,
	    growth_z);

    // remove and insert each of these elements, forcing them to
    // re-calculate themselves
    Stack <WedgeElement *> elist;
    for( long i = 0 ; i < nlist.n ; i ++) {
	Node & n = * nodes[ nlist(i)];
	for( long j = 0 ; j < n.n_refs ; j ++)
	    elist.push_unique( n.refs[j]);
    }
    while( ! elist.is_empty()) {
	WedgeElement * e = elist.pop();
	remove_element( e);
	add_element( new WedgeElement(
			 e-> p[0], e-> p[1], e-> p[2],
			 e-> p[3], e-> p[4], e-> p[5],
			 ym_map, pr_map, ys_map, ft_map));
	delete e;
    }

    return 1;
}

int
    Model::angle_smooth_mesh( void)
// ----------------------------------------------------------------------
// Attempts to improve the mesh using angle-smoothing procedure.
//
// Returns: 1 - if anything was adjusted
//          0 - if nothing was adjusted
//
// Algorithm:
//     - select all nodes that are 'moveable' and adjacent to degenerate
//       wedges
//     - if no nodes selected, return 0
//     - iteratively recalculate positions of the selected nodes using
//       angle smooting (torsion springs)
//     - adjust all elements touching these nodes (delete & add)
//     - return 1
//
// Iterative algorithm:
//     - for each node on nlist, construct a list of double angles
//       affecting the node
//         - each double angle is represented with:
//             - center node
//             - up normal at center node
//             - left node
//             - right node
//     - for each node on nlist, construct a total displacement:
//         - displ. sum
//         - number of displacements
//     - repeat:
//         - for each node <n> on <nlist>
//             - zero out its adjustment structure
//             - for each double angle associated with node <n>
//                 - add displacement due to this double angle to node's
//                   total displacement structure
//         - for each node <n> on <nlist>
//             - set the node's coordinates according to its total
//               displacement structure, remembering the total change
//         - quit when no 'improvement' is 'noticed'
//
// TBD:
//   - what is 'improvement'
//   - what is 'noticed'
// ----------------------------------------------------------------------
{
    //	fprintf( stderr, "ANGLE SMOOTH MESH\n");
    //	return 0;

    // make a list of the nodes to be adjusted
    NodeList nlist;
    for( long i = 0 ; i < n_nodes ; i ++)
	if( nodes[i]-> selected && nodes[i]-> n_refs > 0)
	    nlist.push(i);
    //	fprintf( stderr, "nlist.n=%ld\n", nlist.n);
    if( nlist.n == 0) return 0;

    //	fprintf( stderr, "Smoothing %ld nodes.\n", nlist.n);

    // copy all nodal coordinates into a temporary array
    Vector3d coord[ n_nodes];
    for( long i = 0 ; i < n_nodes ; i ++)
	coord[i].set( nodes[i]-> ox, nodes[i]-> oy, nodes[i]-> oz);

    // for all nodes in <nlist> create a structure containing:
    //     - the node itself
    //     - a list of double edges affecting the node
    //     - the plane in which the node exists
    Stack <SmoothNode> sn_list;
    for( long i = 0 ; i < nlist.n ; i ++) {
	long n0 = nlist.array[i];
	SmoothNode sn;
	sn.n = n0;
	// get a list of all nodes connected to <n0>
	NodeList * lst = get_face_neighbors( n0);
	for( long j = 0 ; j < lst-> n ; j ++) {
	    long nj = (*lst)(j);
	    for( long k = 0 ; k < lst-> n ; k ++) {
		if( j == k) continue;
		long nk = (*lst)(k);
		if( ! face_edge_exists( nj, nk)) continue;
		for( long l = 0 ; l < lst-> n ; l ++) {
		    if( l == j || l == k) continue;
		    long nl = (*lst)(l);
		    if( nl > nk) continue;
		    if( face_edge_exists(nj, nl)) {
			DoubleEdge de;
			de.n0 = nk;
			de.n1 = nj;
			de.n2 = nl;
			sn.de.push( de);
		    }
		}
	    }
	}
	delete lst;
	sn_list.push( sn);
    }

    assert( nlist.n == sn_list.n);

    // calculate planes for every node on the <sn_list>
    Plane planes[n_nodes];
    for( long i = 0 ; i < sn_list.n ; i ++) {
	long n = sn_list(i).n;
	planes[n].norm = get_face_normal( n);
	planes[n].pt = Vector3d( nodes[n]-> ox,
	    nodes[n]-> oy,
	    nodes[n]-> oz);
    }

    // relaxation loop:
    long n_it = 10;
    for( long it = 0 ; it < n_it ; it ++) {
	//		fprintf( stderr, "smooth iteration %ld\n", it);
	// calculate the plane of each node on the list
	for( long i = 0 ; i < sn_list.n ; i ++) {
	    long n = sn_list(i).n;
	    sn_list(i).plane.norm = get_face_normal( n);
	    sn_list(i).plane.pt = Vector3d(
		nodes[ n]-> ox,
		nodes[ n]-> oy,
		nodes[ n]-> oz);
	}
	// for each node on <sn_list>, calculate its optimum
	// position
	for( long i = 0 ; i < sn_list.n ; i ++) {
	    SmoothNode & sn = sn_list(i);
	    sn.calc_opt_pos( coord);
	    coord[sn.n] = sn.opos;
	    assert( ! isnan( coord[sn.n].x));
	    assert( ! isnan( coord[sn.n].y));
	    assert( ! isnan( coord[sn.n].z));
	    //			sn_list.array[i].calc_opt_pos( coord);
	    // 			coord[sn_list.array[i].n] = sn_list.array[i].opos;
	}

	/*		
	// move each node on <sn_list>
	for( long i = 0 ; i < sn_list.n ; i ++)
	coord[sn_list.array[i].n] = sn_list.array[i].opos;
	*/

	// adjust the positions of nodes that are in <sn_list>,
	for( long i = 0 ; i < sn_list.n ; i ++) {
	    long ni = sn_list.array[i].n;
	    Node & n = * nodes[ ni];
	    Vector3d pn = coord[ ni];
	    // project pn onto the original plane
	    planes[ ni].project_onto( pn);
	    n.ox = pn.x; n.oy = pn.y; n.oz = pn.z;
	}
	// figure out the max. diff angle
	// 		double max_angle = 0;
	// 		for( long i = 0 ; i < sn_list.n ; i ++) {
	// 			double angle = sn_list(i).calc_max_angle_diff( coord);
	// 			if( angle > max_angle) max_angle = angle;
	// 		}
	// 		fprintf( stderr, "max. diff angle = %f\n", max_angle);
    }

    // apply growth to every node modified
    for( long i = 0 ; i < sn_list.n ; i ++)
	nodes[sn_list.array[i].n]-> grow( time_curr,
	    growth_x,
	    growth_y,
	    growth_z);

    // remove and insert each of these elements, forcing them to
    // re-calculate themselves
    Stack <WedgeElement *> elist;
    for( long i = 0 ; i < sn_list.n ; i ++) {
	Node & n = * nodes[ sn_list.array[i].n];
	for( long j = 0 ; j < n.n_refs ; j ++)
	    elist.push_unique( n.refs[j]);
    }
    /*
      for( long i = 0 ; i < elist.n ; i ++) {
      WedgeElement * e = elist.array[i];
      remove_element( e);
      add_element( e);
      elist.array[i] = NULL;
      }
    */
    while( ! elist.is_empty()) {
	WedgeElement * e = elist.pop();
	remove_element( e);
	add_element( new WedgeElement(
			 e-> p[0], e-> p[1], e-> p[2],
			 e-> p[3], e-> p[4], e-> p[5],
			 ym_map, pr_map, ys_map, ft_map));
	delete e;
    }

    return 1;
}

Model::NodeList Model::get_node_column( long ind)
// ----------------------------------------------------------------------
// return a list of all nodes above and below <ind>
// ----------------------------------------------------------------------
{
    // create a list of unprocessed nodes
    NodeList ul;
    ul.push( ind);
    // create a result list
    NodeList res;
    // process all unprocessed nodes
    while( ! ul.is_empty()) {
	// get an uprocessed node
	long n = ul.pop();
	// if it was already processed (it is in the result)
	// don't process it again
	if( res.contains( n)) continue;
	// otherwise put this node into the result
	res.push( n);
	// for all elements connected to this node, find
	// all nodes above and below the node
	for( long i = 0 ; i < nodes[n]-> n_refs ; i ++) {
	    WedgeElement * e = nodes[n]-> refs[i];
	    long in = e-> get_internal_index(n);
	    assert( in > -1);
	    if( in < 3)
		ul.push_unique( e-> p[in + 3]);
	    else
		ul.push_unique( e-> p[in - 3]);
	}
    }

    // return the result
    return res;
}

int
    Model::smooth_mesh_around_crack_tips( void)
// ----------------------------------------------------------------------
// attempts to smooth the mesh around crack tips
// ----------------------------------------------------------------------
{
    // unselecte all nodes
    for( long i = 0 ; i < n_nodes ; i ++)
	nodes[i]-> selected = 0;

    // select all crack tip nodes
    for( long i = 0 ; i < ftips.n ; i ++)
	nodes[ ftips.array[i]]-> selected = 1;
    // 	{
    // 		long ns = 0;
    // 		for( long i = 0 ; i < n_nodes ; i ++)
    // 			if( nodes[i]-> selected) ns ++;
    // 		fprintf( stderr, "ns = %ld\n", ns);
    // 	}

    // mark all surrounded nodes as well
    select_nodes( 3);
	
    // 	{
    // 		long ns = 0;
    // 		for( long i = 0 ; i < n_nodes ; i ++)
    // 			if( nodes[i]-> selected) ns ++;
    // 		fprintf( stderr, "ns = %ld\n", ns);
    // 	}

    // for each selected node, also select all nodes above and below it
    char sel[ n_nodes];
    for( long i = 0 ; i < n_nodes ; i ++)
	sel[i] = 0;
    for( long i = 0 ; i < n_nodes ; i ++) {
	if( ! nodes[i]-> selected) continue;
	NodeList lst = get_node_column( i);
	while( ! lst.is_empty())
	    sel[ lst.pop()] = 1;
    }
    for( long i = 0 ; i < n_nodes ; i ++)
	nodes[i]-> selected = sel[i];

    // unmark the nodes which are not moveable
    for( long i = 0 ; i < n_nodes ; i ++)
	if( nodes[i]-> selected)
	    if( ! is_node_moveable(i))
		nodes[i]-> selected = 0;

    // apply angle smooting to the selected nodes
    int res = angle_smooth_mesh();

    // return the result
    return res;
}

int
    Model::smooth_mesh_around_nodes( const NodeList & nlist, long depth)
// ----------------------------------------------------------------------
// attempts to smooth the mesh around nodes on <nlist>
// ----------------------------------------------------------------------
{
    // unselecte all nodes
    for( long i = 0 ; i < n_nodes ; i ++)
	nodes[i]-> selected = 0;

    // select all nodes on the <nlist>
    for( long i = 0 ; i < nlist.n ; i ++)
	nodes[ nlist.array[i]]-> selected = 1;
    // 	{
    // 		long ns = 0;
    // 		for( long i = 0 ; i < n_nodes ; i ++)
    // 			if( nodes[i]-> selected) ns ++;
    // 		fprintf( stderr, "ns = %ld\n", ns);
    // 	}

    // mark all surrounded nodes as well
    select_nodes( depth);
	
    // 	{
    // 		long ns = 0;
    // 		for( long i = 0 ; i < n_nodes ; i ++)
    // 			if( nodes[i]-> selected) ns ++;
    // 		fprintf( stderr, "ns = %ld\n", ns);
    // 	}

    // unmark the nodes which are not moveable
    for( long i = 0 ; i < n_nodes ; i ++)
	if( nodes[i]-> selected)
	    if( ! is_node_moveable(i))
		nodes[i]-> selected = 0;

    // for each selected node, also select all nodes above and below it,
    // but only if they are not moveable
    char sel[ n_nodes];
    for( long i = 0 ; i < n_nodes ; i ++)
	sel[i] = 0;
    for( long i = 0 ; i < n_nodes ; i ++) {
	if( ! nodes[i]-> selected) continue;
	NodeList lst = get_node_column( i);
	while( ! lst.is_empty()) {
	    long ind = lst.pop();
	    if( is_node_moveable( ind))
		sel[ ind] = 1;
	}
    }
    for( long i = 0 ; i < n_nodes ; i ++)
	nodes[i]-> selected = sel[i];

    // unmark the nodes which are not moveable
    for( long i = 0 ; i < n_nodes ; i ++)
	if( nodes[i]-> selected)
	    if( ! is_node_moveable(i))
		nodes[i]-> selected = 0;

    // apply angle smooting to the selected nodes
    int res = angle_smooth_mesh();

    // return the result
    return res;
}

int
    Model::refine_mesh_around_node( long ind, long depth, double req_size)
// ----------------------------------------------------------------------
// makes sure that all nodes in the vicinity of node <ind> are smaller
// than <req_size>. If not, such elements are subdivided.
//
// The neighborhood is defined using parameter 'depth'. For depth = 0,
// no elements are processed. For depth 1, only elements immediately
// touching node <ind> are processed. For depth 2, elements touching
// node <n1> indirectly through another element are selected, etc.
//
// Returns: 0 - if nothing had to be subdivided
//          1 - if something was subdivided
// ----------------------------------------------------------------------
{
    // unselect all nodes
    for( long i = 0 ; i < n_nodes ; i ++)
	nodes[i]-> selected = 0;

    // select the requested node
    nodes[ind]-> selected = 1;

    // mark all surrounded nodes as well
    select_nodes( depth);
	
    // for each selected node, also select all nodes above and below it
    char sel[ n_nodes];
    for( long i = 0 ; i < n_nodes ; i ++)
	sel[i] = 0;
    for( long i = 0 ; i < n_nodes ; i ++) {
	if( ! nodes[i]-> selected) continue;
	NodeList lst = get_node_column( i);
	while( ! lst.is_empty())
	    sel[ lst.pop()] = 1;
    }
    for( long i = 0 ; i < n_nodes ; i ++)
	nodes[i]-> selected = sel[i];

    // select elements which have all nodes selected and are bigger
    // than the requested size
    int res = 0;
    for( long i = 0 ; i < n_elements ; i ++) {
	WedgeElement * e = elements[i];
	e-> selected = 0;
	if( e-> get_size() <= req_size) continue;
	if( nodes[e-> p[0]]-> selected &&
	    nodes[e-> p[1]]-> selected &&
	    nodes[e-> p[2]]-> selected &&
	    nodes[e-> p[3]]-> selected &&
	    nodes[e-> p[4]]-> selected &&
	    nodes[e-> p[5]]-> selected){
	    e-> selected = 1;
	    res = 1;
	}
    }

    // now subdivide all selected elements
    for( long i = 0 ; i < n_elements ; i ++)
	elements[i]-> mark = elements[i]-> selected;
    subdivide_marked_elements();

    // if anything was subdivided, try again
    if( res) refine_mesh_around_node( ind, depth, req_size);

    // return the result
    return res;
	
}

int
    Model::prettify_mesh( void)
// ----------------------------------------------------------------------
// tries to remove degenerate wedges by edge collapsing
//
// Returns: 1 - if there were any changes
//          0 - if no changes were made
//
// Algorithm:
//     - for all triangles:
//           - if a triangle is pretty, continue
//           - try to collapse all three edges of the triangle, starting
//             with the smallest edge
//           - an edge can only be collapsed if one of the two nodes
//             is 'moveable'
//           - if an edge is collapsed, restart the loop
//           - otherwise continue
//     - if no adjustments were made in this loop, break
//     - otherwise repeat the process
// ----------------------------------------------------------------------
{
    fprintf( stderr, "\t- PrettifyMessh()\n");

    //	double too_small_size = min_element_size * 10;
    //	double too_small_size = min_element_size * 1;
    //	double too_small_size = min_refine_size;
    double too_small_size = crack_tip_element_size / 2;
	
    // try to remove elements that are not 'nicely' shaped by
    // subdividing them (actually, their neighbors)
    int change = 0;
    int go_again = 1;
    while( go_again) {
	go_again = 0;
	for( long i = 0 ; i < n_elements ; i ++) {
	    WedgeElement * e = elements[i];

	    // skip elements that are too small
	    if( e-> get_size() < too_small_size) continue;

	    // skip nicely shaped elements
	    if( e-> is_pretty()) continue;

	    // try delaunay edge flipping
	    if( delaunay_edge_flip( e-> p[0], e-> p[1]) ||
		delaunay_edge_flip( e-> p[1], e-> p[2]) ||
		delaunay_edge_flip( e-> p[2], e-> p[0]))
	    {
		change = 1;
		go_again = 1;
		break;
	    }

	    // figure out the edges <e1>, <e2>, <e3> so that
	    // <e1> is the shortest, and <e3> is the longest
	    Edge e1, e2, e3;
	    e1.n1 = e-> p[0]; e1.n2 = e-> p[1];
	    e2.n1 = e-> p[1]; e2.n2 = e-> p[2];
	    e3.n1 = e-> p[2]; e3.n2 = e-> p[0];
	    if( e1.length() > e2.length()) swap( e1, e2);
	    if( e2.length() > e3.length()) swap( e2, e3);
	    if( e1.length() > e2.length()) swap( e1, e2);

	    // find a guy on the opposite side of the
	    // longest edge
	    WedgeElement * ep = find_element(
		e, e3.n1, e3.n2);
	    // 				if( ep == NULL || ep-> get_size()
	    // 				    < min_element_size)
	    if( ep == NULL || ep-> get_size()
		< too_small_size)
		subdivide_element(
		    get_element_index( e));
	    else
		subdivide_element(
		    get_element_index( ep));
	    go_again = 1;
	    change = 1;
	    break;
	}
	if( go_again) continue; else break;
    }

    // try to remove remaining non-pretty elements by nice
    // edge-collapsing, i.e. collapsing the shortest edge, when one of
    // the nodes is moveable
    go_again = 1;
    while( go_again) {
	go_again = 0;
	for( long i = 0 ; i < n_elements ; i ++) {
	    WedgeElement * e = elements[i];

	    // skip non-degenerate elements
	    //			if( ! e-> is_degenerate()) continue;
	    if( e-> is_pretty()) continue;
		
	    // figure out the edges <e1>, <e2>, <e3> so that
	    // <e1> is the shortest, and <e3> is the longest
	    Edge e1, e2, e3;
	    e1.n1 = e-> p[0]; e1.n2 = e-> p[1];
	    e2.n1 = e-> p[1]; e2.n2 = e-> p[2];
	    e3.n1 = e-> p[2]; e3.n2 = e-> p[0];
	    if( e1.length() > e2.length()) swap( e1, e2);
	    if( e2.length() > e3.length()) swap( e2, e3);
	    if( e1.length() > e2.length()) swap( e1, e2);

	    // try collapsing edge <e1>
	    if( is_node_moveable( e1.n1)) {
		collapse_edge( e1.n1, e1.n2);
		change = 1; go_again = 1; break;
	    }
	    if( is_node_moveable( e1.n2)) {
		collapse_edge( e1.n2, e1.n1);
		change = 1; go_again = 1; break;
	    }
	    continue;
	}
	if( go_again) continue; else break;
    }

    // try a more aggressive approach to remove the remaining degenerate
    // elements: remove a degenerate element if on two top edges it does
    // not have a neighbor
    go_again = 1;
    while( go_again) {
	go_again = 0;
	for( long i = 0 ; i < n_elements ; i ++) {
	    WedgeElement * e = elements[i];

	    // skip non-degenerate elements
	    if( ! e-> is_degenerate()) continue;

	    // if this elements has less than 2 neighbors on
	    // the top face, just remove it
	    if( get_neighbor_count( e) > 1) continue;

	    fprintf( stderr, "\t- removing degenerate element\n");
	    remove_element( e);
	    delete e;
	    change = 1; go_again = 1; break;
	}
	if( go_again) continue; else break;
    }

    /*
    // remove all elements with 0 neighbors (islands)
    go_again = 1;
    while( go_again) {
    go_again = 0;
    for( long i = 0 ; i < n_elements ; i ++) {
    WedgeElement * e = elements[i];
    if( get_neighbor_count( e) > 0) continue;

    fprintf( stderr, "\t- removing island\n");
    remove_element( e);
    delete e;
    change = 1; go_again = 1; break;
    }
    if( go_again) continue; else break;
    }
    */

    // return result
    return change;
}

int
    Model::fracture_edge( void)
// ----------------------------------------------------------------------
// - finds an edge with the highest stress/yield stress ratio
// - if the stress of such edge has not exceeded yield stress, return 0
// - bisects triangles sharing this edge and return 1
//
// Assumptions:
//   - only the edges of top faces of each element are considered
// ----------------------------------------------------------------------
{
    // there is actually no explicit list of edges, only lists of
    // elements and nodes. To find an edge whose current stress/yield
    // stress is the highest, process each edge in each triangle
    double max_stress_ratio = -1;
    long max_element_ind = -1;
    long max_edge_j = -1;
    for( long i = 0 ; i < n_elements ; i ++) {
	WedgeElement * e = elements[i];
	// skip degenerate elements
	if( e-> is_degenerate()) continue;
	// consider the three edges of the top face of this element
	for( long j = 0 ; j < 3 ; j ++) {
	    // check edge between n1 and n2
	    Node & n1 = * nodes[ e-> p[j%3]];
	    Node & n2 = * nodes[ e-> p[(j+1)%3]];
	    // calculate the original vector from n1 to n2
	    Vector3d ov(
		n1.ox - n2.ox,
		n1.oy - n2.oy,
		n1.oz - n2.oz);
	    // calculate the new vector from n1 to n2
	    Vector3d nv(
		n1.px - n2.px,
		n1.py - n2.py,
		n1.pz - n2.pz);
	    // calculate the percentual change in length
	    double dl = nv.length() / ov.length() - 1;
	    if( isnan(dl)) {
		fprintf( stderr, "Model::fracture_edge():"
		    " dl = nan!!!\n");
		dl = 0;
	    }
	    // calculate edge stress
	    double es = dl * e-> young_mod;

	    // if the yield stress has not been exceeded,
	    // do not even bother considering this edge
	    if( es < e-> yield_stress) continue;

	    // calculate the ration of the edge stress to
	    // the yield stress
	    double stress_ratio = es / e-> yield_stress;
	    assert( ! isnan(stress_ratio));

	    // check if this is the best edge so far
	    if( max_stress_ratio < stress_ratio ||
		max_element_ind < 0) {
		max_stress_ratio = stress_ratio;
		max_element_ind = i;
		max_edge_j = j;
	    }
	}
    }
	
    // if no edge was found that exceeded yield stress, return 0
    if( max_stress_ratio <= 1.0) return 0;

    fprintf( stderr, "(* %f,%ld,%ld)",
	max_stress_ratio,
	max_element_ind,
	max_edge_j);

    // break the found edge by bisecting the triangles
    // that are neighboring it
    assert( max_element_ind >=0 && max_element_ind < n_elements);
    assert( max_edge_j >= 0 && max_edge_j < 3);

    // set e1 to be the element we already have
    WedgeElement * e1 = elements[ max_element_ind];

    // set tn2 and tn3 the nodes of the edge (the top nodes)
    // set tn1 to be the third not of the top face
    long tn2 = e1-> p[ max_edge_j];
    long bn2 = e1-> p[ max_edge_j + 3];
    long tn3 = e1-> p[ (max_edge_j + 1) % 3];
    long bn3 = e1-> p[ (max_edge_j + 1) % 3 + 3];
    long tn1 = e1-> p[ (max_edge_j + 2) % 3];
    long bn1 = e1-> p[ (max_edge_j + 2) % 3 + 3];

    // find the element e2 which also shares this edge
    WedgeElement * e2 = find_element( e1, tn2, tn3);
	
    // set tn4 and bn4 to be the nodes not on the edge
    long tn4 = -1; long bn4 = -1;
    if( e2 != NULL) {
	if( e2-> p[0] != tn2 && e2-> p[0] != tn3) {
	    tn4 = e2-> p[0]; bn4 = e2-> p[3];
	} else if( e2-> p[1] != tn2 && e2-> p[1] != tn3) {
	    tn4 = e2-> p[1]; bn4 = e2-> p[4];
	} else if( e2-> p[2] != tn2 && e2-> p[2] != tn3) {
	    tn4 = e2-> p[2]; bn4 = e2-> p[5];
	} else {
	    assert( 0);
	}
    }
	
    // delete both elements e1 and e2, they are no longer needed
    remove_element( e1);
    delete e1;
    if( e2 != NULL) {
	remove_element( e2);
	delete e2;
    }

    // create two identical new nodes tn5 and tn5 splitting tn2-tn3,
    // and one node tn56 splitting bn2-bn3
    long tn5 = add_node( new Node(
			     (nodes[tn2]->ox + nodes[tn3]->ox)/2,
			     (nodes[tn2]->oy + nodes[tn3]->oy)/2,
			     (nodes[tn2]->oz + nodes[tn3]->oz)/2,
			     (nodes[tn2]->px + nodes[tn3]->px)/2,
			     (nodes[tn2]->py + nodes[tn3]->py)/2,
			     (nodes[tn2]->pz + nodes[tn3]->pz)/2,
			     0 ));
    long tn6 = add_node( new Node(
			     (nodes[tn2]->ox + nodes[tn3]->ox)/2,
			     (nodes[tn2]->oy + nodes[tn3]->oy)/2,
			     (nodes[tn2]->oz + nodes[tn3]->oz)/2,
			     (nodes[tn2]->px + nodes[tn3]->px)/2,
			     (nodes[tn2]->py + nodes[tn3]->py)/2,
			     (nodes[tn2]->pz + nodes[tn3]->pz)/2,
			     0 ));
    long bn56 = add_node( new Node(
			      (nodes[bn2]->ox + nodes[bn3]->ox)/2,
			      (nodes[bn2]->oy + nodes[bn3]->oy)/2,
			      (nodes[bn2]->oz + nodes[bn3]->oz)/2,
			      (nodes[bn2]->px + nodes[bn3]->px)/2,
			      (nodes[bn2]->py + nodes[bn3]->py)/2,
			      (nodes[bn2]->pz + nodes[bn3]->pz)/2,
			      1 ));

    // interpolate growth vectors for the new nodes
    interpolate_growth_vector( tn2, tn3, tn5, 0.5);
    interpolate_growth_vector( tn2, tn3, tn6, 0.5);
    interpolate_growth_vector( bn2, bn3, bn56, 0.5);
	
    // add up to 4 wedges, replacing e1 and e2
    add_element( new WedgeElement(
		     tn1, tn2, tn6, bn1, bn2, bn56,
		     ym_map, pr_map, ys_map, ft_map));
    add_element( new WedgeElement(
		     tn1, tn5, tn3, bn1, bn56, bn3,
		     ym_map, pr_map, ys_map, ft_map));
    if( e2 != NULL) {
	add_element( new WedgeElement(
			 tn4, tn6, tn2, bn4, bn56, bn2,
			 ym_map, pr_map, ys_map, ft_map));
	add_element( new WedgeElement(
			 tn4, tn3, tn5, bn4, bn3, bn56,
			 ym_map, pr_map, ys_map, ft_map));
    }

    // adjust the triangular mesh to have only pretty triangles
    prettify_mesh();

    // yes, there was breakage
    return 1;
}

long num_edge_neighbors( long n1, long n2)
// ----------------------------------------------------------------------
// returns the number of elements sharing an edge n1-n2 on the top face
//
// return: 0, 1 or 2
// ----------------------------------------------------------------------
{
    assert( n1 >= 0 && n1 < n_nodes);
    assert( n2 >= 0 && n2 < n_nodes);
    WedgeElement * e1 = find_element( NULL, n1, n2);
    if( e1 == NULL) return 0;
    WedgeElement * e2 = find_element( e1, n1, n2);
    if( e2 == NULL) return 1;
    return 2;
}

void
    Model::interpolate_growth_vector(
	long n1, long n2, long n3, double t)
// ----------------------------------------------------------------------
// sets the growth vector for node 'n3' by interpolating the growth
//      vectors from 'n1' and 'n2' at isoparametric coordinate 't'
// ----------------------------------------------------------------------
{
    nodes[n3]-> gv = nodes[n1]-> gv * t + nodes[n2]-> gv * (1-t);
}

int
    Model::_fracture_at_node( long n0, NodeList & new_ftips, double & farea)
// ----------------------------------------------------------------------
// model a fracture at node n0, based on the fracture plane stored in
// the node
//
// - if any of the elements touching <n0) are bigger than 
//   crack_tip_element_size:
//     - refine them
//     - if <n0> is a fracture tip, put it back into <ftips>
//     - return '2'
// - pre-crack:
//     - if any of the elements connected to <n0> is split by the
//       fracture plane:
//         - split it the corresponding edge <n1,n2>, creating a split-node
//           <sn>
//         - collapse the edge <n1,sn> or <n2,sn> if the resulting angle
//           is too small, making sure <n1> or <n2> are moveable, which
//           has the effect of snapping the <n1> or <n2> onto the fracture
//           plane
//         - return 3
// - now all elements are completely above or below the fracture plane
// - do the actual fracture:
//     - make a new copy of node n0, call it 'cn0'
//     - all elements above the stress plane will keep n0
//     - all elements below the stress plane will get cn0
// - figure out what should be the crack tips now:
//     - all nodes connected to both <n0> and <cn0> are fracture tips
//     - add these into <ftips>
// - mark plastic zone around the new fracture tips
//
// - return: 1 - if fracture was modeled
//           2 - if refinement had to be done
//           0 - if nothing was fractured (this should never happend, and
//               if it does, it is an indication of an internal error)
//           - in <new_ftips> returns the new fracture tips
//           - in <farea> returns the area of the fracture
// ----------------------------------------------------------------------
{
    assert( n0 >= 0);
    assert( n0 < n_nodes);
	
    // in case we need to return early, set the fracture length to be 0
    farea = 0;

    // make sure that all elements around the fracture node <n0> are
    // of appropriate size (i.e. smaller than crack_tip_element_size).
    //	int res = refine_mesh_around_node (n0, 2, crack_tip_element_size);
    int res = refine_mesh_around_node (n0, 1, crack_tip_element_size);
    if( res ) {
	fprintf( stderr, "\t- needed refinement first\n");
	//		fprintf( stderr, "\t- needed refinement first, "
	//			 "recalculating equilibrium...\n");
	//		relax();
	//		calculate_stresses();
	//		return 2;
    }

    Node & n = * nodes[ n0];

    // pre-fracture:
    //   - find all elements connected to <n0> which is split by
    //     the fracture plane and split them
    while( 1) {
		
	int go_again = 0;
		
	for( long i = 0 ; i < n.n_refs ; i ++) {
	    WedgeElement & e = * n.refs[i];
	    assert( get_element_index( & e) != -1);
			
	    // into n1 and n2 calculate the indices of the
	    // other 2 nodes that belong to the element on the
	    // same face
	    long in0 = e.get_internal_index( n0);
	    assert( in0 > -1);
	    long n1 = -1; long n2 = -1;
	    if( in0 < 3) {
		// <n0> is on the top face of <e>
		n1 = e.p[(in0+1)%3];
		n2 = e.p[(in0+2)%3];
	    } else {
		n1 = e.p[(in0+1)%3+3];
		n2 = e.p[(in0+2)%3+3];
	    }

	    Vector3d nv( nodes[n0]-> ps_nx,
		nodes[n0]-> ps_ny,
		nodes[n0]-> ps_nz);
	    Vector3d p0( nodes[n0]-> ox,
		nodes[n0]-> oy,
		nodes[n0]-> oz);
	    Vector3d p1( nodes[n1]-> ox,
		nodes[n1]-> oy,
		nodes[n1]-> oz);
	    Vector3d p2( nodes[n2]-> ox,
		nodes[n2]-> oy,
		nodes[n2]-> oz);

	    // calculate 'r' the isoparametric coordinate of the
	    // intersection
	    double r = split_edge_by_plane_raw( nv, p0, p1, p2);
	    // and the point p3 at the split
	    Vector3d p3 = p1 * (1-r) + p2 * r;
			
	    // find out the angles
	    double angle1 = Vector3d::angle( p1-p0,p3-p0);
	    double angle2 = Vector3d::angle( p2-p0,p3-p0);

	    /*
	      {
	      fprintf( stderr, "r = %f\n", r);
	      fprintf( stderr, "n1,n2=%ld,%ld a1,a2=%f,%f\n",
	      n1, n2, angle1, angle2);
	      }
	    */
	    // Determine whether the triangle is more above or
	    // below.
	    double a1 = 90 - Vector3d::angle( p1-p0, nv);
	    double a2 = 90 - Vector3d::angle( p2-p0, nv);
	    double above_angle = 0;
	    if( a1 > above_angle) above_angle = a1;
	    if( a2 > above_angle) above_angle = a2;
	    double below_angle = 0;
	    if( fabs(a1) > below_angle) below_angle = fabs(a1);
	    if( fabs(a2) > below_angle) below_angle = fabs(a2);
	    int more_below = above_angle < below_angle;
	    double p01l = (p0-p1).length();
	    double p02l = (p0-p2).length();
	    double p03l = (p0-p3).length();
	    double p12l = (p1-p2).length();
	    double p13l = (p1-p3).length();
	    double p23l = (p2-p3).length();
	    int too_small =
		(p01l < min_element_size) ||
		(p02l < min_element_size) ||
		(p03l < min_element_size) ||
		(p12l < min_element_size) ||
		(p13l < min_element_size) ||
		(p23l < min_element_size);
				

	    // should the whole triangle be considered above
	    // or below?
	    if( r <= 0 ||
		r >= 1 ||
		! finite( r) ||
		angle1 < min_split_angle ||
		angle2 < min_split_angle ||
		too_small)
	    {
		e.is_below_fp = more_below;
		continue;
	    }

	    // if the fracture plane goes too close to a
	    // boundary node, we are avoiding back-cracking
	    int n1_is_bn = ! is_node_moveable( n1);
	    int n2_is_bn = ! is_node_moveable( n2);
	    int n12_is_be = n1_is_bn && n2_is_bn;
	    int is_too_close_to_bn =
		(n1_is_bn && angle1<min_ftbn_angle) ||
		(n2_is_bn && angle2<min_ftbn_angle);
	    int is_too_close_to_be =
		(n1_is_bn && angle1<min_ftbe_angle) ||
		(n2_is_bn && angle2<min_ftbe_angle);
	    if( (is_too_close_to_bn && ! n12_is_be)||
		(is_too_close_to_be && n12_is_be))
	    {
		// if the element is more above the plane
		if( more_below) {
		    fprintf( stderr,
			"\t- back-crack avoided 1\n");
		    e.is_below_fp = 1;
		    continue;
		} else {
		    fprintf( stderr,
			"\t- back-crack avoided 2\n");
		    e.is_below_fp = 0;
		    continue;
		}
	    }
	    // the element needs to be pre-fractured
	    fprintf( stderr, "\t- pre-fracture (%ld) r=%f "
		"a1=%f a2=%f p01=%f p02=%f\n",
		get_element_index( & e), r,
		angle1, angle2,
		(p1-p0).length(),
		(p2-p0).length());
	    assert( r >= 0);
	    assert( r <= 1);
	    long mn = split_edge( n1, n2, r);
	    fprintf( stderr, "\t- n1=%ld n2=%ld mn=%ld\n",
		n1, n2, mn);
	    // collapse <n1,mn> if the angle1 is not
	    // pretty, and <n1> is moveable
	    if( angle1 < pretty_triangle_min_angle &&
		is_node_moveable( n1)) {
		fprintf( stderr,
		    "\t- node snapped(I) %ld -> %ld\n",
		    n1, mn);
		collapse_edge( n1, mn);
	    }
	    // collapse <n2,mn> if the angle2 is not
	    // pretty, and <n2> is moveable
	    if( angle2 < pretty_triangle_min_angle &&
		is_node_moveable( n2)) {
		fprintf( stderr,
		    "\t- node snapped(II) %ld -> %ld\n",
		    n2, mn);
		collapse_edge( n2, mn);
	    }
	    go_again = 1;
	    break;
	}
		
	if( ! go_again) break;
    }
	
    fprintf( stderr, "\t- splitting: ");

    // create a copy of <n0>, called <cn0>, which will be assigned
    // to all elements below the fracture plane
    long cn0 = add_node( new Node( n));

    // - make a copy of the list of elements connected to node 'n'
    //   (this is necessary, because as we are removing and adding
    //   elements during splitting, and this list can/will change)
    long n_refs = n.n_refs;
    WedgeElement * refs[ n_refs];
    for( long i = 0 ; i < n_refs ; i ++)
	refs[i] = n.refs[i];

    // - for each element 'e' on the list:
    //      - if e is above the plane of stress, do nothing to it
    //      - if e is below the plane of stress, replace its reference to
    //        to node <n0> with a reference to node <cn0>
    //      - if e is cut by the plane, there is an error!!!
    for( long i = 0 ; i < n_refs ; i ++) {
	WedgeElement & e = * refs[i];
	assert( get_element_index( & e) != -1);
	if( e.is_below_fp) {
	    long in0 = e.get_internal_index( n0);
	    fprintf( stderr, "b");
	    long p[6];
	    for( long k = 0 ; k < 6 ; k ++) p[k] = e.p[k];
	    p[in0]=cn0;
	    remove_element( & e); delete & e;
	    add_element( new WedgeElement(
			     p[0], p[1], p[2],
			     p[3], p[4], p[5],
			     ym_map, pr_map, ys_map, ft_map));
	} else {
	    fprintf( stderr, "a");
	}
    }
    fprintf( stderr, "\n");

    // find new crack_tip nodes and update fracture length:
    //  - crack tip nodes are nodes connected to both n0 & cn0, but
    //    not fixed nodes
    for( long i = 0 ; i < n.n_refs ; i ++) {
	WedgeElement & e = * n.refs[i];
	for( long j = 0 ; j < 6 ; j ++) {
	    long fn = e.p[j];
	    // skip fixed nodes
	    if( nodes[fn]-> is_fixed) continue;
	    // skip single point nodes (and split them too)
	    if( split_single_point_node( fn)) continue;
	    if( edge_exists( fn, cn0)) {
		new_ftips.push_unique( fn);
		farea += Node::dist_p( nodes[n0], nodes[fn]) *
		    ( get_mat_thickness_at_node( n0) +
			get_mat_thickness_at_node( fn))/2.0;
	    }
	}
    }

    // mark plastic zone
    //	for( long i = 0 ; i < new_ftips.n ; i ++)
    //		mark_plastic_zone( new_ftips(i));

    // 	// handle single node contacts
    // 	for( long i = 0 ; i < new_ftips.n ; i ++)
    // 		split_single_point_node( new_ftips.array[i]);

    split_single_point_node( n0);
    split_single_point_node( cn0);

    // for statistics reporting, keep the count of fractures
    n_fractures ++;

    return 1;
}

void Model::fracture_extend( void)
// ---------------------------------------------------------------------------
// - assumes that nodal stresses have been calculated (i.e. that
//   calculate_nodal_stresses() was called prior to here)
// - picks one of the fracture tips (the first one in the FIFO ftips)
// - if the area around the fracture tip needs to be refined, refine it
//   and put the fracture tip back onto the list ot crack tips
// - if refinement is needed, refine and put the node into new ftips & return
// - calculate precise stresses at this fracture tip (using multiresolution
//   mesh)
// - calculate the current potential energy
// - fracture the node & remember the fracture tips & length of the fracture
// - relax & compute the new potential energy
// - prettify the mesh, smooth the mesh, handle single point contacts around
//   crack tips
// - calculate the released energy = (pe_old - pe_new)/farea
// - if the released energy is above toughness, put the new fracture tip(s)
//   into the FIFO ftips, otherwise don't
// ---------------------------------------------------------------------------
{
    assert( ! ftips.is_empty());

    // this is where the new fracture tips are accumulated
    NodeList new_ftips;

    // get a fracture tip
    long n0 = ftips.get();
    fprintf( stderr, "\t- continuing fracture at node %ld\n", n0);
    assert( n0 >= 0 && n0 < n_nodes);
    Node & n = * nodes[n0];
    assert( ! n.is_fixed);
    //	if( n.ps_val < n.yield_stress) {
    //		fprintf( stderr, "\t- closing fracture - not enough stress\n");
    //		return;
    //	}

    // calculate precise nodal stresses at the selected fracture tip
    calculate_precise_stress_at_node( n0);

    if( ! is_node_split_valid( n0)) {
	fprintf( stderr, "\t- closing fracture - invalid split\n");
	// calculate_stresses();
	return;
    }

    // calculate the old potential energy
    double old_pe = calculate_potential_energy();

    // fracture the node
    double farea;
    int res = _fracture_at_node( n0, new_ftips, farea);
    // if refinement was needed, put ind back into ftips and return
    if( res == 2) {
	assert( 0);
	ftips.add_unique( n0);
	return;
    }

    // 	fprintf( stderr, "flength=%f\n", flength);
    // 	fprintf( stderr, "1) old_pe=%f new_pe=%f diff=%f\n",
    // 		 old_pe, new_pe, new_pe-old_pe);

    // smooth out the mesh around the new crack tips
    (void) smooth_mesh_around_nodes( new_ftips, 3);

    // refine mesh around the new crack tips
    for( long i = 0 ; i < new_ftips.n ; i ++) {
	//	    refine_mesh_around_node(new_ftips.array[i], 2, crack_tip_element_size);
	refine_mesh_around_node(new_ftips.array[i], 1, crack_tip_element_size);
    }

    // make improvements to the mesh if possible (i.e. get rid of
    // slivers if possible, maybe by edge collapsing)
    prettify_mesh();

    // relax the model & calculate the new potential energy
    relax();
    // calculate_stresses();
    double new_pe = calculate_potential_energy();

    // calculate the released energy
    fprintf( stderr, "\t- farea=%f\n", farea);
    double G = (new_pe - old_pe) / farea;
    fprintf( stderr, "\t- G   = %-20.10f\n", G);
    // calculate the energy needed to close the fracture
    double Kic = n.ftoughness * n.yield_stress;
    double Gic = ((1 - sqr(n.poisson_ratio)) / n.young_modulus) * sqr(Kic);
    fprintf( stderr, "\t- Gic = %-20.10f\n", Gic);

    // append then new fracture tips to the current list of fractures
    // if the released energy is higher than the one required to close
    // the fracture
    if( G > Gic)
	while( ! new_ftips.is_empty())
	    ftips.add_unique( new_ftips.pop());
    else {
	fprintf( stderr, "\t- closing fracture & marking "
	    "plastic zone\n");
	while( ! new_ftips.is_empty())
	    mark_plastic_zone( new_ftips.pop());
    }

    fprintf( stderr, "\t- n_crack_tips after: %ld\n", ftips.n);
}

long Model::fracture_new( int count_mode)
// ---------------------------------------------------------------------------
// - assumes that nodal stresses have been calculated (i.e. that
//   calculate_nodal_stresses() was called prior to here)
// - if count_mode = 1 then it returns the number of nodes where
//   a fracture could occur
// - otherwiser:
//   - skips all nodes in plastic zones
//   - picks a node which has the highest pstress/yield_stress ratio
//   - introduces a fracture at the node
//   - prettify the mesh, smooth the mesh, handle single point contacts around
//     new crack tips
//
// Returns:
//    - if count_mode == 1
//         - returns the number of nodes that could be fractures
//             - if there are fracture tips, return 1
//             - otherwise count nodes that exceed yield stress
//               and their splitting would result is a valid split
//    - if count_mode == 0 return:
//          0 = if nothing was split (tbd: error?)
//          1 = if a node was split
//          2 = if a refinement was needed
// ---------------------------------------------------------------------------
{
    if( count_mode) {
	// return the count of the fracture candidates
	long count = 0;
	for( long i = 0 ; i < n_nodes ; i ++) {
	    Node & n = * nodes[i];
	    // skip fixed nodes
	    if( n.is_fixed != 0) continue;
	    // skip nodes in plastic zones
	    if( n.is_in_pzone) continue;
	    // skip nodes whose principal stress is
	    // less than the material's yield stress
	    if( n.ps_val < n.yield_stress ||
		n.yield_stress < 0)
		continue;
	    // skip nodes that would not split anything
	    if( ! is_node_split_valid( i)) continue;
	    // we have a new candidate for splitting
	    count ++;
	}
	return count;
    }

    //	fprintf( stderr, "FRACTURE NEW FRACTURE NEW FRACTURE NEW FRACTURE\n");

    // this is where the new fracture tips are accumulated
    NodeList new_ftips;

    long n0 = -1;
    // go through the nodes and pick one with highest tensional
    // stress
    // - skip fixed nodes, these are never split (TBD: really?)
    // - test if an attempt at splitting at this node would
    //   actually be succesful (i.e. would the plane of stress
    //   split anything?)
    //   if not, such node is not considered
    // - put the index of the node to be split into <n0>
    double max_ratio = 0;
    for( long i = 0 ; i < n_nodes ; i ++) {
	Node & n = * nodes[i];
	// skip fixed nodes
	if( n.is_fixed != 0) continue;
	// skip nodes in plastic zones
	if( n.is_in_pzone) continue;
	// skip nodes with 0 yield stress
	if( n.yield_stress <= 0) continue;
	// calculate ratio of stress to yield stress
	double r = n.ps_val / n.yield_stress;
	// skip nodes which have not exceeded yield
	// stress
	if( r < 1) continue;
	// skip this node if the node split is not
	// valid
	if( ! is_node_split_valid(i)) continue;
	if( r > max_ratio) {
	    n0 = i;
	    max_ratio = r;
	}
    }

    // if no node has been found where a fracture might occur,
    // return and indicate that nothing has been split
    if( n0 < 0) {
	fprintf( stderr,
	    "****************************************\n"
	    "WARNING WARNING WARNING WARNING WARNING \n"
	    "  - fracture_new() found no candidates!\n"
	    "****************************************\n");
	return 0;
    }

    fprintf( stderr, "\t- new new fracture "
	"at node #%ld ps=%.10f\n",
	n0, nodes[n0]-> ps_val);
    // fracture at the node
    double farea;
    long res = _fracture_at_node( n0, new_ftips, farea);
    if( res == 2) {
	fprintf( stderr, "\t- refinement needed first\n");
	return 2;
    }

    // smooth out the mesh around the new crack tips
    (void) smooth_mesh_around_nodes( new_ftips, 3);

    // refine mesh around the new crack tips
    for( long i = 0 ; i < new_ftips.n ; i ++) {
	//          refine_mesh_around_node(new_ftips.array[i], 2, crack_tip_element_size);
	refine_mesh_around_node(new_ftips.array[i], 1, crack_tip_element_size);
    }

    // append all new fracture tips to the current list of fractures
    while( ! new_ftips.is_empty())
	ftips.add_unique( new_ftips.pop());
    fprintf( stderr, "\t- n_crack_tips after: %ld\n", ftips.n);

    // make improvements to the mesh if possible (i.e. get rid of
    // slivers if possible, maybe by edge collapsing)
    prettify_mesh();

    return 1;
}

void Model::orient_elements( void)
// ---------------------------------------------------------------------------
// TBD: useless routine - can be removed, only used as a debugging hack
// ---------------------------------------------------------------------------
{
    clock_t start_time = clock();
    fprintf( stderr, "Orienting elements:\n");
    fprintf( stderr, "    - this can be removed, it does nothing\n");
    clock_t end_time = clock();
    fprintf( stderr, "\t- time to complete: %.4fsec\n",
	double(end_time - start_time) / CLOCKS_PER_SEC);
}

#include "image_png.hh"

void save_screen_dump( const char * prefix, long count)
{
    char fname[ 4096];
    sprintf( fname, "%s-%010ld.png", prefix, count);
    fprintf( stderr, "Screen dump '%s' ", fname);
    fprintf( stderr, "[%d x %d] ", win_width, win_height);

    fprintf( stderr, "alloc...");
    ImagePNG img( win_width, win_height);

    fprintf( stderr, "glReadPixels()...");

    glReadPixels( 0, 0, win_width, win_height, GL_RGB, GL_UNSIGNED_BYTE,
	img.data);

    fprintf( stderr, "save...");
    if( img.save( fname)) {
	fprintf( stderr, "ERROR saving png\n");
    } else {
	fprintf( stderr, "SUCCESS\n");
    }
}

int Model::simulation_step( void)
// ----------------------------------------------------------------------
// performs a single time step of simulation
// - returns: 0 if simulation is NOT over
//            1 if simulation is over
// ----------------------------------------------------------------------
{
    auto_save();

    // save animation - only if requested and when the model is relaxed
    static int record_frame = 1;
    if( animation_fname_mask != NULL && record_frame)
	save_screen_dump( animation_fname_mask, animation_count ++);
	
    clock_t start_time = clock();

    double dt = min_time_step;
	
    fprintf( stderr,
	"SIMSTEP: t=%.10f dt=%.10f ne=%-6ld nn=%-6ld nf=%ld(%.5fs/f) rt=%.1fs st=%.1fs\n",
	time_curr, dtc * min_time_step,
	n_elements, n_nodes,
	n_fractures,
	(double(clock()) / CLOCKS_PER_SEC)/n_fractures,
	double (relax_time) / CLOCKS_PER_SEC,
	double (stress_time) / CLOCKS_PER_SEC
	     );
    fprintf( stderr, "KE RECALC : %ld\n", ke_recalc_n);

    // finite automaton state:
    static long next_step = 10;
    long last_step = next_step;

    fprintf( stderr, "STEPSTEP : %ld\n", next_step);

    switch( next_step) {
	case 10: // local relaxation
	    {
		fprintf( stderr, "\t- local relaxation\n");
		relax_local();
		(void) flag_nodes_with_error();
		next_step = 15;
		record_frame = 0;
		break;
	    }
	case 15: // global relaxation
	    {
		fprintf( stderr, "\t- global relaxation\n");
		relax_global();
		fprintf( stderr, "\t- calculating stresses\n");
		// calculate_stresses();
		(void) flag_nodes_with_error();
		next_step = 20;
		record_frame = 1;
		break;
	    }
	case 20: // refine the mesh if needed
	    {
		fprintf( stderr, "\t- refining mesh (skipped)\n");
		if (0) {
		    int res = refine_mesh();
		    if( res) next_step = 10;
		    else next_step = 30;
		    record_frame = 0;
		    break;
		} else {
		    next_step = 30;
		    record_frame = 0;
		    break;
		}
	    }
	case 30: // extend existing fracture
	    {
		if( ftips.is_empty()) {
		    fprintf( stderr, "\t- no fractures exist\n");
		    next_step = 40;
		    record_frame = 0;
		    break;
		}

		// try to extend one of the fracture tips
		fprintf( stderr, "\t- extending existing fracture\n");
		fracture_extend();
		next_step = 30;
		record_frame = 1;
		break;
	    }
	case 40: // create a new fracture
	    {
		// force the calculation of all stresses!
		for (long i = 0 ; i < n_nodes ; i ++)
		    nodes[i]-> selected = 1;
		calculate_stresses();
		// count the number of candidates for breaking
		fprintf( stderr, "\t- counting fracture candidates...\n");
		long count = get_number_of_fracture_candidates();
		fprintf( stderr, "\t- fracture candidates: %ld\n", count);
		if( count == 0) {
		    next_step = 50;
		    record_frame = 0;
		    break;
		}
		// create a new fracture
		if( break_method == BreakElement)
		    fracture_element();
		else if( break_method == BreakNode)
		    fracture_new();
		else
		    assert( 0);
		next_step = 10;
		record_frame = 0;
		break;
	    }
	case 50: // advance time (to the closest time when there is one
		 // fracture candidate), and also progressive save is done
		 // here
	case 60: // loop for finding t_max
	    {
		if( next_step == 50)
		{
		    progressive_save();
		    if( time_curr >= time_total) {
			fprintf( stderr, "\t- SIMULATION DONE\n");
			return 1;
		    }
		    fprintf( stderr,
			"\t- Initializing time_advancing...\n");
		    dtc = 1;
		    fc_count = 0;
		    t1 = time_curr;
		}
		fprintf( stderr, "\t- Looking for t_max @ dtc = %ld\n", dtc);
		t1 = time_curr + dt * dtc;
		if( t1 > time_total) {
		    // we have overshot the total simulation time,
		    // there is no need to look there... adjust dtc
		    // to be as small as possible to be just above
		    // simulation end
		    dtc = long( ceil( (time_total - time_curr)/dt));
		    if( dtc <= 0) dtc = 1;
		    t1 = time_curr + dt * dtc;
		    assert( t1 >= time_total);
		}
		if( t1 > time_curr + max_time_step) {
		    // we have overshot the max. time step
		    dtc = long( ceil( max_time_step/dt));
		    if( dtc <= 0) dtc = 1;
		    t1 = time_curr + dt * dtc;
		}
		grow( t1);
		relax();
		calculate_stresses();
		fc_count = get_number_of_fracture_candidates();
		if( fc_count == 0) {
		    if( t1 >= time_total) {
			// this is actually simulation end, but we
			// only want one exit point, so...
			fprintf( stderr, "\t- end is near now :)\n");
			time_curr = t1;
			next_step = 50;
			record_frame = 1;
			break;
		    }
		    if( t1 >= time_curr + max_time_step) {
			fprintf( stderr,
			    "\t- reached max. time step\n");
			time_curr = t1;
			next_step = 50;
			record_frame = 1;
			break;
		    }
		    fprintf( stderr,
			"\t- no fractures @ dtc = %ld, "
			"doubling up\n",
			dtc);
		    dtc = dtc * 2;
		    next_step = 60;
		    record_frame = 1;
		    break;
		}

		// we have found the max. time
		fprintf( stderr, "\t- found t_max @ dtc = %ld\n", dtc);

		// end time but there are still some fracture candidates
		if( t1 >= time_total && 0) {
		    fprintf( stderr,
			"\t- reached end time with fractures\n");
		    time_curr = t1;
		    next_step = 10;
		    record_frame = 0;
		    break;
		}

		// in dtc we have number of time steps to perform to get
		// to time when there is for sure a fracture. Now we have
		// to find the smallest such dtc
		min_dtc = dtc/2;
		max_dtc = dtc;
		next_step = 70;
		record_frame = 1;
		break;
	    }
	case 70: // perform halving
	    {
		fprintf( stderr, "\t- halving t_min and t_max\n");
		assert( max_dtc - min_dtc > 0);
		if( max_dtc - min_dtc == 1) {
		    fprintf( stderr, "\t- found next time step\n");
		    // max_dtc has the number of steps required to
		    // ensure at least one fracture candidate - so
		    // advance the time to that point (if necessary)
		    time_curr = time_curr + dt * max_dtc;
		    if( dtc < max_dtc) {
			grow( time_curr);
			relax();
			calculate_stresses();
		    }
		    next_step = 10;
		    record_frame = 0;
		    break;
		}
		
		// do the halving
		dtc = (max_dtc + min_dtc) / 2;
		t1 = time_curr + dt * dtc;
		grow( t1);
		relax();
		calculate_stresses();
		fc_count = get_number_of_fracture_candidates();
		if( fc_count == 0)
		    min_dtc = dtc;
		else
		    max_dtc = dtc;
		next_step = 70;
		record_frame = 1;
		break;
	    }
	default:
	    die( "Internal error - next_step.");
		
    }
	
    clock_t end_time = clock();
    fprintf( stderr, "\t- time to complete step #%ld: %.4fs\n",
	last_step,
	double( end_time - start_time) / CLOCKS_PER_SEC
	     );
	
    return 0;
}

void Model::_save_header( FILE * fp)
{
    fprintf( fp, "\n");
    fprintf( fp, "sim_time_total         = %.30f\n", time_total);
    fprintf( fp, "min_time_step          = %.30f\n", min_time_step);
    fprintf( fp, "max_time_step          = %.30f\n", max_time_step);
    fprintf( fp, "sim_time_curr          = %.30f\n", time_curr);
    fprintf( fp, "growth_x               = %.30f\n", growth_x);
    fprintf( fp, "growth_y               = %.30f\n", growth_y);
    fprintf( fp, "growth_z               = %.30f\n", growth_z);
    fprintf( fp, "shrink_top_t0          = %.30f\n", shrink_top_t0);
    fprintf( fp, "shrink_top_t1          = %.30f\n", shrink_top_t1);
    fprintf( fp, "shrink_top_val         = %.30f\n", shrink_top_val);
    fprintf( fp, "shrink_bot_t0          = %.30f\n", shrink_bot_t0);
    fprintf( fp, "shrink_bot_t1          = %.30f\n", shrink_bot_t1);
    fprintf( fp, "shrink_bot_val         = %.30f\n", shrink_bot_val);
    fprintf( fp, "shrink_height_t0       = %.30f\n", shrink_height_t0);
    fprintf( fp, "shrink_height_val0     = %.30f\n", shrink_height_val0);
    fprintf( fp, "shrink_height_t1       = %.30f\n", shrink_height_t1);
    fprintf( fp, "shrink_height_val1     = %.30f\n", shrink_height_val1);
    fprintf( fp, "fracture_tip_inertia   = %.30f\n", fracture_tip_inertia);
    fprintf( fp, "gravity_x              = %.30f\n", gravity_x);
    fprintf( fp, "gravity_y              = %.30f\n", gravity_y);
    fprintf( fp, "gravity_z              = %.30f\n", gravity_z);
    fprintf( fp, "max_break_size         = %.30f\n", max_break_size);
    fprintf( fp, "min_refine_size        = %.30f\n", min_refine_size);
    fprintf( fp, "min_element_size       = %.30f\n", min_element_size);
    fprintf( fp, "max_element_size       = %.30f\n", max_element_size);
    fprintf( fp, "plastic_zone_size      = %.30f\n", pzone_size);
    fprintf( fp, "crack_tip_element_size = %.30f\n",
	crack_tip_element_size);
    fprintf( fp, "crack_tip_sub_element_size = %.30f\n",
	crack_tip_sub_element_size);
    fprintf( fp, "precision              = %.30f\n", precision);
    fprintf( fp, "error_radius           = %.30f\n", error_radius);
    fprintf( fp, "min_error_radius       = %.30f\n", min_error_radius);
    fprintf( fp, "\n");
    fprintf( fp, "progressive_save = %s\n",
	progressive_save_on ? "true" : "false");
    if( progressive_save_on) {
	fprintf( fp, "progressive_save_fmask = \"%s\"\n",
	    progressive_save_fmask);
	fprintf( fp, "progressive_save_skip = %.20f\n",
	    progressive_save_skip);
    }
    fprintf( fp, "\n");
    fprintf( fp, "# describe how to randomize material properties\n");
    fprintf( fp, "BEGIN material_properties\n");
    fprintf( fp, "\tyoung_modulus ");
    ym_map-> print( fp);
    fprintf( fp, "\n");
    fprintf( fp, "\tpoisson_ratio ");
    pr_map-> print( fp);
    fprintf( fp, "\n");
    fprintf( fp, "\tyield_stress ");
    ys_map-> print( fp);
    fprintf( fp, "\n");
    fprintf( fp, "\tfracture_toughness ");
    ft_map-> print( fp);
    fprintf( fp, "\n");
    fprintf( fp, "END material_properties\n");
    fprintf( fp, "\n");
}

void Model::save( const char * fname)
{
    FILE * fp = fopen( fname, "w");
    assert( fp != NULL);
	
    _save_header( fp);

    fprintf( fp, "%ld # number of nodes\n", n_nodes);
    for( long i = 0 ; i < n_nodes ; i ++)
	fprintf( fp, "%.10f %.10f %.10f %.10f %.10f %.10f "
	    "%.10f %.10f %.10f %ld\n",
	    nodes[i]-> ox,
	    nodes[i]-> oy,
	    nodes[i]-> oz,
	    nodes[i]-> px - nodes[i]-> ox,
	    nodes[i]-> py - nodes[i]-> oy,
	    nodes[i]-> pz - nodes[i]-> oz,
	    nodes[i]-> gv.x,
	    nodes[i]-> gv.y,
	    nodes[i]-> gv.z,
	    long( nodes[i]-> is_fixed +
		4 * (nodes[i]-> is_in_pzone ? 1 : 0))
		 );

    fprintf( fp, "\n%ld # number of elements\n", n_elements);
    for( long i = 0 ; i < n_elements ; i ++) {
	WedgeElement * e = elements[i];
	fprintf( fp, "fivewall %ld %ld %ld %ld %ld %ld %ld\n",
	    i, e-> p[0], e-> p[1], e-> p[2],
	    e-> p[3], e-> p[4], e-> p[5]);
    }

    fprintf( fp, "\n%ld # number of fracture tips\n", ftips.n);
    for( long i = 0 ; i < ftips.n ; i ++)
	fprintf( fp, "%ld\n", ftips.array[i]);

    fclose( fp);
}

void Model::print_stresses( FILE * fp)
// ---------------------------------------------------------------------------
// output to 'fp' the stresses on all nodes and elements
// ---------------------------------------------------------------------------
{
    fprintf( fp, "\n");
    fprintf( fp, "Time : %.4f\n", time_curr);
    fprintf( fp, "Stresses on ELEMENTS:\n");
    fprintf( fp, "....................................................\n");
    for( long i = 0 ; i < n_elements ; i ++) {
	WedgeElement & e = * elements[i];
	fprintf( fp, "%5ld: %.5f\n",
	    i, e.max_pstress);
    }
    fprintf( fp, "....................................................\n");
    fprintf( fp, "Stresses on nodes:\n");
    fprintf( fp, "....................................................\n");
    for( long i = 0 ; i < n_nodes ; i ++) {
	Node & n = * nodes[i];
	fprintf( fp, "%5ld[%1ld]: %10.5f",
	    i, long(n.is_fixed), n.ps_val);
	if( i % 2 == 1) fprintf( fp, "\n");
    }
    fprintf( fp, "....................................................\n");
}

void Model::auto_save( void)
{
    // if auto-save is on, save the temporary results in the file
    // specified on the command line, every so-many seconds (as
    // specified by autosave_interval)
    if( autosave_interval < 0) return;

    static time_t old_t = time( NULL);
    time_t new_t = time( NULL);
    if( new_t - old_t > autosave_interval) {
	// rotate temporary results
	for( long i = autosave_n_keep - 2 ; i >= 0 ; i --) {
	    char src[4096], dst[4096];
	    sprintf( src, autosave_fname_mask, i);
	    sprintf( dst, autosave_fname_mask, i+1);
	    execute( "rm %s", dst);
	    execute( "mv %s %s", src, dst);
	}
	// save the newest file
	char fname[ 4096];
	sprintf( fname, autosave_fname_mask, 0);
	model-> save( fname);
	// remember the current time
	old_t = new_t;
    }
}

void Model::progressive_save( void)
// - this should be called every time before time is incremented
// - if progressive save is turned the current model is saved into
//   a file determined by the progressive_save_fmask
{
    // if progressive save is not turned on, just return
    if( ! progressive_save_on) return;

    // skip desired number of saves
    if( progressive_save_last + progressive_save_skip > time_curr &&
	time_curr < time_total) return;

    fprintf( stderr, "\t- saving at %.10f (last save at %.10f)\n",
	time_curr, progressive_save_last);

    progressive_save_last = time_curr;
	
    // figure out the name
    char fname[ 4096];
    long src = 0;
    long dst = 0;
    long len = long( strlen( progressive_save_fmask));
    while( 1) {
	assert( src <= len);

	// are we done? if yes, exit the loop
	if( src == len) break;

	// figure out the current and next characters
	int c = progressive_save_fmask[ src];
	int c2 = -1;
	if( src + 1 < len) c2 = progressive_save_fmask[ src+1];

	// if current character is not '%' then just copy it
	// and continue
	if( c != '%') {
	    fname[ dst] = c;
	    src += 1;
	    dst += 1;
	    continue;
	}

	// if the current character is '%' then we have couple of
	// choices
	if( c2 == '%') {
	    // if the sequence is '%%' then insert '%'
	    fname[ dst] = '%';
	    src += 2;
	    dst += 1;
	    continue;
	} else if( tolower( c2) == 'n') {
	    // if the sequence is '%n' then insert
	    // the input filename
	    strcpy( fname + dst, input_fname);
	    // but change all slashes to dashes
	    for( long i = 0; i < long( strlen( input_fname)); i ++)
		if( fname[ dst + i] == '/')
		    fname[ dst + i] = '_';
	    src += 2;
	    dst += strlen( input_fname);
	    continue;
	} else if( tolower( c2) == 't') {
	    // if the sequence is '%t' inser the current time
	    char buff[ 4096];
	    sprintf( buff, "%f", time_curr);
	    strcpy( fname + dst, buff);
	    src += 2;
	    dst += strlen( buff);
	    continue;
	} else {
	    fprintf( stderr,
		"WARNING: Invalid sequence '%%%c' in "
		"progressive_save_fmask.\n", c2);
	    strcpy( fname + dst, "<NULL>");
	    src += 2;
	    dst += strlen( "<NULL>");
	    continue;
	}
    }
    // terminate the fname
    fname[ dst] = '\0';
    //	fprintf( stderr, "Would like to progressive save to '%s'\n", fname);

    // save the model
    save( fname);
}

void idle_func( void)
// ----------------------------------------------------------------------
// called whenever GLUT is idle
// ----------------------------------------------------------------------
{
    if(model-> simulation_step()) {
	glutIdleFunc( NULL);
	return;
    }

    // display the results
    glutPostRedisplay();
}

void menu_func( int menu)
// ----------------------------------------------------------------------
// called when some of the menus are selected
// ----------------------------------------------------------------------
{
    switch( menu) {
	case MENU_START_SIMULATION:
	    glutIdleFunc( idle_func);
	    break;
	case MENU_SIMULATION_STEP:
	    model-> simulation_step();
	    glutPostRedisplay();
	    break;
	case MENU_STOP_SIMULATION:
	    glutIdleFunc( NULL);
	    break;
	case MENU_ORIENT_ELEMENTS:
	    model-> orient_elements();
	    glutPostRedisplay();
	    break;
	case MENU_SAVE_MODEL:
	    model-> save( "/tmp/fem-model.dat" );
	    break;
	case MENU_VIEW_MODEL:
	    model-> save( "/tmp/fem-model.dat" );
	    execute( "%s /tmp/fem-model.dat &", viewer_progname.c_str());
	    break;
    }
}

static long pick_object( int x, int y)
// ----------------------------------------------------------------------
// will return the index of the element under the mouse coordinate x,y
// return values of -1 indicates no element
// ----------------------------------------------------------------------
{
    // toggle selection of an element
    const long BUFSIZE = 51200;
    GLuint buff[BUFSIZE];
    GLint n_hits;
    GLint viewport[4];
	
    glGetIntegerv( GL_VIEWPORT, viewport);
    glSelectBuffer( BUFSIZE, buff);
    (void) glRenderMode( GL_SELECT);
	
    glInitNames( );
    glPushName( 0);

    glMatrixMode( GL_PROJECTION);
    glPushMatrix( );
    glLoadIdentity( );
    //  create 5x5 pixel picking region near cursor location
    gluPickMatrix( (GLdouble) x, (GLdouble) (viewport[3] - y),
	1.0, 1.0, viewport);
    long width = viewport[2];
    long height = viewport[3];
    gluPerspective( 60.0, width/double(height), 0.05, 100.0);
    model-> draw();
    glPopMatrix();
    glFlush();
	
    n_hits = glRenderMode(GL_RENDER);
    //	fprintf( stderr, "n_hits = %d\n", n_hits);

    if( n_hits <= 0) return -1;

    GLuint * ptr = (GLuint *) buff;
    // examine each hit, find the one with the smallest Z
    long min_ind = -1;
    double min_z = 0;
    for( long i = 0 ; i < n_hits; i++) {
	long n_names = ptr[0];
	double z1 = ptr[1];
	// double z2 = ptr[2];
	long name = ptr[3];
	ptr += n_names + 3;
	if( min_ind == -1 || min_z > z1) {
	    min_ind = name;
	    min_z = z1;
	}
    }
    if( min_ind == -1) return -1;

    //	fprintf( stderr, "min_ind = %ld\n", min_ind);

    // return the result
    long pick = min_ind - 1;
    return pick;
}

void debug_info( void)
{
    // just for fun, report number of unreferenced nodes
    long n_un = 0;
    for( long i = 0 ; i < n_nodes ; i ++)
	if( nodes[i]-> n_refs == 0)
	    n_un ++;
    fprintf( stderr, "Unreferenced nodes: %ld (%f%%)\n",
	n_un,
	100*double(n_un)/n_nodes);

    fprintf( stderr, "Kg error = %.30f\n",
	model-> get_Kg_error());
    fprintf( stderr, "consistency = %d\n",
	model-> check_ref_consistency());

    fprintf( stderr, "Enter node index:");
    long ind;
    if( scanf( "%ld", & ind) != 1 ||
	ind < 0 || ind >= n_nodes) return;

    Node * n = nodes[ ind];
    if( n-> n_refs == 0) {
	fprintf( stderr, "Nothing attached to this node!\n");
	return;
    }
	
    // calculate all forces acting on this node
    Vector3d F[ n-> n_refs];
    for( long i = 0 ; i < n-> n_refs ; i ++) {
	WedgeElement * e = n-> refs[i];
	Matrix f(18,1);
	e-> calculate_nodal_forces( f);
	long ii = -1;
	for( long k = 0 ; k < 6 ; k ++)
	    if( e-> p[k] == ind) { ii = k; break; }
	assert( ii != -1);
	F[i].set( f.get( ii * 3 + 0, 0),
	    f.get( ii * 3 + 1, 0),
	    f.get( ii * 3 + 2, 0));
    }

    // print out the nodal forces
    for( long i = 0 ; i < n-> n_refs ; i ++)
	fprintf( stderr, "\t%2ld) %10f %10f %10f\n",
	    i, F[i].x, F[i].y, F[i].z);
    // print out nodal stress and yield force
    printf( "\tMax ps = %f (%f,%f,%f), yf = %f\n",
	n-> ps_val,
	n-> ps_nx, n-> ps_ny, n-> ps_nz,
	n-> yield_stress
	    );

    // report tension given the set of forces
	
	
}

static void keyboard_func( unsigned char key, int x, int y)
{
    key = tolower( key);
    if( key == ' ') {
	// pick an object
	long ind = pick_object( x, y);
	//		fprintf( stderr, "picked ind= %ld\n", ind);
	if( ind == -1) return;
	if( ind < n_elements)
	    elements[ ind]-> selected =
		! elements[ ind]-> selected;
	else {
	    ind = ind - n_elements;
	    assert( ind < n_nodes);
	    nodes[ind]-> selected = ! nodes[ind]-> selected;
	}
	// redraw everything
	glutPostRedisplay();
	return ;
    }
    else if( key == 'u') {
	// unselect all objects
	for( long i = 0 ; i < n_elements ; i ++)
	    elements[i]-> selected = 0;
	for( long i = 0 ; i < n_nodes ; i ++)
	    nodes[i]-> selected = 0;
	glutPostRedisplay();
    } else if( key == 'a') {
	// toggle display of axes
	draw_axes = ! draw_axes;
	glutPostRedisplay();
	return ;
    } else if( key == 'c') {
	// toggle the display of colors
	color_display = (color_display + 1) % 3;
	glutPostRedisplay();
    } else if( key == 'r') {
	// toggle display of surface normals
	draw_growth_vectors = ! draw_growth_vectors;
	glutPostRedisplay();
    } else if( key == 'o') {
	// toggle display of original shapes
	draw_original_shapes = (draw_original_shapes+1) % 3;
	glutPostRedisplay();
    } else if( key == '`') {
	// toggle display of deformed / undeformed model
	draw_deformed_model = ! draw_deformed_model;
	glutPostRedisplay();
    } else if( key == 's') {
	// pick an element
	long ind = pick_object( x, y);
	if( ind < 0 || ind >= n_elements) return;
	model-> subdivide_element( ind);
	model-> relax();
	model-> calculate_stresses();
	glutPostRedisplay();
	return ;
    } else if( key == 'm') {
	// pick an element
	long ind = pick_object( x, y);
	if( ind < 0 || ind < n_elements) return;
	ind = ind - n_elements;
	/*
	  for( long i = 0 ; i < n_nodes ; i ++) nodes[i]-> selected = 0;
	  nodes[ ind]-> selected = 1;
	  model-> select_nodes( 3);
	  glutPostRedisplay();
	*/
	fprintf( stderr, "Splitting single node (%ld)...\n", ind);
	int res = model-> split_single_point_node( ind);
	fprintf( stderr, "res = %d\n", res);
	glutPostRedisplay();
	return;
	fprintf( stderr, "plastic zone computing\n");
	model-> mark_plastic_zone( ind);
	for( long i = 0 ; i < n_nodes ; i ++)
	    nodes[i]-> selected = nodes[i]-> is_in_pzone;
	glutPostRedisplay();
	return ;
    } else if( key == 'j') {
	// pick an element
	long ind = pick_object( x, y);
	if( ind < 0 || ind < n_elements) return;
	ind = ind - n_elements;
	fprintf( stderr, "precise stress at node %ld\n", ind);
	model-> calculate_precise_stress_at_node( ind);
	glutPostRedisplay();
	return ;
    } else if( key == 'k') {
	// pick an element
	long ind = pick_object( x, y);
	if( ind < 0 || ind < n_elements) return;
	ind = ind - n_elements;
	fprintf( stderr, "fracturing node %ld\n", ind);
	double p0 = model-> calculate_potential_energy();
	double len; Model::NodeList ftips;
	int res = model-> _fracture_at_node( ind, ftips, len);
	fprintf( stderr, "\tres = %d len = %f\n", res, len);
	model-> relax();
	model-> calculate_stresses();
	double p1 = model-> calculate_potential_energy();
	fprintf( stderr, "\tp0=%f p1=%f dpe=%f\n", p0, p1, p1-p0);
	glutPostRedisplay();
	return ;
    } else if( key == 13) {
	// pick an element
	glutIdleFunc( NULL);
	model-> simulation_step();
	glutPostRedisplay();
	return ;
    } else if( key == 'i') {
	long ind = pick_object( x, y);
	if( ind < 0) return;
	if( ind < n_elements) {
	    WedgeElement * e = elements[ ind];
	    fprintf( stderr, "Info for element %ld\n", ind);
	    fprintf( stderr, "  ym = %f\n", e->young_mod);
	    fprintf( stderr, "  pr = %f\n", e->poisson_mod);
	    fprintf( stderr, "  ys = %f\n", e->yield_stress);
	    fprintf( stderr, "  ft = %f\n", e->fracture_toughness);
	} else {
	    ind = ind - n_elements;
	    Node * n = nodes[ ind];
	    fprintf( stderr, "Info for node %ld\n", ind);
	    fprintf( stderr, "  ym = %f\n", n-> young_modulus);
	    fprintf( stderr, "  pr = %f\n", n-> poisson_ratio);
	    fprintf( stderr, "  ys = %f\n", n-> yield_stress);
	    fprintf( stderr, "  ft = %f\n", n-> ftoughness);
	    fprintf( stderr, "  ps = %.5f (%.5f,%.5f,%.5f)\n",
		n-> ps_val, n-> ps_nx, n-> ps_ny, n-> ps_nz);
	    fprintf( stderr, "   h = %.10f\n",
		model-> get_mat_thickness_at_node( ind));
	}
	// set the error of all selected nodes to
	// some value, and then reselect all nodes which
	// are affected
	//		model-> mark_local_change_nodes( model-> error_radius);
	//		glutPostRedisplay();
	//		debug_info();
    } else if( key == 'b') {
	// pick an element
	long ind = pick_object( x, y);
	fprintf( stderr, "break ind= %ld\n", ind);
	if( ind < 0 || ind >= n_elements) return;
	WedgeElement * e = elements[ind];
	model-> remove_element( e);
	//		model-> add_element( e);
	delete e;
	model-> relax();
	model-> calculate_stresses();
	glutPostRedisplay();
	return ;
    } else if( key == 'w') {
	draw_wireframe = (draw_wireframe +1) % 3;
	glutPostRedisplay();
	return;
    } else if( key == 'l') {
	draw_solid = (draw_solid + 1) % 3;
	glutPostRedisplay();
	return;
    } else if( key == 'p') {
	draw_nodal_stresses = (draw_nodal_stresses + 1)  % 4;
	glutPostRedisplay();
	return;
    } else if( key == '0') {
	model-> relax();
	model-> calculate_stresses();
	glutPostRedisplay();
	return;
    } else if( key == '1') {
	int res = model-> prettify_mesh();
	fprintf( stderr, "prettify_mesh = %d\n", res);
	glutPostRedisplay();
	return;
    } else if( key == '2') {
	for( long i = 0 ; i < n_nodes ; i ++)
	    nodes[i]-> selected = model-> is_node_moveable( i);
	(void) model-> angle_smooth_mesh();
	glutPostRedisplay();
	return;
	/*
	  long ind = -1;
	  for( long i = 0 ; i < n_nodes ; i ++)
	  if( nodes[i]-> selected) ind = i;
	  if( ind == -1) return;
	  model-> refine_mesh_around_node(
	  ind, 2,
	  model-> crack_tip_element_size * 100);
	  glutPostRedisplay();
	  return;
	*/
	/*
	  for( long i = 0 ; i < n_nodes ; i ++) nodes[i]-> selected = 0;
	  for( long i = 0 ; i < n_nodes ; i ++) {
	  Vector3d norm = model-> get_face_normal( i);
	  nodes[i]-> selected = norm.is_zero();
	  }
	  glutPostRedisplay();
	  return;
	*/
	// 		(void) model-> smooth_mesh_around_crack_tips();
	// 		glutPostRedisplay();
	// 		return;
	for( long i = 0 ; i < n_nodes ; i ++)
	    nodes[i]-> selected = model-> is_node_moveable( i);
	(void) model-> angle_smooth_mesh();
	glutPostRedisplay();
	return;
	// find out the selected node
	long n1 = -1;
	for( long i = 0 ; i < n_nodes ; i ++)
	    if( nodes[i]-> selected) {
		n1 = i;
		break;
	    }
	if( n1 == -1) {
	    fprintf( stderr, "MUST SELECT EXACTLY 1 NODE!!!\n");
	    return;
	}
	Model::NodeList * nlist =
	    model-> get_face_neighbors( n1);
	for( long i = 0 ; i < n_nodes ; i ++)
	    nodes[i]-> selected = 0;
	while( ! nlist-> is_empty())
	    nodes[nlist-> pop()]-> selected = 1;
	delete nlist;
	glutPostRedisplay();
	return;
		
	for( long i = 0 ; i < n_nodes ; i ++)
	    nodes[i]-> selected = model-> is_node_moveable( i);
	glutPostRedisplay();
	return;
    } else if( key == '3') {
	fprintf( stderr, "Peel stress (-1=unpeel): ");
	double peel_stress;
	if( 1 != fscanf( stdin, "%lf", & peel_stress)) {
	    fprintf( stderr, "<INPUT ERROR>\n");
	    return;
	}
	// upeel everything
	for( long i = 0 ; i < n_elements ; i ++) {
	    for( long j = 3 ; j < 6 ; j ++) {
		nodes[elements[i]->p[j]]-> is_fixed = 1;
		// 				nodes[elements[i]->p[j]]-> px = 
		// 					nodes[elements[i]->p[j]]-> ox;
		// 				nodes[elements[i]->p[j]]-> py = 
		// 					nodes[elements[i]->p[j]]-> oy;
		// 				nodes[elements[i]->p[j]]-> pz = 
		// 					nodes[elements[i]->p[j]]-> oz;
		nodes[elements[i]-> p[j]]-> grow(
		    model-> time_curr,
		    model-> growth_x,
		    model-> growth_y,
		    model-> growth_z);
	    }
	}
	model-> relax( 1);
	model-> calculate_stresses();
	if( peel_stress > 0)
	    model-> peel( peel_stress);
	glutPostRedisplay();
	return;
    } else if( key == '9') {
	save_screen_dump( "/tmp/test", 0);
	//		glutPostRedisplay();
    } else if( key == '.') {
	if( calculate_nodal_stress_flag) {
	    FILE * fp = fopen( "/tmp/fem-fine.res", "w");
	    assert( fp != NULL);
	    for( long i = 0 ; i < model-> ftips.n ; i ++) {
		long ind = model-> ftips.array[i];
		fprintf( fp, "%ld ", ind);
		fprintf( fp,
		    "%.10f %.10f %.10f %.10f "
		    "%.10f %.10f\n",
		    nodes[ ind]-> stensor.sxx,
		    nodes[ ind]-> stensor.syy,
		    nodes[ ind]-> stensor.szz,
		    nodes[ ind]-> stensor.syz,
		    nodes[ ind]-> stensor.sxz,
		    nodes[ ind]-> stensor.sxy);
	    }
	    fprintf( fp, "-1\n");
	    fclose( fp);
	    exit( 0);
	}
    } else if( key == 'q') {
	draw_elemental_stresses = (draw_elemental_stresses + 1)  % 2;
	glutPostRedisplay();
	return;
    } else if( key == 'z') {
	model-> print_stresses( stderr);
	return;
    } else if( key == 'n') {
	draw_nodes = (draw_nodes + 1) % 3;
	glutPostRedisplay();
	return;
    } else if( key == '=') {
	draw_node_cube_size *= 1.1;
	glutPostRedisplay();
	return;
    } else if( key == '-') {
	draw_node_cube_size *= 0.5;
	glutPostRedisplay();
	return;
    } else if( key == 'y') {
	// show the plastic zone
	for( long i = 0 ; i < n_nodes ; i ++)
	    nodes[i]-> selected = nodes[i]-> is_in_pzone;
	glutPostRedisplay();
    } else if( key == 'f') {
	draw_font_num = (draw_font_num + 1) % draw_n_fonts;
	glutPostRedisplay();
	return;
    } else if( key == 'g') {
	draw_fog = (draw_fog + 1) % 4;
	glutPostRedisplay();
	return;
    } else if( key == 'h') {
	fprintf( stderr,
	    "=================================================\n"
	    "                H e l p\n"
	    ".................................................\n"
	    " <space> toggles an element\n"
	    "     <0> relax\n"
	    "     <1> prettify mesh\n"
	    "     <2> debug stuff\n"
	    "     <3> peel\n"
	    "     <4> split node\n"
	    "     <k> insert sink\n"
	    "     <s> subdivides an element\n"
	    "     <h> displays this help\n"
	    "     <i> prints some debuggin infor\n"
	    "     <a> toggles display of axes\n"
	    "     <w> toggles wireframe\n"
	    "     <l> toggles wall display\n"
	    "     <c> toggles colors used to draw elements\n"
	    "         normal -> stress -> yield_stree\n"
	    "     <p> show stress planes in the nodes\n"
	    "     <q> show stress planes in the elements\n"
	    "     <z> print all stresses (nodes & elements) on\n"
	    "         stderr\n"
	    "     <n> toggle display of node numbers\n"
	    "     <-> halve the size of nodal cubes\n"
	    "     <=> increase (10%%) the size of nodal cubes\n"
	    "     <f> toggle font\n"
	    "     <g> toggle fog\n"
	    "     <r> toggle disp. of background growth vectors\n"
	    "     <o> draw original shapes toggle\n"
	    "     <`> draw deformed/undeformed model toggle\n"
	    "     <y> show plastic zones\n"
	    "     <j> calculate precise stress at node\n"
	    "     <k> fracture at node\n"
	    "=================================================\n"
		 );
    } else {
	fprintf( stderr, "Unknown key pressed: %d\n", int(key));
	return;
    }
}

void Model::load( const char * fname)
{
    // copy the name of the input file
    input_fname = strdup( fname);

    // some default values
    time_total = 1;
    time_curr = 0;
    growth_x = 0.1;
    growth_y = 0.1;
    growth_z = 0;
    shrink_top_t0 = -1;
    shrink_top_t1 = 0;
    shrink_top_val = 0;
    shrink_bot_t0 = -1;
    shrink_bot_t1 = 0;
    shrink_bot_val = 0;
    shrink_height_t0 = -1;
    shrink_height_val0 = 0;
    shrink_height_t1 = 0;
    shrink_height_val1 = 0;
    fracture_tip_inertia = 1.0;
    gravity_x = 0;
    gravity_y = 0;
    gravity_z = 0;
    max_break_size = 1000;
    min_refine_size = 1000;
    min_element_size = -0.01;
    max_element_size = -0.1;
    pzone_size = -4;
    crack_tip_element_size = -0.1;
    crack_tip_sub_element_size = -0.1;
    precision = 1e-10;
    min_time_step = 1e-6;
    max_time_step = 1;
    error_radius = 1;
    min_error_radius = 0.01;
	
    // open the file for reading
    FILE * fp = fopen( fname, "r");
    Tokenizer tok( fp);

    // read in the  the simulation constants, which are in the
    // format:
    //                constant_name = value
    while( 1) {
	//		char * var_name = tok.read_token();
	string var_name = tok.read_token();
	if( tok.error())
	    die( "Unexpected EOF - while reading variables.\n");
	//		if( strcasecmp( var_name, "BEGIN") == 0) break;
	if( var_name == "BEGIN") break;
	//		var_name = strdup( var_name);
	//		char * eq_sign = tok.read_token();
	//		if( eq_sign == NULL || strcmp( eq_sign, "=") != 0)
	if( tok.read_token() != "=")
	    die( "Expected equal sign after variable '%s'",
		var_name.c_str());
	//		char * value = tok.read_token( fp);
	string value = tok.read_token( );
	if( tok.error())
	    die( "Expected a value for variable '%s'",
		var_name.c_str());
	//		value = strdup( value);

	fprintf( stderr, "variable: %s = %s\n", var_name.c_str(),
	    value.c_str());

	// parse the assignment statement
	if( var_name == "sim_time_total") {
	    if( ! tok.is_double( value))
		die( "Value '%s' for sim_time_total is"
		    " not a real number.", value.c_str());
	    time_total = tok.to_double( value);
	} else if( var_name == "sim_time_curr") {
	    if( ! tok.is_double( value))
		die( "Value '%s' for sim_time_curr is"
		    " not a real number.", value.c_str());
	    time_curr = tok.to_double( value);
	} else if( var_name == "growth_x") {
	    if( ! tok.is_double( value))
		die( "Value '%s' for growth_x is"
		    " not a real number.", value.c_str());
	    growth_x = tok.to_double( value);
	} else if( var_name == "growth_y") {
	    if( ! tok.is_double( value))
		die( "Value '%s' for growth_y is"
		    " not a real number.", value.c_str());
	    growth_y = tok.to_double( value);
	} else if( var_name == "growth_z") {
	    if( ! tok.is_double( value))
		die( "Value '%s' for growth_z is"
		    " not a real number.", value.c_str());
	    growth_z = tok.to_double( value);
	} else if( var_name == "shrink_top_t0") {
	    if( ! tok.is_double( value))
		die( "Value '%s' for shrink_top_t0 is"
		    " not a real number.", value.c_str());
	    shrink_top_t0 = tok.to_double( value);
	} else if( var_name == "shrink_top_t1") {
	    if( ! tok.is_double( value))
		die( "Value '%s' for shrink_top_t1 is"
		    " not a real number.", value.c_str());
	    shrink_top_t1 = tok.to_double( value);
	} else if( var_name == "shrink_top_val") {
	    if( ! tok.is_double( value))
		die( "Value '%s' for shrink_top_val is"
		    " not a real number.", value.c_str());
	    shrink_top_val = tok.to_double( value);
	} else if( var_name == "shrink_bot_t0") {
	    if( ! tok.is_double( value))
		die( "Value '%s' for shrink_bot_t0 is"
		    " not a real number.", value.c_str());
	    shrink_bot_t0 = tok.to_double( value);
	} else if( var_name == "shrink_bot_t1") {
	    if( ! tok.is_double( value))
		die( "Value '%s' for shrink_bot_t1 is"
		    " not a real number.", value.c_str());
	    shrink_bot_t1 = tok.to_double( value);
	} else if( var_name == "shrink_bot_val") {
	    if( ! tok.is_double( value))
		die( "Value '%s' for shrink_bot_val is"
		    " not a real number.", value.c_str());
	    shrink_bot_val = tok.to_double( value);
	} else if( var_name == "shrink_height_t0") {
	    if( ! tok.is_double( value))
		die( "Value '%s' for shrink_height_t0 is"
		    " not a real number.", value.c_str());
	    shrink_height_t0 = tok.to_double( value);
	} else if( var_name == "shrink_height_val0") {
	    if( ! tok.is_double( value))
		die( "Value '%s' for shrink_height_val0 is"
		    " not a real number.", value.c_str());
	    shrink_height_val0 = tok.to_double( value);
	} else if( var_name == "shrink_height_t1") {
	    if( ! tok.is_double( value))
		die( "Value '%s' for shrink_height_t1 is"
		    " not a real number.", value.c_str());
	    shrink_height_t1 = tok.to_double( value);
	} else if( var_name == "shrink_height_val1") {
	    if( ! tok.is_double( value))
		die( "Value '%s' for shrink_height_val1 is"
		    " not a real number.", value.c_str());
	    shrink_height_val1 = tok.to_double( value);
	} else if( var_name == "fracture_tip_inertia") {
	    if( ! tok.is_double( value))
		die( "Value '%s' for fracture_tip_inertia is"
		    " not a real number.", value.c_str());
	    fracture_tip_inertia = tok.to_double( value);
	} else if( var_name == "gravity_x") {
	    if( ! tok.is_double( value))
		die( "Value '%s' for gravity_x is"
		    " not a real number.", value.c_str());
	    gravity_x = tok.to_double( value);
	} else if( var_name == "gravity_y") {
	    if( ! tok.is_double( value))
		die( "Value '%s' for gravity_y is"
		    " not a real number.", value.c_str());
	    gravity_y = tok.to_double( value);
	} else if( var_name == "gravity_z") {
	    if( ! tok.is_double( value))
		die( "Value '%s' for gravity_z is"
		    " not a real number.", value.c_str());
	    gravity_z = tok.to_double( value);
	} else if( var_name == "max_break_size") {
	    if( ! tok.is_double( value))
		die( "Value '%s' for max_break_size is"
		    " not a real number.", value.c_str());
	    max_break_size = tok.to_double( value);
	} else if( var_name == "min_refine_size") {
	    if( ! tok.is_double( value))
		die( "Value '%s' for min_refine_size is"
		    " not a real number.", value.c_str());
	    min_refine_size = tok.to_double( value);
	} else if( var_name == "min_element_size") {
	    if( ! tok.is_double( value))
		die( "Value '%s' for min_element_size is"
		    " not a real number.", value.c_str());
	    min_element_size = tok.to_double( value);
	} else if( var_name == "max_element_size") {
	    if( ! tok.is_double( value))
		die( "Value '%s' for max_element_size is"
		    " not a real number.", value.c_str());
	    max_element_size = tok.to_double( value);
	} else if( var_name == "plastic_zone_size") {
	    if( ! tok.is_double( value))
		die( "Value '%s' for plastic_zone_size is"
		    " not a real number.", value.c_str());
	    pzone_size = tok.to_double( value);
	} else if( var_name == "crack_tip_element_size") {
	    if( ! tok.is_double( value))
		die( "Value '%s' for crack_tip_element_size"
		    "size is not a real number.",
		    value.c_str());
	    crack_tip_element_size =
		tok.to_double( value);
	} else if( var_name == "crack_tip_sub_element_size") {
	    if( ! tok.is_double( value))
		die( "Value '%s' for crack_tip_sub_"
		    "element_size is not a real number.",
		    value.c_str());
	    crack_tip_sub_element_size =
		tok.to_double( value);
	} else if( var_name == "precision") {
	    if( ! tok.is_double( value))
		die( "Value '%s' for precision is"
		    " not a real number.", value.c_str());
	    precision = tok.to_double( value);
	} else if( var_name == "min_time_step") {
	    if( ! tok.is_double( value))
		die( "Value '%s' for min_time_step is"
		    " not a real number.", value.c_str());
	    min_time_step = tok.to_double( value);
	} else if( var_name == "max_time_step") {
	    if( ! tok.is_double( value))
		die( "Value '%s' for max_time_step is"
		    " not a real number.", value.c_str());
	    max_time_step = tok.to_double( value);
	} else if( var_name == "error_radius") {
	    if( ! tok.is_double( value))
		die( "Value '%s' for error_radius is"
		    " not a real number.", value.c_str());
	    error_radius = tok.to_double( value);
	} else if( var_name == "min_error_radius") {
	    if( ! tok.is_double( value))
		die( "Value '%s' for min_error_radius is"
		    " not a real number.", value.c_str());
	    min_error_radius = tok.to_double( value);
	} else if( var_name == "progressive_save") {
	    if( value == "yes" || value == "true")
		progressive_save_on = 1;
	    else if( value == "no" || value == "false")
		progressive_save_on = 0;
	    else
		die( "Expected boolean (true|false) for "
		    "progressive save, instead found '%s'.",
		    value.c_str());
	} else if( var_name == "progressive_save_skip") {
	    if( ! tok.is_double( value))
		die( "Value '%s' for progressive_save_skip is"
		    " not a number.", value.c_str());
	    progressive_save_skip =	tok.to_double( value);
	} else if( var_name == "progressive_save_fmask") {
	    progressive_save_fmask = strdup( value.c_str());
	} else if( var_name == "refinement_map") {
	    tok.push( value);
	    refinement_map = new Map( tok);
	} else {
	    die( "Variable '%s' does not exist.",
		var_name.c_str());
	}
    }

    // min_element_size can be specified relative to max_break_size
    if( crack_tip_element_size < 0)
	crack_tip_element_size = (-crack_tip_element_size) *
	    max_break_size;
    if( crack_tip_sub_element_size < 0)
	crack_tip_sub_element_size = (-crack_tip_sub_element_size) *
	    max_break_size;
    if( min_element_size < 0)
	min_element_size = (-min_element_size) *
	    crack_tip_element_size;
    if( max_element_size < 0)
	max_element_size = (-max_element_size) *
	    crack_tip_element_size;
    if( pzone_size < 0)
	pzone_size = (-pzone_size) *
	    crack_tip_element_size;

    // --------------------------------------------------
    // read in all material properties
    // --------------------------------------------------
    // get the BEGIN token
    // 	char * token = Tokenizer::read_token( fp);
    // 	if( token == NULL) {
    // 		fprintf( stderr, "No material properties found at the end.\n");
    // 		return;
    // 	}
    // 	if( strcasecmp( token, "begin") != 0)
    // 		die( "Invalid keyword '%s' in datafile?\n", token);
    // get the name of the section
    string section = tok.read_string( "section");
    if( section == "material_properties") {
	// one after another, read the list of material
	// properties and how to perturb them
	while( 1) {
	    string token = tok.read_token();
	    if( tok.error()) {
		die( "Unexpected end of file.");
	    } else if( token == "young_modulus") {
		ym_map = new Map( fp);
	    } else if( token == "poisson_ratio") {
		pr_map = new Map( fp);
	    } else if( token == "yield_stress") {
		ys_map = new Map( fp);
	    } else if( token == "fracture_toughness") {
		ft_map = new Map( fp);
	    } else if( token == "END") {
		break;
	    } else {
		die( "unexpected symbol '%s'", token.c_str());
	    }
	}
    } else {
	// complain about unknown section
	die( "Unknown section '%s'.", section.c_str());
    }
    // make sure that all randomness maps have been specified
    if( ym_map == NULL) die( "Young modulus map not specified.");
    if( pr_map == NULL) die( "Poisson ratio map not specified.");
    if( ys_map == NULL) die( "Yield stress map not specified.");
    if( ft_map == NULL) die( "Fracture toughness map not specified.");
    fprintf(stderr,"YM: ");ym_map->print( stderr); fprintf( stderr, "\n");
    fprintf(stderr,"PR: ");pr_map->print( stderr); fprintf( stderr, "\n");
    fprintf(stderr,"YS: ");ys_map->print( stderr); fprintf( stderr, "\n");
    fprintf(stderr,"FT: ");ft_map->print( stderr); fprintf( stderr, "\n");
    string token = tok.read_token( fp);
    if( tok.error() || token != "material_properties")
	die( "Mismatching END of section '%s'.", token.c_str());
    fprintf( stderr, "YS: get_min = %f get_max = %f\n",
	ys_map-> get_min(), ys_map-> get_max());
		

    // --------------------------------------------------
    // read in all nodes
    // --------------------------------------------------
    fprintf( stderr, "Reading in nodes:");
    long nn = tok.read_long( "number of nodes");
    for( long i = 0 ; i < nn ; i ++) {
	if( nn < 100 || i % (nn/60+1) == 0)
	    fprintf( stderr, ".");
	Node * n = new Node();
	n-> read_data( tok);
	add_node( n);
    }
    fprintf( stderr, "\n");

    // --------------------------------------------------
    // read in the number of elements
    // --------------------------------------------------
    long ne = tok.read_long( "n_elements");
    assert( ne > 0);

    // read in all elements
    fprintf( stderr, "Reading in elements:");
    for( long i = 0 ; i < ne ; i ++) {
	if( ne < 100 || i % (ne/60+1) == 0)
	    fprintf( stderr, ".");
	// read in the type of element
	string type = tok.read_string( "element_type");
	if( type == "fivewall") {
	    WedgeElement * e = new WedgeElement();
	    e-> read_data( tok);
	    add_element( new WedgeElement(
			     e->p[0],e->p[1],e->p[2],e->p[3],
			     e->p[4],e->p[5],ym_map,pr_map,ys_map,ft_map));
	    delete e;
	} else {
	    die( "Element %ld has unknown type '%s'",
		i, type.c_str());
	}
    }
    fprintf( stderr, "\n");

    // read in fracture tips:
    long nft = tok.read_long( "number of fracture tips");
    for( long i = 0 ; i < nft ; i ++) {
	long ind = tok.read_long( "fracture tip (i)");
	// special case - for sub model, occasionaly it can happen
	// that a fracture tip node is associated with no elements,
	// and then it is translated to '-1'. So we just ignore it.
	if( ind == -1) continue;
	assert( ind >= 0 && ind < n_nodes);
	ftips.add_unique( ind);
    }

    // we are done witht he file
    fclose( fp);
}

void Model::_calculate_elemental_stresses( void)
// ----------------------------------------------------------------------
// calculate elemental stresses inside elements all elements
// 
// Optimization:
//   - only do this for elements that have at least one node marked
// ----------------------------------------------------------------------
{
    double max_pstress = -1e10;
    for( long i = 0 ; i < n_elements ; i ++) {
	if( elements[i]-> has_marked_node())
	    elements[i]-> calculate_stress_tensor();
	if( max_pstress < elements[i]-> max_pstress)
	    max_pstress = elements[i]-> max_pstress;
    }
    fprintf( stderr, "\t- max. elemental stress = %.7f\n", max_pstress);
}

void Model::_calculate_nodal_stresses_tensor_weighing( void)
// ---------------------------------------------------------------------------
// calculates the principal stress at each node in the model and stores the
// result in:
//
//     the stress value in ps_val
//     the normalized vector in ps_nx, ps_ny, ps_nz
//
// also calculate the yield stress for this node, which is the
// average of the yield stresses of the surrounding elements
//
// Method:
//
//     - for each node find all elements to which the node belongs
//     - calculate the stress tensor for those elements at a gaussian
//         point closest to that node
//     - calculate the average of the stress tensor
//     - find the max. eigenvalue for this average
//     - compute the corresponding eigenvector
//     - normalize the eigenvector
//
// Special consideration:
// 
//	- for nodes that belong to no elements (n_refs = 0):
//	         - ps_val = 0, ps_nx = ps_ny = ps_nz = 0;
//               - yield_stress = -1;
//
// Optimization:
//
//      - we can use the knowledge of nodes that have been modified (they
//        are all marked)
//      - we only need to re-calculate the stresses for nodes that
//        are connected to an element who has a node marked
//
// ---------------------------------------------------------------------------
{
    for( long i = 0 ; i < n_nodes ; i ++) {
	if( nodes[i]-> n_refs <= 0) {
	    nodes[i]-> ps_val = 0.0;
	    nodes[i]-> ps_nx = 0.0;
	    nodes[i]-> ps_ny = 0.0;
	    nodes[i]-> ps_nz = 0.0;
	    nodes[i]-> young_modulus = -1;
	    nodes[i]-> poisson_ratio = -1;
	    nodes[i]-> yield_stress = -1;
	    nodes[i]-> ftoughness = -1;
	    continue;
	}
	// see if this node is connected to a marked node
	int is_connected_to_marked_node = 0;
	for( long j = 0 ; j < nodes[i]->n_refs ; j ++) {
	    WedgeElement * e = nodes[i]->refs[j];
	    if( e-> has_marked_node()) {
		is_connected_to_marked_node = 1;
		break;
	    }
	}
	if( ! is_connected_to_marked_node) continue;

	// the node is connecte to an element with a marked node,
	// so we have to do the calculation
	StressTensor TotalStress;
	double total_young_modulus = 0.0;
	double total_poisson_ratio = 0.0;
	double total_yield_stress = 0.0;
	double total_ftoughness = 0.0;
	double total_weight = 0.0;
	for( long j = 0 ; j < nodes[i]-> n_refs ; j ++) {
	    WedgeElement & e = * (WedgeElement *)
		nodes[i]-> refs[j];
	    StressTensor Stress;
	    if( nodal_stresses_from_elements) {
		Stress = e.get_stress_tensor();
	    } else {
		Stress = e.get_nodal_stress_global_ind(i);
	    }
	    // apply a weight to the tensor
	    double weight = 1;
	    //			double weight = e.get_size();
	    //			double weight = pow( e.get_size(), 2);
	    //			double weight = e.get_top_face_min_angle();
	    //			double weight = pow(e.get_top_face_min_angle(),2);
	    //			double weight = e.get_angle_at_global_ind(i);
	    weight *= weight;
	    //			double weight = e.get_quality();
	    Stress.multiply_by_scalar( weight);
	    // add the weighted stress tensor to the total
	    total_weight += weight;
	    TotalStress.add( Stress);
	    total_young_modulus += e.young_mod;
	    total_poisson_ratio += e.poisson_mod;
	    total_yield_stress += e.yield_stress;
	    total_ftoughness += e.fracture_toughness;
	}

	// calculate the average stress tensor
	TotalStress.multiply_by_scalar( 1.0 / total_weight);
	// calculate and store the average material properties
	nodes[i]-> young_modulus =
	    total_young_modulus / nodes[i]-> n_refs;
	nodes[i]-> poisson_ratio =
	    total_poisson_ratio / nodes[i]-> n_refs;
	nodes[i]-> yield_stress =
	    total_yield_stress / nodes[i]-> n_refs;
	nodes[i]-> ftoughness =
	    total_ftoughness / nodes[i]-> n_refs;
		
	// calculate the eigenvalues
	double e1, e2, e3;
	TotalStress.get_eigenvalues( e1, e2, e3);

	// calculate eigenvectors
	Vector3d ev1, ev2;
	long n = TotalStress.get_eigenvectors( e3, ev1, ev2);

	// store the nodal stress with the node
	nodes[i]-> stensor = TotalStress;

	// assign returned eigenvalue (but eliminate negative
	// principal stresses)
	if( e3 > 0)
	    nodes[i]-> ps_val = e3;
	else
	    nodes[i]-> ps_val = 0;

	// assign returned eigenvector
	if( n == 0) {
	    nodes[i]-> ps_nx = 0;
	    nodes[i]-> ps_ny = 0;
	    nodes[i]-> ps_nz = 0;
	} else if( n == 1) {
	    ev1.normalize();
	    ev1.get( nodes[i]-> ps_nx,
		nodes[i]-> ps_ny,
		nodes[i]-> ps_nz);
	} else if( n == 2) {
	    // TBD
	    // 			fprintf( stderr,
	    // 				 "Warning: Model::calculate_nodal_stresses():"
	    // 				 " two eigenvectors\n"
	    // 				 "         1: ");
	    // 			ev1.print( stderr);
	    // 			fprintf( stderr, " 2:");
	    // 			ev2.print( stderr);
	    // 			fprintf( stderr, "\n");

	    ev1.normalize();
	    ev1.get( nodes[i]-> ps_nx,
		nodes[i]-> ps_ny,
		nodes[i]-> ps_nz);
	} else {
	    // TBD
	    // 			fprintf( stderr,
	    // 				 "Warning: Model::calculate_nodal_stresses():"
	    // 				 " %ld eigenvectors obtained - not "
	    // 				 "handling yet\n", n);
	    // 			fprintf( stderr,
	    // 				 "   - eigenvalues: %.10f %.10f %.10f\n",
	    // 				 e1, e2, e3);
	    //			fprintf( stderr, "   - vector1:");
	    //			ev1.print( stderr);
	    //			fprintf( stderr, "\n");
	    //			fprintf( stderr, "   - vector2:");
	    //			ev2.print( stderr);
	    //			fprintf( stderr, "\n");
			
	    ev1.normalize();
	    ev1.get( nodes[i]-> ps_nx,
		nodes[i]-> ps_ny,
		nodes[i]-> ps_nz);
	}
    }

    // adjust ps_val by the fracture_tip_inertia factor for fracture tips
    for( long i = 0 ; i < ftips.n ; i ++) {
	long n = ftips.array[i];
	nodes[ n]-> ps_val *= fracture_tip_inertia;
    }

    // record maximum nodal stress for reporting
    double max_nodal_stress = -1;
    for( long i = 0 ; i < n_nodes ; i ++)
	if( max_nodal_stress < nodes[i]-> ps_val)
	    max_nodal_stress = nodes[i]-> ps_val;

    // report max. nodal stress
    fprintf( stderr, "\t- max. nodal     stress = %.7f\n",
	max_nodal_stress);
}

void Model::_calculate_nodal_stresses_tensor_LSQF_interpolation( void)
// ---------------------------------------------------------------------------
// calculates the principal stress at each node in the model and stores the
// result in:
//
//     max. principal stress value ins ps_val
//     max. principal stress normal in ps_nx, ps_ny, ps_nz
//
// also calculate the yield stress for this node, which is the
// average of the yield stresses of the surrounding elements
//
// Method:
//
//     - for each node find all elements to which the node belongs
//     - calculate the stress tensor for those elements at gaussian points
//     - use least-square-fit to estimate the value of stress tensor
//       at the node
//     - find the max. eigenvalue for this average
//     - compute the corresponding eigenvector
//     - normalize the eigenvector
//
// Special consideration:
// 
//	- for nodes that belong to no elements:
//	         - ps_val = 0, ps_nx = ps_ny = ps_nz = 0;
//               - yield_stress = -1;
//
// Optimization:
//
//      - we can use the knowledge of nodes that have been modified (they
//        are all marked)
//      - we only need to re-calculate the stresses for nodes that
//        are connected to an element who has a node marked
// ---------------------------------------------------------------------------
{

#ifdef DONT_COMPILE
    for( long i = 0 ; i < n_nodes ; i ++) {
	if( nodes[i]-> n_refs <= 0) {
	    // this node does not belong to any element:
	    nodes[i]-> ps_val = 0.0;
	    nodes[i]-> ps_nx = 0.0;
	    nodes[i]-> ps_ny = 0.0;
	    nodes[i]-> ps_nz = 0.0;
	    nodes[i]-> yield_stress = -1;
	    continue;
	}

	// see if this node is connected to an unmarked node
	int is_connected_to_marked_node = 0;
	for( long j = 0 ; j < nodes[i]->n_refs ; j ++) {
	    WedgeElement * e = nodes[i]->refs[j];
	    if( e-> has_marked_node()) {
		is_connected_to_marked_node = 1;
		break;
	    }
	}
	if( ! is_connected_to_marked_node) continue;

	// prepare the LSQF data
	LSQF fit( nodes[i]-> n_refs * WedgeElement::n_gauss);
	for( long j = 0 ; j < nodes[i]-> n_refs ; j ++) {
	    WedgeElement & e = * (WedgeElement *)
		nodes[i]-> refs[j];
	    // for all gaussian points of the element, find
	    // the corresponding x,y,z coordinates, as well
	    // as the stress tensor at that coordinate, and
	    // record them in the lsq-fit structure
	    for( long k = 0 ; k < e.n_gauss ; k ++) {
		//			for( long k = 0 ; k < e.n_gauss_surface ; k ++) {
		// get the xyz coordinates of the gaussian
		// point
		Vector3d pt = e.iso_to_xyz( e.gps[k]);
		// get the stress tensor at the gaussian point
		StressTensor st;
		e.calculate_stress_tensor_at_gp( st, k);
		// add the xyz coordinate with the stress
		// tensor to LSQF
		fit.add( pt, st);
	    }
	}
	// calculate the least square fit at the node coordinates
	Vector3d pt( nodes[i]-> ox,
	    nodes[i]-> oy,
	    nodes[i]-> oz);
	StressTensor st = fit.get_fit( pt, 1);
		
	// store the nodal stress with the node
	nodes[i]-> stensor = st;

	// calculate the eigenvalues
	double e1, e2, e3;
	st.get_eigenvalues( e1, e2, e3);

	// calculate eigenvectors
	Vector3d ev1, ev2;
	long n = st.get_eigenvectors( e3, ev1, ev2);

	// assign returned eigenvalue (but eliminate negative
	// principal stresses)
	if( e3 > 0)
	    nodes[i]-> ps_val = e3;
	else
	    nodes[i]-> ps_val = 0;

	// assign returned eigenvector
	if( n == 0) {
	    fprintf( stderr,
		"Warning: LSQF: ZERO eigenvectors found\n");
	    nodes[i]-> ps_nx = 0;
	    nodes[i]-> ps_ny = 0;
	    nodes[i]-> ps_nz = 0;
	} else if( n == 1) {
	    ev1.normalize();
	    ev1.get( nodes[i]-> ps_nx,
		nodes[i]-> ps_ny,
		nodes[i]-> ps_nz);
	} else if( n == 2) {
	    // TBD
	    fprintf( stderr,
		"Warning: LSQF: TWO eigenvectors found\n");
	    ev1.print( stderr);
	    fprintf( stderr, " 2:");
	    ev2.print( stderr);
	    fprintf( stderr, "\n");

	    ev1.normalize();
	    ev1.get( nodes[i]-> ps_nx,
		nodes[i]-> ps_ny,
		nodes[i]-> ps_nz);
	} else {
	    // TBD
	    fprintf( stderr,
		"Warning: LSQF: THREE eigenvectors found\n");
	    fprintf( stderr,
		"   - eigenvalues: %.10f %.10f %.10f\n",
		e1, e2, e3);
	    fprintf( stderr, "   - vector1:");
	    ev1.print( stderr);
	    fprintf( stderr, "\n");
	    fprintf( stderr, "   - vector2:");
	    ev2.print( stderr);
	    fprintf( stderr, "\n");
			
	    ev1.normalize();
	    ev1.get( nodes[i]-> ps_nx,
		nodes[i]-> ps_ny,
		nodes[i]-> ps_nz);
	}
	// calculate the yield stress - the average of yield
	// stresses from the elements
	Node * nd = nodes[i];
	nd-> yield_stress = 0;
	for( long j = 0 ; j < nd-> n_refs ; j ++) {
	    WedgeElement * e = nd-> refs[j];
	    nd-> yield_stress +=
		e-> yield_stress / nd-> n_refs;
	}
    }
#endif
}

StressTensor calculate_separation_tensor( Vector3d * flist, long n)
// ----------------------------------------------------------------------
// calculates separation tensor as defines by Hodgins and O'Briend
// ----------------------------------------------------------------------
{
    StressTensor res;

    for( long i = 0 ; i < n ; i ++) {
	StressTensor s;
	s.make_eigen_tensor( flist[i]);
	res.add( s);
    }
    // multiply res by 1/2
    res.multiply_by_scalar( 0.5);
    // return result
    return res;
}

void Model::_calculate_nodal_stresses_JH( void)
// ---------------------------------------------------------------------------
// calculates the principal stress at each node in the model and stores the
// result in:
//
//     the stress value in ps_val
//     the normalized vector in ps_nx, ps_ny, ps_nz
//
// also calculate the yield force for this node
//     - it is very important how the yield force is calculated
//     - yield force depends on the orientation of the fracture plane
//     - yield force will be proportional to yield stress of the elements
//       that the fracture plane intersects
//     - it is also proportional to the area of the intersection between
//       the fracture plane and the elements
//     METHOD:
//     - there will be at most two elements that the fracture plane intersets
//     - calculate the amount of force it would take to fracture each
//       element
//          - this is done by calculating the area of the intersection of the
//            fracture plane and the element, times the yield stress of the
//            element
//     - add the forces from both elements
//
// Method:
//
//     - similar to that described by Jessica Hodgins and O'Brien
//     - calculate all forces acting on a node by the connected elements
//     - from these forces, calculate a new tensor
//     - find the largest eigenvalue of this tensor
//     - find the corresponding eigenvector(s)
//     - normalize the eigenvector(s)
//
// Special consideration:
// 
//	- for nodes that belong to no elements:
//	         - ps_val = 0, ps_nx = ps_ny = ps_nz = 0;
//               - yield_stress = -1;
// ---------------------------------------------------------------------------
{
    // for reporting purposes, keep track of max ps_val
    double max_ps_val = -1e10;
    for( long i = 0 ; i < n_nodes ; i ++) {
	if( nodes[i]-> n_refs <= 0) {
	    nodes[i]-> ps_val = 0.0;
	    nodes[i]-> ps_nx = 0.0;
	    nodes[i]-> ps_ny = 0.0;
	    nodes[i]-> ps_nz = 0.0;
	    nodes[i]-> yield_stress = -1;
	    continue;
	}
	// calculate the array of forces f[]
	Vector3d f[ nodes[i]-> n_refs];
	for( long j = 0 ; j < nodes[i]-> n_refs ; j ++) {
	    WedgeElement & e = * (WedgeElement *)
		nodes[i]-> refs[j];
	    // get the nodal forces F for the element
	    Matrix F(18,1);
	    e.calculate_nodal_forces( F);
	    // find the internal index ii of node i in element e
	    long ii = -1;
	    for( long k = 0 ; k < 6 ; k ++)
		if( e.p[k] == i) { ii = k; break; }
	    assert( ii != -1);
	    // extract the force vector from F
	    f[j].set( F.get( ii * 3 + 0, 0),
		F.get( ii * 3 + 1, 0),
		F.get( ii * 3 + 2, 0));
	    if( 1) {
		double area = e.get_projected_area( f[j]);
		if( area > 1e-10)
		    f[j].scale( 1/area);
		else
		    f[j].set( 0, 0, 0);
	    } else {
		f[j].scale( 1/pow(e.get_size(),2/3.0));
	    }
	}
	// Calculate a tensor from the force vectors. This tensor
	// is something similar to the separation tensor described
	// by Hodgins & O'Brien.
	StressTensor s = calculate_separation_tensor(
	    f, nodes[i]-> n_refs);
	// extract eigenvalues
	double s1, s2, s3;
	s.get_eigenvalues( s1, s2, s3);
	nodes[i]-> ps_val = s3;
	if( nodes[i]-> ps_val > max_ps_val)
	    max_ps_val = nodes[i]-> ps_val;
	// get the eigenvector(s) for the max. eigenvalue 's3'
	Vector3d v1, v2;
	int r = s.get_eigenvectors( s3, v1, v2);
	if( r == 0) {
	    fprintf( stderr, "Jessica: 0 eigenvectors!\n");
	    nodes[i]-> ps_nx = 0;
	    nodes[i]-> ps_ny = 0;
	    nodes[i]-> ps_nz = 0;
	} else if( r == 1) {
	    v1.normalize();
	    v1.get( nodes[i]-> ps_nx,
		nodes[i]-> ps_ny,
		nodes[i]-> ps_nz);
	} else if( r == 2) {
	    fprintf( stderr, "Jessica: 2 eigenvectors (random)\n");
	    double alpha = drand48() * 2 * M_PI;
	    v1.scale( sin( alpha));
	    v2.scale( cos( alpha));
	    v1.add( v2);
	    v1.normalize();
	    v1.get( nodes[i]-> ps_nx,
		nodes[i]-> ps_ny,
		nodes[i]-> ps_nz);
	} else if( r == 3) {
	    fprintf( stderr, "Jessica: 3 eigenvectors (random)\n");
	    double alpha = drand48() * 2 * M_PI;
	    v1.scale( sin( alpha)); v1.normalize();
	    v2.scale( cos( alpha)); v2.normalize();
	    v1.add( v2);
	    v1.normalize();
	    v1.get( nodes[i]-> ps_nx,
		nodes[i]-> ps_ny,
		nodes[i]-> ps_nz);
	}

	// calculate the yield force - depending on the plane
	// of max. stress (its intersection with each element)
	Node * n = nodes[i];
	n-> yield_stress = 0;
	for( long j = 0 ; j < n-> n_refs ; j ++) {
	    WedgeElement * e = n-> refs[j];
	    if( 0)
		nodes[i]-> yield_stress += e-> yield_stress *
		    e-> plane_intersection(
			i, Vector3d( n-> ps_nx,
			    n-> ps_ny,
			    n-> ps_nz));
	    else
		nodes[i]-> yield_stress +=
		    e-> yield_stress / n-> n_refs;
	}
    }
    fprintf( stderr, "Jessica: max_ps_val = %f\n", max_ps_val);
}

void Model::_calculate_nodal_stresses( void)
{
    _calculate_nodal_stresses_tensor_weighing();
    //	_calculate_nodal_stresses_tensor_LSQF_interpolation();
    //	_calculate_nodal_stresses_JH();
}

void Model::calculate_stresses( void)
{
    clock_t start = clock ();

    _calculate_elemental_stresses();
    if( break_method == BreakNode) {
	_calculate_nodal_stresses();
    }

    if( ! calculate_nodal_stress_flag) {
	static double old_pe = 0;
	double pe = calculate_potential_energy();
	fprintf( stderr,
	    "\t- potential energy = pe= %10f dpe= %10f\n",
	    pe, pe-old_pe);
	old_pe = pe;
    }

    clock_t end = clock ();
    stress_time += end - start;

}

void Model::calculate_precise_stress_at_node( long n0)
// ----------------------------------------------------------------------
// attempts to calculate more precise nodal stresses at node n0
//
// Algorithm:
//     -  create a sub-model, with elements superfinely refined around
//        the node
// ----------------------------------------------------------------------
{
    if (! precise_nodal_stress_flag) {
	// just a trick (select only this node, call stress
	// calculation and then reset the selected nodes...
	char sel [n_nodes];
	for (long i = 0 ; i < n_nodes ; i ++) {
	    sel[i] = nodes[i]-> selected;
	    nodes[i]-> selected = 0;
	}
	nodes[n0]-> selected = 1;
	calculate_stresses ();
	for (long i = 0 ; i < n_nodes ; i ++) {
	    nodes[i]-> selected = sel[i];
	}
	return ;
    }
    // select nodes around the node
    // -------------------------------------
    // unselect all nodes
    for( long i = 0 ; i < n_nodes ; i ++) nodes[i]-> selected = 0;
    // select the node
    nodes[ n0]-> selected = 1;
    // select nodes in a radius around selected node
    //	double crack_tip_refinement_radius = 0.1;
    double crack_tip_refinement_radius = crack_tip_element_size * 2;
    mark_local_change_nodes( crack_tip_refinement_radius);

    // select all elements that have at least one node selected
    // ---------------------------------------------------------
    long sub_n_elements = 0;
    for( long i = 0 ; i < n_elements ; i ++) {
	WedgeElement * e = elements[i];
	e-> selected =
	    nodes[e-> p[0]]-> selected ||
	    nodes[e-> p[1]]-> selected ||
	    nodes[e-> p[2]]-> selected ||
	    nodes[e-> p[3]]-> selected ||
	    nodes[e-> p[4]]-> selected ||
	    nodes[e-> p[5]]-> selected;
	if( e-> selected)
	    sub_n_elements ++;
    }
    // it is possible that no elements are actually selected (if n0
    // happens to have n_refs == 0). In this case just return.
    if( sub_n_elements == 0) return;
    // prepare a translation table for the nodes
    long tran[ n_nodes];
    for( long i = 0 ; i < n_nodes ; i ++) tran[i] = -1;
    long sub_n_nodes = 0;
    for( long i = 0 ; i < n_elements ; i ++) {
	WedgeElement * e = elements[i];
	if( ! e-> selected) continue;
	for( long j = 0 ; j < 6 ; j ++)
	    if( tran[ e-> p[j]] == -1) {
		tran[ e-> p[j]] = sub_n_nodes;
		sub_n_nodes ++;
	    }
    }
    // prepare an inverse translation table
    long itran[ sub_n_nodes + 1];
    for( long i = 0 ; i < n_nodes ; i ++) {
	if( tran[i] != -1)
	    itran[ tran[i]] = i;
    }
    // save the sub-model
    FILE * fp = fopen( "/tmp/fem-fine.dat", "w");
    fprintf( fp, "# sub-model\n");
    _save_header( fp);
    fprintf( fp, "\n%ld # number of nodes in sub-model\n",
	sub_n_nodes);
    for( long j = 0 ; j < sub_n_nodes ; j ++) {
	long i = itran[j];
	int fixit = (! nodes[i]-> selected) || nodes[i]-> is_fixed;
	if( fixit)
	    fprintf( fp, "%.10f %.10f %.10f %.10f %.10f %.10f "
		"%.10f %.10f %.10f %ld\n",
		nodes[i]-> ox,
		nodes[i]-> oy,
		nodes[i]-> oz,
		nodes[i]-> px - nodes[i]-> ox,
		nodes[i]-> py - nodes[i]-> oy,
		nodes[i]-> pz - nodes[i]-> oz,
		nodes[i]-> gv.x,
		nodes[i]-> gv.y,
		nodes[i]-> gv.z,
		long( 2));
	else
	    fprintf( fp, "%.10f %.10f %.10f %.10f %.10f %.10f "
		"%.10f %.10f %.10f %ld\n",
		nodes[i]-> ox,
		nodes[i]-> oy,
		nodes[i]-> oz,
		nodes[i]-> px - nodes[i]-> ox,
		nodes[i]-> py - nodes[i]-> oy,
		nodes[i]-> pz - nodes[i]-> oz,
		nodes[i]-> gv.x,
		nodes[i]-> gv.y,
		nodes[i]-> gv.z,
		long( nodes[i]-> is_fixed));
    }

    fprintf( fp, "\n%ld # number of elements in sub-model\n",
	sub_n_elements);
    for( long i = 0 ; i < n_elements ; i ++) {
	WedgeElement * e = elements[i];
	if( ! e-> selected) continue;
	fprintf( fp, "fivewall %ld %ld %ld %ld %ld %ld %ld\n",
	    i,
	    tran[e-> p[0]], tran[e-> p[1]], tran[e-> p[2]],
	    tran[e-> p[3]], tran[e-> p[4]], tran[e-> p[5]]);
    }
    fprintf( fp, "\n1 # number of fracture tips\n");
    fprintf( fp, "%ld\n", tran[ n0]);
    fclose( fp);

    // now call myself on this model
    fprintf (stderr, "\t- executing myself...");
    int res = execute(
	"%s -calculate_nodal_stress /tmp/fem-fine.dat "
	"1>/tmp/fem-fine.log 2>&1", progname.c_str());
    fprintf (stderr, "done\n");
    if( res) die( "Sub-process failed...\n");
    // read in the result
    fp = fopen( "/tmp/fem-fine.res", "r");
    assert( fp != NULL);
    Tokenizer tok( fp);

    // read in the result
    long sub_ind = tok.read_long( "fracture point");
    assert( sub_ind >= 0 && sub_ind < sub_n_nodes + 1);
    long ind = itran[ sub_ind];
    assert( ind == n0);
    StressTensor st;
    st.sxx = tok.read_double( "sxx");
    st.syy = tok.read_double( "syy");
    st.szz = tok.read_double( "szz");
    st.syz = tok.read_double( "syz");
    st.sxz = tok.read_double( "sxz");
    st.sxy = tok.read_double( "sxy");

    // calculate the eigenvalues
    double e1, e2, e3;
    //	ns_debug_flag = 1;
    st.get_eigenvalues( e1, e2, e3);
    // calculate eigenvectors
    Vector3d ev1, ev2;
    long n = st.get_eigenvectors( e3, ev1, ev2);
    // store the nodal stress with the node
    nodes[ind]-> stensor = st;
    // assign returned eigenvalue (but eliminate negative
    // principal stresses)
    if( e3 > 0)
	nodes[ind]-> ps_val = e3 * fracture_tip_inertia;
    else
	nodes[ind]-> ps_val = 0;
    // assign returned eigenvector
    if( n == 0) {
	nodes[ind]-> ps_nx = 0;
	nodes[ind]-> ps_ny = 0;
	nodes[ind]-> ps_nz = 0;
    } else {
	ev1.normalize();
	ev1.get( nodes[ind]-> ps_nx,
	    nodes[ind]-> ps_ny,
	    nodes[ind]-> ps_nz);
    }
    // read in the material properties of the node
    nodes[ind]-> young_modulus = tok.read_double( fp, "ym");
    nodes[ind]-> poisson_ratio = tok.read_double( fp, "pr");
    nodes[ind]-> yield_stress = tok.read_double( fp, "ys");
    nodes[ind]-> ftoughness = tok.read_double( fp, "ft");
    // only one node needed (there should be no other nodes actually)
    ind = tok.read_long( fp, "fracture point (0)");
    assert( ind < 0);
		
    fclose( fp);
}

double
    Model::calculate_potential_energy( void)
{
    // calculate the potential energy
    double pe = 0;
    Matrix f( 18, 1);
    Matrix u( 18, 1);
    for( long i = 0 ; i < n_elements ; i ++) {
	elements[i]-> calculate_nodal_forces( f);
	elements[i]-> calculate_Disp( u);
	double r = 0;
	for( long j = 0 ; j < 18 ; j ++)
	    r += u.data[j][0] * f.data[j][0];
	pe += r;
    }
    return pe;
}

int main( int argc, char ** argv)
{
    // parse command line parameters
    parse_command_line( argc, argv);

    // if memory optimization is requested, set caching of Ke's
    if( optimize_memory_flag)
	WedgeElement::cache_Ke = 0;
    else
	WedgeElement::cache_Ke = 1;

    // set the fast KE removal option
    if( fast_ke_removal_flag)
	WedgeElement::fast_Ke_removal = 1;
    else
	WedgeElement::fast_Ke_removal = 0;

    // initialize a new model and read the datafile
    model = new Model();
    model-> load( fname);
    if( b_method == BreakElementMethod)
	model-> set_break_method( Model::BreakElement);
    else if( b_method == BreakNodeMethod)
	model-> set_break_method( Model::BreakNode);
    else if( b_method == BreakEdgeMethod)
	model-> set_break_method( Model::BreakEdge);

    // report max. size
    double max_size = -1;
    for (long i = 0 ; i < n_elements ; i ++)
	max_size = max (max_size, elements[i]-> get_size ());
    fprintf (stderr, "%ld elements - max. size = %f\n", n_elements, max_size);

    // subdivide elements until all of them are smaller than
    // 'max_element_size' && according to refinement map
    fprintf( stderr, "Applying refinement map...\n");
    long round = 1;
    while( 1) {
	long n_marked = 0;
	for( long i = 0 ; i < n_elements ; i ++) {
	    WedgeElement * e = elements[i];
	    e-> mark = 0;
	    double desired_size = model-> max_element_size;
	    if( model-> refinement_map != NULL)
		desired_size = model-> refinement_map->
		    get_min_value(
			nodes[e-> p[0]]-> ox,
			nodes[e-> p[0]]-> oy,
			nodes[e-> p[0]]-> oz,
			nodes[e-> p[1]]-> ox,
			nodes[e-> p[1]]-> oy,
			nodes[e-> p[1]]-> oz,
			nodes[e-> p[2]]-> ox,
			nodes[e-> p[2]]-> oy,
			nodes[e-> p[2]]-> oz);
	    if( elements[i]-> get_size() > desired_size)
	    {
		elements[i]-> mark = 1;
		n_marked ++;
	    }
	}
	if( n_marked == 0) break;
	fprintf( stderr, "\t- round %ld: subdividing %ld elements\n",
	    round, n_marked);
	model-> subdivide_marked_elements();
	round ++;
    }
    fprintf( stderr, "\t- done.\n");

    if( model-> error_radius < 0) {
	model-> error_radius = - model-> error_radius;
	fprintf( stderr, "\t- computing the multiresolution mesh...");
	for( long i = 0 ; i < model-> ftips.n ; i ++)
	    model-> refine_mesh_around_node(
		model-> ftips.array[i],
		2,
		model-> crack_tip_sub_element_size);
	fprintf( stderr, "done.\n");
    }
	
    // find equilibrium for the model
    // - this is actually only needed to calculate the principal stresses
    fprintf( stderr, "Initial growth:\n");
    model-> grow( model-> time_curr);
    fprintf( stderr, "\t- done.\n");
    fprintf( stderr, "Initial relax:\n");
    model-> relax( 1);
    model-> calculate_stresses();
    fprintf( stderr, "\t- done.\n");

    if( 0) {
	fprintf( stderr, "Angle at 0 = %f\n",
	    elements[0]-> get_angle_at_global_ind( 0));
	fprintf( stderr, "Angle at 1 = %f\n",
	    elements[0]-> get_angle_at_global_ind( 1));
	fprintf( stderr, "Angle at 2 = %f\n",
	    elements[0]-> get_angle_at_global_ind( 2));
	fprintf( stderr, "Angle at 3 = %f\n",
	    elements[0]-> get_angle_at_global_ind( 3));
	fprintf( stderr, "Angle at 4 = %f\n",
	    elements[0]-> get_angle_at_global_ind( 4));
	fprintf( stderr, "Angle at 5 = %f\n",
	    elements[0]-> get_angle_at_global_ind( 5));
	fprintf( stderr, "Quality = %f\n",
	    elements[0]-> get_quality());
	//		return 0;
    }

    if( profile_on) {
	while( ! model-> simulation_step());

	// save the result
	char fname[ 4096];
	sprintf( fname, autosave_fname_mask, -1);
	model-> save( fname);
	return 0;
    }
	
    if( calculate_nodal_stress_flag) {
	// refine the mesh around crack tips
	fprintf( stderr, "\t- computing the multiresolution mesh...");
	for( long i = 0 ; i < model-> ftips.n ; i ++)
	    model-> refine_mesh_around_node(
		model-> ftips.array[i],
		2,
		model-> crack_tip_sub_element_size);
	fprintf( stderr, "done.\n");
	// relax the model (and re-calculate the stresses)
	model-> relax();
	model-> calculate_stresses();
	// put the result into a file
	FILE * fp = fopen( "/tmp/fem-fine.res", "w");
	assert( fp != NULL);
	for( long i = 0 ; i < model-> ftips.n ; i ++) {
	    long ind = model-> ftips.array[i];
	    fprintf( fp, "%ld # node number\n", ind);
	    fprintf( fp, "\t# stress tensor:\n");
	    fprintf( fp,
		"\t%.10f %.10f %.10f %.10f %.10f %.10f\n",
		nodes[ ind]-> stensor.sxx,
		nodes[ ind]-> stensor.syy,
		nodes[ ind]-> stensor.szz,
		nodes[ ind]-> stensor.syz,
		nodes[ ind]-> stensor.sxz,
		nodes[ ind]-> stensor.sxy);
	    fprintf( fp, "\t%.10f # young modulus\n", 
		nodes[ ind]-> young_modulus);
	    fprintf( fp, "\t%.10f # poisson ratio\n",
		nodes[ ind]-> poisson_ratio);
	    fprintf( fp, "\t%.10f # yield stress\n",
		nodes[ ind]-> yield_stress);
	    fprintf( fp, "\t%.10f # fracture toughness\n",
		nodes[ ind]-> ftoughness);
	}
	fprintf( fp, "-1\n");
	fclose( fp);
	// we are done
	exit( 0);
    }

    // open a window
    glutInit( &argc, argv );
    //	glutInitWindowSize( 800, 500);
    glutInitWindowSize( 1024, 768);
    glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    char title[ 4096];
    sprintf( title, "%s %s", argv[0], fname);
    main_window = glutCreateWindow( title);
    assert( main_window);
    glutSetWindow( main_window);

    glutDisplayFunc( disp_func);
    glutReshapeFunc( reshape_func);
    glutMouseFunc( mouse_func);
    glutMotionFunc( mouse_motion_func);
    glutKeyboardFunc( keyboard_func);

    /*
      view.set( Vector3d( 0, 1, 1), // eye
      Vector3d( 0, 0, 0), // center
      Vector3d( 0, 0, 1)  // up-vector
      );
    */
    view.set( Vector3d( 0, 0, 1.1), // eye
	Vector3d( 0, 0, 0), // center
	Vector3d( 0, 1, 0)  // up-vector
	      );
        
    // create a popup menu
    int m = glutCreateMenu( menu_func);
    glutSetMenu( m);
    glutAddMenuEntry( "Simulation step", MENU_SIMULATION_STEP);
    glutAddMenuEntry( "Start simulation", MENU_START_SIMULATION);
    glutAddMenuEntry( "Stop simulation", MENU_STOP_SIMULATION);
    glutAddMenuEntry( "Orient elements", MENU_ORIENT_ELEMENTS);
    glutAddMenuEntry( "Save model", MENU_SAVE_MODEL);
    glutAddMenuEntry( "View model", MENU_VIEW_MODEL);
    glutAttachMenu( GLUT_RIGHT_BUTTON);

    glutMainLoop();
}
