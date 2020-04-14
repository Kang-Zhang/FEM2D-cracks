#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <string.h>
#include <ctype.h>
#include <malloc.h>
#include <vector>
#include <sys/types.h>
#include <unistd.h>
#include <GL/gl.h>
#include <GL/glut.h>
#include <stdarg.h>
#include "die.hh"
#include "Tokenizer.hh"
#include "Vector3d.hh"
#include "glDrawCube.hh"
#include "gl2ps.hh"

char * fname = NULL;
int antialiasing = true;
int draw_axes = 0;
int draw_wireframe = 0;
int draw_solid = 2;
int use_deformed_coordinates = 1;
int draw_type = 4;
const long n_draw_types = 5;
double line_width = 0.1;

long window_width;
long window_height;

void usage( void)
{
    die( "Usage: viewer fname\n");
}

static void glNormal (const Vector3d & v) {
    glNormal3d (v.x, v.y, v.z);
}

static void glVertex (const Vector3d & v) {
    glVertex3d (v.x, v.y, v.z);
}


class Node {
 public:
    Vector3d pos;		// node position
    std::vector <long> els;	// list of element indices that have
				// this node
};
std::vector <Node> nodes;

Vector3d * top_norms = NULL, * side_norms = NULL;

class Element {
 public:
    long p[6];
    bool contains_node (long ind) const
    {
	return ind == p[0] || ind == p[1] || ind == p[2]
	    || ind == p[3] || ind == p[4] || ind == p[5];
    }
};
std::vector <Element> elements;

class Triangle
{
 public:
    long p1, p2, p3;	// node references
    Vector3d n1, n2, n3;	// normals
    long ns1, ns2, ns3;	// number of shared vertices

    Vector3d get_normal () const
    {
	Vector3d & v1 = nodes[p1].pos;
	Vector3d & v2 = nodes[p2].pos;
	Vector3d & v3 = nodes[p3].pos;
	Vector3d norm = cross_product (v2-v1, v3-v1);
	norm.normalize ();
	return norm;
    }
	
    Triangle (long pp1, long pp2, long pp3)
	: p1 (pp1)
	, p2 (pp2)
	, p3 (pp3)
    {;}
};

class Quad
{
 public:
    long p1, p2, p3, p4; // node references

    Vector3d get_normal () const
    {
	Vector3d & v1 = nodes[p1].pos;
	Vector3d & v2 = nodes[p2].pos;
	Vector3d & v3 = nodes[p3].pos;
	Vector3d & v4 = nodes[p4].pos;
	Vector3d norm = (
	    cross_product (v2-v1, v3-v1)
	    + cross_product (v3-v2, v4-v2)
	    + cross_product (v4-v3, v1-v3)
	    + cross_product (v1-v4, v2-v4)
	    ) / 4.0;
	norm.normalize ();
	return norm;
    }
	
    Quad (long pp1, long pp2, long pp3, long pp4)
	: p1 (pp1), p2 (pp2), p3 (pp3), p4 (pp4)
    {;}
};

static void glTriangle (const Triangle & t)
{
	glVertex (nodes[t.p1].pos);
	glVertex (nodes[t.p2].pos);
	glVertex (nodes[t.p3].pos);
}

static void glTriangleNormal (const Triangle & t)
{
    glNormal (t.get_normal ());
    glTriangle (t);
}

static void glQuad (const Quad & t)
{
	glVertex (nodes[t.p1].pos);
	glVertex (nodes[t.p2].pos);
	glVertex (nodes[t.p3].pos);
	glVertex (nodes[t.p4].pos);
}

static void glQuadNormal (const Quad & t)
{
    glNormal (t.get_normal ());
    glQuad (t);
}

typedef std::vector <Triangle> TList;

void calculate_normals (TList & list, double angle_threshold)
{
    // for each node, figure out which triangles share it
    std::vector <std::vector <size_t> > nlist;
    nlist.resize (nodes.size ());
    for (size_t i = 0 ; i < list.size () ; i ++)
    {
	Triangle & t = list [i];
	nlist [t.p1].push_back (i);
	nlist [t.p2].push_back (i);
	nlist [t.p3].push_back (i);
    }

    // for each triangle's vertex, calculate the normal
    for (size_t i = 0 ; i < list.size () ; i ++)
    {
	Triangle & t = list [i];
	// calculate the normal of this triangle
	Vector3d myNorm = t.get_normal ();
	//
	// consider the first vertex
	// ==================================================
	Vector3d sum1 = myNorm;
	long nsum1 = 1;
	// consider normals of all other triangles sharing this node
	for (size_t j = 0 ; j < nlist [t.p1].size () ; j ++)
	{
	    size_t nt = nlist[t.p1][j];
	    // skip 'myself'
	    if (nt == i) continue;
	    // get the triangle's normal
	    Vector3d hisNorm = list [nt].get_normal ();
	    // what is the angle
	    double alpha = Vector3d::angle (myNorm, hisNorm);
	    // if the angle is too big, ignore this neighbor
	    if (alpha > angle_threshold) continue;
	    // otherwise, blend his normal into mine
	    sum1 = sum1 + hisNorm;
	    nsum1 ++;
	}
	// calculate the blended normal
	t.n1 = sum1 / nsum1;
	t.n1.normalize ();
	//
	// consider the second vertex
	// ==================================================
	Vector3d sum2 = myNorm;
	long nsum2 = 1;
	// consider normals of all other triangles sharing this node
	for (size_t j = 0 ; j < nlist [t.p2].size () ; j ++)
	{
	    size_t nt = nlist[t.p2][j];
	    // skip 'myself'
	    if (nt == i) continue;
	    // get the triangle's normal
	    Vector3d hisNorm = list [nt].get_normal ();
	    // what is the angle
	    double alpha = Vector3d::angle (myNorm, hisNorm);
	    // if the angle is too big, ignore this neighbor
	    if (alpha > angle_threshold) continue;
	    // otherwise, blend his normal into mine
	    sum2 = sum2 + hisNorm;
	    nsum2 ++;
	}
	// calculate the blended normal
	t.n2 = sum2 / nsum2;
	t.n2.normalize ();
	//
	// consider the third vertex
	// ==================================================
	Vector3d sum3 = myNorm;
	long nsum3 = 1;
	// consider normals of all other triangles sharing this node
	for (size_t j = 0 ; j < nlist [t.p3].size () ; j ++)
	{
	    size_t nt = nlist[t.p3][j];
	    // skip 'myself'
	    if (nt == i) continue;
	    // get the triangle's normal
	    Vector3d hisNorm = list [nt].get_normal ();
	    // what is the angle
	    double alpha = Vector3d::angle (myNorm, hisNorm);
	    // if the angle is too big, ignore this neighbor
	    if (alpha > angle_threshold) continue;
	    // otherwise, blend his normal into mine
	    sum3 = sum3 + hisNorm;
	    nsum3 ++;
	}
	// calculate the blended normal
	t.n3 = sum3 / nsum3;
	t.n3.normalize ();
    }
}

TList get_top_triangles ()
{
    // add the top face of each element
    TList list;
    for (size_t i = 0 ; i < elements.size () ; i ++)
    {
	Element & e = elements [i];
	list.push_back (Triangle (e.p[0], e.p[1], e.p[2]));
    }
    return list;
}

bool triangle_is_shared (const Triangle & t)
{
    // count how many elements contain all 3 nodes of t
    long count = 0;
    for (size_t i = 0 ; i < nodes [t.p1].els.size () ; i ++)
    {
	Element & e = elements [nodes [t.p1].els [i]];
	if (! e.contains_node (t.p1)) continue;
	if (! e.contains_node (t.p2)) continue;
	if (! e.contains_node (t.p3)) continue;
	count ++;
    }
    return count > 1;
}

bool is_top_edge_shared (long n1, long n2)
{
    long count = 0;
    for (size_t i = 0 ; i < nodes [n1].els.size () ; i ++)
    {
	Element & e = elements [nodes [n1].els [i]];
	if (! e.contains_node (n2)) continue;
	count ++;
    }
    return count > 1;
}

TList get_side_triangles ()
{
    // add the top face of each element
    TList list;
    for (size_t i = 0 ; i < elements.size () ; i ++)
    {
	Element & e = elements [i];
	if (! is_top_edge_shared (e.p[0],e.p[1]))
	{
	    list.push_back (Triangle (e.p[3], e.p[4], e.p[1]));
	    list.push_back (Triangle (e.p[3], e.p[1], e.p[0]));
	}
	if (! is_top_edge_shared (e.p[1],e.p[2]))
	{
	    list.push_back (Triangle (e.p[4], e.p[5], e.p[2]));
	    list.push_back (Triangle (e.p[4], e.p[2], e.p[1]));
	}
	if (! is_top_edge_shared (e.p[2],e.p[0]))
	{
	    list.push_back (Triangle (e.p[5], e.p[3], e.p[0]));
	    list.push_back (Triangle (e.p[5], e.p[0], e.p[2]));
	}
    }
    return list;
}

TList top_triangles, side_triangles;

class EdgeList {
 private:
    class Edge {
	long n_data;
	long * data;
    public:
	Edge() {
	    data = NULL;
	    n_data = 0;
	}
	~Edge() {
	    if( data != NULL)
		free( data);
	}
	void add( long n) {
	    if( n_data == 0)
		data = (long *) malloc( sizeof( long));
	    else
		data = (long *) realloc(
		    data, sizeof( long) * (n_data+1));
	    data[n_data] = n;
	    n_data ++;
	}
	long count( long n) {
	    long res = 0;
	    for( long i = 0 ; i < n_data ; i ++)
		if( data[i] == n) res ++;
	    return res;
	}
    };
    Edge * edges;
    long n_edges;
 public:
    EdgeList( long p_n_edges) {
	n_edges = p_n_edges;
	edges = new Edge[ n_edges];
    }
    ~EdgeList() {
	if( edges != NULL)
	    delete [] edges;
    }
    void add( long n1, long n2) {
	if( n1 < n2) edges[n1].add(n2);
	else edges[n2].add(n1);
    }
    long count( long n1, long n2) {
	if( n1 < n2) return edges[n1].count(n2);
	else return edges[n2].count(n1);
    }
};

class TriangleList {
 private:
    long n_triangles;
    struct Tri {
	long n1, n2, n3;
    };
    Tri * data;
 public:
    TriangleList() {
	data = NULL;
	n_triangles = 0;
    }
    ~TriangleList() {
	if( data != NULL)
	    free( data);
    }
    void add( long n1, long n2, long n3) {
	if( n_triangles == 0)
	    data = (Tri *) malloc( sizeof( Tri));
	else
	    data = (Tri *) realloc( data, sizeof( Tri)
		* (n_triangles + 1));
	data[ n_triangles].n1 = n1;
	data[ n_triangles].n2 = n2;
	data[ n_triangles].n3 = n3;
	n_triangles ++;
    }
    void get( long ind, long & n1, long & n2, long & n3) {
	assert( ind >= 0 && ind < n_triangles);
	n1 = data[ ind].n1;
	n2 = data[ ind].n2;
	n3 = data[ ind].n3;
    }
    long size( void) {
	return n_triangles;
    }
};
TriangleList * top_tlist = NULL, * side_tlist = NULL;

enum Menu {
    MENU_RAYTRACE,
    MENU_RELOAD,
    MENU_ORIG_DEFORMED,
    MENU_SAVE_VIEW,
    MENU_LOAD_VIEW,
    MENU_SAVE_EPS
};

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

    void save( FILE * fp) {
	/*
	  fprintf( fp, "%f %f %f # eye\n",
	  eye_x, eye_y, eye_z);
	*/
	fprintf( fp, "%f %f %f # center\n",
	    center_x, center_y, center_z);
	/*
	  fprintf( fp, "%f %f %f # up\n",
	  up_x, up_y, up_z);
	*/
	fprintf( fp, "%f # alpha\n", alpha);
	fprintf( fp, "%f # beta\n", beta);
	fprintf( fp, "%f # dist\n", dist);
    }

    int load( FILE * fp) {
	center_x = Tokenizer::read_double( fp, "center_x");
	center_y = Tokenizer::read_double( fp, "center_y");
	center_z = Tokenizer::read_double( fp, "center_z");
	double a = Tokenizer::read_double( fp, "alpha");
	double b = Tokenizer::read_double( fp, "beta");
	double d = Tokenizer::read_double( fp, "dist");
	set( a, b, d);
	return 0;
    }

} view;

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
	glutBitmapCharacter( GLUT_BITMAP_TIMES_ROMAN_10, * c);
}

static void setup_material (void)
{
    // setup the wall material
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
}

static void setup_material2 (void)
{
    // setup the wall material
    float amb[] = {0.3, 0.3, 0.3, 1.0 };
    float dif[] = {0.7, 0.7, 0.7, 1.0 };
    float spec[] = {0.02, 0.02, 0.02, 1.0 };
    float emis[] = {0.0, 0.0, 0.0, 1.0 };
    float shin[] = {0.0};
    
    glMaterialfv(GL_FRONT, GL_AMBIENT, amb);
    glMaterialfv(GL_FRONT, GL_DIFFUSE, dif);
    glMaterialfv(GL_FRONT, GL_SPECULAR, spec);
    glMaterialfv(GL_FRONT, GL_EMISSION, emis);
    glMaterialfv(GL_FRONT, GL_SHININESS, shin);
}

static void model_draw_polygons( void)
{
    setup_material ();

    // enable automatic normalization of normals
    glEnable( GL_NORMALIZE);

    // draw the tops
    glBegin( GL_TRIANGLES);
    for (size_t i = 0 ; i < top_triangles.size () ; i ++)
    {
	Triangle & t = top_triangles [i];
	glNormal (t.n1);
	glVertex (nodes[t.p1].pos);
	glNormal (t.n2);
	glVertex (nodes[t.p2].pos);
	glNormal (t.n3);
	glVertex (nodes[t.p3].pos);
    }
    glEnd();

    // draw the walls
    glBegin( GL_TRIANGLES);
    for (size_t i = 0 ; i < side_triangles.size () ; i ++)
    {
	Triangle & t = side_triangles [i];
	glNormal (t.n1);
	glVertex (nodes[t.p1].pos);
	glNormal (t.n2);
	glVertex (nodes[t.p2].pos);
	glNormal (t.n3);
	glVertex (nodes[t.p3].pos);
    }
    glEnd();
}

static void model_draw_flat_polygons( void)
{
    setup_material2 ();

    // enable automatic normalization of normals
    glEnable( GL_NORMALIZE);
    glDisable (GL_NORMALIZE);

    // draw the tops & bottoms
    glBegin (GL_TRIANGLES);
    for (size_t i = 0 ; i < elements.size () ; i ++) {
	Element & e = elements [i];
	glTriangleNormal (Triangle (e.p[0], e.p[1], e.p[2]));
	glTriangleNormal (Triangle (e.p[3], e.p[5], e.p[4]));
    }
    glEnd ();
    // draw the walls (but only non-shared ones)
    glBegin( GL_QUADS);
    for (size_t i = 0 ; i < elements.size () ; i ++) {
	Element & e = elements [i];
	if (! is_top_edge_shared (e.p[1],e.p[0]))
	    glQuadNormal (Quad (e.p[3], e.p[4], e.p[1], e.p[0]));
	if (! is_top_edge_shared (e.p[2],e.p[1]))
	    glQuadNormal (Quad (e.p[4], e.p[5], e.p[2], e.p[1]));
	if (! is_top_edge_shared (e.p[0],e.p[2]))
	    glQuadNormal (Quad (e.p[5], e.p[3], e.p[0], e.p[2]));
    }
    glEnd();
}

static void disp_func( void)
// ----------------------------------------------------------------------
// called whenever the contents of the OGL window need to be redrawn.
// ----------------------------------------------------------------------
{
    // initialize modelview matrix
    glMatrixMode( GL_MODELVIEW );
    glLoadIdentity();

    // setup first light
    {	float light_ambient[] = {0.6, 0.6, 0.6, 0.0 };
	float light_diffuse[] = {0.8, 0.8, 0.8, 0.0 };
	float light_specular[] = {0.8, 0.8, 0.8, 1.0 };
	float light_position[] = {1.0, 0.0, 1.0, 0.0};
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glEnable( GL_LIGHT0 );
    }

    // adjust modelview matrix
    gluLookAt( view.eye_x, view.eye_y, view.eye_z,
	view.center_x, view.center_y, view.center_z,
	view.up_x, view.up_y, view.up_z);

    // setup antialiasing
    if( antialiasing) {
	glEnable( GL_BLEND);
	glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable( GL_LINE_SMOOTH);
	glLineWidth( line_width);
    } else {
	glDisable( GL_BLEND);
	glLineWidth( 1);
    }

    // decide which draw to use
    GLfloat offset_factor = 1;
    GLfloat offset_units = 20;
    if (draw_type == 0) {
	// wireframe with hidden line removal

	// draw the polygons into z buffer only
	glDisable (GL_POLYGON_OFFSET_LINE);
	glClearColor( 0.8, 0.8, 0.8, 1.0);
	glClearDepth (1.0);
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glDepthFunc (GL_LESS);
	glEnable (GL_DEPTH_TEST);
	glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
	glEnable (GL_LIGHTING);
	glDrawBuffer (GL_NONE);
	model_draw_polygons ();

	// draw the visible lines
	glEnable (GL_POLYGON_OFFSET_LINE);
	glPolygonOffset (-offset_factor, - offset_units);
	glDisable (GL_LIGHTING);
	glColor3f (0, 0, 0);
	glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
	glDrawBuffer (GL_BACK);
	model_draw_polygons ();
    }
    else if (draw_type == 1) {
	// just polygons

	// draw the polygons into z buffer only
	glDisable (GL_POLYGON_OFFSET_LINE);
	glClearColor( 0.8, 0.8, 0.8, 1.0);
	glClearDepth (1.0);
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glDepthFunc (GL_LESS);
	glEnable (GL_DEPTH_TEST);
	glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
	glEnable (GL_LIGHTING);
	model_draw_polygons ();
    }
    else if (draw_type == 2) {
	// polygons + wireframe

	// draw the polygons into z buffer only
	glDisable (GL_POLYGON_OFFSET_LINE);
	glClearColor( 0.8, 0.8, 0.8, 1.0);
	glClearDepth (1.0);
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glDepthFunc (GL_LESS);
	glEnable (GL_DEPTH_TEST);
	glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
	glEnable (GL_LIGHTING);
	model_draw_polygons ();

	// draw the visible lines
	glEnable (GL_POLYGON_OFFSET_LINE);
	glPolygonOffset (-offset_factor, - offset_units);
	glDisable (GL_LIGHTING);
	glColor3f (0, 0, 0);
	glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
	glDrawBuffer (GL_BACK);
	model_draw_polygons ();
    }
    else if (draw_type == 3) {
	// wireframe model + faint hidden lines

	// draw all lines (not using the Z buffer
	glClearColor( 0.8, 0.8, 0.8, 1.0);
	glClearDepth (1.0);
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glDisable (GL_LIGHTING);
	glColor3f (0.5, 0.5, 0.5);
//	glColor3f (1, 0,0);
	glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
	glDrawBuffer (GL_BACK);
	glEnable (GL_DEPTH_TEST);
	model_draw_polygons ();

	// draw the polygons into z buffer only
	glDisable (GL_POLYGON_OFFSET_LINE);
	glClearDepth (1.0);
	glClear (GL_DEPTH_BUFFER_BIT);
	glDepthFunc (GL_LESS);
	glEnable (GL_DEPTH_TEST);
	glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
	glEnable (GL_LIGHTING);
	glDrawBuffer (GL_NONE);
	model_draw_polygons ();

	// draw the visible lines
	glEnable (GL_POLYGON_OFFSET_LINE);
	glPolygonOffset (-offset_factor, - offset_units);
	glDisable (GL_LIGHTING);
	glColor3f (0, 0, 0);
	glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
	glDrawBuffer (GL_BACK);
	glDepthMask (GL_FALSE);
	model_draw_polygons ();
	glDepthMask (GL_TRUE);
    }
    else if (draw_type == 4) {
	// polygons + wireframe (non-interpolated normals)

	// draw the polygons into z buffer only
	glDisable (GL_POLYGON_OFFSET_LINE);
	glClearColor( 0.8, 0.8, 0.8, 1.0);
	glClearDepth (1.0);
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glDepthFunc (GL_LESS);
	glEnable (GL_DEPTH_TEST);
	glPolygonMode (GL_FRONT_AND_BACK, GL_FILL);
	glEnable (GL_LIGHTING);
	model_draw_flat_polygons ();

	// draw the visible lines (without updates of the z buffer)
	glEnable (GL_POLYGON_OFFSET_LINE);
	glPolygonOffset (-offset_factor, - offset_units);
	glDisable (GL_LIGHTING);
	glDepthMask (GL_FALSE);
	glColor3f (0, 0, 0);
	glPolygonMode (GL_FRONT_AND_BACK, GL_LINE);
	glDrawBuffer (GL_BACK);
	model_draw_flat_polygons ();
	glDepthMask (GL_TRUE);
    }

    // swap buffers
    glutSwapBuffers();
}

static void reshape_func( int width, int height)
// ----------------------------------------------------------------------
// called whenever the OGL window is reshaped
// ----------------------------------------------------------------------
{
    glViewport( 0, 0, width, height );
    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    gluPerspective( 60.0, width/double(height), 0.1, 10.0);
    window_width = width;
    window_height = height;
}

int mouse_x, mouse_y;
enum MouseState { ROTATE, ZOOM } mouse_state;

void mouse_func( int button, int state, int x, int y)
// ----------------------------------------------------------------------
// called when mouse buttons are pressed
// ----------------------------------------------------------------------
{
    if( button == GLUT_LEFT_BUTTON) mouse_state = ROTATE;
    else if( button == GLUT_MIDDLE_BUTTON) mouse_state = ZOOM;

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
	double alpha = - dx*0.3 + view.alpha;
	double beta = dy*0.3 + view.beta;
	double dist = view.dist;
	view.set( alpha, beta, dist);
	glutPostRedisplay();
    } else if( mouse_state == ZOOM) {
	double alpha = view.alpha;
	double beta = view.beta;
	double dist = view.dist * pow(pow( 2.0, -dy/180.0), 2.0);
	view.set( alpha, beta, dist);
	glutPostRedisplay();
    }
}

void model_raytrace( void)
{
    std::string fname = "/tmp/model.pov";

    // calculate position, direction, up and right vectors
    Vector3d pos( view.eye_x, view.eye_y, view.eye_z);
    Vector3d dir( -view.eye_x, -view.eye_y, -view.eye_z);
    Vector3d up( 0, 0, 1);
    Vector3d right;
    right.assign( dir);
    right.cross_product( up);
    up.assign( dir);
    up.cross_product( right);
    up.normalize();
    right.normalize();
    dir.normalize();
    fprintf( stderr, "pos="); pos.print( stderr); fprintf( stderr, "\n");
    fprintf( stderr, "dir="); dir.print( stderr); fprintf( stderr, "\n");
    fprintf( stderr, "up="); up.print( stderr); fprintf( stderr, "\n");
    fprintf( stderr, "right="); right.print(stderr);fprintf( stderr, "\n");
	
    // save the model into a file
    FILE * fp = fopen( fname.c_str (), "w");
    assert( fp != NULL);
	
    fprintf( fp, "#include \"colors.inc\"\n");
    fprintf( fp, "#include \"woods.inc\"\n");
    fprintf( fp, "background { color White }\n");
    fprintf( fp, "camera {\n");
    fprintf( fp, "\tlocation <%f,%f,%f>\n",
	pos.x, pos.y, pos.z);
    Vector3d cdir = dir; cdir.normalize();
    fprintf( fp, "\tdirection <%.10f,%.10f,%.10f>\n",
	cdir.x, cdir.y, cdir.z);
    Vector3d cup = up; cup.normalize();
    fprintf( fp, "\tup <%.10f,%.10f,%.10f>\n",
	- cup.x, - cup.y, - cup.z);
    Vector3d cright = right; cright.normalize();
    fprintf( fp, "\tright <%.10f,%.10f,%.10f>\n",
	cright.x, cright.y, cright.z);
    fprintf( fp, "\tangle 60.0\n");
    fprintf( fp, "}\n");
    fprintf( fp, "\n");

    fprintf( fp, "light_source { <1, 2, 1> color White }\n");
    fprintf( fp, "\n");

    fprintf( fp, "#declare Tops = union {\n");

    for (size_t i = 0 ; i < top_triangles.size () ; i ++)
    {
	Triangle & t = top_triangles [i];
	fprintf( fp,
	    "\tsmooth_triangle {\n"
	    "\t\t<%f,%f,%f>, <%f,%f,%f>,\n"
	    "\t\t<%f,%f,%f>,<%f,%f,%f>,\n"
	    "\t\t<%f,%f,%f>,<%f,%f,%f> }\n",
	    nodes [t.p1].pos.x, nodes [t.p1].pos.y, nodes [t.p1].pos.z,
	    t.n1.x, t.n1.y, t.n1.z,
	    nodes [t.p2].pos.x, nodes [t.p2].pos.y, nodes [t.p2].pos.z,
	    t.n2.x, t.n2.y, t.n2.z,
	    nodes [t.p3].pos.x, nodes [t.p3].pos.y, nodes [t.p3].pos.z,
	    t.n3.x, t.n3.y, t.n3.z);
    }

    fprintf( fp, "}\n\n");
    fprintf( fp, "#declare Sides = union {\n");

    for (size_t i = 0 ; i < side_triangles.size () ; i ++)
    {
	Triangle & t = side_triangles [i];
	fprintf( fp,
	    "\tsmooth_triangle {\n"
	    "\t\t<%f,%f,%f>, <%f,%f,%f>,\n"
	    "\t\t<%f,%f,%f>,<%f,%f,%f>,\n"
	    "\t\t<%f,%f,%f>,<%f,%f,%f> }\n",
	    nodes [t.p1].pos.x, nodes [t.p1].pos.y, nodes [t.p1].pos.z,
	    t.n1.x, t.n1.y, t.n1.z,
	    nodes [t.p2].pos.x, nodes [t.p2].pos.y, nodes [t.p2].pos.z,
	    t.n2.x, t.n2.y, t.n2.z,
	    nodes [t.p3].pos.x, nodes [t.p3].pos.y, nodes [t.p3].pos.z,
	    t.n3.x, t.n3.y, t.n3.z);
    }
    fprintf( fp, "}\n\n");
    fprintf ( fp, "\n");
    fprintf ( fp, "#declare topsTexture = texture { T_Wood35  scale 0.5 rotate <0,1,0> }\n");
    fprintf ( fp, "#declare sidesTexture = texture { T_Wood35 scale 0.5 rotate <0,1,0> }\n");
    fprintf( fp,
	"object {\n"
	"\tTops \n"
	"\ttexture {\n"
	"\t\ttopsTexture\n"
	// 		 "\t\tfinish { ambient 0.3 }\n"
	// 		 "\t\tnormal { bumps 0.5 scale 0.001}\n"
	"\t}\n"
	"}\n");
    fprintf( fp,
	"object {\n"
	"\tSides \n"
	"\ttexture {\n"
	"\t\tsidesTexture\n"
	// 		 "\t\tfinish { ambient 0.3 }\n"
	// 		 "\t\tnormal { bumps 0.5 scale 0.001}\n"
	"\t}\n"
	"}\n");

    fclose( fp);

    // call raytracer on the result
    char * fname2 = "/tmp/model.png";
    char buff[ 4096];
    sprintf( buff,
	"(povray +P +A0.1 +D +SP16 +EP1 +I%s +O%s +W%ld +H%ld ; "
	"ee %s) &",
	fname.c_str (), fname2, window_width, window_height, fname2);
    system( buff);
}

int load_model( const char * fname);

void menu_func( int menu)
// ----------------------------------------------------------------------
// called when some of the menus are selected
// ----------------------------------------------------------------------
{
    switch( menu) {
	case MENU_RAYTRACE:
	    model_raytrace();
	    break;
	case MENU_RELOAD:
	    load_model( fname);
	    glutPostRedisplay();
	    break;
	case MENU_ORIG_DEFORMED:
	    use_deformed_coordinates = ! use_deformed_coordinates;
	    glutPostRedisplay();
	    break;
	case MENU_SAVE_VIEW:
	    {
		FILE * fp = fopen( "/tmp/viewer.dat", "w");
		if( fp == NULL) {
		    fprintf( stderr, "Could not save the view.\n");
		    break;
		}
		view.save( fp);
		fclose( fp);
		break;
	    }
	case MENU_LOAD_VIEW:
	    {
		FILE * fp = fopen( "/tmp/viewer.dat", "r");
		if( fp == NULL || view.load( fp)) {
		    fprintf( stderr, "Could not load the view.\n");
		    break;
		}
		fclose( fp);
		glutPostRedisplay();
		break;
	    }
	case MENU_SAVE_EPS:
	    {
		char * fname = "out.pdf";
		FILE * fp = fopen( fname, "wb");
		if (fp == NULL) {
		    fprintf( stderr, "Could not open output EPS file "
			"'%s'.\n", fname);
		    break;
		}
		GLint buffsize = 0, state = GL2PS_OVERFLOW;
		GLint viewport[4];
		
		glGetIntegerv(GL_VIEWPORT, viewport);
		gl2psLineWidth(0.1);		
		while( state == GL2PS_OVERFLOW ){ 
		    buffsize += 1024*1024;
		    gl2psBeginPage ( "MyTitle", "MySoftware", viewport,
			//GL2PS_EPS,
			GL2PS_PDF,
			GL2PS_BSP_SORT,
			GL2PS_SIMPLE_LINE_OFFSET | GL2PS_SILENT |
			GL2PS_OCCLUSION_CULL | GL2PS_BEST_ROOT
			| GL2PS_NO_PS3_SHADING
			,
			GL_RGBA, 0, NULL, 0, 0, 0, buffsize,
			fp, NULL );
		    disp_func ();
		    state = gl2psEndPage();
		}

		fclose( fp);
		break;
	    }
    }
}

static void keyboard_func( unsigned char key, int x, int y)
{
    key = tolower( key);
    if( key == 'a') {
	// toggle display of axes
	draw_axes = ! draw_axes;
	glutPostRedisplay();
	return ;
    } else if (key == 't') {
	draw_type = (draw_type + 1) % n_draw_types;
	glutPostRedisplay ();
	return;
    } else if( key == 'w') {
	draw_wireframe = (draw_wireframe +1) % 3;
	glutPostRedisplay();
	return;
    } else if( key == '-') {
	line_width *= 0.9;
	glutPostRedisplay();
	return;
    } else if( key == '+' || key == '=') {
	line_width *= 1.5;
	glutPostRedisplay();
	return;
    } else if( key == 'l') {
	draw_solid = (draw_solid + 1) % 3;
	glutPostRedisplay();
	return;
    } else if( key == 'h') {
	fprintf( stderr,
	    "=================================================\n"
	    "                H e l p\n"
	    ".................................................\n"
	    "     <h> displays this help\n"
	    "     <a> toggles display of axes\n"
	    "     <w> toggles wireframe\n"
	    "     <l> toggles wall display\n"
	    "=================================================\n"
		 );
    } else {
	fprintf( stderr, "Unknown key pressed: %d\n", int(key));
	return;
    }
}

long edge_count( long n1, long n2)
{
    // counts the occurence of edge n1-n2 in the list of elements
    // (assuming n1 & n2 are on the top layer)
    long res = 0;
    for( size_t i = 0 ; i < elements.size () ; i ++) {
	if( elements[i].p[0] == n1 && elements[i].p[1] == n2) res ++;
	if( elements[i].p[1] == n1 && elements[i].p[2] == n2) res ++;
	if( elements[i].p[2] == n1 && elements[i].p[0] == n2) res ++;
	if( elements[i].p[0] == n2 && elements[i].p[1] == n1) res ++;
	if( elements[i].p[1] == n2 && elements[i].p[2] == n1) res ++;
	if( elements[i].p[2] == n2 && elements[i].p[0] == n1) res ++;
    }
    return res;
}

int load_model( const char * fname)
{
    // free up stuff from before
    if( top_tlist != NULL) delete top_tlist;
    if( side_tlist != NULL) delete side_tlist;
    if( top_norms != NULL) delete [] top_norms;
    if( side_norms != NULL) delete [] side_norms;

    // read the file
    FILE * fp = fopen( fname, "r");
    if( fp == NULL) die( "Cannot open file '%s' for reading.", fname);
	
    // skip all variables & material properties
    while( 1) {
	char * token = Tokenizer::read_string( fp, "END");
	if( token == NULL)
	    die( "Unexpected end of file when looking of END.");
	if( strcasecmp( token, "END") == 0) break;
    }
    (void) Tokenizer::read_string( fp, "sec_name");

    // read the number of nodes
    long n_nodes = Tokenizer::read_long( fp, "number of nodes");
    fprintf( stderr, "Number of nodes: %ld\n", n_nodes);

    // allocate room for the nodes and read them all in
    nodes.clear ();
    nodes.resize (n_nodes);
    for( long i = 0 ; i < n_nodes ; i ++) {
	nodes[i].pos.x = Tokenizer::read_double( fp, "ox");
	nodes[i].pos.y = Tokenizer::read_double( fp, "oy");
	nodes[i].pos.z = Tokenizer::read_double( fp, "oz");
	nodes[i].pos.x += Tokenizer::read_double( fp, "dx");
	nodes[i].pos.y += Tokenizer::read_double( fp, "dy");
	nodes[i].pos.z += Tokenizer::read_double( fp, "dz");
	// skip growth normals
	(void) Tokenizer::read_double( fp, "growth x");
	(void) Tokenizer::read_double( fp, "growth y");
	(void) Tokenizer::read_double( fp, "growth z");
	// skip node type
	(void) Tokenizer::read_long( fp, "type");
    }
	
    // read in the number of elements
    long n_elements = Tokenizer::read_long( fp, "n_elements");
    fprintf( stderr, "Number of elements: %ld\n", n_elements);
    elements.clear ();
    elements.resize(n_elements);
    for( long i = 0 ; i < n_elements ; i ++) {
	char * token = Tokenizer::read_string( fp, "FIVEWALL");
	if( token == NULL || strcasecmp( token, "fivewall") != 0)
	    die( "Expected 'fivewall' instead of '%s'.",
		token);
	(void) Tokenizer::read_long( fp, "id");
	elements[ i].p[0] = Tokenizer::read_long( fp, "p0");
	elements[ i].p[1] = Tokenizer::read_long( fp, "p1");
	elements[ i].p[2] = Tokenizer::read_long( fp, "p2");
	elements[ i].p[3] = Tokenizer::read_long( fp, "p3");
	elements[ i].p[4] = Tokenizer::read_long( fp, "p4");
	elements[ i].p[5] = Tokenizer::read_long( fp, "p5");

	// add the element to every node
	for (long j = 0 ; j < 6 ; j ++)
	    nodes [elements [i].p [j]].els.push_back (i);
    }
    fclose( fp);

    // create a list of all triangles
    top_triangles = get_top_triangles ();
    side_triangles = get_side_triangles ();

    // calculate normals for the top list
    calculate_normals (top_triangles, 30);
    calculate_normals (side_triangles, 30);

    return 0;
}

int main( int argc, char ** argv)
{
    // parse the arguments
    if( argc != 2) usage();
    fname = argv[1];

    load_model( fname);

    // initialize GUI
    glutInitWindowSize( 500, 500);
    glutInit( &argc, argv );
    glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    char title[ 4096];
    sprintf( title, "%s %s", argv[0], fname);
    int main_window = glutCreateWindow (title);
    assert( main_window);
    glutSetWindow( main_window);

    glutDisplayFunc( disp_func);
    glutReshapeFunc( reshape_func);
    glutMouseFunc( mouse_func);
    glutMotionFunc( mouse_motion_func);
    glutKeyboardFunc( keyboard_func);

    // attempt to load the view
    view.set( 45, 40, 2);
    {
	FILE * fp = fopen( "/tmp/viewer.dat", "r");
	if( fp == NULL || view.load( fp)) {
	    fprintf( stderr, "Could not load the view.\n");
	}
	if( fp != NULL) fclose( fp);
    }
        
    // create a popup menu
    int m = glutCreateMenu( menu_func);
    glutSetMenu( m);
    glutAddMenuEntry( "Raytrace...", MENU_RAYTRACE);
    glutAddMenuEntry( "Reload", MENU_RELOAD);
    glutAddMenuEntry( "Save view", MENU_SAVE_VIEW);
    glutAddMenuEntry( "Load view", MENU_LOAD_VIEW);
    glutAddMenuEntry( "Save EPS", MENU_SAVE_EPS);
    glutAttachMenu( GLUT_RIGHT_BUTTON);
	
    glutMainLoop();

    return 0;
}

