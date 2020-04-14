#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include <string.h>
#include <ctype.h>
#include <malloc.h>
#include <sys/types.h>
#include <unistd.h>
#include <GL/gl.h>
#include <GL/glut.h>
#include "die.hh"
#include "Tokenizer.hh"
#include "Vector3d.hh"

void usage( void)
{
	die( "Usage: viewer-cyl radius fname\n");
}

class Node {
public:
	double x, y, z;		// position of the node
} * nodes;
long n_nodes;

Vector3d * top_norms, * side_norms;

class Element {
public:
	long p[6];
} * elements;
long n_elements;

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
TriangleList * top_tlist, * side_tlist;

enum Menu {
	MENU_RAYTRACE,
	MENU_ORIG_DEFORMED
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
} view;

int antialiasing = 0;
int draw_axes = 0;
int draw_wireframe = 0;
int draw_solid = 2;
int use_deformed_coordinates = 1;

long window_width;
long window_height;

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
	glNormal3d( nx, ny, nz);
}

void draw_triangle( long pn1, long pn2, long pn3)
{
	Node & n1 = nodes[pn1];
	Node & n2 = nodes[pn2];
	Node & n3 = nodes[pn3];

	draw_triangle_normal( n1.x, n1.y, n1.z,
			      n2.x, n2.y, n2.z,
			      n3.x, n3.y, n3.z);

	glBegin( GL_TRIANGLES);
	glVertex3d( n1.x, n1.y, n1.z);
	glVertex3d( n2.x, n2.y, n2.z);
	glVertex3d( n3.x, n3.y, n3.z);
	glEnd();
}

void model_draw( void)
{
	// draw axes
	if( draw_axes) {
		glDisable( GL_LIGHTING);
		glColor3f( 1, 0, 0);
		glBegin( GL_LINES);
		glVertex3d( 0, 0, 0); glVertex3d( 1, 0, 0);
		glVertex3d( 0, 0, 0); glVertex3d( 0, 1, 0);
		glVertex3d( 0, 0, 0); glVertex3d( 0, 0, 1);
		glEnd();
		glEnable( GL_LIGHTING);
	}

	// draw wireframe top if requested
	if( draw_wireframe > 0) {
//		glDisable( GL_DEPTH_TEST);
		glDisable( GL_LIGHTING);
		glColor3f( 0, 0, 0);
		glLineWidth( 2);
		for( long i = 0 ; i < top_tlist-> size() ; i ++) {
			glBegin( GL_LINE_LOOP);
			long n1, n2, n3;
			top_tlist-> get( i, n1, n2, n3);
			glVertex3d( nodes[n1].x, nodes[n1].y, nodes[n1].z);
			glVertex3d( nodes[n2].x, nodes[n2].y, nodes[n2].z);
			glVertex3d( nodes[n3].x, nodes[n3].y, nodes[n3].z);
			glEnd();
		}
//		glEnable( GL_DEPTH_TEST);
		glEnable( GL_LIGHTING);
	}

	// set the material properties
	{
//		float amb[] = {0.3, 0.13, 0.0, 1.0 };
		float amb[] = {0.15, 0.07, 0.0, 1.0 };
		float dif[] = {0.7, 0.3, 0.0, 1.0 };
//		float dif[] = {0.35, 0.15, 0.0, 1.0 };
		float spec[] = {0.0, 0.0, 0.0, 1.0 };
		float emis[] = {0.0, 0.0, 0.0, 1.0 };
		float shin[] = {0.0};

		glMaterialfv(GL_FRONT, GL_AMBIENT, amb);
		glMaterialfv(GL_FRONT, GL_DIFFUSE, dif);
		glMaterialfv(GL_FRONT, GL_SPECULAR, spec);
		glMaterialfv(GL_FRONT, GL_EMISSION, emis);
		glMaterialfv(GL_FRONT, GL_SHININESS, shin);
	}
	
	// draw solid tops if requested
	if( draw_solid > 0) {
		glBegin( GL_TRIANGLES);
		for( long i = 0 ; i < top_tlist-> size() ; i ++) {
			long n1, n2, n3;
			top_tlist-> get( i, n1, n2, n3);
			glNormal3d( top_norms[n1].x,
				    top_norms[n1].y,
				    top_norms[n1].z);
			glVertex3d( nodes[n1].x, nodes[n1].y, nodes[n1].z);
			glNormal3d( top_norms[n2].x,
				    top_norms[n2].y,
				    top_norms[n2].z);
			glVertex3d( nodes[n2].x, nodes[n2].y, nodes[n2].z);
			glNormal3d( top_norms[n3].x,
				    top_norms[n3].y,
				    top_norms[n3].z);
			glVertex3d( nodes[n3].x, nodes[n3].y, nodes[n3].z);
		}
		glEnd();
	}

	{
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

	// draw solid walls if requested
	if( draw_solid > 1) {
		glBegin( GL_TRIANGLES);
		for( long i = 0 ; i < side_tlist-> size() ; i ++) {
			long n1, n2, n3;
			side_tlist-> get( i, n1, n2, n3);
			glNormal3d( side_norms[n1].x,
				    side_norms[n1].y,
				    side_norms[n1].z);
			glVertex3d( nodes[n1].x, nodes[n1].y, nodes[n1].z);
			glNormal3d( side_norms[n2].x,
				    side_norms[n2].y,
				    side_norms[n2].z);
			glVertex3d( nodes[n2].x, nodes[n2].y, nodes[n2].z);
			glNormal3d( side_norms[n3].x,
				    side_norms[n3].y,
				    side_norms[n3].z);
			glVertex3d( nodes[n3].x, nodes[n3].y, nodes[n3].z);
		}
		glEnd();
	}
}

static void disp_func( void)
// ----------------------------------------------------------------------
// called whenever the contents of the OGL window need to be redrawn.
// ----------------------------------------------------------------------
{
	glMatrixMode( GL_MODELVIEW);
	glLoadIdentity();

/*
	// setup a light
	float light_ambient[] = {0.3, 0.3, 0.3, 0.0 };
	float light_diffuse[] = {0.8, 0.8, 0.8, 0.0 };
	float light_specular[] = {1.0, 1.0, 1.0, 0.0 };
	float light_position[] = {0.0, 2.0, 1.0, 0.0 };
	
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glEnable( GL_LIGHT0 );
*/
	// setup first light
	{
		float light_ambient[] = {0.6, 0.6, 0.6, 0.0 };
		float light_diffuse[] = {0.8, 0.8, 0.8, 0.0 };
		float light_specular[] = {0.8, 0.8, 0.8, 0.0 };
		float light_position[] = {3.0, 3.0, 3.0, 0.0 };
		
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
		
		glEnable( GL_LIGHT1 );
	}


        // set the MODELVIEW matrix
        glMatrixMode( GL_MODELVIEW );
        glLoadIdentity();
        gluLookAt( view.eye_x, view.eye_y, view.eye_z,
                   view.center_x, view.center_y, view.center_z,
                   view.up_x, view.up_y, view.up_z);

        // set the background color to light grey
        glClearColor( 0.8, 0.8, 0.8, 1.0);
//        glClearColor( 0.0, 0.0, 0.0, 1.0);

        
        // clear the background
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

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

        // enable automatic normalization of normals
        glEnable( GL_NORMALIZE);

	// enable Z buffer
	glEnable( GL_DEPTH_TEST);

	// enable lighting
	glEnable( GL_LIGHTING );

	model_draw();

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
        gluPerspective( 60.0, width/double(height), 0.1, 100.0);
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
	char * fname = "/tmp/model.pov";

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
	fprintf( stderr, "right="); right.print( stderr); fprintf( stderr, "\n");
	

	// save the model into a file
	FILE * fp = fopen( fname, "w");
	assert( fp != NULL);
	
	fprintf( fp, "#include \"colors.inc\"\n");
	fprintf( fp, "background { color Black }\n");
	fprintf( fp, "camera {\n");
/*	fprintf( fp, "\tlocation <%f,%f,%f>\n",
		 view.eye_x, view.eye_y, view.eye_z);
	fprintf( fp, "\tsky <0,0,1>\n");
	fprintf( fp, "\tlook_at <%f,%f,%f>\n",
		 view.center_x, view.center_y, view.center_z);
*/
	fprintf( fp, "\tlocation <%f,%f,%f>\n",
		 pos.x, pos.y, pos.z);
	fprintf( fp, "\tdirection <%.10f,%.10f,%.10f>\n",
		 dir.x, dir.y, dir.z);
	fprintf( fp, "\tup <%.10f,%.10f,%.10f>\n",
		 - up.x, - up.y, - up.z);
	fprintf( fp, "\tright <%.10f,%.10f,%.10f>\n",
		 right.x, right.y, right.z);
	fprintf( fp, "\tangle 60.0\n");
	fprintf( fp, "}\n");
	fprintf( fp, "\n");

	fprintf( fp, "light_source { <1, 2, 1> color White }\n");
	fprintf( fp, "\n");

	fprintf( fp, "#declare Tops = union {\n");

	for( long i = 0 ; i < top_tlist-> size() ; i ++) {
		long n1, n2, n3;
		top_tlist-> get( i, n1, n2, n3);
		fprintf( fp,
			 "\tsmooth_triangle {\n"
			 "\t\t<%f,%f,%f>, <%f,%f,%f>,\n"
			 "\t\t<%f,%f,%f>,<%f,%f,%f>,\n"
			 "\t\t<%f,%f,%f>,<%f,%f,%f> }\n",
			 nodes[ n1].x, nodes[ n1].y, nodes[ n1].z,
			 top_norms[ n1].x, top_norms[ n1].y, top_norms[ n1].z,
			 nodes[ n2].x, nodes[ n2].y, nodes[ n2].z,
			 top_norms[ n2].x, top_norms[ n2].y, top_norms[ n2].z,
			 nodes[ n3].x, nodes[ n3].y, nodes[ n3].z,
			 top_norms[ n3].x, top_norms[ n3].y, top_norms[ n3].z
			);
	}

	fprintf( fp, "}\n\n");
	fprintf( fp, "#declare Sides = union {\n");

	for( long i = 0 ; i < side_tlist-> size() ; i ++) {
		long n1, n2, n3;
		side_tlist-> get( i, n1, n2, n3);
		fprintf( fp,
			 "\tsmooth_triangle {\n"
			 "\t<%f,%f,%f>,<%f,%f,%f>,\n"
			 "\t<%f,%f,%f>,<%f,%f,%f>,\n"
			 "\t<%f,%f,%f>,<%f,%f,%f> }\n",
			 nodes[ n1].x, nodes[ n1].y, nodes[ n1].z,
			 side_norms[n1].x, side_norms[n1].y, side_norms[n1].z,
			 nodes[ n2].x, nodes[ n2].y, nodes[ n2].z,
			 side_norms[n2].x, side_norms[n2].y, side_norms[n2].z,
			 nodes[ n3].x, nodes[ n3].y, nodes[ n3].z,
			 side_norms[n3].x, side_norms[n3].y, side_norms[n3].z
			);
	}
	
	fprintf( fp, "}\n\n");
	fprintf( fp,
		 "object {\n"
		 "\tTops \n"
		 "\ttexture {\n"
		 "\t\tpigment { color rgb <0.4,0.4,0.4> }\n"
		 "\t\tfinish { ambient 0.3 }\n"
		 "\t\tnormal { bumps 0.5 scale 0.001}\n"
		 "\n"
		 "\t}\n"
		 "}\n");
	fprintf( fp,
		 "object {\n"
		 "\tSides \n"
		 "\ttexture {\n"
		 "\t\tpigment { color rgb <0.2,0.2,0.2> }\n"
		 "\t\tfinish { ambient 0.3 }\n"
		 "\t\tnormal { bumps 0.5 scale 0.001}\n"
		 "\t}\n"
		 "}\n");

	fclose( fp);

	// call raytracer on the result
	char * fname2 = "/tmp/model.png";
	char buff[ 4096];
	sprintf( buff,
		 "(povray +A0.1 +D +SP16 +EP1 +I%s +O%s +W%ld +H%ld ; "
		 "xv %s) &",
		 fname, fname2, window_width, window_height, fname2);
	system( buff);
}

void menu_func( int menu)
// ----------------------------------------------------------------------
// called when some of the menus are selected
// ----------------------------------------------------------------------
{
	switch( menu) {
	case MENU_RAYTRACE:
		model_raytrace();
		break;
	case MENU_ORIG_DEFORMED:
		use_deformed_coordinates = ! use_deformed_coordinates;
		glutPostRedisplay();
		break;
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
	} else if( key == 'w') {
		draw_wireframe = (draw_wireframe +1) % 3;
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
	for( long i = 0 ; i < n_elements ; i ++) {
		if( elements[i].p[0] == n1 && elements[i].p[1] == n2) res ++;
		if( elements[i].p[1] == n1 && elements[i].p[2] == n2) res ++;
		if( elements[i].p[2] == n1 && elements[i].p[0] == n2) res ++;
		if( elements[i].p[0] == n2 && elements[i].p[1] == n1) res ++;
		if( elements[i].p[1] == n2 && elements[i].p[2] == n1) res ++;
		if( elements[i].p[2] == n2 && elements[i].p[0] == n1) res ++;
	}
	return res;
}

int main( int argc, char ** argv)
{
	// parse the arguments
	char * fname = NULL;
	if( argc != 3) usage();
	fname = argv[2];
	double radius;
	if( 1 != sscanf( argv[1], "%lf", & radius)) usage();

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
	n_nodes = Tokenizer::read_long( fp, "number of nodes");
	fprintf( stderr, "Number of nodes: %ld\n", n_nodes);

	// allocate room for the nodes and read them all in
	nodes = new Node[ n_nodes];
	for( long i = 0 ; i < n_nodes ; i ++) {
		nodes[i].x = Tokenizer::read_double( fp, "ox");
		nodes[i].y = Tokenizer::read_double( fp, "oy");
		nodes[i].z = Tokenizer::read_double( fp, "oz");
		nodes[i].x += Tokenizer::read_double( fp, "dx");
		nodes[i].y += Tokenizer::read_double( fp, "dy");
		nodes[i].z += Tokenizer::read_double( fp, "dz");
		// skip growth normals
		(void) Tokenizer::read_double( fp, "growth x");
		(void) Tokenizer::read_double( fp, "growth y");
		(void) Tokenizer::read_double( fp, "growth z");
		// skip node type
		(void) Tokenizer::read_long( fp, "type");
	}

	// map the nodes onto a cylinder
	if( radius > 0) {
		for( long i = 0 ; i < n_nodes ; i ++) {
			Node & n = nodes[i];
			double alpha = n.x / ( M_PI * radius);
			double h = n.y;
			double r = radius + n.z;
			n.x = sin(alpha) * r;
			n.y = cos(alpha) * r;
			n.z = h;
		}
	}
	
	// read in the number of elements
	n_elements = Tokenizer::read_long( fp, "n_elements");
	fprintf( stderr, "Number of elements: %ld\n", n_elements);
	elements = new Element[ n_elements];
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
	}
	fclose( fp);
	
	fprintf( stderr, "Creating edge list...");
	EdgeList elist( n_nodes);
	for( long i = 0 ; i < n_elements ; i ++) {
		long n0 = elements[i].p[0];
		long n1 = elements[i].p[1];
		long n2 = elements[i].p[2];
		long n3 = elements[i].p[3];
		long n4 = elements[i].p[4];
		long n5 = elements[i].p[5];
		elist.add( n0, n1);
		elist.add( n1, n2);
		elist.add( n2, n0);
/*
		elist.add( n0, n3);
		elist.add( n1, n4);
		elist.add( n2, n5);
		elist.add( n3, n4);
		elist.add( n4, n5);
		elist.add( n5, n3);
*/
	}
	fprintf( stderr, "done.\n");

	fprintf( stderr, "Creating triangle list...");
	// create two lists of triangles - one for the top layer,
	// and one for all the other triangles
	top_tlist = new TriangleList();
	side_tlist = new TriangleList();
	for( long i = 0 ; i < n_elements ; i ++) {
		long p0 = elements[i].p[0];
		long p1 = elements[i].p[1];
		long p2 = elements[i].p[2];
		long p3 = elements[i].p[3];
		long p4 = elements[i].p[4];
		long p5 = elements[i].p[5];
		top_tlist-> add( p0, p1, p2);
		// if edge p0-p1 is not shared by two elements we add
		// to the list p0-p3-p1 and p3-p4-p1, etc
		if( elist.count( p0, p1) < 2) {
			side_tlist-> add( p0, p3, p1);
			side_tlist-> add( p3, p4, p1);
		}
		if( elist.count( p1, p2) < 2) {
			side_tlist-> add( p1, p4, p2);
			side_tlist-> add( p4, p5, p2);
		}
		if( elist.count( p2, p0) < 2) {
			side_tlist-> add( p2, p5, p0);
			side_tlist-> add( p5, p3, p0);
		}
	}
	fprintf( stderr, "\n");

	// create a list of normals for the top layer
	// --------------------------------------------------
	top_norms = new Vector3d[ n_nodes];
	long * counts = new long[ n_nodes];
	for( long i = 0 ; i < n_nodes ; i ++) counts[ i] = 0;
	for( long i = 0 ; i < top_tlist-> size() ; i ++) {
		long n1, n2, n3;
		top_tlist-> get( i, n1, n2, n3);
		Vector3d v1( nodes[n2].x - nodes[n1].x,
			     nodes[n2].y - nodes[n1].y,
			     nodes[n2].z - nodes[n1].z);
		Vector3d v2( nodes[n3].x - nodes[n1].x,
			     nodes[n3].y - nodes[n1].y,
			     nodes[n3].z - nodes[n1].z);
		Vector3d norm;
		norm.assign( v1);
		norm.cross_product( v2);
		norm.normalize();
		top_norms[n1].add( norm); counts[n1] ++;
		top_norms[n2].add( norm); counts[n2] ++;
		top_norms[n3].add( norm); counts[n3] ++;
	}
	for( long i = 0 ; i < n_nodes ; i ++) {
		if( counts[i] > 0) {
			top_norms[i].x /= counts[ i];
			top_norms[i].y /= counts[ i];
			top_norms[i].z /= counts[ i];
		}
	}

	// create a list of normals for the sides
	// --------------------------------------------------
	side_norms = new Vector3d[ n_nodes];
	for( long i = 0 ; i < n_nodes ; i ++) counts[ i] = 0;
	for( long i = 0 ; i < side_tlist-> size() ; i ++) {
		long n1, n2, n3;
		side_tlist-> get( i, n1, n2, n3);
		Vector3d v1( nodes[n2].x - nodes[n1].x,
			     nodes[n2].y - nodes[n1].y,
			     nodes[n2].z - nodes[n1].z);
		Vector3d v2( nodes[n3].x - nodes[n1].x,
			     nodes[n3].y - nodes[n1].y,
			     nodes[n3].z - nodes[n1].z);
		Vector3d norm;
		norm.assign( v1);
		norm.cross_product( v2);
		norm.normalize();
		side_norms[n1].add( norm); counts[n1] ++;
		side_norms[n2].add( norm); counts[n2] ++;
		side_norms[n3].add( norm); counts[n3] ++;
	}
	for( long i = 0 ; i < n_nodes ; i ++) {
		if( counts[i] > 0) {
			side_norms[i].x /= counts[ i];
			side_norms[i].y /= counts[ i];
			side_norms[i].z /= counts[ i];
		}
	}
	delete [] counts;

	// initialize GUI
	glutInitWindowSize( 500, 500);
        glutInit( &argc, argv );
        glutInitDisplayMode( GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
        char title[ 4096];
        sprintf( title, "%s %s", argv[0], fname);
        int main_window = glutCreateWindow( title);
        assert( main_window);
        glutSetWindow( main_window);

        glutDisplayFunc( disp_func);
        glutReshapeFunc( reshape_func);
        glutMouseFunc( mouse_func);
        glutMotionFunc( mouse_motion_func);
	glutKeyboardFunc( keyboard_func);

        view.set( 45, 40, 2);
        
        // create a popup menu
        int m = glutCreateMenu( menu_func);
        glutSetMenu( m);
	glutAddMenuEntry( "Raytrace...", MENU_RAYTRACE);
	glutAddMenuEntry( "Switch original/deformed", MENU_ORIG_DEFORMED);
        glutAttachMenu( GLUT_RIGHT_BUTTON);
	
        glutMainLoop();

	return 0;
}
