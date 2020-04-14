#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "Map.hh"
#include "die.hh"
#include "Tokenizer.hh"
#include "Distribution.hh"

Map::Map( FILE * fp)
{
	// read the type of the map
	char * type = Tokenizer::read_string( fp, "map_type");
	
	// figure out which map type this is, and then read in the
	// appropriate parameters
	if( strcasecmp( type, "texture2D") == 0) {
		map_type = Texture2D;

		// read the name of the texture image file
		fname = Tokenizer::read_string( fp, "img. fname");
		fname = strdup( fname);
		image = new Image();
		if( image-> load( fname))
			die( "Could not load image: %s", fname);
		
		// read in the min and max values
		min_val = Tokenizer::read_double( fp, "min_val");
		max_val = Tokenizer::read_double( fp, "max_val");

		// read in the min & max XY values
		min_x = Tokenizer::read_double( fp, "min_x");
		max_x = Tokenizer::read_double( fp, "max_x");
		min_y = Tokenizer::read_double( fp, "min_y");
		max_y = Tokenizer::read_double( fp, "max_y");
	} else if( strcasecmp( type, "constant") == 0) {
		map_type = Constant;
		min_val = Tokenizer::read_double( fp, "constant map value");
		max_val = min_val;
	} else if( strcasecmp( type, "normal") == 0) {
		map_type = Normal;
		min_x = Tokenizer::read_double( fp, "normal distr. mean");
		min_y = Tokenizer::read_double( fp, "normal distr. variance");
		min_val = min_x - 4 * min_y;
		max_val = min_x + 4 * min_y;
	} else if( strcasecmp( type, "textureCylinder") == 0) {
		map_type = TextureCylinder;

		// read the name of the texture image file
		fname = Tokenizer::read_string( fp, "img. fname");
		fname = strdup( fname);
		image = new Image();
		if( image-> load( fname))
			die( "Could not load image: %s", fname);
		
		// read in the min and max values
		min_val = Tokenizer::read_double( fp, "min_val");
		max_val = Tokenizer::read_double( fp, "max_val");

		// the height of the cylinder
		min_z = Tokenizer::read_double( fp, "min_z");
		max_z = Tokenizer::read_double( fp, "max_z");
	} else {
		die( "Unknown texture type '%s'", type);
	}
}

Map::Map( Tokenizer & tok)
{
	// read the type of the map
	std::string type = tok.read_string( "map_type");
	
	// figure out which map type this is, and then read in the
	// appropriate parameters
	if( type == "texture2D" == 0) {
		map_type = Texture2D;

		// read the name of the texture image file
		fname = strdup( tok.read_string( "img. fname").c_str());
		image = new Image();
		if( image-> load( fname))
			die( "Could not load image: %s", fname);
		
		// read in the min and max values
		min_val = tok.read_double( "min_val");
		max_val = tok.read_double( "max_val");

		// read in the min & max XY values
		min_x = tok.read_double( "min_x");
		max_x = tok.read_double( "max_x");
		min_y = tok.read_double( "min_y");
		max_y = tok.read_double( "max_y");
	} else if( type == "constant") {
		map_type = Constant;
		min_val = tok.read_double( "constant map value");
		max_val = min_val;
	} else if( type == "normal") {
		map_type = Normal;
		min_x = tok.read_double( "normal distr. mean");
		min_y = tok.read_double( "normal distr. variance");
		min_val = min_x - 4 * min_y;
		max_val = min_x + 4 * min_y;
	} else if( type == "textureCylinder") {
		map_type = TextureCylinder;

		// read the name of the texture image file
		fname = strdup( tok.read_string( "img. fname").c_str());
		image = new Image();
		if( image-> load( fname))
			die( "Could not load image: %s", fname);
		
		// read in the min and max values
		min_val = tok.read_double( "min_val");
		max_val = tok.read_double( "max_val");

		// the height of the cylinder
		min_z = tok.read_double( "min_z");
		max_z = tok.read_double( "max_z");
	} else {
		die( "Unknown texture type '%s'", type.c_str());
	}
}

Map::~Map( )
{
	if( fname != NULL) free( fname);
}

void Map::print( FILE * fp)
{
	if( map_type == Constant)
		fprintf( fp, "Constant %f", min_val);
	else if( map_type == Normal)
		fprintf( fp, "Normal %f %f", min_x, min_y);
	else if( map_type == Texture2D)
		fprintf( fp, "Texture2D \"%s\" %.15f %.15f %.15f %.15f"
			 " %.15f %.15f",
			 fname,
			 min_val, max_val, min_x, max_x, min_y, max_y);
	else if( map_type == TextureCylinder)
		fprintf( fp, "TextureCylinder \"%s\" %.15f %.15f %.15f %.15f",
			 fname,
			 min_val, max_val, min_z, max_z);
	else
		fprintf( fp, "Unknown map");
}
void Map::cyl_to_plane( double & x, double & y, double & z)
// assumption:
//    + [x,y,z] is a point on a surface of a cylinder whose base is in
//      the XY plane, and erects in the Z plane
// returns:
//    + [x,y,z] when the surface of the cylinder would be projected on
//      a rectangle of the width: 2*M_PI and of height: the same as that
//      of the cylinder
{
	double d = sqrt( x*x + y*y);
	double a = x / d;
	if( a > 1) a = 1; else if( a < -1) a = -1;
	double alpha = acos( a);
	if( y < 0) alpha = 2*M_PI - alpha;

	// now put the point in the plane
	x = alpha;
	y = z;
	z = 0;
}

double Map::get_value( double x1, double y1, double z1,
		       double x2, double y2, double z2,
		       double x3, double y3, double z3)
// ----------------------------------------------------------------------
// given the coordinates of a triangle, calculate the average pixel value
// in the triangle if it was to be drawn on the texture map
// ----------------------------------------------------------------------
{
	if( map_type == Constant) {

		return min_val;

	} else if( map_type == Normal) {
		
		return Distribution::get_normal_rnd( min_x, min_y);

	} else if( map_type == Texture2D) {
		
		// Iterative process:
		//  - generate random coordinates inside the triangle,
		//    transform them to the picture image, add to the
		//    sum
		//  - from the sum, min_val and max_val and number of
		//    iterations determine the return value

		long n_iterations = 100; // TBD - add interface to the
				// class to set this up, so that it is
				// not a constant

		// figure out the texture coordinates of the triangle
		double tx1 = (x1-min_x) / (max_x-min_x) * image-> width;
		double ty1 = (y1-min_y) / (max_y-min_y) * image-> height;
		double tx2 = (x2-min_x) / (max_x-min_x) * image-> width;
		double ty2 = (y2-min_y) / (max_y-min_y) * image-> height;
		double tx3 = (x3-min_x) / (max_x-min_x) * image-> width;
		double ty3 = (y3-min_y) / (max_y-min_y) * image-> height;

		// do the loop
		double sum = 0;
		for( long i = 0 ; i < n_iterations ; i ++) {
			// generate random isoparametric coordinates p,q=(0,1)
			double p = drand48();
			double q = drand48();
			// make sure p+q <= 1
			if( p+q > 1) { p = 1-p; q = 1-q; }
			// transform p,q into texture coordinates
			long x = long( (tx2-tx1)*p + (tx3-tx1)*q + tx1);
			long y = long( (ty2-ty1)*p + (ty3-ty1)*q + ty1);
			// wrap around to fit into the image
			while( x < 0) x += image-> width;
			while( x >= image-> width) x-= image-> width;
			while( y < 0) y += image-> height;
			while( y >= image-> height) y -= image-> height;
			sum += image-> get( y, x);
		}
		
		// calculate the return value
		double avg = sum / n_iterations;
		double res = avg / 255 * (max_val - min_val) + min_val;
		return res;
	} else if( map_type == TextureCylinder) {

		// project the three points on a rectangle that would be of
		// size 2*PI x height and then proceed as with the 2D
		// mapping
		min_x = 0; max_x = 2*M_PI; min_y = min_z; max_y = max_z;
		cyl_to_plane( x1, y1, z1);
		cyl_to_plane( x2, y2, z2);
		cyl_to_plane( x3, y3, z3);
		
		// Iterative process:
		//  - generate random coordinates inside the triangle,
		//    transform them to the picture image, add to the
		//    sum
		//  - from the sum, min_val and max_val and number of
		//    iterations determine the return value

		long n_iterations = 100; // TBD - add interface to the
				// class to set this up, so that it is
				// not a constant

		// figure out the texture coordinates of the triangle
		double tx1 = (x1-min_x) / (max_x-min_x) * image-> width;
		double ty1 = (y1-min_y) / (max_y-min_y) * image-> height;
		double tx2 = (x2-min_x) / (max_x-min_x) * image-> width;
		double ty2 = (y2-min_y) / (max_y-min_y) * image-> height;
		double tx3 = (x3-min_x) / (max_x-min_x) * image-> width;
		double ty3 = (y3-min_y) / (max_y-min_y) * image-> height;

		// do the loop
		double sum = 0;
		for( long i = 0 ; i < n_iterations ; i ++) {
			// generate random isoparametric coordinates p,q=(0,1)
			double p = drand48();
			double q = drand48();
			// make sure p+q <= 1
			if( p+q > 1) { p = 1-p; q = 1-q; }
			// transform p,q into texture coordinates
			long x = long( (tx2-tx1)*p + (tx3-tx1)*q + tx1);
			long y = long( (ty2-ty1)*p + (ty3-ty1)*q + ty1);
			// wrap around to fit into the image
			while( x < 0) x += image-> width;
			while( x >= image-> width) x-= image-> width;
			while( y < 0) y += image-> height;
			while( y >= image-> height) y -= image-> height;
			sum += image-> get( y, x);
		}
		
		// calculate the return value
		double avg = sum / n_iterations;
		double res = avg / 255 * (max_val - min_val) + min_val;
		return res;
	} else {
		assert( 0); // get_val() invoked on bad map_type
	}
	// just so that the compiler does not complain:
	return 0;
}


double Map::get_min_value(
	double x1, double y1, double z1,
	double x2, double y2, double z2,
	double x3, double y3, double z3)
// ----------------------------------------------------------------------
// given the coordinates of a triangle, calculate the min. pixel value
// in the triangle if it was to be drawn on the texture map
// ----------------------------------------------------------------------
{
	if( map_type == Constant) {

		return min_val;

	} else if( map_type == Normal) {
		
		return min_x;

	} else if( map_type == Texture2D) {
		
		// Iterative process:
		//  - generate random coordinates inside the triangle,
		//    transform them to the picture image, add to the
		//    sum
		//  - from the sum, min_val and max_val and number of
		//    iterations determine the return value

		long n_iterations = 100; // TBD - add interface to the
				// class to set this up, so that it is
				// not a constant

		// figure out the texture coordinates of the triangle
		double tx1 = (x1-min_x) / (max_x-min_x) * image-> width;
		double ty1 = (y1-min_y) / (max_y-min_y) * image-> height;
		double tx2 = (x2-min_x) / (max_x-min_x) * image-> width;
		double ty2 = (y2-min_y) / (max_y-min_y) * image-> height;
		double tx3 = (x3-min_x) / (max_x-min_x) * image-> width;
		double ty3 = (y3-min_y) / (max_y-min_y) * image-> height;

		// do the loop
		long min = -1;
		for( long i = 0 ; i < n_iterations ; i ++) {
			// generate random isoparametric coordinates p,q=(0,1)
			double p = drand48();
			double q = drand48();
			// make sure p+q <= 1
			if( p+q > 1) { p = 1-p; q = 1-q; }
			// transform p,q into texture coordinates
			long x = long( (tx2-tx1)*p + (tx3-tx1)*q + tx1);
			long y = long( (ty2-ty1)*p + (ty3-ty1)*q + ty1);
			// wrap around to fit into the image
			while( x < 0) x += image-> width;
			while( x >= image-> width) x-= image-> width;
			while( y < 0) y += image-> height;
			while( y >= image-> height) y -= image-> height;
			long val = image-> get( y, x);
			if( min == -1 || val < min)
				min = val;
		}
		
		// calculate the return value
		double res = (min / 255.0) * (max_val - min_val) + min_val;
		return res;
	} else if( map_type == TextureCylinder) {

		// project the three points on a rectangle that would be of
		// size 2*PI x height and then proceed as with the 2D
		// mapping
		min_x = 0; max_x = 2*M_PI; min_y = min_z; max_y = max_z;
		cyl_to_plane( x1, y1, z1);
		cyl_to_plane( x2, y2, z2);
		cyl_to_plane( x3, y3, z3);
		
		// Iterative process:
		//  - generate random coordinates inside the triangle,
		//    transform them to the picture image, add to the
		//    sum
		//  - from the sum, min_val and max_val and number of
		//    iterations determine the return value

		long n_iterations = 100; // TBD - add interface to the
				// class to set this up, so that it is
				// not a constant

		// figure out the texture coordinates of the triangle
		double tx1 = (x1-min_x) / (max_x-min_x) * image-> width;
		double ty1 = (y1-min_y) / (max_y-min_y) * image-> height;
		double tx2 = (x2-min_x) / (max_x-min_x) * image-> width;
		double ty2 = (y2-min_y) / (max_y-min_y) * image-> height;
		double tx3 = (x3-min_x) / (max_x-min_x) * image-> width;
		double ty3 = (y3-min_y) / (max_y-min_y) * image-> height;

		// do the loop
		long min = -1;
		for( long i = 0 ; i < n_iterations ; i ++) {
			// generate random isoparametric coordinates p,q=(0,1)
			double p = drand48();
			double q = drand48();
			// make sure p+q <= 1
			if( p+q > 1) { p = 1-p; q = 1-q; }
			// transform p,q into texture coordinates
			long x = long( (tx2-tx1)*p + (tx3-tx1)*q + tx1);
			long y = long( (ty2-ty1)*p + (ty3-ty1)*q + ty1);
			// wrap around to fit into the image
			while( x < 0) x += image-> width;
			while( x >= image-> width) x-= image-> width;
			while( y < 0) y += image-> height;
			while( y >= image-> height) y -= image-> height;
			long val = image-> get( y, x);
			if( min == -1 || val < min)
				min = val;
		}
		
		// calculate the return value
		double res = (min / 255.0) * (max_val - min_val) + min_val;
		return res;
	} else {
		assert( 0); // get_val() invoked on bad map_type
	}
	// just so that the compiler does not complain:
	return 0;
}
