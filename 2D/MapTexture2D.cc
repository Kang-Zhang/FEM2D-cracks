#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "MapTexture2D.hh"
#include "die.hh"
#include "Distribution.hh"

MapTexture2D::MapTexture2D( 
	const char * p_fname,
	double p_min_val, double p_max_val,
	double p_min_x, double p_max_x,
	double p_min_y, double p_max_y)
{
	// read the name of the texture image file
	if( p_fname != NULL)
		fname = strdup( p_fname);
	else
		fname = NULL;
	min_val = p_min_val;
	max_val = p_max_val;
	min_x = p_min_x;
	max_x = p_max_x;
	min_y = p_min_y;
	max_y = p_max_y;

	if( fname != NULL) {
		image = new Image();
		if( image-> load( fname))
			die( "Could not load image: %s", fname);
	} else {
		image = NULL;
	}
}

MapTexture2D::~MapTexture2D( )
{
	if( fname != NULL) free( fname);
	if( image != NULL) delete image;
}

void MapTexture2D::print( FILE * fp)
{
	fprintf( fp, "Texture2D \"%s\" %.15f %.15f %.15f %.15f"
		 " %.15f %.15f",
		 fname == NULL ? "<NULL>" : fname,
		 min_val, max_val, min_x, max_x, min_y, max_y);
}

double MapTexture2D::get_value(
	double x, double y)
{
	if( image == NULL) {
		return min_val;
	}

	// figure out the texture coordinates of the triangle
	long tx = (x-min_x) / (max_x-min_x) * image-> width;
	long ty = (y-min_y) / (max_y-min_y) * image-> height;

	if( tx < 0) tx = 0;
	if( ty < 0) ty = 0;
	if( tx >= image-> width) tx = image-> width-1;
	if( ty >= image-> height) ty = image-> height-1;

	return (image-> get( ty, tx) / 255.0) * (max_val-min_val) + min_val;
	
}

double MapTexture2D::get_value(
	double x1, double y1, double z1,
	double x2, double y2, double z2,
	double x3, double y3, double z3)
// ----------------------------------------------------------------------
// given the coordinates of a triangle, calculate the average pixel value
// in the triangle if it was to be drawn on the texture map
// ----------------------------------------------------------------------
{
	if( image == NULL) {
		return min_val;
	}
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
}
