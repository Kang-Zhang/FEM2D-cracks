#ifndef __MAP_TEXTURE_2D_HH__
#define __MAP_TEXTURE_2D_HH__

#include "Image.hh"

class MapTexture2D {

	double min_val, max_val; // cached min/max for the given distribution
	double min_x, max_x,	//  mapping for 2D textures
		min_y, max_y;
	char * fname;

public:
	Image * image;

	MapTexture2D( const char * p_fname,
		      double p_min_val, double p_max_val,
		      double p_min_x, double p_max_x,
		      double p_min_y, double p_max_y);
	~MapTexture2D( );
	void print( FILE * fp);
	double get_value( double x1, double y1, double z1,
			  double x2, double y2, double z2,
			  double x3, double y3, double z3);
	double get_value( double x, double y);
	double get_min( void) { return min_val; }
	double get_max( void) { return max_val; }
};


#endif
