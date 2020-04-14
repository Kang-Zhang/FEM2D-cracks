#ifndef __MAP_HH__
#define __MAP_HH__

#include "Tokenizer.hh"
#include "Image.hh"

class Map {

	double min_val, max_val; // cached min/max for the given distribution
	double min_x, max_x,	//  mapping for 2D textures
		min_y, max_y;
	double min_z, max_z;	// mapping for a cylinder
	char * fname;

public:
	enum MapType { Constant, Normal, Texture2D, TextureCylinder } map_type;

	Image * image;

	Map( FILE * fp);
	~Map( );
	Map( Tokenizer & tok);
	void print( FILE * fp);
	double get_value(
		double x1, double y1, double z1,
		double x2, double y2, double z2,
		double x3, double y3, double z3);
	double get_min_value(
		double x1, double y1, double z1,
		double x2, double y2, double z2,
		double x3, double y3, double z3);
	double get_min( void) { return min_val; }
	double get_max( void) { return max_val; }
private:
	void cyl_to_plane( double & x, double & y, double & z);
};


#endif
