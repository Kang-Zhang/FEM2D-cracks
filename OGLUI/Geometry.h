#ifndef __OGLUI_GEOMETRY_H__
#define __OGLUI_GEOMETRY_H__

#include <iostream>
#include <Position.h>
#include <Size.h>
#include <algorithm>
#include <ostream>

namespace OGLUI {
#ifdef INDENTATION_HACK
} // to fool emacs into not indenting the rest of the file
#endif

class Geometry {
private:
//	double _x, _y, _width, _height;
	Position _position;
	Size _size;
public:
	// empty geometry constructor
	Geometry ()
		: _position (0,0)
		, _size (100,100)
		{ ; }
	// constructor using position & size
	Geometry ( const Position & pos
		   , const Size & size )
		: _position ( pos )
		, _size ( size )
		{ ; }
	// constructor using two positions
	Geometry ( const Position & pos1, const Position & pos2)
		{
			double x1 = std::min (pos1.x, pos2.x);
			double y1 = std::min (pos1.y, pos2.y);
			double x2 = std::max (pos1.x, pos2.x);
			double y2 = std::max (pos1.y, pos2.y);
			_position.x = x1;
			_position.y = y1;
			_size.width = x2-x1;
			_size.height = y2-y1;
		}
	/// convenience functions = to get the bounding box
	double minx () const { return _position.x ; }
	double maxx () const { return _position.x + _size.width; }
	double miny () const { return _position.y; }
	double maxy () const { return _position.y + _size.height; }
	/// calculate the horiz/vert center
	double midx () const { return (minx () + maxx ()) / 2.0; }
	double midy () const { return (miny () + maxy ()) / 2.0; }
	/// easy access to width and height
	double width () const { return _size.width; }
	double height () const { return _size.height; }
	/// get position
	const Position & getPosition () const { return _position; }
	/// set position
	void setPosition (const Position & pos) { _position = pos; }
	/// get size
	const Size & getSize () const { return _size; }
	/// set size
	void setSize (const Size & size) { _size = size; }
	/// comparison operator
	bool operator == (const Geometry & g) const;
	bool operator != (const Geometry & g) const;
	// access 
	// get/set the position
//	const double & x () const { return _x; }
//	double & x () { return _x; }
//	void x (double px) { _x = px; }
//	const double & y () const { return _y; }
//	double & y () { return _y; }
//	void y (double py) { _y = py; }
	/// get/set width
//	const double & width () const { return _width; }
//	double & width () { return _width; }
	/// get/set height
//	const double & height () const { return _height; }
//	double & height () { return _height; }
	/// get center
//	double cx () const { return _x + _width/2 ; }
//	double cy () const { return _y + _height/2 ; }
};

std::ostream & operator<< (std::ostream & os, const Geometry & g);

/// calculates the intersection of two rectangles specified
Geometry intersection (const Geometry & g1, const Geometry & g2);

#ifdef INDENTATION_HACK
{ // to fool emacs into not indenting the rest of the file
#endif
} // namespace OGLUI


#endif

