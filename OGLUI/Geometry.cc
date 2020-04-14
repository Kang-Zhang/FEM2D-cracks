#include <algorithm>
#include "Geometry.h"

OGLUI_BEGIN_NAMESPACE

/*
ostream & operator<< (ostream & os, const Geometry & g)
{
	os << "Geometry:(" << g.minx() << "," << g.miny() << " "
	   << g.width() << "x" << g.height() << ")";
	return os;
}
*/

static void interval_intersection (double x1, double x2
				  , double x3, double x4
				  , double & rx1, double & rx2)
// ======================================================================
// Calculates inteval intersection
//  - first interval is specified by x1, x2
//  - second interval is specified by x3, x4
//  - result is stored in rx1, rx2
//  - if no intersection exists, rx1 = rx2
// ......................................................................
{
	// order x1..4
	if (x1 > x2) std::swap<double> (x1, x2);
	if (x3 > x4) std::swap<double> (x3, x4);

	// assume the result is x1, x2
	rx1 = x1; rx2 = x2;

	// clip by x3...
	rx1 = std::max (rx1, x3);
	rx2 = std::max (rx2, x3);
	// clip by ...x4
	rx1 = std::min (rx1, x4);
	rx2 = std::min (rx2, x4);

	// return the result
	return;
	
}

Geometry intersection (const Geometry & g1, const Geometry & g2)
// ======================================================================
// calculates the intersection of two rectangles specified if no
// intersection exists, the result is a geometry with Size (0,0)
// ......................................................................
{
	// calculate the intersection for the X-dimension
	double x1, x2;
	interval_intersection (g1.minx (), g1.maxx ()
			       , g2.minx (), g2.maxx ()
			       , x1, x2);
	// calculate the intersection for the Y-dimension
	double y1, y2;
	interval_intersection (g1.miny (), g1.maxy ()
			       , g2.miny (), g2.maxy ()
			       , y1, y2);
	// combine the result
	if (x1 == x2 || y1 == y2)
		return Geometry (Position (0,0), Size (0,0));
	else
		return Geometry (Position (x1, y1), Size (x2-x1, y2-y1));
}

bool Geometry::operator == (const Geometry & g) const
// ======================================================================
// compares itself to another geometry
// ......................................................................
{
	return (_position == g._position && _size == g._size);
}

bool Geometry::operator != (const Geometry & g) const
// ======================================================================
// compares itself to another geometry
// ......................................................................
{
	return ! (_position == g._position && _size == g._size);
}
OGLUI_END_NAMESPACE
