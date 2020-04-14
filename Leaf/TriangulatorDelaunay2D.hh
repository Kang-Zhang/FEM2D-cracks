#ifndef __TRIANGULATOR_DELAUNAY_2D_HH__
#define __TRIANGULATOR_DELAUNAY_2D_HH__

#include <vector>

class TriangulatorDelaunay2D
{
private:
	class Point {
	public: double x, y;
		Point( double px, double py) {
			x = px; y = py;
		}
	};
	std::vector<Point> pts;
public:
	class Triangle { public: long p1, p2, p3; };

	void add_point( double x, double y);
	std::vector<Triangle> triangulate( void);
//	typedef std::vector<Triangle> tlist;
//	tlist triangulate( void);
};

#endif
