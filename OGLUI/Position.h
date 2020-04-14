#ifndef __OGLUI_POSITION_H__
#define __OGLUI_POSITION_H__

#include <OGLUInamespace.h>

OGLUI_BEGIN_NAMESPACE

class Position {
public:
	double x, y;
	Position (double px, double py) : x (px), y (py) {;}
	Position () : x (0), y (0) {;}
};

Position operator + (const Position & p1, const Position & p2);
Position operator - (const Position & p1, const Position & p2);
bool operator == (const Position & p1, const Position & p2);
bool operator != (const Position & p1, const Position & p2);

OGLUI_END_NAMESPACE

#endif

