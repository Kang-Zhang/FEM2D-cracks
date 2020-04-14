#include <Position.h>

OGLUI_BEGIN_NAMESPACE

Position operator + (const Position & p1, const Position & p2)
{
	return Position (p1.x + p2.x, p1.y + p2.y);
}

Position operator - (const Position & p1, const Position & p2)
{
	return Position (p1.x - p2.x, p1.y - p2.y);
}

bool operator == (const Position & p1, const Position & p2)
{
	return p1.x == p2.x && p1.y == p2.y;
}
bool operator != (const Position & p1, const Position & p2)
{
	return ! (p1.x == p2.x && p1.y == p2.y);
}

OGLUI_END_NAMESPACE

