#include <Size.h>

OGLUI_BEGIN_NAMESPACE

bool operator == (const Size & p1, const Size & p2)
{
	return p1.width == p2.width && p1.height == p2.height;
}
bool operator != (const Size & p1, const Size & p2)
{
	return ! (p1.width == p2.width && p1.height == p2.height);
}

OGLUI_END_NAMESPACE

