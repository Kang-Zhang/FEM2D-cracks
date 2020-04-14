#include "MouseEvent.h"

ostream & operator<< (ostream & os, const MouseEvent & ev)
{
	os << "MouseEvent:(" << ev.x << "," << ev.y << " ";
	if (ev.reason == MouseEvent::Push)
		os << "Push button " << ev.button;
	else if (ev.reason == MouseEvent::Release)
		os << "Release button " << ev.button;
	else if (ev.reason == MouseEvent::Motion)
		os << "Motion";

	return os;
}
