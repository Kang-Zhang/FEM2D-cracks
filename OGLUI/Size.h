#ifndef __OGLUI_SIZE_H__
#define __OGLUI_SIZE_H__

#include <OGLUInamespace.h>

namespace OGLUI {
#ifdef INDENTATION_HACK
} // to fool emacs into not indenting the rest of the file
#endif

class Size {
public:
	double width, height;
	Size () : width (0), height (0) {;}
	Size (double w, double h) : width (w), height (h) {;}
};

bool operator == (const Size & s1, const Size & s2);
bool operator != (const Size & s1, const Size & s2);

#ifdef INDENTATION_HACK
{ // to fool emacs into not indenting the rest of the file
#endif
} // namespace OGLUI


#endif

