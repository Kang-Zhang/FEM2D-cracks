#ifndef __OGLUI_HH__
#define __OGLUI_HH__

#include <vector>
#include <string>
#include <cassert>
#include <iostream>
#include <GL/glut.h>
#include <algorithm>

#include <OGLUInamespace.h>
#include <Widget.h>
#include <Color.h>
#include <SmartPointer.h>
#include <Font.h>

OGLUI_BEGIN_NAMESPACE

#define NI { cerr << "Not implemented:" << __FILE__ << ":" << __LINE__ << "\n"; assert (0); }

// default theme
SmartPointer <Theme> getDefaultTheme ();

// some convenience functions
void glBox (double x1, double y1, double x2, double y2);
void glBeveledBox (
	double minx, double miny, // bottom left corner
	double maxx, double maxy, // top right corner
	double bt,		// bevel size top
	double bb,		// bevel size bottom
	double bl,		// bevel size left
	double br,		// bevel size right
	Color::Pointer col[3],
	double colorWeights[12]
	);

// error reporting function
void error (std::ostream & os, const Widget * w, const std::string & str);
void error (const Widget * w, const std::string & str);
void error (const std::string & str);

OGLUI_END_NAMESPACE

#endif
