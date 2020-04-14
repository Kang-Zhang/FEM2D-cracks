#include <cstdlib>
#include <cstdio>
#include "OGLUI.h"

namespace OGLUI {
#ifdef INDENTATION_HACK
} // to fool emacs into not indenting the rest of the file
#endif

// ======================================================================
// Default theme
// ......................................................................
static SmartPointer <Theme> _defaultTheme = NULL;

SmartPointer <Theme> getDefaultTheme ()
{
	if (_defaultTheme.isNull ())
	{
		_defaultTheme = new Theme ("Default OGLUI theme");
		_defaultTheme-> set ("*background", new Color(0,0,0,0));
	}
	return _defaultTheme;
}

void glBox (double x1, double y1, double x2, double y2)
{
	glBegin (GL_QUADS);
	glVertex2d (x1, y1);
	glVertex2d (x2, y1);
	glVertex2d (x2, y2);
	glVertex2d (x1, y2);
	glEnd ();
}

void glBeveledBox (
	double minx, double miny, // bottom left corner
	double maxx, double maxy, // top right corner
	double bt,		// bevel size top
	double bb,		// bevel size bottom
	double bl,		// bevel size left
	double br,		// bevel size right
	Color::Pointer col[3],
	double colorWeights[12]
	)
{
	// name the sizes for convenience
	double x1 = minx;
	double x2 = minx + bl;
	double x3 = maxx - br;
	double x4 = maxx;
	double y1 = miny;
	double y2 = miny + bb;
	double y3 = maxy - bt;
	double y4 = maxy;
	double & tw1 = colorWeights [0];
	double & tw2 = colorWeights [1];
	double & tw3 = colorWeights [2];
	double & bw1 = colorWeights [3];
	double & bw2 = colorWeights [4];
	double & bw3 = colorWeights [5];
	double & lw1 = colorWeights [6];
	double & lw2 = colorWeights [7];
	double & lw3 = colorWeights [8];
	double & rw1 = colorWeights [9];
	double & rw2 = colorWeights [10];
	double & rw3 = colorWeights [11];
	Color::Pointer _col1 = col[0];
	Color::Pointer _col2 = col[1];
	Color::Pointer _col3 = col[2];
	// top edge
	Color tc (_col1-> r * tw1 + _col2-> r * tw2 + _col3-> r * tw3,
		  _col1-> g * tw1 + _col2-> g * tw2 + _col3-> g * tw3,
		  _col1-> b * tw1 + _col2-> b * tw2 + _col3-> b * tw3,
		  _col1-> a * tw1 + _col2-> a * tw2 + _col3-> a * tw3);
	tc.glColor ();
	glBegin (GL_QUADS);
 	glVertex2d (x2,y3);
 	glVertex2d (x3,y3);
 	glVertex2d (x4,y4);
 	glVertex2d (x1,y4);
	glEnd ();
	// bottom edge
	Color bc (_col1-> r * bw1 + _col2-> r * bw2 + _col3-> r * bw3,
		  _col1-> g * bw1 + _col2-> g * bw2 + _col3-> g * bw3,
		  _col1-> b * bw1 + _col2-> b * bw2 + _col3-> b * bw3,
		  _col1-> a * bw1 + _col2-> a * bw2 + _col3-> a * bw3);
	bc.glColor ();
	glBegin (GL_QUADS);
 	glVertex2d (x1,y1);
 	glVertex2d (x4,y1);
 	glVertex2d (x3,y2);
 	glVertex2d (x2,y2);
	glEnd ();
	// left edge
	Color lc (_col1-> r * lw1 + _col2-> r * lw2 + _col3-> r * lw3,
		  _col1-> g * lw1 + _col2-> g * lw2 + _col3-> g * lw3,
		  _col1-> b * lw1 + _col2-> b * lw2 + _col3-> b * lw3,
		  _col1-> a * lw1 + _col2-> a * lw2 + _col3-> a * lw3);
	lc.glColor ();
	glBegin (GL_QUADS);
	glVertex2d (x1,y1);
	glVertex2d (x2,y2);
	glVertex2d (x2,y3);
	glVertex2d (x1,y4);
	glEnd ();
	// right edge
	Color rc ( _col1-> r * rw1 + _col2-> r * rw2 + _col3-> r * rw3,
		   _col1-> g * rw1 + _col2-> g * rw2 + _col3-> g * rw3,
		   _col1-> b * rw1 + _col2-> b * rw2 + _col3-> b * rw3,
		   _col1-> a * rw1 + _col2-> a * rw2 + _col3-> a * rw3);
	rc.glColor ();
	glBegin (GL_QUADS);
	glVertex2d (x4,y1);
	glVertex2d (x4,y4);
	glVertex2d (x3,y3);
	glVertex2d (x3,y2);
	glEnd ();
}


void error (std::ostream & os, const Widget * w, const std::string & str)
{
	os << "OGLUI::";
	if (w) os << w-> className () << "(" << w-> instanceName () << "):";
	os << "Error:";
	os << str;
	os << "\n";
}

void error (const Widget * w, const std::string & str)
{
	error (std::cerr, w, str);
}

void error (const std::string & str)
{
	error (NULL, str);
}


#ifdef INDENTATION_HACK
} // to fool emacs into not indenting the rest of the file
#endif
}
