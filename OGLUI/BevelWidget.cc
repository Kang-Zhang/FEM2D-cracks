#include "BevelWidget.h"

#include <OGLUI.h>

namespace OGLUI {
#ifdef INDENTATION_HACK
} // to fool emacs into not indenting the rest of the file
#endif

// constructor
BevelWidget::BevelWidget (
	Widget * p_parent,
	const std::string & p_name)
	: Widget (p_parent, p_name)
	, _bevelSizeTop (2)
	, _bevelSizeBottom (2)
	, _bevelSizeLeft (2)
	, _bevelSizeRight (2)
	, _col2 (Color::White)
	, _col3 (Color::Black)
	, _lightWeightsTop (0.5,0.5,0)
	, _lightWeightsBottom (0.5,0,0.5)
	, _lightWeightsLeft (0.5,0.5,0)
	, _lightWeightsRight (0.5,0,0.5)
{;}

void
BevelWidget::redrawBorder ()
{
	// draw the Widget's border
	Widget::redrawBorder ();
	// draw the beveled border
	// name the sizes for convenience
	double bt = _bevelSizeTop;
	double bb = _bevelSizeBottom;
	double bl = _bevelSizeLeft;
	double br = _bevelSizeRight;
	double t = _borderThickness;
	Geometry g = getOutGeom ();
	double w = g.width ();
	double h = g.height ();
	double x1 = t; double x4 = w-t;
	double y1 = t; double y4 = h-t;

	Color::Pointer cols [3] = {
		getBgColor(),
		_col2,
		_col3};
	double colWeights [12] =
		{
			_lightWeightsTop.w1 (),
			_lightWeightsTop.w2 (),
			_lightWeightsTop.w3 (),
			_lightWeightsBottom.w1 (),
			_lightWeightsBottom.w2 (),
			_lightWeightsBottom.w3 (),
			_lightWeightsLeft.w1 (),
			_lightWeightsLeft.w2 (),
			_lightWeightsLeft.w3 (),
			_lightWeightsRight.w1 (),
			_lightWeightsRight.w2 (),
			_lightWeightsRight.w3 ()
		};

	glBeveledBox (x1,y1,x4,y4,
		      bt,bb,bl,br,
		      cols,
		      colWeights);
}

void
BevelWidget::redraw ()
{
	// draw the border
	redrawBorder ();
	// name the sizes for convenience
	double bt = _bevelSizeTop;
	double bb = _bevelSizeBottom;
	double bl = _bevelSizeLeft;
	double br = _bevelSizeRight;
	double t = _borderThickness;
	Geometry g = getOutGeom ();
	double w = g.width ();
	double h = g.height ();
	double x2 = t + bl; double x3 = w-t-br;
	double y2 = t + bb; double y3 = h-t-bt;
	// redraw inside
	getBgColor () -> glColor ();
	glBox (x2,y2, x3, y3);
}

#ifdef INDENTATION_HACK
{ // to fool emacs into not indenting the rest of the file
#endif
};

