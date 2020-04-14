#include <OGLUI.h>
#include <LabelWidget.h>

namespace OGLUI {
#ifdef INDENTATION_HACK
} // to fool emacs into not indenting the rest of the file
#endif

void
LabelWidget::redraw ()
{
	BevelWidget::redraw ();
	// if there is no label, we are done
	if (label.empty ()) return;
	// setup OGL inside
	BevelWidget::setupInOGLclipped ();
	// apply scissors - so that the text is not outside of the widget
	Geometry g = BevelWidget::getInGeom ();
//	glEnable (GL_SCISSOR_TEST);
	// set the color of the font
	fontColor-> glColor ();
	// get the dimensions of the label
	Size size = font-> getSize (label);
	// figure out the position
	Position pos;
	if (_justHoriz == LEFT)
		pos.x = 0;
	else if (_justHoriz == RIGHT)
		pos.x = g.width () - size.width;
	else
		pos.x = int(g.width()/2-size.width/2.0);
	if (_justVert == TOP)
		pos.y = 0;
	else if (_justVert == BOTTOM)
		pos.y = g.height () - size.height;
	else
		pos.y = int(g.height ()/2-size.height/2);
	// apply offset
	pos.x += _justOffset.x;
	pos.y += _justOffset.y;
	// render the text
	font-> render (pos.x, pos.y, label);

//	glDisable (GL_SCISSOR_TEST);
}

#ifdef INDENTATION_HACK
{ // to fool emacs into not indenting the rest of the file
#endif
};

