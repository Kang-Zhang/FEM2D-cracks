#include "OGLUInamespace.h"
#include "WindowWidget.h"
#include "UIWidget.h"

OGLUI_BEGIN_NAMESPACE

/// constructor
WindowWidget::WindowWidget (
		Widget * p_parent,
		const std::string & p_name)
	: Widget ( p_parent, p_name)
	, _glutDispMode (GLUT_RGB|GLUT_DOUBLE)
	, _glutWinId (-1)
{
	borderThickness (0);
}

bool
WindowWidget::init ()
{
	// create a window for ourselves
	glutInitDisplayMode (_glutDispMode);
	Size size = getPrefSize ();
	glutInitWindowSize (GLint (size.width), GLint (size.height));
	_glutWinId = glutCreateWindow ("WindowWidget");
	glutSetWindow (_glutWinId);
	// make sure the parent is a UI
	UI * ui = dynamic_cast <UI *> (parent);
	if (ui == NULL)
	{
		std::cerr << "OGLUI::Error:WindowWidget can only have UI "
			"as parent.\n";
		exit (-1);
	}
	// register our callbacks
	glutSetWindow (_glutWinId);
	glutDisplayFunc (UI::dispFunc);
	glutReshapeFunc (UI::reshapeFunc);
	glutMouseFunc (UI::mouseFunc);
	glutMotionFunc (UI::mouseMotionFunc);
	glutPassiveMotionFunc (UI::mouseMotionFunc);

        // enable anti-aliasing of lines
	glEnable( GL_BLEND);
	glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	return true;
}

void
WindowWidget::redraw ()
{
	Widget::redraw ();
}

void
WindowWidget::redraw2 ()
{
	glutSwapBuffers();
}

void
WindowWidget::reshape (int width, int height)
{
	onGeometryChanged (Geometry (Position(0,0), Size (width,height)));
//	_geom.setSize (Size (width, height));
}

void WindowWidget::postRedisplay ()
{
	if (_glutWinId != -1)
	{
//		glutSetWindow (_glutWinId);
//		glutPostRedisplay ();
		glutPostWindowRedisplay (_glutWinId);
	}
}

Geometry WindowWidget::getOutGeom () const
// ======================================================================
// my geometry is the one I know
// ......................................................................
{
	return _geom;
}

Geometry WindowWidget::getOutGeomClipped () const
// ======================================================================
// my geometry is the one I know
// ......................................................................
{
	return _geom;
}

Geometry WindowWidget::getInGeom () const
// ======================================================================
// returns the inside geometry
// ......................................................................
{
	return Widget::getInGeom ();
}

Geometry WindowWidget::getInGeomClipped () const
// ======================================================================
// returns the inside geometry
// ......................................................................
{
	return Widget::getInGeom ();
}

OGLUI_END_NAMESPACE
