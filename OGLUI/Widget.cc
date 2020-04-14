#include <GL/glu.h>
#include <Widget.h>
#include <OGLUI.h>
#include <sstream>
#include <OGLUInamespace.h>

OGLUI_BEGIN_NAMESPACE

Widget::Widget ( Widget * p_parent
		 , const std::string & p_name)
// ======================================================================
// Constructor
// ......................................................................
	: _bgColor (new Color (0.7,0.7,0.7,1.0))
	, _instanceName (p_name)
	, _borderColor (new Color (0,0,0,1))
	, _minSize (10, 10)
	, _maxSize (-1, -1)
	, _prefSize (100, 25)
	, _layMan (NULL)
	, _borderThickness (1)
	, parent (p_parent)
	, _invisible (false)
{
	if (parent)
		parent -> addChild (this);
}

Widget::~Widget ()
// ======================================================================
// destructor
// ......................................................................
{
	// empty
}

void Widget::_redraw ()
{
	// if this widget is invisible, don't do anything at all
	if (_invisible) return;

	// setup OpenGL to draw inside this widget, clipped properly
	setupOutOGLclipped ();

	// run my own redraw
	redraw ();
	
	// now run children's _redraw()
	for (WidgetList::iterator w = children.begin ();
	     w != children.end ();
	     w ++)
	{
		(*w)-> _redraw ();
	}

	// and, run my own redraw2
	redraw2 ();
}

bool Widget::mouse (const MouseEvent & ev)
{
	// if the mouse is not in our bounding box, return unprocessed event
	if (! ev.inBox (getOutGeom ())) return false;
	// ask our children to process the event
	propagateMouseEvent (ev);
	// in any case, return true
	return true;
}

bool Widget::propagateMouseEvent (const MouseEvent & ev)
{
	// check the children in reverse order!
//	for (WidgetList::iterator w = children.begin ();
//	     w != children.end();
//	     w ++)
	for (WidgetList::reverse_iterator w = children.rbegin ();
	     w != children.rend ();
	     w ++)
	{
		if ((*w)-> mouse (ev)) return true;
	}
	// none of the children accepted the event, so return false
	return false;
}

void Widget::redraw ()
{
	redrawBorder ();
	// draw the background
	_bgColor-> glColor();
	double t = _borderThickness;
	Geometry g = getOutGeom ();
	glBox (t, t, g.width()-t, g.height()-t);
}

void Widget::redrawBorder ()
{
	// draw the border (only if it is larger than 0 pixels)
	double t = _borderThickness;
	if (t > 0)
	{
		_borderColor -> glColor();
		glBegin (GL_QUADS);
		// bottom border
		Geometry g = getOutGeom ();
		glBox (0, 0, g.width(), t);
		// top border
		glBox (0, g.height()-t,
		       g.width(), g.height());
		// left border
		glBox (0, t, t, g.height()-t);
		// right border
		glBox (g.width()-t, t,
		       g.width(), g.height ()-t);
		glEnd ();
	}
}

void Widget::redraw2 ()
{
	// again, do nothing
}

bool Widget::_init ()
{
	// run my own method, and then all children's run methods
	if( ! init ()) return false;

	for (WidgetList::iterator w = children.begin();
	     w != children.end();
	     w ++)
	{
		bool status = (*w)-> _init();
		if( ! status) return false;
	}

	return true;
}

bool Widget::init ()
{
	// do nothing
	return true;
}

void Widget::postRedisplay ()
{
	// default behaviour is to ask the parent to redisplay itself
	if (parent != NULL) parent-> postRedisplay();
}

void Widget::transparency (double t)
{
	_bgColor-> a = t;
	_borderColor-> a = t;
	postRedisplay ();
}

/// widget can call this if it needs to have exclusive access to
/// the mouse (the returned value indicates success / failure)
bool Widget::captureMouse (Widget * w)
{
	// who wants the mouse
	Widget * who = (w == NULL) ? this : w;

	// delegate the responsibility to the parent
	if (parent)
	{
		return parent-> captureMouse (who);
	}
	else
	{
		OGLUI::error (
			this,
			"captureMouse(): unhandled in top level widget.");
		return false;
	}
}

/// releases the mouse
void Widget::releaseMouse (Widget * w)
{
	// who is releasing the mouse
	Widget * who = (w == NULL) ? this : w;
	// ask the parent to release the mouse
	if (parent)
	{
		return parent-> releaseMouse (who);
	}
	else
	{
		OGLUI::error (
			this,
			"releaseMouse(): unhandled in top level widget.");
	}
}

/// returns the full identifier of the widget
std::string Widget::id()
{
	std::ostringstream os;
	os << className() << "(" << _instanceName << "," << this << ")";
	return os.str();
}

const Size & Widget::getMinSize () const
// ======================================================================
// sets the min. size of the widget
// ......................................................................
{
	return _minSize;
}

void Widget::setMinSize (const Size & size)
// ======================================================================
// gets the min. size of the widget
// ......................................................................
{
	_minSize = size;
}

const Size & Widget::getMaxSize () const
// ======================================================================
// sets the max. size of the widget
// ......................................................................
{
	return _maxSize;
}

void Widget::setMaxSize (const Size & size)
// ======================================================================
// gets the max. size of the widget
// ......................................................................
{
	_maxSize = size;
}

const Size & Widget::getPrefSize () const
// ======================================================================
// sets the preferred size of the widget
// ......................................................................
{
	return _prefSize;
}

void Widget::setPrefSize (const Size & size)
// ======================================================================
// gets the preferred size of the widget
// ......................................................................
{
	_prefSize = size;
}

Geometry Widget::getInGeom () const
// ======================================================================
// returns the inside geometry of the widget
//   - this method is used usually by the derived classes, so that they
//     can use the super-class's drawOutside() method, and then know
//     where the 'inside' is
// ......................................................................
{
	Geometry g = getOutGeom ();
	Position pos = g.getPosition ();
	pos.x += _borderThickness;
	pos.y += _borderThickness;
	Size size = g.getSize ();
	size.width -= 2 * _borderThickness;
	size.height -= 2 * _borderThickness;
	return Geometry (pos, size);
}

Geometry Widget::getInGeomClipped () const
// ======================================================================
// get inside geometry - but clipped w.r.t. to all parents
// ......................................................................
{
	// get my own inside geometry
	Geometry g = getInGeom ();
	// if I don't have a parent, return my unclipped geometry
	if (! parent)
		return g;
	// get parent's clipped inside geometry
	Geometry gp = parent-> getInGeomClipped ();
	// return the intersection of these two geometries
	return intersection (g, gp);
}

Geometry Widget::getOutGeom () const
// ======================================================================
// returns the geometry of the whole widget - the geometry is specified
// in global terms, i.e. in the coordinate system of the closest Window
// ......................................................................
{
	return _geom;
/*
	if (parent)
	{
		// ask the parent for my geometry
		return parent-> getChildGeometry (this);
	}
	else
	{
		// nobody to ask for geometry -- weird - report error
		OGLUI::error (
			this
			, "getOutGeom(): unhandled in top level widget.");
		return Geometry (Position (0,0), _prefSize);
	}
*/
}

Geometry Widget::getOutGeomClipped () const
// ======================================================================
// get inside geometry - but clipped w.r.t. to all parents
// ......................................................................
{
	// get my own outside geometry
	Geometry g = getOutGeom ();
	// if I don't have a parent, my entire outside geometry is the result
	if (! parent)
		return getOutGeom ();
	// get the clipped geometry of the parent
	Geometry gp = parent-> getInGeomClipped ();
	// return the intersection
	return intersection (g, gp);
}

void Widget::addChild (Widget * w)
// ======================================================================
// adds a child to this widget
// ......................................................................
{
	children.push_back (w); w-> parent = this;
}

LayoutManager::Pointer Widget::setLayoutManager (LayoutManager::Pointer lm)
// ======================================================================
// sets the layout manager for this widget
// WARNING - destroys the old layout manager
// ......................................................................
{
	_layMan = lm;
	return lm;
}

LayoutManager::Pointer Widget::getLayoutManager ()
// ======================================================================
// returns the layout manager for this widget
// ......................................................................
{
	return _layMan;
}

/*
Geometry Widget::getChildGeometry (const Widget * w)
// ======================================================================
// returns the geometry of the child - in 'this' widget's local coordinate
// system
// ......................................................................
{
	if (_layMan != NULL)
	{
		// if there is a layout manager - ask it for geometry
		return _layMan -> getGeometry (w, getInGeom ());
	}
	else
	{
		// otherwise return trivial geometry
		return Geometry (Position (0,0), w-> getPrefSize ());
	}
}
*/

void Widget::setupInOGL () const
// ======================================================================
// sets up opengl state to be able to draw directly inside the widget
// ......................................................................
{
	// setup the viewport for the inside of the widget
	Geometry g = getInGeom ();
	glViewport (GLint (g.getPosition ().x)
		    , GLint (g.getPosition ().y)
		    , GLint (g.getSize ().width)
		    , GLint (g.getSize ().height));
	// setup the matrices
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	Size s = g.getSize ();
	gluOrtho2D (0, s.width, 0, s.height);
	glMatrixMode (GL_MODELVIEW);
	glLoadIdentity ();
}

void Widget::setupOutOGLclipped () const
// ======================================================================
// sets up opengl state to draw into the whole area of the widget
// ......................................................................
{
	Geometry g = getOutGeom ();
	glViewport (GLint (g.getPosition ().x)
		    , GLint (g.getPosition ().y)
		    , GLint (g.getSize ().width)
		    , GLint (g.getSize ().height));
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	Size s = g.getSize ();
	gluOrtho2D (0, s.width, 0, s.height);
	glMatrixMode (GL_MODELVIEW);
	glLoadIdentity ();
	// get the outside geometry clipped
	g = getOutGeomClipped ();
	// setup the scissor box
	glScissor (GLint(g.minx()), GLint(g.miny()),
		   GLint(g.width()), GLint(g.height()));
	glEnable (GL_SCISSOR_TEST);

	// reset OGL...
	glDisable (GL_DEPTH_TEST);
	glDisable (GL_LIGHTING);
	glEnable (GL_BLEND);
}

void Widget::setupInOGLclipped () const
// ======================================================================
// sets up opengl state to be able to draw directly inside the widget,
// plus it sets up glScissor() to clip anything outside
// ......................................................................
{
	// setup unclipped OGL
	setupInOGL ();
	// get the inside geometry - clipped
	Geometry g = getInGeomClipped ();
	// setup the scissor box
	glScissor (GLint(g.minx()), GLint(g.miny()),
		   GLint(g.width()), GLint(g.height()));
}

const Color::Pointer Widget::getBgColor () const
// ======================================================================
// get the background color
// ......................................................................
{
	return _bgColor;
}

void Widget::setBgColor (const Color::Pointer col)
// ======================================================================
// sets the background color
// ......................................................................
{
	_bgColor = col;
	postRedisplay ();
}

const std::string & Widget::instanceName() const
// ======================================================================
// name of the instance
// ......................................................................
{
	return _instanceName;
}

const std::string Widget::className () const
// ======================================================================
// returns name of the class
// ......................................................................
{
	return "Widget";
}

bool Widget::invisible () const
// ======================================================================
// gets invisibility status
// ......................................................................
{
	return _invisible;
}

void Widget::invisible (bool i)
// ======================================================================
// sets invisibility status
// ......................................................................
{
	if (_invisible != i)
	{
		postRedisplay ();
		_invisible = i;
	}
}

/*
void Widget::setGeometry (const Geometry & p_geom)
// ======================================================================
// sets the geometry
// ......................................................................
{
//	_geom = p_geom;
}
*/

void Widget::borderThickness (double t)
// ======================================================================
// sets the border thickness
// ......................................................................
{
	_borderThickness = t;
	postRedisplay();
}

double Widget::borderThickness () const 
// ======================================================================
// gets the border thickness
// ......................................................................
{
	return _borderThickness;
}

void Widget::setBorderColor (const Color::Pointer c)
// ======================================================================
// sets the border color
// ......................................................................
{
	_borderColor = c;
	postRedisplay ();
}

void Widget::printTree (const std::string prefix)
// ======================================================================
// debuging method - to list the children
// ......................................................................
{
	std::cerr << prefix
		  << className() << ":" << _instanceName
		  << " (" << this << "):\n";
	for (WidgetList::iterator i = children.begin ();
	     i != children.end();
	     i ++)
	{
		(*i) -> printTree (prefix + "    ");
	}
}

void Widget::onGeometryChanged (const Geometry & g)
// ======================================================================
// This method is called by the parent whenever the geometry of the
// child has been modified.
// ......................................................................
{
	// cache the geometry given to us
	_geom = g;
	// signal redisplay
	postRedisplay ();
	// ask the layout manager to recompute kids
	if (_layMan != NULL)
		_layMan -> recompute ( getInGeom ());
	return ;
	// ask the layout manager if any of the children should change
	// their geometry
	if (_layMan == NULL) return;
	for (WidgetList::iterator i = children.begin ()
		     ; i < children.end ()
		     ; i ++)
	{
		Geometry cg = _layMan -> getGeometry (*i, getInGeom ());
		if ((*i) -> _geom != cg)
			(*i) -> onGeometryChanged (cg);
	}
	
}

OGLUI_END_NAMESPACE

