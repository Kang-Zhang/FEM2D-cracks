#ifndef __OGLUI_SLIDER_H__
#define __OGLUI_SLIDER_H__

#include <OGLUInamespace.h>
#include <OGLUI.h>
#include <Slot.h>
#include <BevelWidget.h>
#include <cmath>
#include <algorithm>

OGLUI_BEGIN_NAMESPACE

//typedef BevelWidget SliderWidgetSuper;
typedef BevelWidget SliderWidgetSuper;

template <class T = int>
class SliderWidget : public SliderWidgetSuper {
public:
	typedef SmartPointer <SliderWidget> Pointer;
private:
	T _minVal;		// minimum value
	T _maxVal;		// max. value
	T _currVal;		// current value
	Color::Pointer _col1;	// color of the indicator on the left
	Color::Pointer _col2;	// color of the indicator on the right
	Color::Pointer _col3;	// color below the indicator
	double _marginWidth;	// margins around the slider
	double _marginHeight;
	double _grooveDepth;	// embossing size
	bool mouseDown;		// is the mouse down
	MouseEvent::Button mouseButton;	// which button is pushed
	int lastx, lasty;	// last mouse coordinates
public:
	Signal1<T> valueChangedSignal;
	// constructor
	SliderWidget (
		Widget * p_parent,
		const std::string & p_name = "slider"
		) :
		SliderWidgetSuper (p_parent, p_name),
		_minVal (0),
		_maxVal (100),
		_currVal (25),
		_col1 (Color::Orange),
		_col2 (Color::Red3),
//		_col3 (Color::Blue),
		_col3 (new Color(0,0,1,0.5)),
		_marginWidth (3),
		_marginHeight (3),
		_grooveDepth (1),
		mouseDown (false)
		{
			setBgColor (_col2);
		}
	virtual void redraw ();
	virtual bool mouse (const MouseEvent & ev);
	void minVal (double val) { _minVal = val; postRedisplay (); }
	void maxVal (double val) { _maxVal = val; postRedisplay (); }
	void currVal (double val) { _currVal = val; postRedisplay (); }
	
	void foo ();
};

template <class T>
void SliderWidget <T> ::redraw ()
{
	// draw the widget underneath
	SliderWidgetSuper::redrawBorder ();
	// setup opengl
	SliderWidgetSuper::setupInOGL ();
	// get the geometry of the area where I can draw
	Geometry g = SliderWidgetSuper::getInGeom ();
	// get the size of the area where I can draw
//		Size s = SliderWidgetSuper::getInSize ();
	double x1 = _marginWidth;
	double x2 = g.width () - _marginWidth;
	double y1 = _marginHeight;
	double y2 = g.height () - _marginHeight;
	// draw the top of the widget (with a hole)
	// --------------------------------------------------
	getBgColor () -> glColor ();
	glBegin (GL_TRIANGLE_STRIP);
	glVertex2d (0,0); // 0
	glVertex2d (x1, y1); // 1
	glVertex2d (g.width (), 0); // 2
	glVertex2d (x2, y1); // 3
	glVertex2d (g.width (), g.height ()); // 4
	glVertex2d (x2, y2); // 5
	glVertex2d (0,g.height()); // 6
	glVertex2d (x1,y2); // 7
	glVertex2d (0, 0); // 8
	glVertex2d (x1,y1); // 9
	glEnd ();
	// draw a groove
	// --------------------------------------------------
	double & gd = _grooveDepth;
	if (gd > 0)
	{
		Color::Pointer cols [3] = {
			getBgColor (),
			SliderWidgetSuper::_col3,
			SliderWidgetSuper::_col2};
		double weights [12] = {
			0.5, 0.5, 0.0,
			0.5, 0.0, 0.5,
			0.5, 0.5, 0.0,
			0.5, 0.0, 0.5
		};
		glBeveledBox (x1,y1,x2,y2,
			      gd, gd, gd, gd,
			      cols,
			      weights);
	}
	
	// draw the tick marks
	// --------------------------------------------------
	Color::Pointer _col4 = Color::Black;
	long _tickMarks = 10;
	long nt = _tickMarks;
	if (nt > 0)
	{
		glViewport (GLint (g.minx () + _marginWidth + gd),
			    GLint (g.miny ()),
			    GLint (g.width ()
				   - 2 * _marginWidth - 2 * gd),
			    GLint (_marginHeight));
		glMatrixMode (GL_PROJECTION);
		glLoadIdentity ();
		gluOrtho2D (0, 1.01, 0, 1);
		_col4-> glColor ();
		glBegin (GL_LINES);
		for (long i = 0 ; i < nt ; i ++)
		{
			glVertex2d (double (i)/(nt-1),1);
			glVertex2d (double (i)/(nt-1),0.5);
		}
		glEnd ();
	}
	
	
	// draw the indicator
	// --------------------------------------------------
	glViewport (GLint (g.minx () + _marginWidth + gd),
		    GLint (g.miny () + _marginHeight + gd),
		    GLint (g.width () - 2 * _marginWidth - 2 * gd),
		    GLint (g.height () - 2 * _marginHeight - 2 * gd));
	glMatrixMode (GL_PROJECTION);
	glLoadIdentity ();
	gluOrtho2D (0, 1, 0, 1);
	// what is current value normalized:
	double xpos = (_currVal - _minVal) / (_maxVal - _minVal);
	// background
	_col3-> glColor ();
	glBox (xpos,0,1,1);
	// foreground
	glBegin (GL_QUADS);
	_col1 -> glColor ();
	glVertex2d (0, 1);
	glVertex2d (0, 0);
//	_col2 -> glColor ();
	Color cend = (1-xpos) * (*_col1) + xpos * (*_col2);
	cend.glColor();
	glVertex2d (xpos, 0);
	glVertex2d (xpos, 1);
	glEnd ();
}

template <class T>
bool SliderWidget <T> ::mouse (const MouseEvent & ev)
{
	// ask for the inside geometry
	Geometry inGeom = SliderWidgetSuper::getInGeom ();
	Geometry outGeom = getOutGeom ();
	if (! mouseDown)
	{
		// mouse is not down (at least not yet), so if the event
		// happened outside of the slider, don't handle it
		if (! ev.inBox (outGeom)) return false;

		// was the left button was pushed ?
		if (ev.reason == MouseEvent::Push)
		{
			// if so, try to capture the mouse
			if (! captureMouse ()) return false;
		
			// set the new value
			mouseDown = true;
			mouseButton = ev.button;

			// depending on which button was pushed, modify
			// the value
			if (ev.button == MouseEvent::ButtonLeft)
			{
				double r = (ev.x - inGeom.minx () -
					    _marginWidth) /
					(inGeom.width () - 2 * _marginWidth);
				if (r < 0) r = 0;
				if (r > 1) r = 1;
				T newVal = (1-r) * _minVal + r * _maxVal;
				if (newVal != _currVal)
				{
					_currVal = newVal;
					valueChangedSignal.emit (_currVal);
					postRedisplay ();
				}
			} else {
				lastx = ev.x;
				lasty = ev.y;
			}
			return true;
		}
		else
		{
			// slider doesn't care about other events
			return true;
		}
	}

	// mouse is down
	// --------------------------------------------------
	if (ev.reason == MouseEvent::Release &&
	    ev.button == mouseButton)
	{
		releaseMouse ();
		mouseDown = false;
		return true;
	}
	
	if (ev.reason == MouseEvent::Motion &&
	    mouseButton == MouseEvent::ButtonLeft)
	{
		// get the relative position
		double r = (ev.x - inGeom.minx () - _marginWidth) /
			(inGeom.width () - 2 * _marginWidth);
		if (r < 0) r = 0;
		if (r > 1) r = 1;
		T newVal = (1-r) * _minVal + r * _maxVal;
		if (newVal != _currVal)
		{
			_currVal = newVal;
			valueChangedSignal.emit (_currVal);
			postRedisplay ();
		}
		return true;
	}
	if (ev.reason == MouseEvent::Motion &&
	    mouseButton == MouseEvent::ButtonRight)
	{
		// get the difference in positions
		double dx = ev.x - lastx;
		double dy = std::fabs (ev.y - outGeom.midy ());
		lastx = ev.x; lasty = ev.y;
		// calculate how much the value should
		// change, directly proportional to dx and
		// inversely proportional to dy
		double weight = (1-tanh(dy/30-2))/2;
		T valDiff = dx * weight * (_maxVal - _minVal) /
			(inGeom.width () - 2 * _marginWidth);
		T newVal = _currVal + valDiff;
		newVal = std::max (_minVal, newVal);
		newVal = std::min (_maxVal, newVal);
		if (newVal != _currVal)
		{
			_currVal = newVal;
			valueChangedSignal.emit (_currVal);
			postRedisplay ();
		}
	}
	
	// otherwise ignore the event
	return true;
}

OGLUI_END_NAMESPACE

#endif

