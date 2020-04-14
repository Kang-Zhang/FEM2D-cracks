#ifndef __OGLUI_WIDGET_H__
#define __OGLUI_WIDGET_H__

#include <vector>
#include <memory>
#include <Geometry.h>
#include <MouseEvent.h>
#include <Color.h>
#include <SmartPointer.h>
#include <Theme.h>
#include <Size.h>
#include <GL/glu.h>
#include <OGLUInamespace.h>
#include <LayoutManager.h>

OGLUI_BEGIN_NAMESPACE

class Widget;
typedef std::vector<Widget *> WidgetList;

class Widget
{
public:
	typedef SmartPointer <Widget> Pointer;
	/// returns name of the class
	virtual const std::string className () const;
private:
	/// background color
	Color::Pointer _bgColor;
	/// stores the name of the instance (different from class name)
	const std::string _instanceName;
	/// border color
	Color::Pointer _borderColor;
	/// min, max & preferred size
	Size _minSize;
	Size _maxSize;
	Size _prefSize;
	/// the layout manager
	LayoutManager::Pointer _layMan;
protected:
	/// border thickness
	double _borderThickness;
	/// cached geometry
	Geometry _geom;
protected:
	/// points to the parent
	Widget * parent;
	/// whether the widget is visible or not
	bool _invisible;
	/// contains a list of widgets
	WidgetList children;
	/// recursively calls all init() in the widget tree
	bool _init ();
	/// recursively calls all redraw() and redraw2() in the widget tree
	void _redraw ();
	/// recursively propagates the event to the children widgets
	bool propagateMouseEvent (const MouseEvent & e);
	/// widget can call this if it needs to have exclusive access to
	/// the mouse (the returned value indicates success / failure)
	virtual bool captureMouse (Widget * w = NULL);
	/// releases the mouse
	virtual void releaseMouse (Widget * w = NULL);
	/// convenience function - to redraw only border
	virtual void redrawBorder ();
public:
	/// sets the layout manager
	LayoutManager::Pointer setLayoutManager (LayoutManager::Pointer lm);
	/// returns the layout manager
	LayoutManager::Pointer getLayoutManager ();
	/// returns the geometry of the child - through the layout manager
	Geometry getChildGeometry (const Widget * w);
	/// sets the min. size of the widget
	virtual const Size & getMinSize () const;
	/// gets the min. size of the widget
	virtual void setMinSize (const Size & size);
	/// sets the max. size of the widget
	virtual const Size & getMaxSize () const;
	/// gets the max. size of the widget
	virtual void setMaxSize (const Size & size);
	/// sets the preferred size of the widget
	virtual const Size & getPrefSize () const;
	/// gets the preferred size of the widget
	virtual void setPrefSize (const Size & size);
	/// when geometry is changed, this should be called
	void onGeometryChanged (const Geometry & g);
	/// returns the inside geometry
	virtual Geometry getInGeom () const;
	/// return inside geometry, clipped
	virtual Geometry getInGeomClipped () const;
	/// returns the cached outside geometry
	virtual Geometry getOutGeom () const;
	/// returns the outside geometry, clipped
	virtual Geometry getOutGeomClipped () const;
	/// sets up opengl state to be able to draw directly inside the widget
	virtual void setupInOGL () const;
	/// as above but with clipping
	virtual void setupInOGLclipped () const;
	/// setup OGL for the entire area, with appropriate clipping
	virtual void setupOutOGLclipped () const;
	/// get the background color
	const Color::Pointer getBgColor () const;
	/// sets the background color
	virtual void setBgColor (const Color::Pointer col);
	/// name of the instance
	const std::string & instanceName() const;
	/// full id of the widget: <classname>(<instancename>,<ptr>
	std::string id ();
	/// gets invisibility status
	bool invisible () const;
	/// sets invisibility status
	void invisible (bool i);
//	/// sets the geometry
//	void setGeometry (const Geometry & p_geom);
	/// sets the border thickness
	void borderThickness (double t);
	/// gets the border thickness
	double borderThickness () const;
	/// sets the border color
	void setBorderColor (const Color::Pointer c);
	/// adds a child to the widget
	virtual void addChild (Widget * w);
	/// called when the widget is first initialized (before its
	/// childre)
	virtual bool init ();
	/// destructor
	virtual ~Widget ();
	/// called when the widget needs to be redrawn (before the children)
	virtual void redraw ();
	/// called when the widget needs to be redrawn (after its children)
	virtual void redraw2 ();
	/// called when the widget gets the mouse
	virtual bool mouse (const MouseEvent & e);
	/// calling this schedules a redraw event for this widget
	virtual void postRedisplay ();
	/// constructor
	Widget (Widget * p_parent, const std::string & p_name = "widget");
	/// debuging method - to list the children
	void printTree ( const std::string prefix = "");
	// set the transparency of the widget
	virtual void transparency (double t);
};

OGLUI_END_NAMESPACE

#endif

