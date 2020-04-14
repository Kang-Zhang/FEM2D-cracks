#ifndef __OGLUI_WINDOW_WIDGET_H__
#define __OGLUI_WINDOW_WIDGET_H__

#include <GL/glut.h>
#include "Widget.h"

namespace OGLUI {
#ifdef INDENTATION_HACK
} // to fool emacs into not indenting the rest of the file
#endif

class WindowWidget : public Widget {
private:
//	Geometry _geom;
protected:
	friend class UI;
	unsigned int _glutDispMode;
	int _glutWinId;
	virtual void reshape (int width, int height);
	virtual const std::string className () const { return "WindowWidget"; }
public:
	/// constructor
	WindowWidget ( Widget * p_parent,
		       const std::string & p_name = "window");
	void setDisplayMode (unsigned int pmode)
		{
			_glutDispMode = pmode;
		}
	virtual void postRedisplay ();
	virtual bool init ();
	virtual void redraw ();
	virtual void redraw2 ();
	/// overloaded geometry methods
	virtual Geometry getInGeom () const;
	virtual Geometry getInGeomClipped () const;
	virtual Geometry getOutGeom () const;
	virtual Geometry getOutGeomClipped () const;
	
};

typedef std::vector<WindowWidget *> WindowWidgetList;

#ifdef INDENTATION_HACK
{ // to fool emacs into not indenting the rest of the file
#endif
} // namespace OGLUI


#endif

