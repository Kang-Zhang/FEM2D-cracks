#ifndef __OGLUI_UI_H__
#define __OGLUI_UI_H__

#include "WindowWidget.h"

namespace OGLUI {
#ifdef INDENTATION_HACK
} // to fool emacs into not indenting the rest of the file
#endif

class UI : public Widget {
protected:
	friend class WindowWidget;
	static void dispFunc ();
	static void reshapeFunc (int w, int h);
	static void mouseFunc (int button, int state, int x, int y);
	static void mouseMotionFunc (int x, int y);
	void redraw ();
	void reshape (int w, int h);
	virtual const std::string className () const { return "UIWidget"; }
	/// captures the mouse for some widget
	virtual bool captureMouse (Widget * w = NULL);
	/// releases the mouse for some widget
	virtual void releaseMouse (Widget * w = NULL);
	/// this is where we hold the pointer to the widget that has
	/// captured the mouse
	Widget * mouseCaptor;
	/// handle mouse events
	virtual bool mouse (const MouseEvent & e);
public:
	// constructor
	UI (int & argc, char ** & argv) ;
	virtual void theme (SmartPointer <Theme> th) {;}
	void run () ;
};

#ifdef INDENTATION_HACK
{ // to fool emacs into not indenting the rest of the file
#endif
} // namespace OGLUI


#endif

