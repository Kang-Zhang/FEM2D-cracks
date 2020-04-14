#include "UIWidget.h"

#include <OGLUI.h>

namespace OGLUI {
#ifdef INDENTATION_HACK
} // to fool emacs into not indenting the rest of the file
#endif

static UI * gui = NULL;

UI::UI (int & argc, char ** & argv) :
	Widget (NULL, "ui"),
	mouseCaptor (NULL)
{
	// make sure we are the only interface
	if( gui != NULL)
	{
		std::cerr << "You can only create one UI.\n";
		std::exit (-1);
	}
	gui = this;
	glutInit (& argc, argv);
}

void UI::dispFunc ()
{
	gui-> redraw ();
}

void UI::reshapeFunc (int w, int h)
{
	gui-> reshape (w,h);
}

void UI::mouseFunc (int button, int state, int x, int y)
{
// 	cerr << "UI::mouseFunc ("
// 	     << button << ","
// 	     << state << ","
// 	     << x << ","
// 	     << y << ")\n";
	MouseEvent ev;
	ev.x = x; ev.y = glutGet (GLUT_WINDOW_HEIGHT) - y;
	MouseEvent::Button buttons [10] = {
		MouseEvent::Button1,
		MouseEvent::Button2,
		MouseEvent::Button3,
		MouseEvent::Button4,
		MouseEvent::Button5,
		MouseEvent::Button6,
		MouseEvent::Button7,
		MouseEvent::Button8,
		MouseEvent::Button9,
	};
	if (button >= 0 && button < 10)
		ev.button = buttons [button];
	else
		ev.button = MouseEvent::ButtonUnknown;
	if (state == GLUT_DOWN)
		ev.reason = MouseEvent::Push;
	else
		ev.reason = MouseEvent::Release;
	gui-> mouse (ev);
}

void UI::mouseMotionFunc (int x, int y)
{
// 	cerr << "UI::mouseMotionFunc ("
// 	     << x << ","
// 	     << y << ")\n";
	MouseEvent ev;
	ev.x = x; ev.y = glutGet (GLUT_WINDOW_HEIGHT) - y;
	ev.reason = MouseEvent::Motion;
	gui-> mouse (ev);
}

void UI::redraw ()
{
	// find out which window widget does this redraw belong to, and
	// call its _redraw()
	int id = glutGetWindow();
	for (WidgetList::iterator w = children.begin();
	     w != children.end();
	     w ++)
	{
		WindowWidget * ww = dynamic_cast<WindowWidget *> (*w);
		if (ww == NULL)
		{
			std::cerr << "OGLUI::Error: UI can only have "
				"Windows for children\n";
			exit (-1);
		}
		if (ww-> _glutWinId == id)
		{
			ww-> _redraw();
			return;
		}
	}
	std::cerr << "OGLUI::UI::display event for unknown widget\n";
	std::exit (-1);
}

void UI::reshape (int width, int height)
{
	// find out which window widget does this reshape belong to, and
	// call its reshape()
	int id = glutGetWindow();
	for (WidgetList::iterator i = children.begin();
	     i != children.end();
	     i ++)
	{
		WindowWidget * ww = dynamic_cast<WindowWidget *> (* i);
		if (ww == NULL)
		{
			std::cerr << "OGLUI:ERROR - UI can only have "
				"Windows for children\n";
			exit (-1);
		}
		if (ww-> _glutWinId == id)
		{
			ww-> reshape (width,height);
			return;
		}
	}
	std::cerr << "OGLUI::UI::reshape event for unknown widget\n";
	std::exit (-1);
}

void
UI::run ()
{
	std::cerr << "UI::_init()...\n";
	// call everybody's init() functions
	_init();
	std::cerr << "UI::_init() done\n";

	// now enter the main loop
	glutMainLoop ();
}

/// handle mouse events
bool UI::mouse (const MouseEvent & ev)
{
	// if the mouse is captured, send the event straight to that widget
	if (mouseCaptor)
		return mouseCaptor-> mouse (ev);
	// ask our children to process the event
	return propagateMouseEvent (ev);
}

/// widget can call this if it needs to have exclusive access to
/// the mouse (the returned value indicates success / failure)
bool UI::captureMouse (Widget * w)
{
	// we cannot capture the mouse for ourselves!!!
	if (w == NULL)
	{
		OGLUI::error (
			this,
			"captureMouse(): cannot capture for myself.");
		return false;
	}
	else
	{
		// if already captured by someone else, report warning
		if (mouseCaptor != NULL)
		{
			OGLUI::error (
				this,
				"captureMouse(): mouse already captured by "
				+ mouseCaptor-> id());
			return false;
		}
		else
		{
			mouseCaptor = w;
			return true;
		}
	}
}

/// releases the mouse
void UI::releaseMouse (Widget * w)
{
	// if NULL is the argument, just release the mouse
	if (w == NULL)
	{
		mouseCaptor = NULL;
		return;
	}
	if (mouseCaptor == NULL)
	{
		OGLUI::error(
			this,
			"releaseMouse(): mouse was not captured, but "
			+ w-> id() + "wants to release it.");
		return;
	}
	// report mismatch
	if( w != mouseCaptor)
	{
		OGLUI::error(
			this,
			"releaseMouse(): mouse was captured by " +
			mouseCaptor-> id() + "but " + w-> id() +
			" released it.");
		return;
	}
	mouseCaptor = NULL;
	return;
}

#ifdef INDENTATION_HACK
{ // to fool emacs into not indenting the rest of the file
#endif
};

