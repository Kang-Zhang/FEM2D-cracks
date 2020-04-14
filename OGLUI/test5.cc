#include <cstdlib>
#include <cstdio>
#include <memory>
#include <UIWidget.h>
#include <WindowWidget.h>
#include <BevelWidget.h>
#include <SmartLayoutManager.h>
#include <PushButtonWidget.h>

OGLUI_USE_NAMESPACE

static int real_main (int & argc, char ** & argv)
{
	// initialize the user interface library
	UI * ui = new UI (argc, argv);
	// create a new window
	WindowWidget * win = new WindowWidget (ui);
	win -> setPrefSize (Size (500, 500));
	win -> setBgColor (Color::DarkGreen);
	win -> setDisplayMode (GLUT_RGB | GLUT_DOUBLE);
	// create a layout manager
	SLM::Pointer lm = win -> setLayoutManager (new SLM ());
	// create a push button
	PushButtonWidget * pb = new PushButtonWidget (win);
	pb -> setLabel ("nicely centered");
	pb -> setBevelSize (1);
	pb -> setPrefSize (Size (100, 30));
	pb -> setBgColor (Color::Red);
	lm -> attach (SLM::LeftEdge (pb), SLM::Relative (33));
	lm -> attach (SLM::RightEdge (pb), SLM::Relative (67));
	lm -> attach (SLM::BottomEdge (pb), SLM::Relative (33));
	lm -> attach (SLM::TopEdge (pb), SLM::Relative (67));
	// run the application
	ui-> run ();

	return 0;
}

int main (int argc, char ** argv)
{
	try
	{
		real_main (argc, argv);
	}
	catch (...)
	{
	}
}
