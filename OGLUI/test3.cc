#include <cstdlib>
#include <cstdio>
#include <memory>
#include <UIWidget.h>
#include <WindowWidget.h>
#include <BevelWidget.h>
#include <SmartLayoutManager.h>
#include <PushButtonWidget.h>

using namespace OGLUI;

static BevelWidget * createPanel (Widget * parent)
{
	// create the panel
	BevelWidget * pan = new BevelWidget (parent, "panel");
//	pan -> setBgColor (Color::Red);
	pan -> setBgColor (new Color(1,0,drand48()*0.7));
	pan -> setBevelSize (8);
	pan -> setPrefSize ( Size (900, 900));
	// create first button
	PushButtonWidget * pb1 = new PushButtonWidget (pan);
	pb1 -> setBgColor (Color::Yellow);
	pb1 -> setLabel ("Yellow");
	pb1 -> setPrefSize (Size (100, 30));
	// create second button
	PushButtonWidget * pb2 = new PushButtonWidget (pan);
	pb2 -> setBgColor (Color::Green);
	pb2 -> setLabel ("Green");
	pb2 -> setPrefSize (Size (100, 30));
	// set the layouts
	SLM::Pointer lm = pan -> setLayoutManager (new SLM ());
	lm -> attach (SLM::LeftEdge (pb1), SLM::Absolute (10));
	lm -> attach (SLM::BottomEdge (pb1), SLM::Absolute (5));
	lm -> attach (SLM::LeftEdge (pb2), SLM::Absolute (120));
	lm -> attach (SLM::BottomEdge (pb2), SLM::Absolute (5));

	return pan;
}

static int real_main (int & argc, char ** & argv)
{
	// initialize the user interface library
	UI * ui = new UI (argc, argv);
	// create a new window
	WindowWidget * win = new WindowWidget (ui);
	win -> setPrefSize (Size (1000, 1000));
	win -> setBgColor (Color::DarkGreen);
	win -> setDisplayMode (GLUT_RGB | GLUT_DOUBLE);
	// create a panel (with 2 push-buttons)
	Widget * parent = win;
	for (long i = 0 ; i < 20 ; i ++)
	{
		Widget * panel = createPanel (parent);
		// set the layout of the panel
		SLM::Pointer lm = parent-> getLayoutManager ();
		if (lm == NULL)
			lm = parent -> setLayoutManager (new SLM);
		lm -> attach (SLM::LeftEdge (panel), SLM::Absolute (10));
		lm -> attach (SLM::BottomEdge (panel), SLM::Absolute (10));
		parent = panel;
	}
	std::cerr << "Timing getOutGeom () on the last widget\n";
	parent-> getOutGeom ();
	std::cerr << "...timing done\n";
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
	catch (const std::string & str)
	{
		std::cerr << "main():Error:" << str << "\n";
	}
	catch (const char * str)
	{
		std::cerr << "main():Error:" << str << "\n";
	}
	catch (...)
	{
		std::cerr << "main():Error-Unknown\n";
	}
}
