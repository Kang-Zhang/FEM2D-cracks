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
	BevelWidget * pan = new BevelWidget (parent);
	pan -> setBgColor (Color::Red);
	pan -> setBevelSize (10);
	pan -> setPrefSize ( Size (70, 300));
	// create first button
	PushButtonWidget * pb1 = new PushButtonWidget (pan);
	pb1 -> setBgColor (Color::Yellow);
	pb1 -> setLabel ("PB1hahahahahaha");
	pb1 -> setPrefSize ( Size (100, 100));
	// create second button
	PushButtonWidget * pb2 = new PushButtonWidget (pan);
	pb2 -> setBgColor (Color::Green);
	pb2 -> setLabel ("PB2");
	// set the layouts
	SmartLayoutManager::Pointer lm = pan -> setLayoutManager (
		new SmartLayoutManager);
//	lm -> attach (SLM::LeftEdge (pb1), SLM::Absolute (10));
	lm -> attach (SLM::LeftEdge (pb1), SLM::Relative (10));
	lm -> attach (SLM::BottomEdge (pb1), SLM::Absolute (10));
	lm -> attach (SLM::LeftEdge (pb2), SLM::Absolute (15));
	lm -> attach (SLM::BottomEdge (pb2), SLM::Absolute (50));
	
	return pan;
}

static int real_main (int & argc, char ** & argv)
{
	// initialize the user interface library
	UI * ui = new UI (argc, argv);
	// create a layout manager
	SmartLayoutManager * lm = new SmartLayoutManager;
	// create a new window
	WindowWidget * win = new WindowWidget (ui);
	win -> setPrefSize (Size (500, 500));
	win -> setBgColor (Color::DarkGreen);
	win -> setDisplayMode (GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	win -> setLayoutManager (lm);
	// create a panel (with 2 push-buttons)
	Widget * w1 = createPanel (win);
	lm -> attach (SLM::LeftEdge (w1), SLM::Absolute (20));
	lm -> attach (SLM::BottomEdge (w1), SLM::Absolute (50));
	// another panel (with 2 push-buttons)
	Widget * w2 = createPanel (win);
	lm -> attach (SLM::LeftEdge (w2), SLM::Absolute (300));
	lm -> attach (SLM::BottomEdge (w2), SLM::Absolute (10));
	// run the application
	ui-> run ();

	return 0;
}

int main (int argc, char ** argv)
{
/*
	Geometry g1 (Position (2,1), Position (3,5));
	Geometry g2 (Position (1,3), Position (6,4));
	Geometry g3 = intersection (g1, g2);
	std::cerr << "result = " << g3.minx () << "," << g3.miny () << " "
		  << g3.maxx () << "," << g3.maxy () << "\n";
	return 0;
*/
	try
	{
		real_main (argc, argv);
	}
	catch (...)
	{
	}
}
