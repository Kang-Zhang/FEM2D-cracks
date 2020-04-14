#include <cstdlib>
#include <cstdio>
#include <memory>
#include <UIWidget.h>
#include <WindowWidget.h>
#include <BevelWidget.h>
#include <SmartLayoutManager.h>
#include <PushButtonWidget.h>

OGLUI_USE_NAMESPACE

PushButtonWidget * makeButton (Widget * parent, const std::string name)
{
	PushButtonWidget * pb = new PushButtonWidget (parent, name);
	pb -> setLabel (name);
	pb -> setBevelSize (2);
	pb -> setPrefSize (Size (100, 30));
	pb -> setBgColor (new Color(drand48(), drand48(), drand48 ()));
	return pb;
}
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
	SmartLayoutManager::Pointer lm = 
		win -> setLayoutManager (new SmartLayoutManager ());
	// First button - top and left are attached
	PushButtonWidget * pb1 = makeButton ( win, "Button 1");
	lm -> attach (SLM::TopEdge (pb1), SLM::Relative (90, 0));
	lm -> attach (SLM::BottomEdge (pb1), SLM::None ());
	lm -> attach (SLM::LeftEdge (pb1), SLM::Relative (10,0));
	lm -> attach (SLM::RightEdge (pb1), SLM::None ());
	// second button - top is attached to previous button's bottom
	//                 left is attached to relative position
	PushButtonWidget * pb2 = makeButton ( win, "Button2");
	lm -> attach (SLM::TopEdge (pb2), SLM::OtherWidget (pb1));
	lm -> attach (SLM::LeftEdge (pb2), SLM::Relative (5));
	// third button - top is attached to previous button's top
	//                 left is attached to previous button's right
	PushButtonWidget * pb3 = makeButton ( win, "Button3");
	lm -> attach (SLM::TopEdge (pb3), SLM::OppositeWidget (pb2));
	lm -> attach (SLM::LeftEdge (pb3), SLM::OtherWidget (pb2));
	// fourth button - top is attached to previous button's bottom
	//                 left is attached to previous button's right
	PushButtonWidget * pb4 = makeButton ( win, "Button4");
	lm -> attach (SLM::TopEdge (pb4), SLM::OtherWidget (pb3, 20));
	lm -> attach (SLM::LeftEdge (pb4), SLM::OtherWidget (pb3));
	lm -> attach (SLM::BottomEdge (pb4), SLM::Relative (0,10));
	// 
	PushButtonWidget * pb5 = makeButton ( win, "Button5");
	lm -> attach (SLM::BottomEdge (pb5), SLM::OtherWidget (pb4, 20));
	lm -> attach (SLM::LeftEdge (pb5), SLM::OtherWidget (pb4));
	
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
