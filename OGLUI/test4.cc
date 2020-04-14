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
	double width = 1000;
	double height = 1000;
	long n_cols = 50;
	long n_rows = 20;
	double bwidth = width / n_cols;
	double bheight = height / n_rows;
	WindowWidget * win = new WindowWidget (ui);
	win -> setPrefSize (Size (width, height));
	win -> setBgColor (Color::DarkGreen);
	win -> setDisplayMode (GLUT_RGB | GLUT_DOUBLE);
	// create a layout manager
	SLM::Pointer lm = win -> setLayoutManager ( new SLM);
	// create a grid of buttons
	for (long r = 0 ; r < n_rows ; r ++)
	{
		for (long c = 0 ; c < n_cols ; c ++)
		{
			PushButtonWidget * pb = new PushButtonWidget (win);
			pb -> setLabel ("a");
			pb -> setBevelSize (1);
			pb -> setPrefSize (Size (bwidth, bheight));
			pb -> setBgColor (new Color (drand48(),
						     drand48(),
						     drand48()));
			lm -> attach (SLM::LeftEdge (pb)
				      , SLM::Absolute (c*bwidth));
			lm -> attach (SLM::BottomEdge (pb)
				      , SLM::Absolute (r*bheight));
//			lm -> setPosition (pb, Position (c*bwidth, r*bheight));
		}
	}
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
