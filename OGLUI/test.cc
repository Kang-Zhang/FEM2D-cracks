#include <cstdlib>
#include <cstdio>
#include <PushButtonWidget.h>
#include <ToggleButtonWidget.h>
#include <UIWidget.h>
#include <SliderWidget.h>
#include <SmartLayoutManager.h>

// using OGLUI::Geometry;
// using OGLUI::Widget;

using namespace OGLUI;

//using OGLUI::Signal;

class MyWidget : public PushButtonWidget
{
public:
	MyWidget (
		Widget * p_parent,
		const std::string & p_name = "mywidget"
		) :
		PushButtonWidget( p_parent, p_name)
		{
			setLabel ("My Widget");
		}
	virtual void redraw () {
		PushButtonWidget::redraw ();
		glBegin (GL_LINES);
		glColor4f (1, 0, 0, 0.5);
		glVertex2d (0,0);
		glColor4f (1, 1, 0, 0.5);
		Geometry g = getOutGeom ();
		glVertex2d (g.width (),g.height ());
		glColor4f (0, 0, 1, 0.5);
		glVertex2d (g.width (),0);
		glColor4f (0, 1, 0, 0.5);
		glVertex2d (0,g.height ());
		glEnd ();
	}
};

class MyPBListener : public PushButtonWidget::Listener
{
private:
	const std::string _str;
public:
	MyPBListener (const std::string & str) : _str(str) {;}
	void pushed ( const PushButtonWidget & pb) {
		std::cerr << "pressed button - name = " << _str << "\n";
	}
};

/*
class MyToggleListener : public ToggleButtonWidget::Listener
{
private:
	const std::string _str;
public:
	MyToggleListener (const std::string & str) : _str(str) {;}
	void stateChanged ( const ToggleButtonWidget & tb) {
		cerr << "Toggle(" << _str << ") changed state: " <<
			tb.state() << "\n";
	}
//	virtual ~MyToggleListener () {;}
};
*/

static void toggleSlot (bool state)
{
	std::cerr << "Toggle state = " << state << "\n";
}

class MyTheme : public Theme
{
public:
//	MyTheme () :

	MyTheme (const std::string & id) :
		Theme (id)
		{
			std::cerr << "MyTheme()\n";
			set ("*.foregroundColor", Color::Red);
			set ("*.backgroundColor", Color::Pink);
		}
	virtual ~MyTheme ()
		{
			std::cerr << "~MyTheme()\n";
		}
};

void sliderCB (double d)
{
	std::cerr << "SliderCB (" << d << ")\n";
}

int real_main (int & argc, char ** & argv)
{
	UI * ui = new UI (argc, argv);
//	ui -> theme (new MyTheme ("MyTheme"));
	WindowWidget * win = new WindowWidget (ui);
	win-> setPrefSize (Size (500, 500));
	win-> setBgColor (Color::DarkGreen);
	win-> setDisplayMode (GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	SLM::Pointer lm = win -> setLayoutManager (new SLM);

	LabelWidget * lab = new LabelWidget (win);
	lm -> setGeometry (lab, Geometry (Position (10, 20), Size (200, 50)));
	lab-> setLabel ("I am pretty!");
	lab-> setBgColor (new Color (1,0,1));
//	win-> add (lab);
	LabelWidget * lab2 = new LabelWidget (win);
	lm -> setGeometry (lab2, Geometry (Position (200, 100), Size (80, 40)));
	lab2-> setLabel ("Semi visible");
	lab2-> setBgColor (new Color (0.6,0.6,0,1));
	lab2-> borderThickness (5);
	lab2-> setBevelSize (2);
	lab2-> lightWeightsBottom (BevelWidget::LightWeights(1,0,0));
	lab2-> lightWeightsRight (BevelWidget::LightWeights(1,0,0));
	lab2-> setBorderColor (new Color (1,1,1,0.7));
	lab2-> setFontColor (new Color (1,0,0,0.5));
	lab2-> setFont (BitmapFont::TimesRoman24);
//	win-> add (lab2);
	MyWidget * myw = new MyWidget (win);
	lm -> setGeometry (myw, Geometry (Position (50, 10), Size (30, 80)));
	myw-> setBgColor (new Color (0.7,0.7,0.0,0.5));
	myw-> addListener (new MyPBListener ("myPB2"));
//	win-> add (myw);
	PushButtonWidget * pb = new PushButtonWidget (win);
//	pb-> setFont (BitmapFont::TimesRoman10);
	pb-> setFont (BitmapFont::TimesRoman24);
	pb-> setLabel ("ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz");
	pb-> setBgColor (Color::Ivory);
	lm -> setGeometry (pb, Geometry (Position (5,250), Size (490,30)));
	pb-> addListener (new MyPBListener ("PB1"));
//	win-> add (pb);
	BevelWidget * bw = new BevelWidget (win);
	lm -> setGeometry (bw, Geometry (Position (5,200), Size (100,30)));
	bw-> setBgColor (new Color (0.5,0.6,0.3));
	bw-> setBevelSize (1);
//	bw-> setColor (Color::Green3);
//	win-> add (bw);

	SliderWidget<double> * slider = new SliderWidget<double> (win);
	lm -> setGeometry (slider, Geometry (Position (5,150), Size (200,25)));
	slider-> setBgColor (Color::PeachPuff);
	slider-> valueChangedSignal.connect (slot(sliderCB));

	Color::Pointer bg = new Color (0,0,0,0.5);
	Color::Pointer fg = new Color (1,1,1,0.5);
	ToggleButtonWidget * tog1 = new ToggleButtonWidget (win);
	tog1-> setBgColor (bg);
	tog1-> setFontColor (fg);
	lm -> setGeometry (tog1, Geometry (Position (105,5), Size (300,80)));
	tog1-> setFont (BitmapFont::Helvetica12);
	tog1-> setLabel ("I am a toggle");
//	tog1-> addListener (new MyToggleListener ("I am a toggle."));
	tog1-> valueChangedSignal.connect (slot (& ::toggleSlot));
	tog1-> setToggleSize ( Size (6,6));
//	win-> add (tog1);

//	ui-> add (win);

	ui-> printTree ( "#.#.#.#.#        ");
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
		std::cerr << "Exception occured!\n";
	}
}
