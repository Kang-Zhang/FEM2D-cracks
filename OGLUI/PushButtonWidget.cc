#include "PushButtonWidget.h"

namespace OGLUI {
#ifdef INDENTATION_HACK
} // to fool emacs into not indenting the rest of the file
#endif

// constructor
PushButtonWidget::PushButtonWidget (
	Widget * p_parent,
	const std::string & p_name
	) :
	LabelWidget (p_parent, p_name),
	mouseDown (false),
	depressed (false)
{
	// empty constructor
}

// this gets called whenever the button gets pushed
void PushButtonWidget::pushed ()
{
	// call onPushed() in case it has been overloaded
	onPushed ();
	// call all listeners
	for (std::vector<Listener *>::iterator l = listeners.begin ();
	     l != listeners.end ();
	     l ++)
	{
		(*l)-> pushed (* this);
	}
	// emit a signal
	pushedSignal.emit ();
}

// this gets called whenever the button gets pushed
void PushButtonWidget::onPushed ()
{
}

bool PushButtonWidget::mouse (const MouseEvent & ev)
{
	if (mouseDown)
	{
		if (ev.reason == MouseEvent::Release &&
		    ev.button == MouseEvent::ButtonLeft)
		{
			if (depressed)
			{
				depressed = false;
				postRedisplay ();
				pushed ();
			}
			mouseDown = false;
			releaseMouse ();
			return true;
		}
		if (ev.reason == MouseEvent::Motion)
		{
			if (ev.inBox (getOutGeom ()))
			{
				if (! depressed)
				{
					depressed = true;
					postRedisplay ();
				}
				return true;
			}
			else
			{
				if (depressed)
				{
					depressed = false;
					postRedisplay ();
				}
				return true;
			}
		}
		// otherwise ignore events
		return true;
	}
	// deal with this event if it is for us:
	if (ev.inBox (getOutGeom ()))
	{
		if (ev.reason == MouseEvent::Push &&
		    ev.button == MouseEvent::ButtonLeft)
		{
			if (captureMouse ())
			{
				mouseDown = true;
				depressed = true;
				postRedisplay ();
				return true;
			}
		}
		// otherwise ignore events
		return false;
	}
	// it was not a button or it was outside of the button,
	// so propagate the event to our kids
	return propagateMouseEvent (ev);
}

void
PushButtonWidget::redraw ()
{
	if (depressed)
	{
		// flip the color weights
		std::swap (_lightWeightsTop, _lightWeightsBottom);
		std::swap (_lightWeightsLeft, _lightWeightsRight);
	}
	LabelWidget::redraw ();
	if (depressed)
	{
		// flip the color weights
		std::swap (_lightWeightsTop, _lightWeightsBottom);
		std::swap (_lightWeightsLeft, _lightWeightsRight);
	}
}

#ifdef INDENTATION_HACK
{ // to fool emacs into not indenting the rest of the file
#endif
};

