#include "ToggleButtonWidget.h"

#include "OGLUI.h"

namespace OGLUI {
#ifdef INDENTATION_HACK
} // to fool emacs into not indenting the rest of the file
#endif

/// constructor
ToggleButtonWidget::ToggleButtonWidget (
	Widget * p_parent,
	const std::string & p_name
	) :
//	PushButtonWidget (p_parent, p_name),
	ToggleButtonWidgetSuper (p_parent, p_name),
	_state (false),
	_toggleColor (Color::Yellow),
	_toggleSize (Size (6,10))
{
	_justHoriz = LEFT;
	_justOffset = Position (5+6+10, 0);
}

void ToggleButtonWidget:: onPushed ()
{
	// flip the state & post redisplay
	_state = ! _state;
	postRedisplay ();
	// call onValueChanged() // in case somebody overloaded it
	onValueChanged ();
	// call all listeners - if any
	for (std::vector<Listener *>::iterator l = listeners.begin ();
	     l != listeners.end ();
	     l ++)
	{
		(*l)-> stateChanged (* this);
	}
	// emit the value changed signal
	valueChangedSignal.emit (_state);
};

void
ToggleButtonWidget::redraw ()
{
	ToggleButtonWidgetSuper::redraw ();
	if (_state)
	{
		_toggleColor-> glColor ();
	}
	else
	{
		Color shadedColor = * getBgColor ();
		shadedColor.tint( -0.1);
		shadedColor.glColor ();
	}
	// ask our superclass what is the inside area
	Geometry g = ToggleButtonWidgetSuper::getInGeom ();
	ToggleButtonWidgetSuper::setupInOGL ();
	double cx = 5 + _toggleSize.width / 2;
	double cy = g.getSize ().height / 2;
	glBox (cx - _toggleSize.width / 2, cy - _toggleSize.height / 2,
	       cx + _toggleSize.width / 2, cy + _toggleSize.height / 2);
}

#ifdef INDENTATION_HACK
{ // to fool emacs into not indenting the rest of the file
#endif
};

