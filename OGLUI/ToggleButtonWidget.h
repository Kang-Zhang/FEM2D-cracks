#ifndef __OGLUI_TOGGLE_BUTTON_H__
#define __OGLUI_TOGGLE_BUTTON_H__

#include <GL/gl.h>
#include <PushButtonWidget.h>
#include <Slot.h>

namespace OGLUI {
#ifdef INDENTATION_HACK
} // to fool emacs into not indenting the rest of the file
#endif

// define the superclass
typedef PushButtonWidget ToggleButtonWidgetSuper;

//class ToggleButtonWidget : public PushButtonWidget {
class ToggleButtonWidget : public ToggleButtonWidgetSuper {
public:
	class Listener {
	public:
		virtual void stateChanged (const ToggleButtonWidget & pb) = 0;
		virtual ~Listener () {;}
	};
private:
	/// represents the current state of the button (on/off)
	bool _state;
	/// list of listeners
	std::vector<Listener *> listeners;
	/// the color of the toggle when it is on
	Color::Pointer _toggleColor;
	/// the size of the toggle in pixels
	Size _toggleSize;
	/// overloaded pushbutton method: flip state, post redisplay and
	/// call all listeners
	virtual void onPushed ();
protected:
	/// redraws the widget
	virtual void redraw ();
	/// this function gets called on value changed - so it can be
	/// overloaded
	virtual void onValueChanged () {;}
public:
	/// the class name
	virtual const std::string className () const
		{
			return "ToggleButtonWidget";
		}
	/// the signal which toggle emits when the value is changed
	Signal1 <bool> valueChangedSignal;
	/// constructor
	ToggleButtonWidget (
		Widget * p_parent,
		const std::string & p_name = "togglebutton");
	virtual bool getState () const
		{
			return _state;
		}
	virtual void setState (bool s)
		{
			if (_state != s)
			{
				postRedisplay ();
				_state = s;
			}
		}
	/// adds a listener to the list of listeners
	virtual void addListener (Listener * lis)
		{
			listeners.push_back (lis);
		}
	/// changes the toggle size
	void setToggleSize (const Size & s)
		{
			_toggleSize = s;
			_justOffset = Position (5+10+s.width,0);
			postRedisplay();
		}
	/// reports the toggle size
	const Size & getToggleSize () const
		{
			return _toggleSize;
		}
	/// interface to setting / getting toggle color
	Color::Pointer getToggleColor () const
		{ 
			return _toggleColor;
		}
	void setToggleColor (Color::Pointer c)
		{
			_toggleColor = c;
			postRedisplay ();
		}
};

#ifdef INDENTATION_HACK
{ // to fool emacs into not indenting the rest of the file
#endif
} // namespace OGLUI


#endif

