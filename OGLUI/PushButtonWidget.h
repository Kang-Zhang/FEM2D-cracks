#ifndef __OGLUI_PUSH_BUTTON_H__
#define __OGLUI_PUSH_BUTTON_H__

#include <GL/gl.h>
#include <LabelWidget.h>
#include <Slot.h>

namespace OGLUI {
#ifdef INDENTATION_HACK
} // to fool emacs into not indenting the rest of the file
#endif

class PushButtonWidget : public LabelWidget {
public:
	class Listener {
	public:
		virtual void pushed (const PushButtonWidget & pb) = 0;
		virtual ~Listener () {;}
	};
private:
	/// list of listeners
	std::vector<Listener *> listeners;
	/// indicates whether the mouse is down or not
	bool mouseDown;
	/// if the mouse is down, this basically indicates whether the
	/// mouse is on top of the button
	bool depressed;
	/// called when the mouse is pushed
	void pushed ();
protected:
	/// this gets called when the button is pushed (can be overloaded)
	virtual void onPushed ();
	/// overloaded mouse event
	virtual bool mouse (const MouseEvent & ev);
	/// redraws the button
	virtual void redraw ();
public:
	/// the signal which gets emitted whenever the button is pushed
	Signal0 pushedSignal;
	/// the class name
	virtual const std::string className () const
		{
			return "PushButtonWidget";
		}
	/// constructor
	PushButtonWidget ( Widget * p_parent,
			   const std::string & p_name = "pushbutton"
		);
	// adds a listener
	void addListener (Listener * lis) {
		listeners.push_back (lis);
	}
};

#ifdef INDENTATION_HACK
{ // to fool emacs into not indenting the rest of the file
#endif
} // namespace OGLUI


#endif

