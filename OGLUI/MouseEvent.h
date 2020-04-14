#ifndef __OGLUI_MOUSE_EVENT_H__
#define __OGLUI_MOUSE_EVENT_H__

namespace OGLUI {
#ifdef INDENTATION_HACK
} // to fool emacs into not indenting the rest of the file
#endif

class MouseEvent {
public:
	enum Button {
		Button1=0,Button2,Button3,Button4,Button5,
		Button6,Button7,Button8,Button9,
		ButtonLeft = Button1,
		ButtonMiddle = Button2,
		ButtonRight = Button3,
		ButtonUnknown
	};
	/// mouse coordinates
	int x, y;
	/// which button was pushed (if Pushed event)
	Button button;
	/// reason for generating this event
	enum Reason { Push, Release, Motion, Entry, Exit } reason;
	/// is the mouse in the given box?
	bool inBox (const Geometry & geom) const
		{
			return x >= geom.minx () && x <= geom.maxx () &&
				y >= geom.miny () && y <= geom.maxy ();
		}
};

std::ostream & operator<< (std::ostream & os, const MouseEvent & g);

#ifdef INDENTATION_HACK
{ // to fool emacs into not indenting the rest of the file
#endif
} // namespace OGLUI


#endif

