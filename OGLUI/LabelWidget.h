#ifndef __OGLUI_LABEL_H__
#define __OGLUI_LABEL_H__

#include "BevelWidget.h"
#include "Font.h"
#include "Position.h"

namespace OGLUI {
#ifdef INDENTATION_HACK
} // to fool emacs into not indenting the rest of the file
#endif

class LabelWidget : public BevelWidget {
protected:
	std::string label;
	Font::Pointer font;
	Color::Pointer fontColor;
	enum JustificationHoriz { LEFT, RIGHT, CENTEREDH } _justHoriz;
	enum JustificationVert { TOP, BOTTOM, CENTEREDV } _justVert;
	Position _justOffset;
	virtual const std::string className () const { return "LabelWidget"; }
public:
	void setLabel (const std::string & str) { label = str; }
	void setFontColor (const Color::Pointer c) { fontColor = c; }
	void setFont (Font::Pointer f) { font = f; }
	virtual void redraw ();
	// constructor
	LabelWidget (
		Widget * p_parent,
		const std::string & p_name = "label"
		) :
		BevelWidget (p_parent, p_name),
		label ("Label"),
		font (BitmapFont::TimesRoman24),
		fontColor (new Color (0,0,0,1)),
		_justHoriz (CENTEREDH),
		_justVert (CENTEREDV),
		_justOffset (0,0)
		{
			setBevelSize (2);
		}
};

#ifdef INDENTATION_HACK
{ // to fool emacs into not indenting the rest of the file
#endif
} // namespace OGLUI


#endif

