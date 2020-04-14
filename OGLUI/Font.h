#ifndef __OGLUI_FONT_H__
#define __OGLUI_FONT_H__

#include <GL/glut.h>
#include "Theme.h"
#include "Size.h"
#include "SmartPointer.h"

namespace OGLUI {
#ifdef INDENTATION_HACK
} // to fool emacs into not indenting the rest of the file
#endif

class Font :
	public Theme::Themable
{
// interface class
public:
	typedef SmartPointer <Font> Pointer;
	virtual Size getSize (const std::string & str) const = 0;
	virtual void render (double x, double y, const std::string & str)
		const = 0;
};

class BitmapFont :
	public Font
{
private:
	const void * font;	// GLUT font
	const int height;	// the height of the font
public:
	BitmapFont (const void * p_font, int p_height);
	virtual Size getSize (const std::string & str) const;
	void render (double x, double y, const std::string & str) const;

	// some predefined bitmap fonts
	static SmartPointer <Font> Reg9by15;
	static SmartPointer <Font> Reg8by13;
	static SmartPointer <Font> Helvetica10;
	static SmartPointer <Font> Helvetica12;
	static SmartPointer <Font> Helvetica18;
	static SmartPointer <Font> TimesRoman10;
	static SmartPointer <Font> TimesRoman24;
};

#ifdef INDENTATION_HACK
{ // to fool emacs into not indenting the rest of the file
#endif
} // namespace OGLUI


#endif

