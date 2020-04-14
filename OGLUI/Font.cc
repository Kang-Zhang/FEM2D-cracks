#include "Font.h"

namespace OGLUI {
#ifdef INDENTATION_HACK
} // to fool emacs into not indenting the rest of the file
#endif

// constructor
// ----------------------------------------------------------------------
BitmapFont::BitmapFont (
	const void * p_font
	, int p_height
	)
	: font (p_font)
	, height (p_height)
{
};


// returns the size of the text rendered using this font
// ----------------------------------------------------------------------
Size BitmapFont::getSize( const std::string & str) const
{
	Size res;
	res.height = height;
	res.width = glutBitmapLength (
		(void *) font
		, (const unsigned char *) str.c_str ());
	return res;
}

// renders a given text
// ----------------------------------------------------------------------
void BitmapFont::render (
	double x
	, double y
	, const std::string & str
	) const
{
	// Trickery - to get valid raster position even if outside
	// of the viewport, first set some valid position
	glRasterPos2d (0,0);
	// then move the raster
	glBitmap (0, 0, 0, 0, x, y, NULL);
	// finally, render the string
	for (const char * c = str.c_str () ; * c != '\0' ; c ++)
		glutBitmapCharacter ((void *) font, *c);
}



// ======================================================================
// Default fonts
// ......................................................................
SmartPointer <Font> BitmapFont::Reg9by15
= new BitmapFont (GLUT_BITMAP_9_BY_15,10);
SmartPointer <Font> BitmapFont::Reg8by13
= new BitmapFont (GLUT_BITMAP_8_BY_13,9);
SmartPointer <Font> BitmapFont::Helvetica10
= new BitmapFont (GLUT_BITMAP_HELVETICA_10,8);
SmartPointer <Font> BitmapFont::Helvetica12
= new BitmapFont (GLUT_BITMAP_HELVETICA_12,9);
SmartPointer <Font> BitmapFont::Helvetica18
= new BitmapFont (GLUT_BITMAP_HELVETICA_18,14);
SmartPointer <Font> BitmapFont::TimesRoman10
= new BitmapFont (GLUT_BITMAP_TIMES_ROMAN_10,7);
SmartPointer <Font> BitmapFont::TimesRoman24
= new BitmapFont (GLUT_BITMAP_TIMES_ROMAN_24,17);
	
#ifdef INDENTATION_HACK
{ // to fool emacs into not indenting the rest of the file
#endif
};
