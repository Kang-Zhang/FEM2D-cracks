#ifndef __OGLUI_BEVEL_H__
#define __OGLUI_BEVEL_H__

#include <GL/gl.h>
#include "Widget.h"

namespace OGLUI {
#ifdef INDENTATION_HACK
} // to fool emacs into not indenting the rest of the file
#endif

class BevelWidget : public Widget {
public:
	class LightWeights : public std::vector<double> {
	public:
		LightWeights () :
			std::vector<double> (3)
			{
				(*this) [0] = 1;
				(*this) [1] = 0;
				(*this) [2] = 0;
			}
		LightWeights (double pw1, double pw2, double pw3) :
			std::vector<double> (3)
			{
				(*this) [0] = pw1;
				(*this) [1] = pw2;
				(*this) [2] = pw3;
			}
		double w1 () { return (*this) [0]; }
		double w2 () { return (*this) [1]; }
		double w3 () { return (*this) [2]; }
	};

private:
	double _bevelSizeTop;
	double _bevelSizeBottom;
	double _bevelSizeLeft;
	double _bevelSizeRight;
protected:
	Color::Pointer _col2;
	Color::Pointer _col3;
	LightWeights _lightWeightsTop;
	LightWeights _lightWeightsBottom;
	LightWeights _lightWeightsLeft;
	LightWeights _lightWeightsRight;
	/// name of the class
	virtual const std::string className () const { return "BevelWidget"; }
	/// convenience function - to redraw only border
	virtual void redrawBorder ();
public:
	// constructor
	BevelWidget ( Widget * p_parent,
		      const std::string & p_name = "bevel");
	virtual void redraw ();
	void setBevelSize (double bs)
		{
			_bevelSizeTop = bs;
			_bevelSizeBottom = bs;
			_bevelSizeLeft = bs;
			_bevelSizeRight = bs;
		}
	double getBevelSizeTop () const { return _bevelSizeTop; }
	double getBevelSizeBottom () const { return _bevelSizeBottom; }
	double getBevelSizeLeft () const { return _bevelSizeLeft; }
	double getBevelSizeRight () const { return _bevelSizeRight; }
	void setBevelSizeTop (double s) {
		_bevelSizeTop = s; postRedisplay ();
	}
	void setBevelSizeBottom (double s) {
		_bevelSizeBottom = s; postRedisplay ();
	}
	void setBevelSizeLeft (double s) {
		_bevelSizeLeft = s; postRedisplay ();
	}
	void setBevelSizeRight (double s) {
		_bevelSizeRight = s; postRedisplay ();
	}
	void lightWeightsTop (const LightWeights & l) {
		_lightWeightsTop = l; postRedisplay ();
	}
	void lightWeightsBottom (const LightWeights & l) {
		_lightWeightsBottom = l; postRedisplay ();
	}
	void lightWeightsLeft (const LightWeights & l) {
		_lightWeightsLeft = l; postRedisplay ();
	}
	void lightWeightsRight (const LightWeights & l) {
		_lightWeightsRight = l; postRedisplay ();
	}
	/// returns the current size of the widget (minus borders, etc)
	virtual Geometry getInGeom () const
		{
			Geometry g = Widget::getInGeom ();
			Position pos = g.getPosition ();
			pos.x += _bevelSizeRight;
			pos.y += _bevelSizeBottom;
			Size size = g.getSize ();
			size.width -= _bevelSizeLeft + _bevelSizeRight;
			size.height -= _bevelSizeTop + _bevelSizeBottom;
			return Geometry (pos, size);
		}
};

#ifdef INDENTATION_HACK
{ // to fool emacs into not indenting the rest of the file
#endif
} // namespace OGLUI


#endif

