#ifndef __LAYOUT_MANAGER_H__
#define __LAYOUT_MANAGER_H__

#include <OGLUInamespace.h>
#include <Geometry.h>
#include <SmartPointer.h>
//#include <Widget.h>

OGLUI_BEGIN_NAMESPACE

class Widget;

class LayoutManager
{
public:
	/// Handle
	typedef SmartPointer <LayoutManager> Pointer;
	/// returns the geometry of this child, given the geometry of the
	/// parent
	virtual Geometry getGeometry (Widget * child
				      , const Geometry & geom) const = 0;
	/// virtual destructor - so that the compiler does not complain
	virtual ~LayoutManager () {;}
	/// this is called by the client widget to ask the layout manager
	/// to recompute the geometry of all kids
	virtual void recompute (const Geometry & geom) = 0;
};

OGLUI_END_NAMESPACE

#endif
