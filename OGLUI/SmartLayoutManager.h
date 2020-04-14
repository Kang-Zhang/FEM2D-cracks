#ifndef __SMART_LAYOUT_MANAGER__
#define __SMART_LAYOUT_MANAGER__

#include <LayoutManager.h>
#include <map>

OGLUI_BEGIN_NAMESPACE
OGLUI_USE_NAMESPACE

class SmartLayoutManager
	: public LayoutManager
{
public:
	class Edge;
	class TopEdge;
	class BottomEdge;
	class LeftEdge;
	class RightEdge;
	class Attachment;
	class Relative;
	class Absolute;
	class None;
	class OtherWidget;
	class OppositeWidget;
	typedef SmartPointer <SmartLayoutManager> Pointer;
	/// default constructor
	SmartLayoutManager ();
	/// attaches a given edge to the given attachment
	void attach (const Edge & e, const Attachment & a);
	/// returns the geometry of the widget
	virtual Geometry getGeometry (Widget * child
				      , const Geometry & geom) const;
	/// recomputes all geometries
	virtual void recompute (const Geometry & geom);
	/// convenience function - for setting an overall geometry
	virtual void setGeometry (Widget * w, const Geometry & g);
private:
	typedef	std::map <Edge, Attachment> ListType;
	ListType _list;
};

class SmartLayoutManager::Edge
{
public:
	enum Which {TOP=0, LEFT, RIGHT, BOTTOM, NONE};
/*
	Edge ()
		: _widget (NULL)
		, _which (NONE)
		{;}
*/
	Edge (Widget * w, Which t)
		: _widget (w)
		, _which (t)
		{;}
	bool operator < (const Edge & e) const
		{
			if (_widget < e._widget) return true;
			if (_widget > e._widget) return false;
			return _which < e._which;
		}

	Widget * _widget;
	Which _which;
};

class SmartLayoutManager::TopEdge
	: public SmartLayoutManager::Edge
{
public:
	TopEdge (Widget * w)
		: Edge (w, TOP)
		{;}
};

class SmartLayoutManager::BottomEdge
	: public SmartLayoutManager::Edge
{
public:
	BottomEdge (Widget * w)
		: Edge (w, BOTTOM)
		{;}
};

class SmartLayoutManager::LeftEdge
	: public SmartLayoutManager::Edge
{
public:
	LeftEdge (Widget * w)
		: Edge (w, LEFT)
		{;}
};

class SmartLayoutManager::RightEdge
	: public SmartLayoutManager::Edge
{
public:
	RightEdge (Widget * w)
		: Edge (w, RIGHT)
		{;}
};

class SmartLayoutManager::Attachment
{
public:
	enum Type {RELATIVE, ABSOLUTE, WIDGET, OPPWIDGET, NONE};
	double _val;
	double _offset;
	Type _type;
	Widget * _widget;
	Attachment (double val, double offset, Type type)
		: _val (val)
		, _offset (offset)
		, _type (type)
		, _widget (NULL)
		{;}
	Attachment (Widget * w, double offset, Type type)
		: _val (0)
		, _offset (offset)
		, _type (type)
		, _widget (w)
		{;}
	Attachment ()
		: _val (0)
		, _offset (0)
		, _type (NONE)
		, _widget (NULL)
		{;}
};

class SmartLayoutManager::Relative
	: public SmartLayoutManager::Attachment
{
public:
	Relative (double val, double offset = 0)
		: Attachment (val, offset, RELATIVE)
		{;}
};

class SmartLayoutManager::Absolute
	: public SmartLayoutManager::Attachment
{
public:
	Absolute (double val, double offset = 0)
		: Attachment (val, offset, ABSOLUTE)
		{;}
};

class SmartLayoutManager::OtherWidget
	: public SmartLayoutManager::Attachment
{
public:
	OtherWidget (Widget * w, double offset = 0)
		: Attachment (w, offset, WIDGET)
		{;}
};

class SmartLayoutManager::OppositeWidget
	: public SmartLayoutManager::Attachment
{
public:
	OppositeWidget (Widget * w, double offset = 0)
		: Attachment (w, offset, OPPWIDGET)
		{;}
};

class SmartLayoutManager::None
	: public SmartLayoutManager::Attachment
{
public:
	None ()
		: Attachment (0.0, 0.0, NONE)
		{;}
};

/// shortcut
typedef SmartLayoutManager SLM;

OGLUI_END_NAMESPACE


#endif

