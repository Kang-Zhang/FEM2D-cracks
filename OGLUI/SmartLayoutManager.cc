#include <SmartLayoutManager.h>
#include <Widget.h>

OGLUI_BEGIN_NAMESPACE

Geometry SmartLayoutManager::getGeometry (
	Widget * w
	, const Geometry & pgeom
	) const
// ======================================================================
// Reports the geometry of this widget.
// Method:
//    - lookup the stored position of this widget
//    - if not found, set to the default position: 0,0
//    - combine the extracted position with widget's prefered size
//    - return as a geometry
// ----------------------------------------------------------------------
{
	// we need to figure out x1, x2, y1, y2 of the widget
	double x1, x2, y1, y2;
	
	// query left attachment
	Attachment left;
	ListType::const_iterator i = _list.find (LeftEdge (w));
	if (i != _list.end ())
		left = i-> second;
	// query right attachment
	Attachment right;
	i = _list.find (RightEdge (w));
	if (i != _list.end ())
		right = i-> second;
	// query top attachment
	Attachment top;
	i = _list.find (TopEdge (w));
	if (i != _list.end ())
		top = i-> second;
	// query bottom attachment
	Attachment bottom;
	i = _list.find (BottomEdge (w));
	if (i != _list.end ())
		bottom = i-> second;

	// get the prefered size of the widget
	Size prefSize = w -> getPrefSize ();
	// get the real size of the parent
	Size psize = pgeom.getSize ();

	// do the horizontal portion of geometry
	// --------------------------------------------------
	// figure out x1 - if possible
	if (left._type == Attachment::RELATIVE)
		x1 = left._val/100 * psize.width + left._offset;
	else if (left._type == Attachment::ABSOLUTE)
		x1 = left._val + left._offset;
	else
		x1 = 0;
	// figure out x2 - if possible
	if (right._type == Attachment::RELATIVE)
		x2 = right._val/100 * psize.width - right._offset;
	else
		x2 = right._val - right._offset;
	// special cases
	if (right._type == Attachment::NONE)
		x2 = x1 + prefSize.width;
	if (left._type == Attachment::NONE)
		x1 = x2 - prefSize.width;

	// do the vertical portion of geometry
	// --------------------------------------------------
	// figure out y1 - if possible
	if (bottom._type == Attachment::RELATIVE)
		y1 = bottom._val/100 * psize.height + bottom._offset;
	else if (bottom._type == Attachment::ABSOLUTE)
		y1 = bottom._val + bottom._offset;
	else
		y1 = 0;
	// figure out y2 - if possible
	if (top._type == Attachment::RELATIVE)
		y2 = top._val/100 * psize.height - top._offset;
	else
		y2 = top._val - top._offset;
	// special cases
	if (top._type == Attachment::NONE)
		y2 = y1 + prefSize.height;
	if (bottom._type == Attachment::NONE)
		y1 = y2 - prefSize.height;

	return Geometry (
		Position (x1,y1) + pgeom.getPosition ()
		, Position (x2, y2) + pgeom.getPosition ());
}

SmartLayoutManager::SmartLayoutManager ()
// ======================================================================
// default constructor
// ......................................................................
{
	// empty
}

void SmartLayoutManager::attach (const Edge & e, const Attachment & a)
// ======================================================================
// attaches a given edge to the given attachment
// ......................................................................
{
	// special case for 'Attachment::NONE' - simply remove the entry
	if (a._type == Attachment::NONE)
		_list.erase (e);
	else
		_list [e] = a;
}

#include <algorithm>
bool memberOf (const std::vector <Widget *> & list, const Widget * w)
{
	return find (list.begin (), list.end (), w) != list.end ();
}

struct RHS {
	RHS ()
		: e(NULL,SLM::Edge::NONE)
		, hasEdge (false)
		, inProgress (false)
		{;}
	RHS (double pk)
		: e(NULL,SLM::Edge::NONE)
		, k (pk)
		, hasEdge (false)
		, inProgress (false)
		{;};
	RHS (const SLM::Edge & pedge, double pk)
		: e (pedge)
		, k (pk)
		, hasEdge (true)
		, inProgress (false)
		{;}
	RHS & operator = (const RHS & rhs)
		{
			e = rhs.e;
			k = rhs.k;
			hasEdge = rhs.hasEdge;
			inProgress = rhs.inProgress;
			return * this;
		}
	SLM::Edge e;
	double k;
	bool hasEdge;
	bool inProgress;
};

static double figureEdge (
	const SLM::Edge & e
	, std::map <SLM::Edge, RHS> & eq
	, bool & cyclic)
{
	cyclic = false;
	RHS & rhs = eq [e];
	// figure out whether this is a cycle
	if (rhs.inProgress)
	{
		cyclic = true;
		return 0;
	}
	// if RHS is known (i.e. it is not associated with another
	// edge), we are done
	if (! rhs.hasEdge)
	{
		cyclic = false;
		return rhs.k;
	}
	// otherwise we have to ask that edge to figure itself out,
	// which is a recursive call. To prevent cycles, we set
	// 'inProgress' of this edge to true before descending into
	// recursion
	rhs.inProgress = true;
	// and do the recursion
	double res = figureEdge (rhs.e, eq, cyclic);
	// if there is a cycle, report it
	if (cyclic)
	{
		return 0;
	}
	// otherwise add the result to the constant from the rhs
	double total = res + rhs.k;
	// overwrite the rhs
	eq [e] = total;
	// we are happy - no cycles
	cyclic = false;
	return total;
}

#include <iostream>

std::ostream & operator<< (std::ostream & os, const SLM::Edge & e)
{
	// id first
	if (e._widget != NULL)
		os << e._widget -> instanceName () << ".";
	else
		os << "NULL.";
	if (e._which == SLM::Edge::LEFT)
		os << "x1";
	else if (e._which == SLM::Edge::RIGHT)
		os << "x2";
	else if (e._which == SLM::Edge::BOTTOM)
		os << "y1";
	else if (e._which == SLM::Edge::TOP)
		os << "y2";
	else
		os << "???";
	return os;
}

std::ostream & operator << (std::ostream & os, const RHS & rhs)
{
	if (rhs.hasEdge)
		os << rhs.e << " + ";
	os << rhs.k;
	return os;
}

static void printEquations ( const std::string & name
			     , const std::vector <Widget *> & wList
			     , std::map <SLM::Edge,RHS> & eq)
{
	std::cerr << "Default equations:\n";
	for (size_t i = 0 ; i < wList.size () ; i ++)
	{
		std::cerr << "\t" << (SLM::LeftEdge (wList [i])) << " = "
			  << (eq [SLM::LeftEdge (wList [i])]) << "\n";
		std::cerr << "\t" << (SLM::RightEdge (wList [i])) << " = "
			  << (eq [SLM::RightEdge (wList [i])]) << "\n";
		std::cerr << "\t" << (SLM::BottomEdge (wList [i])) << " = "
			  << (eq [SLM::BottomEdge (wList [i])]) << "\n";
		std::cerr << "\t" << (SLM::TopEdge (wList [i])) << " = "
			  << (eq [SLM::TopEdge (wList [i])]) << "\n";
	}
}

void SmartLayoutManager::recompute (const Geometry & geom)
// ======================================================================
// Recomputes the geometries of all clients
// ......................................................................
{
	// create a list of all widgets I know about
	std::vector <Widget *> wList;
	for (ListType::const_iterator i = _list.begin ();
	     i != _list.end ();
	     i ++)
	{
		if (! memberOf (wList, i-> first._widget))
			wList.push_back (i-> first._widget);
	}
	std::cerr << "SLM::recompute: managing " << wList.size ()
		  << "widgets\n";
	// if there are no kids we are managing, there is nothing to do
	if (wList.size() == 0) return;
	// create a set of default equations - 4 for each managed widget
	typedef std::map <Edge, RHS> EqType;
	EqType eq;
	for (size_t i = 0 ; i < wList.size () ; i ++)
	{
		Size ps = wList [i] -> getPrefSize ();
		// left edge = right-edge - prefWidth
		eq [LeftEdge (wList [i])]
			= RHS (RightEdge(wList [i]), -ps.width);
		// right edge = left-edge + prefWidth
		eq [RightEdge (wList [i])]
			= RHS (LeftEdge(wList [i]), ps.width);
		// bottom edge = top-edge - prefHeight
		eq [BottomEdge (wList [i])]
			= RHS (TopEdge(wList [i]), -ps.height);
		// top edge = bottom-edge + prefHeight
		eq [TopEdge (wList [i])]
			= RHS (BottomEdge(wList [i]), ps.height);
	}
//	printEquations ("Default equations:", wList, eq);
	// modify equations according to every attachment
	for (ListType::const_iterator i = _list.begin ();
	     i != _list.end ();
	     i ++)
	{
		const Edge & e = i -> first;
		const Attachment & att = i -> second;
		if (att._type == Attachment::RELATIVE)
		{
			if (e._which == Edge::LEFT)
				eq [e] = att._val / 100 * geom.width ()
					+ att._offset;
			else if (e._which == Edge::RIGHT)
				eq [e] = att._val / 100 * geom.width ()
					- att._offset;
			else if (e._which == Edge::BOTTOM)
				eq [e] = att._val / 100 * geom.height ()
					+ att._offset;
			else // TOP
				eq [e] = att._val / 100 * geom.height ()
					- att._offset;
		}
		else if (att._type == Attachment::ABSOLUTE)
		{
			eq [e] = att._val + att._offset;
		}
		else if (att._type == Attachment::WIDGET)
		{
			if (e._which == Edge::TOP)
				eq [e] = RHS (BottomEdge (att._widget)
					      , - att._offset);
			else if (e._which == Edge::BOTTOM)
				eq [e] = RHS (TopEdge (att._widget)
					      , att._offset);
			else if (e._which == Edge::LEFT)
				eq [e] = RHS (RightEdge (att._widget)
					      , att._offset);
			else if (e._which == Edge::RIGHT)
				eq [e] = RHS (LeftEdge (att._widget)
					      , - att._offset);
		}
		else if (att._type == Attachment::OPPWIDGET)
		{
			if (e._which == Edge::TOP)
				eq [e] = RHS (TopEdge (att._widget)
					      , - att._offset);
			else if (e._which == Edge::BOTTOM)
				eq [e] = RHS (BottomEdge (att._widget)
					      , att._offset);
			else if (e._which == Edge::LEFT)
				eq [e] = RHS (LeftEdge (att._widget)
					      , att._offset);
			else if (e._which == Edge::RIGHT)
				eq [e] = RHS (RightEdge (att._widget)
					      , - att._offset);
		}
		else // must be 'no' attachment, which in theory we should
		     // not even allow...
		{
			std::cerr << "SML::recompute():NONE attachment !!!\n";
		}
	}
//	printEquations ("Updated equations:", wList, eq);
	// proceed to solve all unknown equations
	for (EqType::iterator i = eq.begin (); i != eq.end (); i ++)
	{
		bool cyclic;
		figureEdge (i-> first, eq, cyclic);
		// if there is a cycle, report it
		if (cyclic)
		{
			std::cerr << "SLM::recompute(): Encountered cycle.\n";
		}
	}
//	printEquations ("After solving:", wList, eq);

	// call 'onGeometryChanged' for all widget that we manage
	for (size_t i = 0 ; i < wList.size () ; i ++)
	{
		Widget * w = wList [i];
		// get the widget's current geometry
		Geometry currG = w-> getOutGeom ();
		// get the widget's new geometry
		double x1 = eq [LeftEdge (w)].k;
		double x2 = eq [RightEdge (w)].k;
		double y1 = eq [BottomEdge (w)].k;
		double y2 = eq [TopEdge (w)].k;
		Geometry newG ( Position (x1, y1) + geom.getPosition ()
				, Position (x2, y2));
		// if the new geometry is different, notify the widget
		if (currG != newG)
			w -> onGeometryChanged (newG);
	}

}

void SmartLayoutManager::setGeometry (Widget * w, const Geometry & g)
{
	attach (LeftEdge (w), Absolute (g.minx ()));
	attach (RightEdge (w), Absolute (g.maxx ()));
	attach (BottomEdge (w), Absolute (g.miny ()));
	attach (TopEdge (w), Absolute (g.maxy ()));
}
OGLUI_END_NAMESPACE
