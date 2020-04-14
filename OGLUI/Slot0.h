#ifndef __OGLUI_SLOT0_H__
#define __OGLUI_SLOT0_H__

#include <iostream>
#include <vector>
#include <SmartPointer.h>

namespace OGLUI {
#ifdef INDENTATION_HACK
} // to fool emacs into not indenting the rest of the file
#endif

template <class RetType>
class Slot0
{
private:
	class _Callback0
	{
	public:
		virtual RetType call () = 0;
		virtual ~_Callback0 () {;};
	};
	
	template <class Any>
	class _Callback0MethodPtr : public _Callback0
	// this is a helper callback - for member functions
	{
	private:
		Any & objref;
		RetType (Any::*funcref) ();
	public:
		_Callback0MethodPtr (
			Any & p_objref,
			RetType (Any::*p_funcref) ()
			) :
			objref (p_objref),
			funcref (p_funcref)
			{;}
		
		virtual RetType call ()
			{
				return (objref.*funcref) ();
			}
	};
	
	class _Callback0StaticFunc : public _Callback0
	// this ia a helper callback - for global functions
	{
	private:
		RetType (*funcref) ();
	public:
		_Callback0StaticFunc (
			RetType (*p_funcref) ()
			) :
			funcref (p_funcref)
			{;}
		
		virtual RetType call ()
			{
//				cerr << "_Callback0StaticFunc::call()\n";
				return (*funcref) ();
			}
	};
	
	SmartPointer <_Callback0> cbptr;
	// no default constructors
	Slot0 ();
public:
	template <class Any>
	Slot0 (Any & objref, RetType (Any::*fptr) ())
	{
		cbptr = new _Callback0MethodPtr <Any>
			(objref, fptr);
	}
	Slot0 (RetType (*fptr) ()) {
		cbptr = new _Callback0StaticFunc
			(fptr);
	}		
	RetType call ()
	{
		return cbptr-> call ();
	}

	RetType operator () ()
	{
		return call ();
	}

	void operator = (const Slot0 <RetType> & sl)
	{
		cbptr = sl.cbptr;
	}

	~Slot0 ()
	{
	}
};

template <class RetType>
Slot0 <RetType>
slot (RetType (* funcptr) ())
{
	return Slot0<RetType> (funcptr);
}

template <class Any, class RetType>
Slot0 <RetType>
slot (Any & objref, RetType (Any::*methptr) ())
{
	return Slot0<RetType> (objref, methptr);
}

class Signal0
{
private:
	typedef Slot0<void> SlotType;
	std::vector<SlotType> slots;
public:
	void connect (const SlotType & slot)
	{
		slots.push_back (slot);
	}

	void emit ()
	{
		for (std::vector<SlotType>::iterator sl = slots.begin ();
		     sl != slots.end ();
		     sl ++)
		{
			sl-> call ();
		}
	}
};

#ifdef INDENTATION_HACK
{ // to fool emacs into not indenting the rest of the file
#endif
} // namespace OGLUI


#endif

