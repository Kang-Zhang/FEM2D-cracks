#ifndef __OGLUI_SLOT1_H__
#define __OGLUI_SLOT1_H__

#include <iostream>
#include <vector>
#include <SmartPointer.h>

namespace OGLUI {
#ifdef INDENTATION_HACK
} // to fool emacs into not indenting the rest of the file
#endif

template <class RetType, class Arg1Type>
class Slot1
{
private:
	class _Callback1
	{
	public:
		virtual RetType call (Arg1Type arg1) = 0;
		virtual ~_Callback1 () {;}
	};
	
	template <class Any>
	class _Callback1MethodPtr : public _Callback1
	// this is a helper callback - for member functions
	{
	private:
		Any & objref;
		RetType (Any::*funcref) (Arg1Type arg1);
	public:
		_Callback1MethodPtr (
			Any & p_objref,
			RetType (Any::*p_funcref) (Arg1Type arg1)
			) :
			objref (p_objref),
			funcref (p_funcref)
			{;}
		
		virtual RetType call (Arg1Type arg1)
			{
//				cerr << "_Callback1MethodPtr::call()\n";
				return (objref.*funcref) (arg1);
			}
	};
	
	class _Callback1StaticFunc : public _Callback1
	// this ia a helper callback - for global functions
	{
	private:
		RetType (*funcref) (Arg1Type arg1);
	public:
		_Callback1StaticFunc (
			RetType (*p_funcref) (Arg1Type arg1)
			) :
			funcref (p_funcref)
			{;}
		
		virtual RetType call (Arg1Type arg1)
			{
//				cerr << "_Callback1StaticFunc::call()\n";
				return (*funcref) (arg1);
			}
	};
	
	SmartPointer <_Callback1> cbptr;
	// no default constructors
	Slot1 ();
public:
	template <class Any>
	Slot1 (Any & objref, RetType (Any::*fptr) (Arg1Type arg1))
	{
		cbptr = new _Callback1MethodPtr <Any>
			(objref, fptr);
	}
	Slot1 (RetType (*fptr) (Arg1Type arg1)) {
		cbptr = new _Callback1StaticFunc
			(fptr);
	}		
	RetType call (Arg1Type arg1)
	{
		return cbptr-> call (arg1);
	}

	RetType operator () (Arg1Type arg1)
	{
		return call (arg1);
	}

	void operator = (const Slot1 <RetType,Arg1Type> & sl)
	{
		cbptr = sl.cbptr;
	}

	~Slot1 ()
	{
	}
};

template <class RetType, class Arg1Type>
Slot1 <RetType,Arg1Type>
slot (RetType (* funcptr) (Arg1Type arg1))
{
	return Slot1<RetType,Arg1Type> (funcptr);
}

template <class Any, class RetType, class Arg1Type>
Slot1 <RetType,Arg1Type>
slot (Any & objref, RetType (Any::*methptr) (Arg1Type arg1))
{
	return Slot1<RetType,Arg1Type> (objref, methptr);
}

template <class Arg1Type>
class Signal1
{
private:
	typedef Slot1<void,Arg1Type> SlotType;
	std::vector<SlotType> slots;
public:
	void connect (const SlotType & slot)
	{
		slots.push_back (slot);
	}

	void emit (Arg1Type arg1)
	{
		for (typename std::vector<SlotType>::iterator sl = slots.begin ();
		     sl != slots.end ();
		     sl ++)
		{
			sl-> call (arg1);
		}
	}
};

#ifdef INDENTATION_HACK
{ // to fool emacs into not indenting the rest of the file
#endif
} // namespace OGLUI


#endif

