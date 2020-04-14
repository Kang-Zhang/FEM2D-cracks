#ifndef __POINTER_H__
#define __POINTER_H__

#include <cstdlib>
#include <iostream>

template <class T>
class SmartPointer
{
private:
	class PointerOwner
	{
		friend class SmartPointer <T>;
	private:
		T * _objPtr;	// the actual pointer to the data
		long count;	// how many things are pointing to me
	public:
		explicit PointerOwner (T * p_ptr) :
			_objPtr (p_ptr),
			count (0)
			{;}
		void incCount ()
			{
				count ++;
			}
		void decCount ()
			{
				count --;
				if (count == 0)
					release ();
				if (count < 0)
					throw "PointerOwner::count negative";
			}
		bool empty () const
			{
				return count == 0;
			}
		void release ()
			{
				if (count > 0)
					throw "PointerOwner: release on "
						"positive count.";
				delete _objPtr;
				_objPtr = NULL;
			}
		~PointerOwner ()
			{
				release ();
			}
	private:
		PointerOwner (const PointerOwner &);
		PointerOwner & operator= (const PointerOwner &);
	};

public:
	PointerOwner * _ptr;
public:
	SmartPointer () :
		_ptr (NULL)
		{;}
	SmartPointer (T * ptr)
		{
			if (ptr)
			{
				_ptr = new PointerOwner (ptr);
				_ptr-> incCount ();
			}
			else
			{
				_ptr = NULL;
			}
		}
	// copy constructor
	SmartPointer (const SmartPointer & p)
		{
			_ptr = p._ptr;
			if (_ptr)
			{
				_ptr-> incCount ();
			}
		}
	// copy constructor on a different type
	template <class U>
	SmartPointer (const SmartPointer <U> & p)
		{
			// do the dynamic cast thingy, so that the compiler
			// can complain when it is impossible
			T * dummy = dynamic_cast <T*> (p.getPtr ());
			// check the run-time result
			if (dummy == NULL)
			{
				_ptr = NULL;
			}
			else
			{
				_ptr = (PointerOwner *) p._ptr;
				if (_ptr)
				{
					_ptr-> incCount ();
				}
				else
				{
					std::cerr << "WARNING " << __FILE__
						  << ":" << __LINE__ << "\n";
				}
			}
		}
	bool isNull () const
		{
			return _ptr == NULL;
		}
	// assignment operator - for a regular pointer
	void operator = (T * ptr)
		{
			release ();
			if (ptr == NULL)
			{
				// special case
				_ptr = NULL;
			}
			else
			{
				_ptr = new PointerOwner (ptr);
				_ptr-> incCount ();
			}
		}
	// assignment operator - for a smart pointer of different type
	template <class U>
	void operator = (const SmartPointer <U> & p)
		{
			// make sure this is not a self-assignment test
			if ((void *) this == (void *) & p) return;
			// release whatever we pointed to
			release ();
			// do the dynamic cast thingy, so that the compiler
			// can complain when it is impossible
			T * dummy = dynamic_cast <T*> (p.getPtr ());
			// check the run-time result
			if (dummy == NULL)
			{
				_ptr = NULL;
			}
			else
			{
				_ptr = (PointerOwner *) p._ptr;
				if (_ptr)
				{
					_ptr-> incCount ();
				}
			}
		}
	// regular assignment operator - on the pointer of the same type
	void operator = (const SmartPointer & p)
		{
			// make sure this is not a self-assignment test
			if ((void *) this == (void *) & p) return;
			// release whatever we pointed to
			release ();
			_ptr = p._ptr;
			if (_ptr)
			{
				_ptr-> incCount ();
			}
		}
	// releases the data pointed to by this smart pointer
	void release ()
		{
			if (_ptr) {
				_ptr-> decCount();
				if( _ptr-> empty ())
					delete _ptr;
				_ptr = NULL;
			}
		}
	// member access operator
	T * operator -> ()
		{
			if (_ptr == NULL)
				throw "SmartPointer: NULL -> invoked\n";
			return _ptr-> _objPtr;
		}
	// member access operator
	const T * operator -> () const
		{
			if (_ptr == NULL)
				throw "SmartPointer: NULL -> invoked\n";
			return _ptr-> _objPtr;
		}
	// equal sign operator
	bool operator == (const SmartPointer & p)
		{
			if (_ptr == NULL) return p._ptr == NULL;
			if (p._ptr == NULL) return false;
			return _ptr == p._ptr;
		}
	bool operator == (const T * p)
		{
			if (_ptr == NULL) return p == NULL;
			return _ptr-> _objPtr == p;
		}
	bool operator != (const T * p)
		{
			return ! operator == (p);
		}
	// dereference operator
	T & operator * ()
		{
			if (_ptr == NULL)
				throw "SmartPointer: NULL . invoked\n";
			return * _ptr-> _objPtr;
		}
	// dereference operator
	const T & operator * () const
		{
			if (_ptr == NULL)
				throw "SmartPointer: NULL . invoked\n";
			return * _ptr-> _objPtr;
		}
	// destructor
	~SmartPointer ()
		{
			release ();
		}
	// get the actual pointer
	T * getPtr () const
		{
			if (_ptr == NULL)
				return NULL;
			return _ptr-> _objPtr;
		}

	// get the current count
	long getCount () const
		{
			if (isNull ()) return 0;
			return _ptr-> count;
		}
/*
	operator T * ()
		{
			if (_ptr == NULL)
				throw "SmartPointer: NULL T * invoked\n";
			return * _ptr-> _objPtr;
		}

	template <class U>
	operator U * ()
		{
			return dynamic_cast <U *> (_ptr-> _obrPtr);
		}
*/
/*
	/// Cast operation
	template <class U>
	SmartPointer RefCountPtr& cast( const RefCountPtr<U>& ptr )
		{
			if (__ptr != ptr.toPtr())
			{
				if (__ptr) __ptr->removeReference();
				__ptr = dynamic_cast<T *>(ptr.toPtr());
				if (__ptr) __ptr->addReference();
			};
			return *this;
		}
*/
	
};

// convenience function
template <class T>
SmartPointer <T> sptr (T * ptr)
{
	return SmartPointer <T> (ptr);
}

#endif
