#ifndef __OGLUI_THEME_H__
#define __OGLUI_THEME_H__

#include <string>
#include <map>
#include <iostream>
#include "SmartPointer.h"

namespace OGLUI {
#ifdef INDENTATION_HACK
} // to fool emacs into not indenting the rest of the file
#endif

class Theme
{
public:
	typedef SmartPointer <Theme> Pointer;
	class Themable
	{
	public:
		typedef SmartPointer <Themable> Pointer;
		virtual ~Themable () {;}
	};
protected:
	/// list of attributes
	typedef  std::map <std::string,Themable::Pointer> _List;
	_List _list;
	/// name of the theme
	std::string _name;
public:
	/// empty constructor
	Theme (const std::string & name) :
		_name (name)
		{
		}
	/// put a new attribute into the list, or rewrite an existing
	/// attribute
	void set (const std::string name, Themable::Pointer ptr)
		{
			_list [name] = ptr;
		}
	/// get an attribute out of the list. Returns NULL if the attribute
	/// cannot be found.
	Themable::Pointer get (const std::string name)
		{
			_List::iterator i = _list.find (name);
			if (i == _list.end ())
				return Themable::Pointer ();
			else
				return i-> second;
		}
	/// prints a list of attributes to stderr (for debugging purposes)
	void printList ()
		{
			std::cerr << "Theme (" << _name << "):\n";
			for (_List::iterator i = _list.begin ();
			     i != _list.end ();
			     i ++)
			{
				std::cerr << "\t" << i-> first << "\n";
			}
		}
	
};

#ifdef INDENTATION_HACK
{ // to fool emacs into not indenting the rest of the file
#endif
} // namespace OGLUI


#endif
