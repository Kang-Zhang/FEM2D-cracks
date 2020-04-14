#ifndef __SERIALIZABLE_H__
#define __SERIALIZABLE_H__

#include <iostream>
#include <string>

namespace OGLUI {
#ifdef INDENTATION_HACK
} // to fool emacs into not indenting the rest of the file
#endif

class Serializable
{
public:
	virtual void writeYourself (std::ostream & is) const = 0;
	virtual void readYourself (std::istream & is) = 0;
};

std::istream & operator >> (std::istream & is, Serializable & obj);
std::ostream & operator << (std::ostream & os, const Serializable & obj);
std::string serialize (const std::string & str);

#ifdef INDENTATION_HACK
{ // to fool emacs into not indenting the rest of the file
#endif
} // namespace OGLUI

#endif
