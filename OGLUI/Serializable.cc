#include <sstream>
#include "Serializable.h"

namespace OGLUI {
#ifdef INDENTATION_HACK
} // to fool emacs into not indenting the rest of the file
#endif

std::istream & operator >> (std::istream & is, Serializable & obj)
{
	obj.readYourself (is);
	return is;
}

std::ostream & operator << (std::ostream & os, const Serializable & obj)
{
	obj.writeYourself (os);
	return os;
}

static std::string toHex (const void * ptr, long n_bytes)
{
	static char table[] = "01234567890abcdefgh";
	std::string res;
	const unsigned char * array = (unsigned char *) ptr;
	for (long i = 0 ; i < n_bytes ; i ++)
	{
		unsigned char c = array[i];
		res.push_back (table[c/16]);
		res.push_back (table[c%16]);
	}
	return res;
}

std::string serialize (const std::string & str)
{
	std::ostringstream res;
	res << "\"";
	for (size_t i = 0 ; i < str.length () ; i ++)
	{
		unsigned char c = str[i];
		if (c == '"')
			res << "\\\"";
		else if (c == '\\')
			res << "\\\\";
		else if (c == '\n')
			res << "\\n";
		else if (c < 32)
			res << "\\0x" << toHex (&c, sizeof (c));
		else if (c > 126)
			res << "\\0x" << toHex (&c, sizeof (c));
		else
			res << (char) c;
	}
	res << "\"";
	return res.str();
}

#ifdef INDENTATION_HACK
{ // to fool emacs into not indenting the rest of the file
#endif
} // namespace OGLUI
