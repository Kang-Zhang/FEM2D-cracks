#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include "die.hh"

void die( const char * fmtstr, ...)
{
	va_list ap;
	va_start( ap, fmtstr);
	fprintf( stderr, "Error: ");
	vfprintf( stderr, fmtstr, ap);
	fprintf( stderr, "\n");
	va_end( ap);
	exit( -1);
}

