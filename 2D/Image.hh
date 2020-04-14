#ifndef __IMAGE_HH__
#define __IMAGE_HH__

class Image {

private:
	unsigned char ** data;

public:
	Image();
	int load( const char * fname);
	unsigned char get( long row, long col);

	long width, height;
};

#endif
