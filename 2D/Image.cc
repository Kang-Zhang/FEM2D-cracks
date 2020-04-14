#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <malloc.h>
#include "Image.hh"
#include "readrgb.hh"

Image::Image()
{
	data = NULL;
	width = 0;
	height = 0;
}

int Image::load( const char * fname)
{
	fprintf( stderr, "Image::load( %s)\n", fname);

	int tmp_width, tmp_height, tmp_components;
	unsigned * tmp_data = read_texture( fname, & tmp_width, & tmp_height,
					    & tmp_components);

	fprintf( stderr, "  - size: %dx%d\n", tmp_width, tmp_height);
	fprintf( stderr, "  - components: %d\n", tmp_components);

	// make sure there were no errors
	if( tmp_data == NULL) {
		fprintf( stderr, "Image::load(%s): could not load.\n", fname);
		return -1;
	}

	// make sure this was a GREYSCALE image
	if( tmp_components != 1) {
		fprintf( stderr, "Image::load(%s): not a grey-scale image.\n",
			 fname);
		free( tmp_data);
		return -1;
	}

	// transform the data to our own
	height = tmp_height;
	width = tmp_width;
	data = new unsigned char * [height];
	for( long i = 0 ; i < height ; i ++) {
		data[i] = new unsigned char [width];
		for( long j = 0 ; j < width ; j ++)
			data[i][j] = ((unsigned char *) tmp_data)
				[i*width*4 + j*4 + 0];
	}
	free( tmp_data);

	return 0;
}

unsigned char Image::get( long row, long col)
{
	assert( row >= 0 && row < height);
	assert( col >= 0 && col < width);

	return data[row][col];
}

