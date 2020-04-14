#ifndef __IMAGE_PNG_HH__
#define __IMAGE_PNG_HH__

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <malloc.h>
#include <png.h>

class ImagePNG {
public:
	long width, height;
	unsigned char * data;

	ImagePNG( long p_w, long p_h) {
		width = p_w;
		height = p_h;
		data = (unsigned char *) malloc( width * height * 3 + 1024);
		assert( data != NULL);
	}

	~ImagePNG() {
		if( data != NULL) free( data);
	}

	int save_ppm( const char * fname) {
		FILE * fp = fopen( fname, "wb");
		if( fp == NULL) return -1;

		fprintf( fp,
			 "P6\n"
			 "%ld %ld\n"
			 "255\n",
			 width, height);

		long row_width = width * 3;
		for( long i = height - 1 ; i >= 0 ; i --) {
			fwrite( data + row_width * i, row_width, 1, fp);
		}
		fclose( fp);
		return 0;
	}

	int save( const char * fname) {
		
		png_structp png_ptr = png_create_write_struct(
			PNG_LIBPNG_VER_STRING, NULL,
			NULL, NULL);
		if( png_ptr == NULL) {
			fprintf( stderr, "Cannot write PNG to file %s.\n",
				 fname);
			return -1;
		}

		png_infop info_ptr = png_create_info_struct(png_ptr);
		if (!info_ptr)
		{
			png_destroy_write_struct(&png_ptr,
						 (png_infopp)NULL);
			return -1;
		}

		FILE * fp = fopen( fname, "wb");
		if( fp == NULL) {
			fprintf( stderr, "Cannot write PNG to file %s.\n",
				 fname);
			return -1;
		}

		png_init_io( png_ptr, fp);

		png_set_compression_level( png_ptr, Z_BEST_COMPRESSION);

		png_set_IHDR( png_ptr, info_ptr, width, height, 8,
			      PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
			      PNG_COMPRESSION_TYPE_DEFAULT,
			      PNG_FILTER_TYPE_DEFAULT);
		
		png_color_8 sig_bit;
		sig_bit.red = 8;
		sig_bit.green = 8;
		sig_bit.blue = 8;
		png_set_sBIT( png_ptr, info_ptr, & sig_bit);

		png_write_info( png_ptr, info_ptr);

		png_byte ** row_ptrs = (png_byte **) malloc(
			height * sizeof( png_byte *));
		for( long i = 0 ; i < height ; i ++)
			row_ptrs[height-1-i] =
				(png_byte *) (data + i*width*3);

		png_write_image( png_ptr, row_ptrs);

		png_write_end( png_ptr, info_ptr);

		png_destroy_write_struct( & png_ptr, & info_ptr);

		free( row_ptrs);

		fclose( fp);

		return 0;

	}

};

#endif
