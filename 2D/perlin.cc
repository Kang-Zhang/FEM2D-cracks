#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>
#include <string.h>

/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0df   /* constant vector a */
#define UPPER_MASK 0x80000000 /* most significant w-r bits */
#define LOWER_MASK 0x7fffffff /* least significant r bits */

/* Tempering parameters */   
#define TEMPERING_MASK_B 0x9d2c5680
#define TEMPERING_MASK_C 0xefc60000
#define TEMPERING_SHIFT_U(y)  (y >> 11)
#define TEMPERING_SHIFT_S(y)  (y << 7)
#define TEMPERING_SHIFT_T(y)  (y << 15)
#define TEMPERING_SHIFT_L(y)  (y >> 18)

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializing the array with a NONZERO seed */
void
sgenrand( unsigned long seed)
{
	/* setting initial seeds to mt[N] using         */
	/* the generator Line 25 of Table 1 in          */
	/* [KNUTH 1981, The Art of Computer Programming */
	/*    Vol. 2 (2nd Ed.), pp102]                  */
	mt[0]= seed & 0xffffffff;
	for (mti=1; mti<N; mti++)
		mt[mti] = (69069 * mt[mti-1]) & 0xffffffff;
}

double genrand( void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if sgenrand() has not been called, */
            sgenrand(4357); /* a default initial seed is used   */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1];

        mti = 0;
    }
  
    y = mt[mti++];
    y ^= TEMPERING_SHIFT_U(y);
    y ^= TEMPERING_SHIFT_S(y) & TEMPERING_MASK_B;
    y ^= TEMPERING_SHIFT_T(y) & TEMPERING_MASK_C;
    y ^= TEMPERING_SHIFT_L(y);

    return ( (double)y / (unsigned long)0xffffffff ); /* reals */
    /* return y; */ /* for integer generation */
}

void my_setrand( unsigned long seed)
{
//	srand48( seed);
	sgenrand( seed);
}

double my_getrand( void)
{
//	return drand48();
	return genrand();
}

double interpolate( double v1, double v2, double x)
{
	assert( x >= 0.0 && x <= 1.0);

	double ft = x * M_PI;
	double f = (1 - cos( ft)) * 0.5;

	return  v1*(1-f) + v2*f;

	// linear interpolation:
	// return v1 * (1-x) + v2 * x;
}

void usage( void)
{
	fprintf( stderr,
		 "Usage: perlin width height persistence octaves\n"
		 "              init_size size_ratio [rgb|gray]\n");
	exit( -1);
}

void octave_get_2d(
	long width,
	long height,
	double dx,
	double dy,
	double ** & img)
{
	// prepare a map of random points
	long map_width = long( width / dx + 2);
	long map_height = long( height / dy + 2);
	fprintf( stderr, "map = %ld x %ld points\n", map_width, map_height);
	assert( map_width > 0);
	assert( map_height > 0);
//	double map[map_height][map_width];
	double ** map = (double **) malloc( sizeof( double *) * map_height);
	assert( map != NULL);
	for( long y = 0 ; y < map_height ; y ++) {
	        map[y] = (double *) malloc( sizeof( double) * map_width);
		assert( map[y] != NULL);
		for( long x = 0 ; x < map_width ; x ++)
		        map[y][x] = my_getrand() * 2 - 1;
	}
	
	// fill the picture with interpolated values
	for( long y = 0 ; y < height ; y ++)
		for( long x = 0 ; x < width ; x ++)
		{
			// calculate the adjacent lower left grid point
			long gx = long( x / dx);
			long gy = long( y / dy);
			assert( gx+1 < map_width);
			assert( gy+1 < map_height);

			// get the values in the surrounding corners
			double v1 = map[gy+0][gx+0];
			double v2 = map[gy+0][gx+1];
			double v3 = map[gy+1][gx+0];
			double v4 = map[gy+1][gx+1];

			// what are the isoparameteric coordinates inside
			// the grid square
			double px = (x - gx*dx) / dx;
			double py = (y - gy*dy) / dy;

			// interpolate on the edges of the grid square
			double i1 = interpolate( v1, v2, px);
			double i2 = interpolate( v3, v4, px);
			double i3 = interpolate( i1, i2, py);

			// put the result in the image
			img[y][x] = i3;
		}
}

void perlin_generate_2d(
	long width,
	long height,
	double persistence,
	long n_octaves,
	double init_size,
	double size_ratio,
	double ** & img)
{
	// allocate room for the result image
	img = (double **) malloc( sizeof( double *) * height);
	for( long i = 0 ; i < height ; i ++)
		img[i] = (double *) malloc( sizeof( double) * width);

	// allocate room for the temporary (octave) image
	double ** oct = (double **) malloc( sizeof( double *) * height);
	for( long i = 0 ; i < height ; i ++)
		oct[i] = (double *) malloc( sizeof( double) * width);

	// zero out the resulting image
	for( long y = 0 ; y < height ; y ++)
		for( long x = 0 ; x < width ; x ++)
			img[y][x] = 0;

	// add all octaves
	double size = init_size;
	double weight = persistence;
	double total_weight = 0.0;
	for( long i = 0 ; i < n_octaves ; i ++) {
		fprintf( stderr, "Processing octave %ld...\n", i);
		fprintf( stderr, "\tsize = %f\n", size);
		fprintf( stderr, "\tweight = %f\n", weight);
		fprintf( stderr, "\ttotal_weight = %f\n", total_weight);
		octave_get_2d( width, height, size, size, oct);
		// add the octave to the result
		for( long y = 0 ; y < height ; y ++)
			for( long x = 0 ; x < width ; x ++)
				img[y][x] += weight * oct[y][x];
		total_weight += weight;
		weight *= persistence;
		size = size * size_ratio;
		fprintf( stderr, "done.\n");
	}

	// divide the image by the total weight
	for( long y = 0 ; y < height ; y ++)
		for( long x = 0 ; x < width ; x ++)
			img[y][x] /= total_weight;
}

int main( int argc, char ** argv)
{
	long width, height, octaves;
	double persistence, init_size, size_ratio;
	int rgb_mode = 0;

	// get the parameters:
	//   - image dimensions (width/height)
	//   - persistence & octaves
	if( argc != 8) usage();
	if( 1 != sscanf( argv[1], "%ld", & width)) usage();
	if( 1 != sscanf( argv[2], "%ld", & height)) usage();
	if( 1 != sscanf( argv[3], "%lf", & persistence)) usage();
	if( 1 != sscanf( argv[4], "%ld", & octaves)) usage();
	if( 1 != sscanf( argv[5], "%lf", & init_size)) usage();
	if( 1 != sscanf( argv[6], "%lf", & size_ratio)) usage();
	if( strcasecmp( argv[7], "rgb") == 0)
		rgb_mode = 1;
	else if( strcasecmp( argv[7], "gray") == 0)
		rgb_mode = 0;
	else usage();
	if( width <= 0) usage();
	if( height <= 0) usage();
	if( persistence < 0) usage();
	if( octaves <= 0) usage();

	// get the noise
	double ** imgr = NULL;
	double ** imgg = NULL;
	double ** imgb = NULL;
	perlin_generate_2d( width, height,
			    persistence, octaves,
			    init_size, size_ratio,
			    imgr);

	if( rgb_mode) {
		perlin_generate_2d( width, height,
				    persistence, octaves,
				    init_size, size_ratio,
				    imgg);
		
		perlin_generate_2d( width, height,
				    persistence, octaves,
				    init_size, size_ratio,
				    imgb);
	} else {
		imgg = imgr;
		imgb = imgr;
	}

	// print out the data
	printf( "P3\n");
	printf( "# PERLIN NOISE creator\n");
	printf( "%ld %ld\n", width, height);
	printf( "255\n");
	long count = 0;
	for( long row = 0 ; row < height ; row ++) {
		for( long col = 0 ; col < width ; col ++) {
			long c = long(256*(imgr[row][col]+1)/2);
			printf( "%3ld ", c);
			c = long(256*(imgg[row][col]+1)/2);
			printf( "%3ld ", c);
			c = long(256*(imgb[row][col]+1)/2);
			printf( "%3ld", c);
			count += 3;
			if( count > 15) {
				printf( "\n");
				count = 0;
			} else {
				printf( " ");
			}
		}
	}
	printf( "\n");

#ifdef DONT_COMPILE
	// get the noise
	double ** img;
	perlin_generate_2d( width, height,
			    persistence, octaves,
			    init_size, size_ratio,
			    img);

	// print out the data
	printf( "P3\n");
	printf( "# PERLIN NOISE creator\n");
	printf( "%ld %ld\n", width, height);
	printf( "255\n");
	long count = 0;
	for( long row = 0 ; row < height ; row ++) {
		for( long col = 0 ; col < width ; col ++) {
			long c = long(256*(img[row][col]+1)/2);
			printf( "%3ld %3ld %3ld", c, c, c);
			count += 3;
			if( count > 15) {
				printf( "\n");
				count = 0;
			} else {
				printf( " ");
			}
		}
	}
	printf( "\n");
#endif
}
