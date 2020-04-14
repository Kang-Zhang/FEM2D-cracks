#ifndef __VECTOR_HH__
#define __VECTOR_HH__

#include <assert.h>

template <class Type>
class Vector {
public:
	Type * array;
	long size;
private:
	long size_alloc;
public:
	
	Vector( long p_size) {
		assert( p_size > 0);
		size = p_size;
		size_alloc = p_size;
		array = new Type[ size_alloc];
		assert( array != NULL);
		for( long i = 0 ; i < size ; i ++)
			array[i] = 0;
	}
	
	Vector( const Vector & v) {
		size = v.size;
		size_alloc = v.size_alloc;
		array = new Type[ size_alloc];
		assert( array != NULL);
		for( long i = 0 ; i < size ; i ++)
			array[i] = v.array[i];
	}
	
	Vector() {
		size = 0;
		array = NULL;
	}

	int empty( void) {
		return array == NULL;
	}

	~Vector() {
		if( ! empty())
			delete [] array;
	}

	void resize( long new_size) {
		assert( new_size >= size);
		assert( array != NULL);
		// only resize if necessary
		if( new_size > size_alloc) {
			size_alloc = long( new_size * 1.5) + 10;
			Type * new_array = new Type[ size_alloc];
			assert( new_array != NULL);
			if( array != NULL)
			{
				// move the old stuff into the new array
				for( long i = 0 ; i < size ; i ++)
					new_array[i] = array[i];
				// delete the old array and replace with
				// the new one
				delete [] array;
			}
			array = new_array;
		}
		// set to zero's the new entries
		for( long i = size ; i < new_size ; i ++)
			array[i] = 0;
		size = new_size;
	}

	Type operator() (long ind) const {
		assert( ind >= 0);
		assert( ind < size);
		return array[ind];
	}

	Type & operator() (long ind) {
		assert( ind >= 0);
		assert( ind < size);
		return array[ind];
	}
	
	Vector & operator= ( const Vector & v) {
		assert( size == v.size);
		for( long i = 0 ; i < size ; i ++)
			array[i] = v.array[i];
		return * this;
	}

	void copy( const Vector & a) {
		if( size != a.size) {
			Type * tmp = new Type[size];
			assert( tmp != NULL);
			for( long i = 0 ; i < a.size ; i ++)
				tmp[i] = a.array[i];
			if( array != NULL) delete [] array;
			size = a.size;
			array = tmp;
		} else {
			for( long i = 0 ; i < size ; i ++)
				array[i] = a.array[i];
		}
	}

	void sub( const Vector & a) {
		// computes this = this - a
		assert( size == a.size);
		for( long i = 0 ; i < size ; i ++)
			array[i] -= a.array[i];
	}

	void add( const Vector & a) {
		// computes this = this + a
		assert( size == a.size);
		for( long i = 0 ; i < size ; i ++)
			array[i] += a.array[i];
	}

	Type max_abs( void) const {
		assert( size > 0);
		Type max = array[0];
		for( long i = 1 ; i < size ; i ++)
			if( fabs( max) < fabs( array[i]))
				max = array[i];
		return max;
	}

	Type max_sq( void) const {
		assert( size > 0);
		Type max = array[0]*array[0];
		for( long i = 1 ; i < size ; i ++) {
			double val = array[i] * array[i];
			if( max < val) max = val;
		}
		return max;
	}

	void multiply_by_scalar( const Type a) {
		for( long i = 0 ; i < size ; i ++)
			array[i] *= a;
	}

	Type norm( void) const {
		Type sum = 0;
		for( long i = 0 ; i < size ; i ++)
			sum += array[i] * array[i];
		return sum;
	}

	void square( void) {
		for( long i = 0 ; i < size ; i ++)
			array[i] *= array[i];
	}

	Type dot_product( const Vector & b) const {
		assert( b.size == size);
		Type sum = 0;
		for( long i = 0 ; i < size ; i ++)
			sum += array[i] * b.array[i];
		return sum;
	}

	void print( char * str, char * fmtstr = NULL) {

		if( fmtstr == NULL) fmtstr = " %8.4f";
		printf( "%s[", str);
		for( long i = 0 ; i < size ; i ++)
			printf( fmtstr, array[i]);
		printf( "]\n");
	}
};

typedef Vector<double> Vector_double;

#endif
