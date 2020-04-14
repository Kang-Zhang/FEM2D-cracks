#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <assert.h>

template <class Type> class Fifo {
private:
public:
	Type * array;
	long n;
	long n_alloc;

	Fifo( void) {
		n_alloc = 10;
		array = new Type[ n_alloc];
		assert( array != NULL);
		n = 0;
	}

	~Fifo() {
		delete [] array;
	}
	
	void add_unique( Type val) {
		// if this value is already in, ignore it
		for( long i = 0 ; i < n ; i ++)
			if( val == array[i]) return;
		// reallocate space if needed
		if( n == n_alloc) {
			n_alloc = long (1.1*n_alloc) + 10;
			Type * old_array = array;
			array = new Type[ n_alloc];
			for( long i = 0 ; i < n ; i ++)
				array[i] = old_array[i];
			delete [] old_array;
		}
		array[ n] = val;
		n ++;
	}

	int is_empty( void) {
		return n == 0;
	}

	Type get( void) {
		assert( n > 0);
		Type ret = array[0];
		for( long i = 1 ; i < n ; i ++)
			array[i-1] = array[i];
		n --;
		return ret;
	}
		
};
