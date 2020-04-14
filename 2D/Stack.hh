#ifndef __STACK_HH_
#define __STACK_HH_

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

template <class Type> class Stack {
private:
public:
	Type * array;
	long n;
	long n_alloc;

	Stack( void) {
		n_alloc = 10;
		array = new Type[ n_alloc];
		assert( array != NULL);
		n = 0;
	}

	Stack( const Stack & s) {
		n_alloc = s.n_alloc;
		n = s.n;
		Type * tmp = new Type[ n_alloc];
		assert( tmp != NULL);
		for( long i = 0 ; i < n ; i ++)
			tmp[i] = s.array[i];
		array = tmp;
	}

	~Stack() {
		assert( array != NULL);
		delete [] array;
	}
	
	void push( const Type & val) {
		assert( n <= n_alloc);
		// reallocate space if needed
		if( n == n_alloc) {
			n_alloc = long (1.5*n_alloc) + 10;
			Type * old_array = array;
			array = new Type[ n_alloc];
			for( long i = 0 ; i < n ; i ++)
				array[i] = old_array[i];
			delete [] old_array;
		}
		array[ n] = val;
		n ++;
	}

	void push_unique( const Type & val) {
		// if entry is already on the stack, don't add it
		for( long i = 0 ; i < n ; i ++)
			if( array[i] == val) return;
		
		// reallocate space if needed
		if( n == n_alloc) {
			n_alloc = long (1.5*n_alloc) + 10;
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

	int contains( const Type & val) {
		for( long i = 0 ; i < n ; i ++)
			if( val == array[i])
				return 1;
		return 0;
	}

	int remove_val( Type & val) {
		long ind = -1;
		for( long i = 0 ; i < n ; i ++)
			if( val == array[i]) {
				ind = i;
				break;
			}
		if( ind == -1) return 0;
		array[ind] = array[n-1];
		n --;
		return 1;
	}

	Type & pop( void) {
		assert( n > 0);
		n --;
		return array[n];
	}

	Type & pop_max( void) {
		assert( n > 0);
		Type & max = array[0];
		long max_ind = 0;
		for( long i = 1 ; i < n ; i ++)
			if( max < array[i]) {
				max = array[i];
				max_ind = i;
			}
		// replace max_ind with the last entry
		array[max_ind] = array[n-1];
		// return result
		n --;
		return max;
	}

	Stack<Type> & operator= ( const Stack<Type> & s) {
		n_alloc = s.n_alloc;
		n = s.n;
		Type * tmp = new Type[ n_alloc];
		assert( tmp != NULL);
		for( long i = 0 ; i < n ; i ++)
			tmp[i] = s.array[i];
		delete [] array;
		array = tmp;
		return * this;
	}

	Type operator() (long ind) const {
		assert( ind >= 0);
		assert( ind < n);
		return array[ind];
	}

	Type & operator() (long ind) {
		assert( ind >= 0);
		assert( ind < n);
		return array[ind];
	}
		
};

#endif
