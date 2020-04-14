#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "SparseMatrix.hh"

SparseMatrix::SparseMatrix( long p_n_rows, long p_n_cols)
{
	n_rows = p_n_rows;
	n_cols = p_n_cols;

	// allocate rows
	rows = (Row **) malloc( n_rows * sizeof( Row *));
	assert( rows != NULL);
	for( long i = 0 ; i < n_rows ; i ++)
		rows[i] = new Row();
}

SparseMatrix::SparseMatrix( const SparseMatrix & m)
{
	// copy the matrix
	n_rows = m.n_rows;
	n_cols = m.n_cols;
	rows = (Row **) malloc( n_rows * sizeof( Row *));
	assert( rows != NULL);
	for( long r = 0 ; r < n_rows ; r ++) {
		rows[r] = new Row();
		rows[r]-> copy( * m.rows[r]);
	}
}

SparseMatrix::~SparseMatrix()
{
//	delete [] rows;
	for( long i = 0 ; i < n_rows ; i ++)
		delete rows[i];
	free( rows);
}

void SparseMatrix::expand( long new_n_rows, long new_n_cols)
{
	rows = (Row **) realloc( rows, new_n_rows * sizeof( Row *));
	for( long i = n_rows ; i < new_n_rows ; i ++)
		rows[i] = new Row();
	n_rows = new_n_rows;
	n_cols = new_n_cols;
}

double SparseMatrix::get( long row, long col) const 
{
	assert( row >= 0 && row < n_rows);
	assert( col >= 0 && col < n_cols);
	return rows[row]-> get( col);
}

void SparseMatrix::set( long row, long col, double val)
{
	assert( row >= 0 && row < n_rows);
	assert( col >= 0 && col < n_cols);
	assert( finite( val));

	rows[row]-> set( col, val);
}

void SparseMatrix::add( long row, long col, double val)
{
	assert( row >= 0 && row < n_rows);
	assert( col >= 0 && col < n_cols);
	assert( finite( val));

	rows[row]-> add( col, val);
}

void SparseMatrix::print( char * str, char * fmtstr)
{
	long len = strlen( str);

	// print the matrix row by row
	for( long r = 0 ; r < n_rows ; r ++) {
		if( r == n_rows / 2)
			printf( "%s", str);
		else
			for( long i = 0 ; i < len ; i ++)
				printf( " ");
		if( r == 0) printf( " /");
		else if( r == n_rows-1) printf( " \\");
		else printf( "| ");
		for( long c =0 ; c < n_cols ; c ++)
			printf( fmtstr, get(r,c));
		if( r == 0) printf( " \\");
		else if( r == n_rows-1) printf( " /");
		else printf( "  |");
		printf( "\n");
	}
}

void SparseMatrix::reset_row( long row)
{
	assert( row >= 0 && row < n_rows);

	rows[row]-> reset();
}


void SparseMatrix::reset( void)
{
	for( long r = 0 ; r < n_rows ; r ++)
		rows[r]-> reset();
}

void SparseMatrix::del_entry( long row, long col)
{
	assert( row >= 0 && row < n_rows);
	assert( col >= 0 && col < n_cols);
	rows[row]-> del_entry( col);
}

void SparseMatrix::copy( SparseMatrix & m)
{
	// erase everything
	for( long i = 0 ; i < n_rows ; i ++)
		delete rows[i];
	free( rows);

	// copy the matrix
	n_rows = m.n_rows;
	n_cols = m.n_cols;
	rows = (Row **) malloc( n_rows * sizeof( Row *));
	assert( rows != NULL);
	for( long r = 0 ; r < n_rows ; r ++) {
		rows[r] = new Row();
		rows[r]-> copy( * m.rows[r]);
	}
}

void SparseMatrix::mult( const Vector_double & a, Vector_double & res) const
// ----------------------------------------------------------------------
// Calculates [res] = [this].[a]
// ----------------------------------------------------------------------
{
	assert( a.size == n_rows);
	assert( res.size = a.size);

	for( long i = 0 ; i < n_rows ; i ++) {
		Row * r = rows[i];
		Row::Entry * entries = rows[i]-> entries;
		double sum = 0;
		for( long j = 0 ; j < r-> n_entries ; j ++)
//			sum += entries[j].val * a( entries[j].col);
			sum += entries[j].val * a.array[entries[j].col];
		res.array[i] = sum;
	}
	
}
