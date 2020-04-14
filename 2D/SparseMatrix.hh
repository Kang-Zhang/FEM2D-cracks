#ifndef __SPARSE_MATRIX_H__
#define __SPARSE_MATRIX_H__

#include <string.h>
#include <malloc.h>
#include "Vector.hh"

class SparseMatrix {
	
public:
	
	class Row {
	public:
		struct Entry {
			long col;
			double val;
		};
		long n_entries;
		Entry * entries;
		Row() { n_entries = 0; entries = NULL; }
		~Row() { if( entries != NULL) free( entries); }
		double get( long col) const {
			for( long i = 0 ; i < n_entries ; i ++)
				if( entries[i].col == col)
					return entries[i].val;
			return 0.0;
		}
		void set( long col, double val) {
			int found = 0;
			for( long i = 0 ; i < n_entries ; i ++)
				if( entries[i].col == col) {
					entries[i].val = val;
					found = 1;
					break;
				}
			if( ! found) {
				entries = (Entry *)
					realloc( entries, sizeof( Entry)
						 * (n_entries + 1));
				assert( entries != NULL);
				entries[ n_entries].col = col;
				entries[ n_entries].val = val;
				n_entries ++;
			}
		}
		void add( long col, double val) {
			int found = 0;
			for( long i = 0 ; i < n_entries ; i ++)
				if( entries[i].col == col) {
					entries[i].val += val;
					found = 1;
					break;
				}
			if( ! found) {
				entries = (Entry *)
					realloc( entries, sizeof( Entry)
						 * (n_entries + 1));
				assert( entries != NULL);
				entries[ n_entries].col = col;
				entries[ n_entries].val = val;
				n_entries ++;
			}
		}
		void reset( void) {
			n_entries = 0;
			if( entries != NULL) { free( entries); entries=NULL;}
		}
		void del_entry( long col)
		// delete an entry who's .col is 'col'
		{
			long ind = -1;
			for( long i = 0 ; i < n_entries ; i ++)
				if( entries[i].col == col) {
					ind = i;
					break;
				}
			if( ind == -1) {
				// nothing to do, the entry is not found
				return;
			}

			// replace the entry at ind with the last entry
			if( ind < n_entries-1)
				entries[ind] = entries[n_entries-1];
			n_entries --;
		}
		void delete_marked_entries( void)
		// delete all marked entries (marked entry is one with
		// negative .col
		// Method: iteratively move all non-marked entries
		//         into places that were marked (this is sort of
		//         a simple defragmentation)
		//
		//         repeat:
		//            find the left-most marked entry
		//            if no such is found, we are done
                //            find the right-most unmarked entry
		{
			// handle the special case
			if( n_entries == 0) return;

			// determine the number of entries after
			// degragmentation is done
			long new_n_entries = n_entries;
			for( long i = 0 ; i < n_entries ; i ++)
				if( entries[i].col < 0)
					new_n_entries --;

			long m_ind = 0;
			long u_ind = n_entries -1;
			while( 1) {
				while( entries[m_ind].col >=0 &&
				       m_ind < u_ind) m_ind ++;
				while( entries[u_ind].col < 0 &&
				       u_ind > m_ind) u_ind --;
				if( m_ind == u_ind) {
					break;
				}
				assert( m_ind < u_ind);
				assert( entries[m_ind].col < 0);
				assert( entries[u_ind].col >= 0);
				entries[ m_ind] = entries[ u_ind];
				u_ind --;
			}
			n_entries = new_n_entries;
		}
			
		void copy( Row & row) {
			n_entries = row.n_entries;
			entries = (Entry *) realloc( entries, sizeof( Entry) *
						     n_entries);
			memcpy( entries, row.entries, sizeof( Entry) *
				n_entries);
		}
	};

	double get( long row, long col) const;
	void set( long row, long col, double value);
	void add( long row, long col, double value);
	void reset_row( long row);
	void reset( void);
	void del_entry( long row, long col);
	void copy( SparseMatrix & m);
	void mult( const Vector_double & a, Vector_double & res) const;

	long n_rows, n_cols;
	Row ** rows;

	SparseMatrix( long p_n_rows, long p_n_cols);
	SparseMatrix( const SparseMatrix & m);
	~SparseMatrix();

	void print( char * str, char * fmtstr);
	void expand( long new_n_rows, long new_n_cols);
};

#endif
