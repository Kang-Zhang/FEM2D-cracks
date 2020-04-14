#include <cstdio>
#include <cstdlib>
#include <cassert>
#include <malloc.h>
#include <math.h>

double FMAX( double a, double b)
{
	if( a > b) return a ; else return b;
}

long IMIN( long a, long b)
{
	if( a < b) return a; else return b;
}

#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

void nrerror( const char * str)
{
	fprintf( stderr, "nrerror: %s\n", str);
	exit( 0);
}

double ** dmatrix( long nrl, long nrh, long ncl, long nch)
{
	double ** m = (double **) malloc( (nrh-nrl+1) * sizeof( double *));
	assert( m != NULL);
	m -= nrl;

	// allocate rows and set pointers to them
	for( long i = nrl ; i <= nrh ; i ++) {
		m[i] = (double *) malloc( (nch-ncl+1)*sizeof(double));
		assert( m[i] != NULL);
		m[i] -= ncl;
	}

	return m;
}

void free_dmatrix( double ** m, long nrl, long nrh, long ncl, long nch)
{
	for( long i = nrh ; i >= nrl ; i--)
		free( m[i] + ncl);
	free( m+nrl);
}

double * dvector( long nl, long nh)
{
	double * v = (double *) malloc( (nh-nl+1) * sizeof( double));
	assert( v != NULL);
	return v-nl;
}

void free_dvector( double * v, long nl, long nh)
{
	free( v+nl);
}

double SQR( double a) { return a*a; }

double pythag(double a, double b)
// Computes (a^2 + b^2)^(1/2) without destructive underflow or overflow.
{
	double absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+SQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
}

void svdcmp( double **a, int m, int n, double w[], double **v)
// Given a matrix a[1..m][1..n], this routine computes its singular value
// decomposition, A = U W V T. ThematrixU replaces a on output. The
// diagonal matrix of singular values W is output as a vector
// w[1..n]. ThematrixV (not the transpose V T ) is output as v[1..n][1..n].
{
	double pythag(double a, double b);
	int flag,i,its,j,jj,k,l,nm;
	double anorm,c,f,g,h,s,scale,x,y,z,*rv1;

	rv1=dvector(1,n);
	g=scale=anorm=0.0; // Householder reduction to bidiagonal form.
	for (i=1;i<=n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= m) {
			for (k=i;k<=m;k++) scale += fabs(a[k][i]);
			if (scale) {
				for (k=i;k<=m;k++) {
					a[k][i] /= scale;
					s += a[k][i]*a[k][i];
				}
				f=a[i][i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][i]=f-g;
				for (j=l;j<=n;j++) {
					for (s=0.0,k=i;k<=m;k++)
						s += a[k][i]*a[k][j];
					f=s/h;
					for (k=i;k<=m;k++)
						a[k][j] += f*a[k][i];
				}
				for (k=i;k<=m;k++) a[k][i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i <= m && i != n) {
			for (k=l;k<=n;k++) scale += fabs(a[i][k]);
			if (scale) {
				for (k=l;k<=n;k++) {
					a[i][k] /= scale;
					s += a[i][k]*a[i][k];
				}
				f=a[i][l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i][l]=f-g;
				for (k=l; k<=n; k++) rv1[k]=a[i][k]/h;
				for (j=l; j<=m; j++) {
					for (s=0.0,k=l; k<=n; k++)
						s += a[j][k]*a[i][k];
					for (k=l; k<=n; k++)
						a[j][k] += s*rv1[k];
				}
				for (k=l; k<=n; k++) a[i][k] *= scale;
			}
		}
		anorm=FMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n; i>=1; i--) { // Accumulation of right-hand transformations.
		if (i < n) {
			if (g) {
				for (j=l; j<=n; j++) // Double division to
						     // avoid possible
						     // underflow.
					v[j][i]=(a[i][j]/a[i][l])/g;
				for (j=l; j<=n; j++) { 
					for (s=0.0,k=l;k<=n;k++)
						s += a[i][k]*v[k][j];
					for (k=l; k<=n; k++)
						v[k][j] += s*v[k][i];
				}
			}
			for (j=l;j<=n;j++) v[i][j]=v[j][i]=0.0;
		}
		v[i][i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=IMIN(m,n); i>=1; i--) { // Accumulation of left-hand
				       // transformations.
		l=i+1;
		g=w[i];
		for (j=l; j<=n; j++) a[i][j]=0.0;
		if (g) {
			g=1.0/g;
			for (j=l; j<=n; j++) {
				for (s=0.0,k=l;	k<=m; k++)
					s += a[k][i]*a[k][j];
				f=(s/a[i][i])*g;
				for (k=i; k<=m; k++) a[k][j] += f*a[k][i];
			} for (j=i; j<=m; j++) a[j][i] *= g;
		} else for (j=i; j<=m; j++) a[j][i]=0.0;
		++a[i][i];
	} for (k=n; k>=1;k--) { // Diagonalization of the bidiagonal form:
				// Loop over singular values, and over
				// allowed iterations.
		for (its=1; its<=30; its++) {
			flag=1;
			for (l=k;l>=1; l--) { // Test for splitting.
				nm=l-1; // Note that rv1[1] is always zero.
				if ((double)(fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((double)(fabs(w[nm])+anorm) == anorm)
					break;
			} if (flag) {
				c=0.0; // Cancellation of rv1[l], if l > 1.
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((double)(fabs(f)+anorm) == anorm)
						break;
					g=w[i];
					h=pythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=1;j<=m;j++) {
						y=a[j][nm];
						z=a[j][i];
						a[j][nm]=y*c+z*s;
						a[j][i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) { // Convergence.
				if (z < 0.0) { // Singular value is made
					       // nonnegative.
					w[k] = -z;
					for (j=1; j<=n; j++)
						v[j][k] = -v[j][k];
				}
				break;
			} if (its == 30)
				nrerror( "no convergence in 30 "
					 "svdcmp iterations");
			x=w[l];	// Shift from bottom 2-by-2 minor.
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=pythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0; // Next QR transformation:
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=pythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=1;jj<=n;jj++) {
					x=v[jj][j];
					z=v[jj][i];
					v[jj][j]=x*c+z*s;
					v[jj][i]=z*c-x*s;
				}
				z=pythag(f,h);
				w[j]=z;	// Rotation can be arbitrary if z = 0.
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=1; jj<=m; jj++) {
					y=a[jj][j];
					z=a[jj][i];
					a[jj][j]=y*c+z*s;
					a[jj][i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	free_dvector(rv1,1,n);
}

int pinvT(
	double **A, 	/* matrix[0..Nrows-1][0..Ncolumns-1] */
	int Nrows, 	/* number of rows    */
	int Ncolumns	/* number of columns */
	)
	
{
	double **U, **V, *s, max;
	int i, j, k;
	
	U = dmatrix(1,Nrows,1,Ncolumns);
	V = dmatrix(1,Ncolumns,1,Ncolumns);
	s = dvector(1,Ncolumns);
	
	for (i=1;i<=Nrows;i++) {
		for( j = 1 ; j <= Ncolumns ; j ++)
			U[i][j] = A[i-1][j-1];
	}
	
	svdcmp(U, Nrows, Ncolumns, s, V);
	
	/* calculate the pseudoinverse matrix                */
	/* remove the near singular values from the vector s */
	
	max = s[1];
	for( i = 2 ; i <= Ncolumns ; i ++)
		if( max < s[i]) max = s[i];
	
	for (i=1;i<=Ncolumns;i++)
		for (j=1;j<=Nrows;j++) {
			A[j-1][i-1] = 0;
			for (k=1;k<=Ncolumns;k++) {
				if( s[k] >= (max * 1e-6))
					A[j-1][i-1] +=
						V[i][k] * U[j][k] / s[k];
			}
		}
	
	free_dmatrix(U,1,Nrows,1,Ncolumns);
	free_dmatrix(V,1,Ncolumns,1,Ncolumns);
	free_dvector(s,1,Ncolumns);
	return(0);
}

void matrix_print( char * str, double ** m,
		   long nrl, long nrh, long ncl, long nch)
{
	fprintf( stdout, "%s %ldx%ld\n", str, nrh-nrl+1, nch-ncl+1);
	for( long i = nrl ; i <= nrh ; i ++) {
		for( long j = ncl ; j <= nch ; j ++) {
			fprintf( stdout, "%8.3f", m[i][j]);
		}
		fprintf( stdout, "\n");
	}
}


#ifdef _DONT_COMPILE_

int main( void)
{
	// read in the matrix
	long m, n;
	assert( 2 == fscanf( stdin, "%ld %ld", & m, & n));
	assert( m > 0 && n > 0);
	
	// allocate the matrix
	double ** A = dmatrix( 0, m-1, 0, n-1);
	// populate the matrix
	for( long i = 0 ; i < m ; i ++)
		for( long j = 0 ; j < n ; j ++)
			assert( 1 == fscanf( stdin, "%lf", & A[i][j]));

	// print the matrix
	matrix_print( "A=", A, 0, m-1, 0, n-1);

	// calculate the pseudo inverse
	pinvT( A, m, n);

	// print the result
	matrix_print( "A+=", A, 0, m-1, 0, n-1);
}

#endif
