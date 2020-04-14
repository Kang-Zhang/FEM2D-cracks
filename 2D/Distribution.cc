#include <stdio.h>
#include <stdlib.h>

#include "Distribution.hh"

double Distribution::get_normal_rnd( double mean, double var) {
	static int has_extra = 0;
	static double gset;
	double fac,rsq,v1,v2;
	double result;
	if ( has_extra == 0) {
		// We don't have an extra deviate handy, so compute
		// two deviates
		do {
			// pick two uniform numbers in the square
			// extending from -1 to +1 in each
			// direction,
			v1=2.0*drand48()-1.0;
			v2=2.0*drand48()-1.0;
			rsq=v1*v1+v2*v2;
			// see if they are in the unit circle, and
			// if they are not, try again.
		} while (rsq >= 1.0 || rsq < 1E-30);
		// Now make the Box-Muller transformation to get
		// two normal deviates. Return one and save the
		// other for next time.
		fac=sqrt(-2.0*log(rsq)/rsq);
		gset=v1*fac;
		has_extra = 1;
		result = (v2*fac)*var + mean;
	} else {
		// We have an extra deviate handy,
		has_extra = 0;
		result = (gset)*var + mean;
	}
	return result;
}


double Distribution::get_poisson_rnd( double xm) {
	// Returns as a oating-point number an integer value that
	// is a random deviate drawn from a Poisson distribution of
	// mean xm, usingran1(idum) as a source of uniform random
	// deviates.
	static double sq = 0;
	static double alxm = 0;
	static double g = 0;
	static double oldm = -1;
	// oldm is a ag for whether xm has changed since
	// last call.
	double em,t,y;
	if (xm < 12.0) {
		// Use direct method.
		if (xm != oldm) {
			oldm=xm;
			// If xm is new, compute the exponential.
			g=exp(-xm); 
		}
		em = -1;
		t=1.0;
		do {
			// Instead of adding exponential deviates
			// it is equivalent to multiply uniform
			// deviates. We never actually have to take
			// the log, merely compare to the
			// pre-computed exponential.
			++em;
			t *= drand48();
		} while (t > g);
	} else { 
		// Use rejection method.
		if (xm != oldm) {
			// If xm has changed since the last call,
			// then precompute some functions that
			// occur below.
			oldm=xm;
			sq=sqrt(2.0*xm);
			alxm=log(xm);
			g=xm*alxm-log(xm+1.0);
			// The function gammln is the natural log
			// of the gamma function, as given in x
			// 6:1.
		}
		do {
			do {
				// y is a deviate from a Lorentzian
				// comparison function.
				y=tan(M_PI*drand48());
				em=sq*y+xm;
				// em is y, shifted and scaled.
			} while (em < 0.0);
			// Reject if in regime of zero probability.
			em=floor(em);
			// The trick for integer-valued
			// distributions.
			t=0.9*(1.0+y*y)*
				exp(em*alxm-log(em+1.0)-g);
			// The ratio of the desired distribution to
			// the comparison function; we accept or
			// reject by comparing it to another
			// uniform deviate. The factor 0.9 is
			// chosen so that t never exceeds 1.
		} while ( drand48() > t);
	}
	return em;
}
