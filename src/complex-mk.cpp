/* mkfilter -- given n, compute recurrence relation
   to implement Butterworth, Bessel or Chebyshev filter of order n
   A.J. Fisher, University of York   <fisher@minster.york.ac.uk>
   September 1992 */

/* Routines for complex arithmetic */

#include <math.h>

//#include "mkfilter.h"
#include "complex.h"

static complex eval(complex[], int, complex);
static lsmp_t Xsqrt(lsmp_t);


complex evaluate(complex topco[], int nz, complex botco[], int np, complex z)
  { /* evaluate response, substituting for z */
    return eval(topco, nz, z) / eval(botco, np, z);
  }

static complex eval(complex coeffs[], int npz, complex z)
  { /* evaluate polynomial in z, substituting for z */
    complex sum = complex(0.0);
    for (int i = npz; i >= 0; i--) sum = (sum * z) + coeffs[i];
    return sum;
  }

complex csqrt(complex x)
  { lsmp_t r = hypot(x);
    complex z = complex(Xsqrt(0.5 * (r + x.re)),
			Xsqrt(0.5 * (r - x.re)));
    if (x.im < 0.0) z.im = -z.im;
    return z;
  }

static lsmp_t Xsqrt(lsmp_t x)
  { /* because of deficiencies in hypot on Sparc, it's possible for arg of Xsqrt to be small and -ve,
       which logically it can't be (since r >= |x.re|).	 Take it as 0. */
    return (x >= 0.0) ? sqrt(x) : 0.0;
  }

complex cexp(complex z)
  { return exp(z.re) * expj(z.im);
  }

complex expj(lsmp_t theta)
  { return complex(cos(theta), sin(theta));
  }

complex operator * (complex z1, complex z2)
  { return complex(z1.re*z2.re - z1.im*z2.im,
		   z1.re*z2.im + z1.im*z2.re);
  }

complex operator / (complex z1, complex z2)
  { lsmp_t mag = (z2.re * z2.re) + (z2.im * z2.im);
    return complex (((z1.re * z2.re) + (z1.im * z2.im)) / mag,
		    ((z1.im * z2.re) - (z1.re * z2.im)) / mag);
  }

