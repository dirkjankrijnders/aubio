typedef long double  lsmp_t;

struct c_complex
  { lsmp_t re, im;
  };

struct complex
  { lsmp_t re, im;
    complex(lsmp_t r, lsmp_t i = 0.0) { re = r; im = i; }
    complex() { }					/* uninitialized complex */
    complex(c_complex z) { re = z.re; im = z.im; }	/* init from denotation */
  };

extern complex csqrt(complex), cexp(complex), expj(lsmp_t);	    /* from complex.C */
extern complex evaluate(complex[], int, complex[], int, complex);   /* from complex.C */

inline lsmp_t hypot(complex z) { return ::hypot(z.im, z.re); }
inline lsmp_t atan2(complex z) { return ::atan2(z.im, z.re); }

inline complex cconj(complex z)
  { z.im = -z.im;
    return z;
  }

inline complex operator * (lsmp_t a, complex z)
  { z.re *= a; z.im *= a;
    return z;
  }

inline complex operator / (complex z, lsmp_t a)
  { z.re /= a; z.im /= a;
    return z;
  }

inline void operator /= (complex &z, lsmp_t a)
  { z = z / a;
  }

extern complex operator * (complex, complex);
extern complex operator / (complex, complex);

inline complex operator + (complex z1, complex z2)
  { z1.re += z2.re;
    z1.im += z2.im;
    return z1;
  }

inline complex operator - (complex z1, complex z2)
  { z1.re -= z2.re;
    z1.im -= z2.im;
    return z1;
  }

inline complex operator - (complex z)
  { return 0.0 - z;
  }

inline bool operator == (complex z1, complex z2)
  { return (z1.re == z2.re) && (z1.im == z2.im);
  }

inline complex sqr(complex z)
  { return z*z;
  }

