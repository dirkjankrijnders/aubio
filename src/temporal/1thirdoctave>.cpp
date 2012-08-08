/*
  Copyright (C) 2003-2009 Paul Brossier <piem@aubio.org>

  This file is part of aubio.

  aubio is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  aubio is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with aubio.  If not, see <http://www.gnu.org/licenses/>.

*/

#include <math.h>
#include "complex.h"

#include "aubio_priv.h"
#include "types.h"
#include "fvec.h"
#include "lvec.h"
#include "temporal/filter.h"
#include "temporal/1thirdoctave.h"

#define TWOPI	    (2.0 * PI)
#define MAXPZ   10
#define EPS	    1e-10
#define _DEBUG 1

struct pzrep
{ complex poles[MAXPZ], zeros[MAXPZ];
    int numpoles, numzeros;
};



void multin(complex w, int npz, complex coeffs[])
{ /* multiply factor (z-w) into coeffs */
    complex nw = -w;
    for (int i = npz; i >= 1; i--) coeffs[i] = (nw * coeffs[i]) + coeffs[i-1];
    coeffs[0] = nw * coeffs[0];
}

static void expand(complex pz[], int npz, complex coeffs[])
{ /* compute product of poles or zeros as a polynomial of z */
    int i;
    coeffs[0] = 1.0;
    for (i=0; i < npz; i++) coeffs[i+1] = 0.0;
    for (i=0; i < npz; i++) multin(pz[i], npz, coeffs);
    /* check computed coeffs of z^k are all real */
    for (i=0; i < npz+1; i++)
    { if (fabs(coeffs[i].im) > EPS)
    { fprintf(stderr, "mkfilter: coeff of z^%d is not real; poles/zeros are not complex conjugates\n", i);
	    exit(1);
    }
    }
}

static void prcomplex(complex z)
{ printf("%14.10f + j %14.10f", z.re, z.im);
}

static void printpz(complex *pzvec, int num)
{ int n1 = 0;
    while (n1 < num)
    { putchar('\t'); prcomplex(pzvec[n1]);
        int n2 = n1+1;
        while (n2 < num && pzvec[n2] == pzvec[n1]) n2++;
        if (n2-n1 > 1) printf("\t%d times", n2-n1);
        putchar('\n');
        n1 = n2;
    }
    putchar('\n');
}

static void printgain(const char *str, complex gain)
{ 
    double r = hypot(gain);
    printf("gain at %s:   mag = %15.9e", str, r);
    if (r > EPS) printf("   phase = %14.10f pi", atan2(gain) / PI);
    putchar('\n');
}

static complex blt(complex pz)
{ return (2.0 + pz) / (2.0 - pz);
}

uint_t
aubio_filter_set_1thirdoctave (aubio_filter_t * f, uint_t samplerate, double centerfrequency)
{
  aubio_filter_set_samplerate (f, samplerate);
  lvec_t *bs = aubio_filter_get_feedforward (f);
  lvec_t *as = aubio_filter_get_feedback (f);
  lsmp_t *b = bs->data, *a = as->data;
  uint_t order = aubio_filter_get_order (f);
    
    if (order != 7) {
        AUBIO_ERROR ("order of 1/3 octave filter must be 7, not %d\n", order);
        return 1;
    }

    if (centerfrequency > 0.88*((double)samplerate/2.))
        AUBIO_ERROR("Design not possible. Check frequencies.");
    
    pzrep splane, zplane;
    complex dc_gain, fc_gain, hf_gain;
    u_int32_t polemask;
  
    double *xcoeffs = b, *ycoeffs = a;
    order = 3; // 3rd order Butterworth
    
    
    double f1 = centerfrequency/(exp2f(1./6.))/2; // Cut-off frequencies
    double f2 = centerfrequency*(exp2f(1./6.))/2; 

    /* select coefficients according to sampling frequency */
    double fs= samplerate/2;
    
    double raw_alpha1 = f1/fs;
    double raw_alpha2 = f2/fs;
#ifdef _DEBUG
    printf("samplerate    =   %f\n", fs);
    printf("raw alpha1    =   %f (%f)\n", raw_alpha1, f1);
    printf("raw alpha2    =   %f (%f)\n", raw_alpha2, f2);
#endif    
    splane.numpoles = 0;
    polemask = ~0;
    int i;
    complex z;
    for (i = 0; i < 2*order; i++)
    { 
        double theta = (order & 1) ? (i*PI) / order : ((i+0.5)*PI) / order;
        
       z = (expj(theta));
        if (z.re < 0.0)
        { 
            if (polemask & 1) splane.poles[splane.numpoles++] = z;
            polemask >>= 1;
        }
#ifdef _DEBUG
       printf("i : %i; splane.numpoles : %i, theta: %f, z: ", i, splane.numpoles, theta);
        prcomplex(z);
        printf("\n");
#endif    
    }
    
    //prewarp alpha
    double warped_alpha1 = tan(PI * raw_alpha1) / PI;
    double warped_alpha2 = tan(PI * raw_alpha2) / PI;
#ifdef _DEBUG
    printf("warped alpha1    =   %f\n", warped_alpha1);
    printf("warped alpha2    =   %f\n", warped_alpha2);
#endif    
    
    double w1 = TWOPI * warped_alpha1;
    double w2 = TWOPI * warped_alpha2;

    // Bandpass filter
    double w0 = sqrt(w1*w2), bw = w2-w1;
    for (i=0; i < splane.numpoles; i++)
    { 
        complex hba = 0.5 * (splane.poles[i] * bw);
        complex temp = csqrt(1.0 - sqr(w0 / hba));
        splane.poles[i] = hba * (1.0 + temp);
        splane.poles[splane.numpoles+i] = hba * (1.0 - temp);
    }
    
    for (i=0; i < splane.numpoles; i++) 
        splane.zeros[i] = 0.0;	 /* also N zeros at (0,0) */

    splane.numzeros = splane.numpoles;
    splane.numpoles *= 2;
    
    // Bi-linear transform
    zplane.numpoles = splane.numpoles;
    zplane.numzeros = splane.numzeros;
    for (i=0; i < zplane.numpoles; i++) zplane.poles[i] = blt(splane.poles[i]);
    for (i=0; i < zplane.numzeros; i++) zplane.zeros[i] = blt(splane.zeros[i]);
    while (zplane.numzeros < zplane.numpoles) zplane.zeros[zplane.numzeros++] = -1.0;
    
    complex topcoeffs[MAXPZ+1], botcoeffs[MAXPZ+1]; 
    
    expand(zplane.zeros, zplane.numzeros, topcoeffs);
    expand(zplane.poles, zplane.numpoles, botcoeffs);
    
    dc_gain = evaluate(topcoeffs, zplane.numzeros, botcoeffs, zplane.numpoles, 1.0);

    double theta = TWOPI * 0.5 * (raw_alpha1 + raw_alpha2); /* "jwT" for centre freq. */
    
    fc_gain = evaluate(topcoeffs, zplane.numzeros, botcoeffs, zplane.numpoles, expj(theta));
    hf_gain = evaluate(topcoeffs, zplane.numzeros, botcoeffs, zplane.numpoles, -1.0);
    
    for (i = 0; i <= zplane.numzeros; i++) xcoeffs[zplane.numzeros-i] = +(topcoeffs[i].re / botcoeffs[zplane.numpoles].re);
    for (i = 0; i <= zplane.numpoles; i++) ycoeffs[zplane.numpoles-i] = +(botcoeffs[i].re / botcoeffs[zplane.numpoles].re);
    
    //for (i = 0; i <= zplane.numzeros; i++) xcoeffs[zplane.numzeros-i] = xcoeffs[zplane.numzeros-i] / fc_gain.re;
//    aubio_filter_set_gain(f, 1);
    aubio_filter_set_gain(f, fc_gain.re);
    //xcoeffs[0] = xcoeffs[0] /fc_gain.re;
#ifdef _DEBUG

    printgain("dc    ", dc_gain);
    printgain("centre", fc_gain);
    printgain("hf    ", hf_gain);
    
    printf("S-plane zeros:\n");
    printpz(splane.zeros, splane.numzeros);
    printf("S-plane poles:\n");
    printpz(splane.poles, splane.numpoles);
    
    printf("Z-plane zeros:\n");
    printpz(zplane.zeros, zplane.numzeros);
    printf("Z-plane poles:\n");
    printpz(zplane.poles, zplane.numpoles);
    
    printf("Recurrence relation:\n");
    printf("y[n] = ");
    for (i = 0; i < zplane.numzeros+1; i++)
    { 
        if (i > 0) printf("     + ");
        double x = xcoeffs[i];
        double f = fmod(fabs(x), 1.0);
        const char *fmt = (f < EPS || f > 1.0-EPS) ? "%3g" : "%14.10f";
        putchar('('); printf(fmt, x); printf(" * x[n-%2d])\n", zplane.numzeros-i);
    }
    putchar('\n');
    for (i = 0; i < zplane.numpoles; i++)
    { 
        printf("     + (%14.10f * y[n-%2d])\n", ycoeffs[i], zplane.numpoles-i);
    }
    putchar('\n');
#endif    
    
  return 0;
}


aubio_filter_t *
new_aubio_filter_1thirdoctave (uint_t samplerate, double centerfrequency)
{
  aubio_filter_t *f = new_aubio_filter (7);
  aubio_filter_set_1thirdoctave (f, samplerate, centerfrequency);
    lvec_t *bs = aubio_filter_get_feedforward (f);
    lvec_t *as = aubio_filter_get_feedback (f);
    lsmp_t *b = bs->data, *a = as->data;
    printf("b[2] = %f", b[2]);
  return f;
}
