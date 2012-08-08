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

static u_int32_t polemask;

struct pzrep
{ complex poles[MAXPZ], zeros[MAXPZ];
    int numpoles, numzeros;
};

static pzrep splane, zplane;
static complex dc_gain, fc_gain, hf_gain;
double xcoeffs[MAXPZ+1], ycoeffs[MAXPZ+1];

void multin(complex w, int npz, complex coeffs[])
{ /* multiply factor (z-w) into coeffs */
    complex nw = -w;
    for (int i = npz; i >= 1; i--) coeffs[i] = (nw * coeffs[i]) + coeffs[i-1];
    coeffs[0] = nw * coeffs[0];
}

static void choosepole(complex z)
{ 
    if (z.re < 0.0)
    { 
        if (polemask & 1) splane.poles[splane.numpoles++] = z;
            polemask >>= 1;
    }
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

complex eval(complex coeffs[], int npz, complex z)
{ /* evaluate polynomial in z, substituting for z */
    complex sum = complex(0.0);
    for (int i = npz; i >= 0; i--) sum = (sum * z) + coeffs[i];
    return sum;
}

complex evaluate(complex topco[], int nz, complex botco[], int np, complex z)
{ /* evaluate response, substituting for z */
    return eval(topco, nz, z) / eval(botco, np, z);
}

uint_t
aubio_filter_set_1thirdoctave (aubio_filter_t * f, uint_t samplerate, double centerfrequency)
{
  aubio_filter_set_samplerate (f, samplerate);
  lvec_t *bs = aubio_filter_get_feedforward (f);
  lvec_t *as = aubio_filter_get_feedback (f);
  lsmp_t *b = bs->data, *a = as->data;
  uint_t order = aubio_filter_get_order (f);

//    uint_t order = 3;
    
  if (order != 3) {
    AUBIO_ERROR ("order of A-weighting filter must be 7, not %d\n", order);
    return 1;
  }

  if (centerfrequency > 0.88*((double)samplerate/2.))
    AUBIO_ERROR("Design not possible. Check frequencies.");
  
    //pi = 3.14159265358979;
    double f1 = centerfrequency/(2^(1/6)); // Cut-off frequencies
    double f2 = centerfrequency*(2^(1/6)); 
    /*double Qr = centerfrequency/(f2-f1);   // Bandwidth
    double Qd = (pi/2/N)/(sin(pi/2/N))*Qr; 
    double alpha = (1 + sqrt(1+4*Qd^2))/2/Qd; 
    double W1 = centerfrequency/(samplerate/2)/alpha; // Warped alpha
    double W2 = centerfrequency/(samplerate/2)*alpha; // Warped alpha
    [B,A] = butter(N,[W1,W2]);*/
  /* select coefficients according to sampling frequency */
    double fs= samplerate;
    
    double raw_alpha1 = f1/fs;
    double raw_alpha2 = f2/fs;
    printf("raw alpha1    =   %f", raw_alpha1);
    printf("raw alpha2    =   %f", raw_alpha2);
    
    splane.numpoles = 0;
    int i;
    for (i = 0; i < 2*order; i++)
    { double theta = (order & 1) ? (i*PI) / order : ((i+0.5)*PI) / order;
        choosepole(expj(theta));
    }
    
    //prewarp alpha
    double warped_alpha1 = tan(PI * raw_alpha1) / PI;
    double warped_alpha2 = tan(PI * raw_alpha2) / PI;
    
    double w1 = TWOPI * warped_alpha1;
    double w2 = TWOPI * warped_alpha2;

    double w0 = sqrt(w1*w2), bw = w2-w1;
    for (i=0; i < splane.numpoles; i++)
    { complex hba = 0.5 * (splane.poles[i] * bw);
        complex temp = csqrt(1.0 - sqr(w0 / hba));
        splane.poles[i] = hba * (1.0 + temp);
        splane.poles[splane.numpoles+i] = hba * (1.0 - temp);
    }
    
    for (i=0; i < splane.numpoles; i++) 
        splane.zeros[i] = 0.0;	 /* also N zeros at (0,0) */

    splane.numzeros = splane.numpoles;
    splane.numpoles *= 2;
    
    complex topcoeffs[MAXPZ+1], botcoeffs[MAXPZ+1]; 
    
    expand(zplane.zeros, zplane.numzeros, topcoeffs);
    expand(zplane.poles, zplane.numpoles, botcoeffs);
    
    dc_gain = evaluate(topcoeffs, zplane.numzeros, botcoeffs, zplane.numpoles, 1.0);

    double theta = TWOPI * 0.5 * (raw_alpha1 + raw_alpha2); /* "jwT" for centre freq. */
    
    fc_gain = evaluate(topcoeffs, zplane.numzeros, botcoeffs, zplane.numpoles, expj(theta));
    hf_gain = evaluate(topcoeffs, zplane.numzeros, botcoeffs, zplane.numpoles, -1.0);
    
    for (i = 0; i <= zplane.numzeros; i++) xcoeffs[i] = +(topcoeffs[i].re / botcoeffs[zplane.numpoles].re);
    for (i = 0; i <= zplane.numpoles; i++) ycoeffs[i] = -(botcoeffs[i].re / botcoeffs[zplane.numpoles].re);

    
  return 0;
}


aubio_filter_t *
new_aubio_filter_1thirdoctave (uint_t samplerate, double centerfrequency)
{
  aubio_filter_t *f = new_aubio_filter (7);
  aubio_filter_set_1thirdoctave (f, samplerate, centerfrequency);
  return f;
}
