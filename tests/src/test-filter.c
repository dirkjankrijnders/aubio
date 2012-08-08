
#define HAVE_AUBIO_DOUBLE 1

#include <aubio.h>
#include <stdio.h>

int
main (void)
{
  /* allocate some memory */
  uint_t win_s = 32;            /* window size */
  fvec_t *in = new_fvec (win_s);      /* input buffer */
  fvec_t *out = new_fvec (win_s);     /* input buffer */
    smpl_t fc = (smpl_t)13500.;
    aubio_filter_t *o = new_aubio_filter_1thirdoctave(48000, fc);
  print_aubio_filter( o);
  in->data[2] = 1;
  aubio_filter_do_outplace(o, in, out);
  int i;
  for (i = 0; i < win_s; i++) {
    printf(" %f\n", out->data[i]);
  }
  
  del_fvec (in);
  del_fvec (out);
  aubio_cleanup ();


  return 0;
}
