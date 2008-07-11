/**
 * @file   rta_ifft_mex.c
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Thu May 15 21:12:57 2008
 * 
 * @brief  rta_ifft mex function
 * 
 * Copyright (C) 2008 by IRCAM-Centre Georges Pompidou, Paris, France.
 *
 */

#include "mex.h"
#include "rta_fft.h"
#include "rta_fft_setup_mex.h"

/** 
 * Entry point to C MEX-file
 * 
 * @param nlhs number of expected outputs
 * @param plhs array of pointers to the expected output 
 * @param nrhs number of inputs
 * @param prhs array of pointers to the inputs
 */
void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
  /* matlab inputs */
  rta_ptr_t setup_address;
  double * input;  
  unsigned int input_size;
  unsigned int input_m;
  rta_real_t scale;

  /* rta inputs */
  rta_fft_setup_mex_t * fft_setup_mex;
  rta_real_t * real_fft_complex; /* for conversion */
  rta_real_t * real_input; 

  /* matlab outputs */
  rta_real_t * fft_re; /* for conversion */
  rta_real_t * fft_im; /* for conversion */
  unsigned int fft_m;
  unsigned int fft_n;

  /* other */
  int i,r;

  /* check proper input and output */
  if(nrhs < 2)
  {
    mexErrMsgTxt("Two input required.");
  }
  else if(nlhs > 1)
  {
    mexErrMsgTxt("Too many output arguments.");
  }

  input = mxGetData(prhs[0]); 
  input_size = mxGetNumberOfElements(prhs[0]);    
  input_m=mxGetM(prhs[0]);

  setup_address = mxGetScalar(prhs[1]);
  fft_setup_mex = (rta_fft_setup_mex_t *) setup_address;

  if(input_size > fft_setup_mex->fft_size)
  {
    mexErrMsgTxt("FFT size must be at least input size.");
  }

#if (RTA_REAL_TYPE == RTA_FLOAT_TYPE)
  if(mxGetClassID(prhs[0]) != mxSINGLE_CLASS)
  {
    /* input float precision conversion */
    real_input = (rta_real_t *) mxMalloc(input_size * sizeof(rta_real_t)); 
    for (i=0; i<input_size ;i++)
    {
      real_input[i] = (rta_real_t) input[i]; 
    }
  }
  else
  {
    /* no conversion for matlab single */
    real_input = (rta_real_t *) input;
#else
    /* no conversion for rta double */
    real_input = (rta_real_t *) input;
#endif

#if (RTA_REAL_TYPE == RTA_FLOAT_TYPE)
  }
#endif

  rta_fft_execute(fft_setup_mex->fft, real_input, 
                  input_size, fft_setup_mex->fft_setup);

  /* output results */
  if(input_m == 1)
  {
    fft_m = 1;
    fft_n = fft_setup_mex->fft_size;
  }
  else
  {
    fft_m = fft_setup_mex->fft_size;
    fft_n = 1;
  }

  plhs[0] = mxCreateNumericMatrix(fft_m, fft_n, RTA_MEX_REAL_TYPE, mxCOMPLEX);
  fft_re = (rta_real_t *) mxGetPr(plhs[0]);
  fft_im = (rta_real_t *) mxGetPi(plhs[0]);
  real_fft_complex = (rta_real_t *) fft_setup_mex->fft;

  /* only first half spectrum is computed for a real FFT */
  for(i=0, r=0; i<fft_setup_mex->fft_size/2; i++, r+=2)
  {
    fft_re[i] = real_fft_complex[r];
    fft_im[i] = real_fft_complex[r+1];
  }
  
  /* nyquist */
  fft_re[i] = fft_setup_mex->nyquist;
  fft_im[i] = 0.;

  /* build the half upper spectrum */
  for(i++, r-=2; i<fft_setup_mex->fft_size; i++, r-=2)
  {
    fft_re[i] = real_fft_complex[r];
    fft_im[i] = -real_fft_complex[r+1];
  }
  
#if (RTA_REAL_TYPE != RTA_DOUBLE_TYPE)
  if(mxGetClassID(prhs[0]) != mxSINGLE_CLASS)
  {
    /* free mem of tmp vec for float precision conversion */
    mxFree(real_input);
  }
#endif

  return;
}
