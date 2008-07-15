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
  unsigned int fft_size;
  unsigned int fft_m;
  double * fft_re; 
  double * fft_im; /* for conversion */
  mxArray * mx_input_im = NULL;

  /* rta inputs */
  rta_fft_setup_mex_t * fft_setup_mex;
  rta_real_t * real_fft_complex; /* for conversion */
  rta_real_t * real_fft_re;
  rta_real_t * real_fft_im;

  /* matlab outputs */

  rta_real_t * output;
  unsigned int output_m;
  unsigned int output_n;

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

  fft_re = mxGetPr(prhs[0]);
  fft_size = mxGetNumberOfElements(prhs[0]);
  fft_m=mxGetM(prhs[0]);

  /* input arguments */
  if(mxIsComplex(prhs[0]))
  {
    fft_im = mxGetPi(prhs[0]);
  }
  else
  {
    /* create vector of zeroes for the missing complex part */
    mx_input_im = mxCreateNumericMatrix(
      1, fft_size, mxGetClassID(prhs[0]), mxREAL);
    fft_im = mxGetData(mx_input_im);
  }

  setup_address = mxGetScalar(prhs[1]);
  fft_setup_mex = (rta_fft_setup_mex_t *) setup_address;

  if(fft_size > fft_setup_mex->fft_size)
  {
    mexErrMsgTxt("Input FFT size must no be greater than that of the setup.");
  }

  real_fft_complex = (rta_real_t *) fft_setup_mex->fft;

/* process only the half lower part of the spectrum */
#if (RTA_REAL_TYPE == RTA_FLOAT_TYPE)
  if(mxGetClassID(prhs[0]) != mxSINGLE_CLASS)
  {
    for (i=0, r=0; i<fft_size/2 ;i++, r+=2)
    {
      real_fft_complex[r] = (rta_real_t) fft_re[i];
      real_fft_complex[r+1] = (rta_real_t) fft_im[i];
    }
    fft_setup_mex->nyquist = (rta_real_t) fft_re[fft_size/2];
  }
  else
  {
    /* no conversion for matlab single */
    real_fft_re = (rta_real_t *) fft_re;
    real_fft_im = (rta_real_t *) fft_im;
    for (i=0, r=0; i<fft_size/2 ;i++, r+=2)
    {
      real_fft_complex[r] = real_fft_re[i];
      real_fft_complex[r+1] = real_fft_im[i];
    }
    fft_setup_mex->nyquist = real_fft_re[fft_size/2];
#else
    /* no conversion for rta double */
    real_fft_re = (rta_real_t *) fft_re;
    real_fft_im = (rta_real_t *) fft_im;
    for (i=0, r=0; i<fft_size/2 ;i++, r+=2)
    {
      real_fft_complex[r] = real_fft_re[i];
      real_fft_complex[r+1] = real_fft_im[i];
    }
    fft_setup_mex->nyquist = real_fft_re[fft_size/2];
#endif

#if (RTA_REAL_TYPE == RTA_FLOAT_TYPE)
  }
#endif

  /* output results */
  if(fft_m == 1)
  {
    output_m = 1;
    output_n = fft_setup_mex->fft_size;
  }
  else
  {
    output_m = fft_setup_mex->fft_size;
    output_n = 1;
  }

  plhs[0] = mxCreateNumericMatrix(output_m, output_n,
                                  RTA_MEX_REAL_TYPE, mxREAL);
  output = mxGetData(plhs[0]);

  rta_fft_execute(output, fft_setup_mex->fft, 
                  fft_size/2, fft_setup_mex->fft_setup);

  if(mx_input_im != NULL)
  {
    mxDestroyArray(mx_input_im);
  }

  return;
}
