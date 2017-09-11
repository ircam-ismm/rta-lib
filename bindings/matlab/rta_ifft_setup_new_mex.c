/**
 * @file   rta_fft_setup_new_mex.c
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Thu May 15 21:12:57 2008
 * 
 * @brief  rta_fft mex initialisation function (for real fft)
 * 
 * Copyright (C) 2008 by IRCAM-Centre Georges Pompidou, Paris, France.
 *
 */

#include "mex.h"
#include "rta_fft.h"
#include "rta_fft_setup_mex.h"
#include "rta_stdlib.h"


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
  double * fft_re;
  double * fft_im;
  unsigned int fft_size;
  rta_real_t scale = 0.; /* default for 1/fft_size normalisation */
  mxArray * mx_output;
  mxArray * mx_input_im = NULL;
  
  /* rta inputs */
  rta_fft_setup_mex_t * fft_setup_mex;
  rta_real_t * real_fft_re;
  rta_real_t * real_fft_im;

  rta_real_t * real_fft_complex; /* for conversion */
  rta_real_t * output; 

  /* matlab outputs */
  rta_ptr_t * setup_address;

  /* other */
  int i,r;

  /* check proper input and output */
  if(nlhs > 1)
  {
    mexErrMsgTxt("Too many output arguments.");
  }

  if(nrhs < 1)
  {
    mexErrMsgTxt("One input required.");
  }

  fft_re = mxGetPr(prhs[0]);
  fft_size = mxGetNumberOfElements(prhs[0]);

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
  
  fft_setup_mex = (rta_fft_setup_mex_t *)
    rta_malloc(sizeof(rta_fft_setup_mex_t));

  if(nrhs > 2)
  {
    scale = mxGetScalar(prhs[2]);
  }

  if(scale == 0.) /* default value */
  {
    fft_setup_mex->scale = 1. / fft_size;
  }
  else
  {
    fft_setup_mex->scale = scale;
  }

  fft_setup_mex->fft_size = fft_size;

  fft_setup_mex->fft = (rta_complex_t *) rta_malloc(
    fft_size/2 * sizeof(rta_complex_t));

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
#else
    /* no conversion for rta double */
    real_fft_re = (rta_real_t *) fft_re;
    real_fft_im = (rta_real_t *) fft_im;
    for (i=0, r=0; i<fft_size/2 ;i++, r+=2)
    {
      real_fft_complex[r] = real_fft_re[i];
      real_fft_complex[r+1] = real_fft_im[i];
    }
#endif

#if (RTA_REAL_TYPE == RTA_FLOAT_TYPE)
  }
#endif
  
/*   for(r=0; r<fft_size; r+=2) */
/*   { */
/*     printf("fft[%d] = %lf + %lf i\n", */
/*            r/2, real_fft_complex[r], real_fft_complex[r+1]); */
/*   } */

  /* add the nyquist middle value */
  fft_setup_mex->nyquist = real_fft_re[fft_size/2];

  mx_output = mxCreateNumericMatrix(1, fft_size, RTA_MEX_REAL_TYPE, mxREAL);
  output = mxGetData(mx_output);

  if(rta_fft_real_setup_new(&(fft_setup_mex->fft_setup), 
                            rta_fft_complex_to_real_1d,
                            &(fft_setup_mex->scale),
                            fft_setup_mex->fft, fft_size/2,
                            output, fft_size,
                            &(fft_setup_mex->nyquist) ) )
  {
    plhs[0] = mxCreateNumericMatrix(1, 1, RTA_MEX_PTR_TYPE, mxREAL);
    setup_address = mxGetData(plhs[0]);
    *setup_address = (rta_ptr_t) fft_setup_mex;
  }
  else
  {
    printf("rta_fft setup failed.\n");
  }
  
  if(mx_input_im != NULL)
  {
    mxDestroyArray(mx_input_im);
  }

  mxDestroyArray(mx_output);

  return;
}
