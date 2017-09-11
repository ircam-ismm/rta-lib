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
  double * input;  
  unsigned int input_size;

  unsigned int fft_size;
  
  /* rta inputs */
  rta_fft_setup_mex_t * fft_setup_mex;
  rta_real_t * real_input; 

  /* matlab outputs */
  rta_ptr_t * setup_address;

  /* other */
  int i;

  /* check proper input and output */
  if(nlhs > 1)
  {
    mexErrMsgTxt("Too many output arguments.");
  }

  if(nrhs < 2)
  {
    mexErrMsgTxt("Two inputs required.");
  }

  /* input arguments */
  input = mxGetData(prhs[0]); 
  input_size = mxGetNumberOfElements(prhs[0]);
  fft_size = mxGetScalar(prhs[1]);

  if(input_size > fft_size)
  {
    mexErrMsgTxt("FFT size must be at least input size.");
  }
  
  fft_setup_mex = (rta_fft_setup_mex_t *) rta_malloc(
    sizeof(rta_fft_setup_mex_t));

  if(nrhs > 2)
  {
    fft_setup_mex->scale = mxGetScalar(prhs[2]);
  }
  else
  {
    fft_setup_mex->scale = 1.;
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

  fft_setup_mex->fft_size = fft_size;
  fft_setup_mex->fft = (rta_complex_t *) rta_malloc(
    fft_size/2 * sizeof(rta_complex_t));
  
  if(rta_fft_real_setup_new(&(fft_setup_mex->fft_setup), 
                            rta_fft_real_to_complex_1d,
                            &(fft_setup_mex->scale),
                            real_input, input_size,
                            fft_setup_mex->fft, fft_size, 
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

#if (RTA_REAL_TYPE != RTA_DOUBLE_TYPE)
  if(mxGetClassID(prhs[0]) != mxSINGLE_CLASS)
  {
    /* free mem of tmp vec for float precision conversion */
    mxFree(real_input);
  }
#endif

  return;
}
