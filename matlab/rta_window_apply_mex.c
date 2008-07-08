/**
 * @file   rta_window_apply_mex.c
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Thu May 15 21:12:57 2008
 * 
 * @brief  rta_window_apply mex function
 * 
 * Copyright (C) 2008 by IRCAM-Centre Georges Pompidou, Paris, France.
 *
 */

#include "mex.h" 
#include "rta_window.h"


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
  unsigned int input_m;

  double * window;
  unsigned int window_size;

  /* rta inputs */
  rta_real_t * real_input;
  rta_real_t * real_window;

  /* rta outputs */

  /* matlab outputs */
  rta_real_t * output; /* output size is input_size */
  unsigned int output_m;
  unsigned int output_n;

  /* other */
  int i;

  /* check proper input and output */
  if(nrhs < 2)
  {
    mexErrMsgTxt("Two inputs required.");
  }
  else if(nlhs > 1)
  {
    mexErrMsgTxt("Too many output arguments.");
  }
  
  input = mxGetData(prhs[0]); 
  input_size = mxGetNumberOfElements(prhs[0]);    
  input_m=mxGetM(prhs[0]);

  if(input_m == 1)
  {
    output_m = 1;
    output_n = input_size;
  }
  else
  {
    output_m = input_size;
    output_n = 1;
  }

  window = mxGetData(prhs[1]); 
  window_size = mxGetNumberOfElements(prhs[1]);

#if (RTA_REAL_TYPE == RTA_FLOAT_TYPE)
  if(mxGetClassID(prhs[0]) != mxSINGLE_CLASS)
  {
    /* input float precision conversion */
    real_input = mxMalloc( input_size * sizeof(rta_real_t)); 
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

#if (RTA_REAL_TYPE == RTA_FLOAT_TYPE)
  if(mxGetClassID(prhs[1]) != mxSINGLE_CLASS)
  {
    /* window float precision conversion */
    real_window = mxMalloc( window_size * sizeof(rta_real_t)); 
    for (i=0; i<window_size ;i++)
    {
      real_window[i] = (rta_real_t) window[i]; 
    }
  }
  else
  {
    /* no conversion for matlab single */
    real_window = (rta_real_t *) window;
#else
    /* no conversion for rta double */
    real_window = (rta_real_t *) window;
#endif

#if (RTA_REAL_TYPE == RTA_FLOAT_TYPE)
  }
#endif

  plhs[0] = mxCreateNumericMatrix(output_m, output_n, RTA_MEX_REAL_TYPE, mxREAL);
  output = mxGetData(plhs[0]);

  if(input_size == window_size)
  {
    rta_window_apply(output, input_size, real_input, real_window);
  }
  else
  {
    rta_window_rounded_apply(output, input_size, real_input, real_window, 
                             window_size);
  }

  return;
}
