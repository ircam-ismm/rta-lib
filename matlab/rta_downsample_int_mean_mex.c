/**
 * @file   rta_mex_downsample_int_mean.c
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Thu May 15 21:12:57 2008
 * 
 * @brief  rta_yin mex function for simple floating-point precision
 * 
 * Copyright (C) 2008 by IRCAM-Centre Georges Pompidou, Paris, France.
 *
 */

#include "mex.h"
#include "rta_resample.h"


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
  unsigned int factor;

  /* rta inputs */
  rta_real_t * real_input; 
  unsigned int input_size;

  /* rta outputs */
  rta_real_t * output;

  /* matlab outputs */

  /* other */
  int i;

  /* input arguments */
  
  input = mxGetData(prhs[0]); 
  input_size = mxGetNumberOfElements(prhs[0]);    

  factor = (unsigned int) mxGetScalar(prhs[1]);

  /* output results */
  plhs[0] = mxCreateNumericMatrix(1, input_size / factor,
                                  RTA_MEX_REAL_TYPE, mxREAL);
  output = mxGetData(plhs[0]);
  
#if (RTA_REAL_TYPE != RTA_DOUBLE_TYPE)
  /* input float precision conversion */
  real_input = mxMalloc( input_size * sizeof(rta_real_t)); 
  for (i=0; i<input_size ;i++)
  {
    real_input[i] = (rta_real_t) input[i]; 
  }
#else
  real_input = input;
#endif

  /* computation */
  rta_downsample_int_mean(output, real_input, input_size, factor);

#if (RTA_REAL_TYPE != RTA_DOUBLE_TYPE)
  /* free mem of tmp vec for float precision conversion */
  mxFree(real_input);
#endif

  return;
}
