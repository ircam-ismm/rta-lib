/**
 * @file   rta_dct_apply_mex.c
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Thu May 15 21:12:57 2008
 * 
 * @brief  rta_window_apply mex function
 * 
 * Copyright (C) 2008 by IRCAM-Centre Georges Pompidou, Paris, France.
 *
 */

#include "mex.h" 
#include "rta_dct.h"


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

  double * weights;
  unsigned int weights_size;
  unsigned int weights_m;

  /* rta inputs */
  rta_real_t * real_input;
  rta_real_t * real_weights;

  /* rta outputs */

  /* matlab outputs */
  rta_real_t * output; 
  unsigned int output_size; /* weights_n */
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
  input_m = mxGetM(prhs[0]);

  weights = mxGetData(prhs[1]); 
  weights_size = mxGetNumberOfElements(prhs[1]);    
  weights_m = mxGetM(prhs[1]);
  output_size = mxGetN(prhs[1]);

  if(weights_m != input_m)
  {
    mexErrMsgTxt("Input size does not fit the weights size.");
  }

  if(input_m == 1)
  {
    output_m = 1;
    output_n = output_size;
  }
  else
  {
    output_m = output_size;
    output_n = 1;
  }


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
    /* weights float precision conversion */
    real_weights = mxMalloc( weights_size * sizeof(rta_real_t)); 
    for (i=0; i<weights_size ;i++)
    {
      real_weights[i] = (rta_real_t) weights[i]; 
    }
  }
  else
  {
    /* no conversion for matlab single */
    real_weights = (rta_real_t *) weights;
#else
    /* no conversion for rta double */
    real_weights = (rta_real_t *) weights;
#endif

#if (RTA_REAL_TYPE == RTA_FLOAT_TYPE)
  }
#endif

  plhs[0] = mxCreateNumericMatrix(output_m, output_n, RTA_MEX_REAL_TYPE, mxREAL);
  output = mxGetData(plhs[0]);

  rta_dct(output, real_input, real_weights, input_size, output_size);

#if (RTA_REAL_TYPE != RTA_DOUBLE_TYPE)
  if(mxGetClassID(prhs[0]) != mxSINGLE_CLASS)
  {
    /* free mem of tmp vec for float precision conversion */
    mxFree(real_input);
  }

  if(mxGetClassID(prhs[1]) != mxSINGLE_CLASS)
  {
    /* free mem of tmp vec for float precision conversion */
    mxFree(real_weights);
  }
#endif

  return;
}
