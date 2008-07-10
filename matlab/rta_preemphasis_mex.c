/**
 * @file   rta_preemphasis_mex.c
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Thu May 15 21:12:57 2008
 * 
 * @brief  rta_preemphasis mex function 
 * 
 * Copyright (C) 2008 by IRCAM-Centre Georges Pompidou, Paris, France.
 *
 */

#include "mex.h"
#include "rta_preemphasis.h"


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
  rta_real_t factor;
  rta_real_t * previous_sample;

  /* rta inputs */
  rta_real_t * real_input; 
  unsigned int input_size;
  unsigned int input_m;


  /* rta outputs */
  rta_real_t * output;

  /* matlab outputs */
  unsigned int output_m;
  unsigned int output_n;

  /* other */
  int i;

  /* input arguments */

  /* check proper input and output */
  if(nrhs<2)
  {
    mexErrMsgTxt("Two inputs required (third is optional).");
  }
  else if(nlhs > 1)
  {
    mexErrMsgTxt("Too many output arguments.");
  }

  input = mxGetData(prhs[0]); 
  input_size = mxGetNumberOfElements(prhs[0]);
  input_m=mxGetM(prhs[0]);

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

  factor = (rta_real_t) mxGetScalar(prhs[1]);

  /* output results */
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

  plhs[0] = mxCreateNumericMatrix(output_m, output_n,
                                  RTA_MEX_REAL_TYPE, mxREAL);
  output = mxGetData(plhs[0]);
  
  plhs[1] = mxCreateNumericMatrix(1, 1, RTA_MEX_REAL_TYPE, mxREAL);
  previous_sample = mxGetData(plhs[1]);

  if(nrhs > 2)
  {
    *previous_sample = mxGetScalar(prhs[2]);
  }
  else
  {
    *previous_sample = 0.;
  }

  /* computation */
  rta_preemphasis_signal(output, real_input, input_size, previous_sample, factor);

#if (RTA_REAL_TYPE != RTA_DOUBLE_TYPE)
  if(mxGetClassID(prhs[0]) != mxSINGLE_CLASS)
  {
    /* free mem of tmp vec for float precision conversion */
    mxFree(real_input);
  }
#endif

  return;
}
