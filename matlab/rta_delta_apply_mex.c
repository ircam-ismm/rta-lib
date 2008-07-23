/**
 * @file   rta_delta_apply_mex.c
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Thu May 15 21:12:57 2008
 * 
 * @brief  rta_delta_apply mex function
 * 
 * Copyright (C) 2008 by IRCAM-Centre Georges Pompidou, Paris, France.
 *
 */

#include "mex.h" 
#include "rta_delta.h"


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
  unsigned int input_n;

  double * weights;
  unsigned int weights_size;

  int vector_mode; /* rta_delta_vector() or rta_delta() */

  /* rta inputs */
  rta_real_t * real_input;
  rta_real_t * real_weights;

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
  input_n=mxGetN(prhs[0]);


  weights = mxGetData(prhs[1]); 
  weights_size = mxGetNumberOfElements(prhs[1]);

  printf("m = %d, n = %d, w = %d\n", input_m, input_n, weights_size);

  if(input_size == weights_size)
  {
    vector_mode = 0;
    output_m = 1;
    output_n = 1;
  }
  else if(input_n == weights_size)
  {
    vector_mode = 1;
    output_m = input_m;
    output_n = 1;
  }
  else
  {
    mexErrMsgTxt("Input size (or number of columns) must be the same as weights size.");
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

  if(vector_mode == 1)
  {
    rta_delta_vector(output, real_input, input_m, 
                     real_weights, weights_size);
  }
  else
  {
    rta_delta(output, real_input, real_weights, weights_size);
  }
  
  /* scale output */
  if(nrhs > 2)
  {
    rta_real_t scale = mxGetScalar(prhs[2]);
    for(i = 0; i<output_m; i++)
    {
      output[i] *= scale;
    }
  }

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
