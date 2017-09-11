/**
 * @file   rta_onepole_mex.c
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Mon Sep  1 19:08:52 2008
 * 
 * @brief  RTA onepole filter.
 * 
 * Copyright (C) 2008 by IRCAM-Centre Georges Pompidou, Paris, France.
 * 
 */

#include "mex.h" 
#include "rta_onepole.h"
#include "rta_filter.h"

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
  char * filter_name;
  rta_filter_t filter_type;
  rta_real_t f0;
  rta_real_t * state;
  double * input;
  unsigned int input_size;
  unsigned int input_m;
  unsigned int input_n;
  int onepole_dim = 0;

  /* rta inputs */
  rta_real_t * real_input;
  
  /* rta outputs */
  rta_real_t * output;
  
  /* other */
  int i;

  /* input arguments */

  /* check proper input and output */
  if(nrhs < 3)
  {
    mexErrMsgTxt("Three inputs required.");
  }
  else if(nlhs > 2)
  {
    mexErrMsgTxt("Too many output arguments.");
  }
  
  if(mxIsChar(prhs[0]) != 1)
  {
    mexErrMsgTxt("Filter type must be a string.");
  }
  else
  {
    filter_name = mxArrayToString(prhs[0]);
    if(0 == strcmp(filter_name, "lowpass"))
    {
      filter_type = rta_lowpass;
    }
    else if(0 == strcmp(filter_name, "highpass"))
    {
      filter_type = rta_highpass;
    }
    else
    {
      mexErrMsgTxt("Bad filter type.");
    }
  }

  f0 = mxGetScalar(prhs[1]);

  input = mxGetData(prhs[2]); 
  input_size = mxGetNumberOfElements(prhs[2]);    
  input_m=mxGetM(prhs[2]);
  input_n=mxGetN(prhs[2]);

#if (RTA_REAL_TYPE == RTA_FLOAT_TYPE)
  if(mxGetClassID(prhs[2]) != mxSINGLE_CLASS)
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

  plhs[1] = mxCreateNumericMatrix(1, 1, RTA_MEX_REAL_TYPE, mxREAL);
  state = mxGetData(plhs[1]);

  /* copy state as it will change (and may be output) */
  if(nrhs > 3 && mxGetNumberOfElements(prhs[3]) >= 1)
  {
    *state = mxGetScalar(prhs[3]);
  }

  if(nrhs > 4)
  {
    onepole_dim = mxGetScalar(prhs[4]);
  }

  if(onepole_dim == 0)
  {
    if(input_m > 1)
    {
      onepole_dim = 1;
    }
    else
    {
      onepole_dim = 2;
    }    
  }
  
  plhs[0] = mxCreateNumericMatrix(input_m, input_n, RTA_MEX_REAL_TYPE, mxREAL);
  output = mxGetData(plhs[0]);

  switch(filter_type)
  {
    case rta_lowpass:
      if(onepole_dim == 1)
      {
        const rta_real_t initial_state = *state;
        for(i=0; i<input_n*input_m; i+=input_m)
        {
          *state = initial_state;
          rta_onepole_lowpass_vector(output + i, real_input + i, input_m, 
                                     f0, state);
        }
      }
      else
      {
        const rta_real_t initial_state = *state;
        for(i=0; i<input_m; i++)
        {
          *state = initial_state;
          rta_onepole_lowpass_vector_stride(output + i, input_m,
                                            real_input + i, input_m, input_n,
                                            f0, state);
        }
      }
      break;

    case rta_highpass:
      if(onepole_dim == 1)
      {
        const rta_real_t initial_state = *state;
        for(i=0; i<input_n*input_m; i+=input_m)
        {
          *state = initial_state;
          rta_onepole_highpass_vector(output + i, real_input + i, input_m, 
                                     f0, state);
        }
      }
      else
      {
        const rta_real_t initial_state = *state;
        for(i=0; i<input_m; i++)
        {
          *state = initial_state;
          rta_onepole_highpass_vector_stride(output + i, input_m,
                                            real_input + i, input_m, input_n,
                                            f0, state);
        }
      }
      break;
  }

#if (RTA_REAL_TYPE == RTA_FLOAT_TYPE)  
  if(mxGetClassID(prhs[2]) != mxSINGLE_CLASS)
  {
    mxFree(real_input);
  }
#endif

  return;
}

