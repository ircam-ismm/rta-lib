/**
 * @file   rta_biquad_coefs_mex.c
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Mon Sep  1 19:08:52 2008
 * 
 * @brief  RTA biquad filter coefficients.
 * 
 * Copyright (C) 2008 by IRCAM-Centre Georges Pompidou, Paris, France.
 * 
 */

#include "mex.h" 
#include "rta_biquad.h"
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
  rta_real_t f0;
  rta_real_t q;
  rta_real_t gain = 1.; /* default is unitary gain */

  /* rta inputs */
  rta_filter_t filter_type;
  

  /* rta outputs */
  rta_real_t * b;
  rta_real_t * a; /* [1, a1, a2] */

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
    else if(0 == strcmp(filter_name, "bandpass_cst_skirt"))
    {
      filter_type = rta_bandpass_constant_skirt;
    }
    else if(0 == strcmp(filter_name, "bandpass_cst_peak"))
    {
      filter_type = rta_bandpass_constant_peak;
    }
    else if(0 == strcmp(filter_name, "notch"))
    {
      filter_type = rta_notch;
    }
    else if(0 == strcmp(filter_name, "allpass"))
    {
      filter_type = rta_allpass;
    }
    else if(0 == strcmp(filter_name, "peaking"))
    {
      filter_type = rta_peaking;
    }
    else if(0 == strcmp(filter_name, "lowshelf"))
    {
      filter_type = rta_lowshelf;
    }
    else if(0 == strcmp(filter_name, "highshelf"))
    {
      filter_type = rta_highshelf;
    }
    else
    {
      mexErrMsgTxt("Bad filter type.");
    }
  }

  f0 = mxGetScalar(prhs[1]);
  q = mxGetScalar(prhs[2]);

  /* third argument is gain */
  if(nrhs > 3)
  {
    gain = mxGetScalar(prhs[3]);
  }


  plhs[0] = mxCreateNumericMatrix(3, 1, RTA_MEX_REAL_TYPE, mxREAL);
  b = mxGetData(plhs[0]);

  plhs[1] = mxCreateNumericMatrix(3, 1, RTA_MEX_REAL_TYPE, mxREAL);
  a = mxGetData(plhs[1]);

  a[0] = 1.;
  rta_biquad_coefs(b, a + 1, filter_type, f0, q, gain);

  return;
}

