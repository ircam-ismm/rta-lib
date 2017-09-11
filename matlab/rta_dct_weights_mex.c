/**
 * @file   rta_dct_weights_mex.c
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Thu May 15 21:12:57 2008
 * 
 * @brief  rta_window mex initialisation function
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
  unsigned int input_size;
  unsigned int output_size;
  char * dct_name;
  
  /* rta inputs */
  rta_dct_t dct_type = rta_dct_slaney;

  /* rta outputs */
  int ret = 0;
  rta_real_t * dct_weights;

  /* matlab outputs */
  

  /* input arguments */

  /* check proper input and output */
  if(nrhs < 2)
  {
    mexErrMsgTxt("Two inputs required.");
  }
  else if(nlhs > 1)
  {
    mexErrMsgTxt("Too many output arguments.");
  }
  
  input_size = mxGetScalar(prhs[0]); 

  output_size = mxGetScalar(prhs[1]); 

  if(nrhs > 2)
  {
    if(mxIsChar(prhs[2]) != 1)
    {
      mexErrMsgTxt("First input must be a string.");
    }

    dct_name = mxArrayToString(prhs[2]);
    if(0 == strcmp(dct_name, "slaney"))
    {
      dct_type = rta_dct_slaney;
    }
    else if(0 == strcmp(dct_name, "htk"))
    {
      dct_type = rta_dct_htk;
    }
    else
    {
      mexErrMsgTxt("Bad DCT type.");
    }
  }

  plhs[0] = mxCreateNumericMatrix(input_size, output_size, 
                                  RTA_MEX_REAL_TYPE, mxREAL);
  dct_weights = mxGetData(plhs[0]);

  ret = rta_dct_weights(dct_weights, input_size, output_size, dct_type);

  if(ret == 0)
  {
    mexErrMsgTxt("Error while creating the DCT weigths");
  }

  return;
}
