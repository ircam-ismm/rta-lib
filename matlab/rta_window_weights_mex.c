/**
 * @file   rta_window_weights_mex.c
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Thu May 15 21:12:57 2008
 * 
 * @brief  rta_window mex initialisation function
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
  unsigned int window_size;
  char * window_name;
  rta_real_t window_coef = 0.08 ; /* for hamming */
  
  /* rta inputs */

  /* rta outputs */
  int ret = 0;
  rta_real_t * window_weights;

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
  
  /* second arguments is window size */
  window_size = mxGetScalar(prhs[1]); 

  plhs[0] = mxCreateNumericMatrix(window_size, 1, RTA_MEX_REAL_TYPE, mxREAL);
  window_weights = mxGetData(plhs[0]);

  /* first argument is window type*/
  if(mxIsChar(prhs[0]) != 1)
  {
    mexErrMsgTxt("First input must be a string.");
  }

  window_name = mxArrayToString(prhs[0]);
  
  if(0 == strcmp(window_name,"hann"))
  {
    ret = rta_window_hann_weights(window_weights, window_size);
  }
  else if(0 == strcmp(window_name,"hamming"))
  {
    if(nrhs >= 3)
    {
      window_coef = mxGetScalar(prhs[2]);
    }
    ret = rta_window_hamming_weights(window_weights, window_size, window_coef);
  }
  else
  {
    mexErrMsgTxt("Bad window type.");
  }

  if(ret == 0)
  {
    mexErrMsgTxt("Error while creating the window weigths");
  }


  return;
}
