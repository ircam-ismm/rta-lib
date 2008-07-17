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
  enum {hann, hamming} window_type = hann;
  rta_real_t window_coef = 0.08 ; /* for hamming */
  
  /* rta inputs */

  /* rta outputs */
  int ret = 0;
  rta_real_t * window_weights;

  /* matlab outputs */
  

  /* input arguments */

  /* check proper input and output */
  if(nrhs < 1)
  {
    mexErrMsgTxt("One input required.");
  }
  else if(nlhs > 1)
  {
    mexErrMsgTxt("Too many output arguments.");
  }
  
  /* first arguments is window size */
  window_size = mxGetScalar(prhs[0]); 

  plhs[0] = mxCreateNumericMatrix(window_size, 1, RTA_MEX_REAL_TYPE, mxREAL);
  window_weights = mxGetData(plhs[0]);

  /* second argument is window type*/
  if(nrhs > 1)
  {
    if(mxIsChar(prhs[1]) != 1)
    {
      mexErrMsgTxt("First input must be a string.");
    }
    
    window_name = mxArrayToString(prhs[1]);
  
    if(0 == strcmp(window_name,"hann"))
    {
      window_type = hann;
    }
    else if(0 == strcmp(window_name,"hamming"))
    {
      if(nrhs > 2)
      {
        window_coef = mxGetScalar(prhs[2]);
      }
    }
    else
    {
      mexErrMsgTxt("Bad window type.");
    }
  }

  if(window_type == hann)
  {
    ret = rta_window_hann_weights(window_weights, window_size);
  }
  else if(window_type == hamming)
  {
    ret = rta_window_hamming_weights(window_weights, window_size, window_coef);
  }

  if(ret == 0)
  {
    mexErrMsgTxt("Error while creating the window weigths");
  }


  return;
}
