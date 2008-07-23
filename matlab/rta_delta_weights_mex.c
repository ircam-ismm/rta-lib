/**
 * @file   rta_delta_weights_mex.c
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Thu May 15 21:12:57 2008
 * 
 * @brief  rta_delta mex initialisation function
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
  unsigned int delta_size;
  char * delta_name;
  rta_real_t delta_coef = 5 ;
  
  /* rta inputs */

  /* rta outputs */
  int ret = 0;
  rta_real_t * delta_weights;
  rta_real_t * delta_normalisation_factor;

  /* matlab outputs */
  

  /* input arguments */

  /* check proper input and output */
  if(nrhs < 1)
  {
    mexErrMsgTxt("One input required.");
  }
  else if(nlhs > 2)
  {
    mexErrMsgTxt("Too many output arguments.");
  }
  
  /* first arguments is delta size */
  delta_size = mxGetScalar(prhs[0]); 
  if(delta_size < 1)
  {
    mexErrMsgTxt("Delta filter size must be > 0");
  }
  else if(!(delta_size & 1)) /* even */
  {
    mexErrMsgTxt("Delta filter size must be odd.");
  }

  plhs[0] = mxCreateNumericMatrix(delta_size, 1, RTA_MEX_REAL_TYPE, mxREAL);
  delta_weights = mxGetData(plhs[0]);

  ret = rta_delta_weights(delta_weights, delta_size);

  plhs[1] = mxCreateNumericMatrix(1, 1, RTA_MEX_REAL_TYPE, mxREAL);
  delta_normalisation_factor = mxGetData(plhs[1]);
    
  *delta_normalisation_factor = rta_delta_normalization_factor(delta_size);

  if(ret == 0)
  {
    mexErrMsgTxt("Error while creating the delta weigths");
  }


  return;
}
