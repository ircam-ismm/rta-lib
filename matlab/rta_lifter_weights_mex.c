/**
 * @file   rta_lifter_weights_mex.c
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Thu May 15 21:12:57 2008
 * 
 * @brief  rta_lifter mex initialisation function
 * 
 * Copyright (C) 2008 by IRCAM-Centre Georges Pompidou, Paris, France.
 *
 */

#include "mex.h" 
#include "rta_lifter.h"


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
  unsigned int lifter_size;
  char * lifter_name;
  rta_lifter_mode_t lifter_mode = rta_lifter_mode_normal;
  rta_lifter_t lifter_type = rta_lifter_exponential;
  rta_real_t lifter_coef = 0.6 ; /* default for exponential */
  
  /* rta inputs */

  /* rta outputs */
  int ret = 0;
  rta_real_t * lifter_weights;

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
  
  /* first arguments is lifter size */
  lifter_size = mxGetScalar(prhs[0]); 

  plhs[0] = mxCreateNumericMatrix(lifter_size, 1, RTA_MEX_REAL_TYPE, mxREAL);
  lifter_weights = mxGetData(plhs[0]);

  /* second argument is coefficient */
  if(nrhs > 1)
  {
    lifter_coef = mxGetScalar(prhs[1]); 
  }

  /* third argument is lifter type*/
  if(nrhs > 2)
  {
    if(mxIsChar(prhs[2]) != 1)
    {
      mexErrMsgTxt("First input must be a string.");
    }
    
    lifter_name = mxArrayToString(prhs[2]);
  
    if(0 == strcmp(lifter_name,"exp") || 0 == strcmp(lifter_name,"slaney"))
    {
      lifter_type = rta_lifter_exponential;
    }
    else if(0 == strcmp(lifter_name,"sin") || 0 == strcmp(lifter_name,"htk"))
    {
      if(lifter_coef <= 0.)
      {
        mexErrMsgTxt("Coefficient must be > 0 for 'sin' type.");
      }
      lifter_type = rta_lifter_sinusoidal;
    }
    else
    {
      mexErrMsgTxt("Bad lifter type.");
    }
  }

  /* fourth argument is mode */
  if(nrhs > 3)
  {
    if(mxIsChar(prhs[3]) != 1)
    {
      mexErrMsgTxt("First input must be a string.");
    }
    
    lifter_name = mxArrayToString(prhs[3]);
  
    if(0 == strcmp(lifter_name,"normal"))
    {
      lifter_mode = rta_lifter_mode_normal;
    }
    else if(0 == strcmp(lifter_name,"inverse"))
    {
      lifter_mode = rta_lifter_mode_inverse;
    }
    else
    {
      mexErrMsgTxt("Bad lifter mode.");
    }
  }

  ret = rta_lifter_weights(lifter_weights, lifter_size, 
                           lifter_coef, lifter_type, lifter_mode);

  if(ret == 0)
  {
    mexErrMsgTxt("Error while creating the lifter weigths");
  }


  return;
}
