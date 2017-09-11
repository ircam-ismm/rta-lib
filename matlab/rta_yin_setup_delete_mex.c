/**
 * @file   rta_yin_setup_delete_mex.c
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Thu May 15 21:12:57 2008
 * 
 * @brief  rta_yin mex delete function
 * 
 * Copyright (C) 2008 by IRCAM-Centre Georges Pompidou, Paris, France.
 *
 */

#include "mex.h"
#include "rta_yin.h"


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

  /* check proper input and output */
  if(nrhs!=1)
  {
    mexErrMsgTxt("One input required.");
  }
  else if(nlhs > 0)
  {
    mexErrMsgTxt("Too many output arguments.");
  }


  /* matlab inputs */
  rta_ptr_t setup_address = mxGetScalar(prhs[0]);

  /* rta inputs */
  rta_yin_setup_t * yin_setup = (rta_yin_setup_t *) setup_address;

  rta_yin_setup_delete(yin_setup);

  return;
}
