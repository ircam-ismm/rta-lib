/**
 * @file   rta_yin_setup_new_mex.c
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Thu May 15 21:12:57 2008
 * 
 * @brief  rta_yin mex initialisation function
 * 
 * Copyright (C) 2008 by IRCAM-Centre Georges Pompidou, Paris, France.
 *
 */

#include "mex.h"
#include "rta_yin.h"

static const unsigned int yin_max_mins = 128;

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
  unsigned int max_mins = yin_max_mins;
  /* rta inputs */
  rta_yin_setup_t * yin_setup;

  /* matlab outputs */
  rta_ptr_t * setup_address;

  /* input arguments */
  if(nrhs > 0)
  {
    max_mins = mxGetScalar(prhs[0]); 
  }

  /* yin setup */
  if(rta_yin_setup_new(&(yin_setup), yin_max_mins))
  {
    plhs[0] = mxCreateNumericMatrix(1, 1, RTA_MEX_PTR_TYPE, mxREAL);
    setup_address = mxGetData(plhs[0]);
    *setup_address = (rta_ptr_t) yin_setup;
  }
  else
  {
    printf("rta_yin setup failed.\n");
  }

  return;
}
