/**
 * @file   rta_fft_setup_delete_mex.c
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Thu May 15 21:12:57 2008
 * 
 * @brief  rta_fft mex delete function
 * 
 * Copyright (C) 2008 by IRCAM-Centre Georges Pompidou, Paris, France.
 *
 */

#include "mex.h"
#include "rta_fft.h"
#include "rta_fft_setup_mex.h"
#include "rta_stdlib.h"

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
  rta_fft_setup_mex_t * fft_setup_mex = (rta_fft_setup_mex_t *) setup_address;

  if(fft_setup_mex != NULL)
  {
    if(fft_setup_mex->fft_setup != NULL)
    {
      rta_fft_setup_delete(fft_setup_mex->fft_setup);
    }
    if(fft_setup_mex->fft != NULL)
    {
      rta_free(fft_setup_mex->fft);
    }
    rta_free(fft_setup_mex);
  }

  return;
}
