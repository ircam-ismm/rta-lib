/**
 * @file   rta_svd_mex.c
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Mon Aug 18 17:21:26 2008
 * 
 * @brief  rta_window mex Singular Value Decomposition
 * 
 * Copyright (C) 2008 by IRCAM-Centre Georges Pompidou, Paris, France.
 *
 */

#include "mex.h" 
#include "rta_svd.h"
#include "rta_int.h" /* rta_imin */

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
  double * input;
  unsigned int m;
  unsigned int n;
  int min_mn;
  
  /* rta inputs */
  rta_svd_setup_t * svd_setup;
  rta_real_t * U = NULL;
  rta_real_t * S;
  rta_real_t * V = NULL;
  rta_real_t * real_input;
  
  /* rta outputs */
  int ret;

  /* other */
  int i,j;

  /* input arguments */

  /* check proper input and output */
  if(nrhs < 1)
  {
    mexErrMsgTxt("One input required.");
  }
  else if(nlhs > 3)
  {
    mexErrMsgTxt("Maximum three [U,S,V] output arguments.");
  }
  
  m = mxGetM(prhs[0]); 
  n = mxGetN(prhs[0]); 
  min_mn = rta_imin(m,n);
  input = mxGetData(prhs[0]);

  /* Output only S */
  if(nlhs == 1)
  {
    plhs[0] = mxCreateNumericMatrix(min_mn, 1, RTA_MEX_REAL_TYPE, mxREAL);
    S = mxGetData(plhs[0]);

  }
  else /* [U, S, V] */
  {
    plhs[0] = mxCreateNumericMatrix(m, min_mn, RTA_MEX_REAL_TYPE, mxREAL);
    U = mxGetData(plhs[0]);

    plhs[1] = mxCreateNumericMatrix(min_mn, 1, RTA_MEX_REAL_TYPE, mxREAL);
    S = mxGetData(plhs[1]);

    if(nlhs > 2) /* V requested as output result */
    {
      plhs[2] =  mxCreateNumericMatrix(n, min_mn, RTA_MEX_REAL_TYPE, mxREAL);
      V = mxGetData(plhs[2]);
    }
  }

  /* input copy, as it may change  
     (and float precision conversion as a result) */
  real_input = mxMalloc( m * n * sizeof(rta_real_t)); 
  j = m*n;
  for (i=0; i<j ;i++)
  {
    real_input[i] = (rta_real_t) input[i]; 
  }

  /* use an in place calculation as the input may be anyway copied for 
     precision adaptation: it must be always copied, then. */
  ret = rta_svd_setup_new(&svd_setup, rta_svd_in_place, 
                          U, S, V, real_input, m, n);

  rta_svd(U, S, V, real_input, svd_setup);

  rta_svd_setup_delete(svd_setup);

  /* free mem of tmp vectors */
  mxFree(real_input);

  return;
}
