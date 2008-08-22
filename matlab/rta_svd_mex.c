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

  
  /* copy input for float precision conversion (and use in-place
     calculations form now to avoid further copies */
  real_input = mxMalloc( m * n * sizeof(rta_real_t)); 
  j = m*n;
  for (i=0; i<j ;i++)
  {
    real_input[i] = (rta_real_t) input[i]; 
  }

  /* rta 2D matrices are in row-major order while matlab 2D matrices
   * are in column-major order:
   * pseudo-transpose input (by swapping m and n) and swap U and V (and
   * later reorder U and V) */
  ret = rta_svd_setup_new(&svd_setup, rta_svd_in_place, 
                          V, S, U, real_input, n, m);

  rta_svd(V, S, U, real_input, svd_setup);

  rta_svd_setup_delete(svd_setup);

  if(nlhs > 1) /* U reqested, reorder (transpose-like) */
  {
    j = m * min_mn;
    for(i=0; i<j; i++)
    {
      real_input[i] = U[i]; /* tmp */
    }

    for (i=0; i<m ;i++)
    {
      for(j=0; j<min_mn; j++)
      {
        U[j*m + i] = real_input[i*min_mn + j];
      }
    }

    if(nlhs > 2) /* V requested, reorder (transpose-like) */
    {
      j = min_mn * n;
      for(i=0; i<j; i++)
      {
        real_input[i] = V[i]; /* tmp */
      }

      for (i=0; i<min_mn ;i++)
      {
        for(j=0; j<n; j++)
        {
          V[i*n + j] = real_input[j*min_mn + i];
        }
      }
    }
  }


  /* free mem of tmp vectors */
  mxFree(real_input);

  return;
}
