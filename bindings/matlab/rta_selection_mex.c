/**
 * @file   rta_var_mex.c
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Mon Aug 25 16:13:42 2008
 * 
 * @brief  variance mex function for RTA 
 * 
 * Copyright (C) 2008 by IRCAM-Centre Georges Pompidou, Paris, France.
 *
 */

#include "mex.h"
#include "rta_selection.h"


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
  float * float_input; /* matlab single-precision input */
  int selection_dim = 0; /* automatic dimension */
  rta_real_t selection_index;

  /* rta inputs */
  rta_real_t * real_input; 
  unsigned int input_m;
  unsigned int input_n;

  /* rta outputs */
  rta_real_t * selection;

  /* matlab outputs */
  unsigned int output_m;
  unsigned int output_n;

  /* other */
  int i, j;

  /* input arguments */

  /* check proper input and output */
  if(nrhs<1)
  {
    mexErrMsgTxt("One input required.");
  }
  else if(nlhs > 2)
  {
    mexErrMsgTxt("Too many output arguments.");
  }

  if(nrhs>1)
  {
    /* matlab indexes start at 1 */
    selection_index = mxGetScalar(prhs[1]) - 1.;
  }

  if(nrhs>2)
  {
    selection_dim = mxGetScalar(prhs[2]);
  }

  input = mxGetData(prhs[0]); 
  input_m=mxGetM(prhs[0]);
  input_n=mxGetN(prhs[0]);
  
  /* copy input for float precision conversion (and use in-place
     calculations form now to avoid further copies */
  real_input = mxMalloc( input_m * input_n * sizeof(rta_real_t)); 
  j = input_m*input_n;

  if(mxGetClassID(prhs[0]) != mxSINGLE_CLASS)
  {
    for (i=0; i<j ;i++)
    {
      real_input[i] = (rta_real_t) input[i]; 
    }
  }
  else
  {
    float_input = (float *) input;
    for (i=0; i<j ;i++)
    {
      real_input[i] = (rta_real_t) float_input[i]; 
    }
  }

  /* automatic dimension */
  if(selection_dim == 0)
  {
    if(input_m > 1)
    {
      selection_dim = 1;
    }
    else
    {
      selection_dim = 2;
    }
  }

  if(selection_dim == 1)
  {
    output_m = 1;
    output_n = input_n;
  }
  else
  {
    output_m = input_m;
    output_n = 1;
  }

  plhs[0] = mxCreateNumericMatrix(output_m, output_n,
                                  RTA_MEX_REAL_TYPE, mxREAL);
  selection = mxGetData(plhs[0]);
  
  if(selection_dim == 1)
  {
    for(i=0; i<input_n; i++)
    {
      selection[i] = rta_selection(real_input+i*input_m, input_m, selection_index);
    }
  }
  else /* strides */
  {
    for(i=0; i<input_m; i++)
    {
      selection[i] = rta_selection_stride(
        real_input+i, input_m, input_n, selection_index);
    }    
  }

  /* free mem of tmp vec for float precision conversion */
  mxFree(real_input);

  return;
}
