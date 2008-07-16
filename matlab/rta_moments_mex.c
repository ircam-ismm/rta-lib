/**
 * @file   rta_moments_mex.c
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Thu May 15 21:12:57 2008
 * 
 * @brief  rta_moments mex function
 * 
 * Copyright (C) 2008 by IRCAM-Centre Georges Pompidou, Paris, France.
 *
 */

#include "mex.h"
#include "rta_moments.h"
#include "rta_math.h" /* rta_sqrt */


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
  unsigned int max_order;
  char * moments_name;
  int std = 1; /* default is to compute the standardised moments */

  /* rta inputs */
  rta_real_t * real_input; 
  unsigned int input_size;
  unsigned int input_m;
  rta_real_t deviation; /* square root of variance, which is
                         * output[1], and needed by standardised
                         * moments */
  /* rta outputs */
  rta_real_t * output;
  rta_real_t * sum;

  /* matlab outputs */
  unsigned int output_m;
  unsigned int output_n;

  /* other */
  int i;

  /* input arguments */

  /* check proper input and output */
  if(nrhs<2)
  {
    mexErrMsgTxt("Two inputs required.");
  }
  else if(nlhs > 2)
  {
    mexErrMsgTxt("Too many output arguments.");
  }

  input = mxGetData(prhs[0]); 
  input_size = mxGetNumberOfElements(prhs[0]);
  input_m=mxGetM(prhs[0]);

  max_order = mxGetScalar(prhs[1]);
  if(max_order < 1)
  {
    mexErrMsgTxt("Moments maximum order must be > 0");
  }

  if(nrhs > 2)
  {
    if(mxIsChar(prhs[2]) != 1)
    {
      mexErrMsgTxt("Moments type must be a string.");
    }
    else
    {
      moments_name = mxArrayToString(prhs[2]);
      if(0 == strcmp(moments_name, "nostd"))
      {
        std = 0;
      }
      else if(0 == strcmp(moments_name, "std"))
      {
        /* default is std = 1 */
      }
      else
      {
        mexErrMsgTxt("Bad bands type.");
      }
    }
  }

#if (RTA_REAL_TYPE == RTA_FLOAT_TYPE)
  if(mxGetClassID(prhs[0]) != mxSINGLE_CLASS)
    {
      /* input float precision conversion */
      real_input = mxMalloc( input_size * sizeof(rta_real_t)); 
      for (i=0; i<input_size ;i++)
      {
        real_input[i] = (rta_real_t) input[i]; 
      }
    }
    else
    {
      /* no conversion for matlab single */
      real_input = (rta_real_t *) input;
#else
      /* no conversion for rta double */
      real_input = (rta_real_t *) input;
#endif

#if (RTA_REAL_TYPE == RTA_FLOAT_TYPE)
    }
#endif


  /* output results */
  if(input_m == 1)
  {
    output_m = max_order;
    output_n = 1;
  }
  else
  {
    output_m = 1;
    output_n = max_order;
  }

  plhs[0] = mxCreateNumericMatrix(output_m, output_n,
                                  RTA_MEX_REAL_TYPE, mxREAL);
  output = mxGetData(plhs[0]);

  plhs[1] = mxCreateNumericMatrix(1, 1, RTA_MEX_REAL_TYPE, mxREAL);
  sum = mxGetData(plhs[1]);

  output[0] = rta_weighted_moment_1_indexes(sum, real_input, input_size);
  
  if(max_order >= 2)
  {
    if(output[0] != 0.)
    {
      output[1] = rta_weighted_moment_2_indexes(
        real_input, input_size, output[0], *sum);
    }
    else
    {
      /* not enough to approximate noise, but we do not want to
         calculate anything here */
      output[1] = input_size;
    }
  }

  if(std == 1)
  {
    if(max_order >= 3)
    {
      deviation = rta_sqrt(output[1]);
      if(output[0] != 0. && deviation != 0.)
      {
        output[2] = rta_std_weighted_moment_3_indexes(
          real_input, input_size, output[0], *sum, deviation);
      }
      else
      {
        output[2] = 0.;
      }
    }

    if(max_order >= 4)
    {
      if(output[0] != 0. && deviation != 0.)
      {
        output[3] = rta_std_weighted_moment_4_indexes(
          real_input, input_size, output[0], *sum, deviation);
      }
      else
      {
        output[3] = 2.;
      }
    }

    for(i = 5; i <= max_order; i++)
    {
      if(output[0] != 0. && deviation != 0.)
      {
        output[i-1] = rta_std_weighted_moment_indexes(
          real_input, input_size, output[0], *sum, deviation, i);
      }
      else
      {
        /* pseudo-similar to noise until order <= 10 */
        if(i & 1) /* even */
        {
          output[i-1] = 0.;
        }
        else
        {
          output[i-1] = i;
        }
      }
    }
  }
  else /* moments not standardised */
  {
    if(max_order >= 3)
    {
      if(output[0] != 0.)
      {
        output[2] = rta_weighted_moment_3_indexes(
          real_input, input_size, output[0], *sum);
      }
      else
      {
        output[2] = 0.;
      }
    }

    if(max_order >= 4)
    {
      if(output[0] != 0.)
      {
        output[3] = rta_weighted_moment_4_indexes(
          real_input, input_size, output[0], *sum);
      }
      else
      {
        output[3] = 2.;
      }
    }

    for(i = 5; i <= max_order; i++)
    {
      if(output[0] != 0.)
      {
        output[i-1] = rta_weighted_moment_indexes(
          real_input, input_size, output[0], *sum, i);
      }
      else
      {
        /* pseudo-similar to noise until order <= 10 */
        if(i & 1) /* even */
        {
          output[i-1] = 0.;
        }
        else
        {
          output[i-1] = output[1];
        }
      }
    }
  }

#if (RTA_REAL_TYPE != RTA_DOUBLE_TYPE)
  if(mxGetClassID(prhs[0]) != mxSINGLE_CLASS)
  {
    /* free mem of tmp vec for float precision conversion */
    mxFree(real_input);
  }
#endif

  return;
}
