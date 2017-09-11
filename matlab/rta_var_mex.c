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
#include "rta_mean_variance.h"


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
  int variance_bias = 0; /* same as matlab var argument */
  enum {fast, accurate} variance_precision = accurate;
  char * variance_name;
  int variance_dim = 0; /* automatic by default, 1 for columns, 2 for rows */

  /* rta inputs */
  rta_real_t * real_input; 
  unsigned int input_size;
  unsigned int input_m;
  unsigned int input_n;

  /* rta outputs */
  rta_real_t * variance;
  rta_real_t * mean;

  /* matlab outputs */
  unsigned int output_m;
  unsigned int output_n;

  /* other */
  int i;

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
    variance_bias = mxGetScalar(prhs[1]);
  }

  if(nrhs>2)
  {
    variance_dim = mxGetScalar(prhs[2]);
  }

  if(nrhs > 3)
  {
    if(mxIsChar(prhs[3]) != 1)
    {
      mexErrMsgTxt("Fourth input must be a string.");
    }

    variance_name = mxArrayToString(prhs[3]);
    if(0 == strcmp(variance_name, "fast"))
    {
      variance_precision = fast;
    }
    else if(0 == strcmp(variance_name, "accurate"))
    {
      variance_precision = accurate;
    }
    else
    {
      mexErrMsgTxt("Bad variance precision type.");
    }
  }

  

  input = mxGetData(prhs[0]); 
  input_size = mxGetNumberOfElements(prhs[0]);
  input_m=mxGetM(prhs[0]);
  input_n=mxGetN(prhs[0]);
  

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

  if(variance_dim == 0)
  {
    if(input_m > 1)
    {
      variance_dim = 1;
    }
    else
    {
      variance_dim = 2;
    }
  }
 

  if(variance_dim == 1)
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
  variance = mxGetData(plhs[0]);
  

  plhs[1] = mxCreateNumericMatrix(output_m, output_n, 
                                    RTA_MEX_REAL_TYPE, mxREAL);
  mean =  mxGetData(plhs[1]);

  if(variance_dim == 1)
  {
    switch(variance_precision)
    {
      case accurate:
        if(variance_bias == 0)
        {
          for(i=0; i<input_n; i++)
          {
            mean[i] = rta_mean(real_input+i*input_m, input_m);
            variance[i] = rta_variance_unbiased(real_input+i*input_m, 
                                                input_m, mean[i]);
          }
        }
        else
        {
          for(i=0; i<input_n; i++)
          {
            mean[i] = rta_mean(real_input+i*input_m, input_m);
            variance[i] = rta_variance(real_input+i*input_m, input_m, mean[i]);
          }
        }
        break;

      case fast:
        if(variance_bias == 0)
        {
          for(i=0; i<input_n; i++)
          {
            rta_mean_variance_unbiased(mean+i, variance+i,
                                       real_input+i*input_m, input_m);
          }
        }
        else
        {
          for(i=0; i<input_n; i++)
          {
            rta_mean_variance(mean+i, variance+i, real_input+i*input_m, input_m);
          }
        }
        break;
    }
  }
  else /* strides */
  {
    switch(variance_precision)
    {
      case accurate:
        if(variance_bias == 0)
        {
          for(i=0; i<input_m; i++)
          {
            mean[i] = rta_mean_stride(real_input+i, input_m, input_n);
            variance[i] = rta_variance_unbiased_stride(
              real_input+i, input_m, input_n, mean[i]);
          }
        }
        else
        {
          for(i=0; i<input_m; i++)
          {
            mean[i] = rta_mean_stride(real_input+i, input_m, input_n);
            variance[i] = rta_variance_stride(
              real_input+i, input_m, input_n, mean[i]);
          }
        }
        break;

      case fast:
        if(variance_bias == 0)
        {
          for(i=0; i<input_m; i++)
          {
            rta_mean_variance_unbiased_stride(
              mean+i, variance+i,real_input+i, input_m, input_n);
          }
        }
        else
        {
          for(i=0; i<input_m; i++)
          {
            rta_mean_variance_stride(
              mean+i, variance+i, real_input+i, input_m, input_n);
          }
        }
        break;
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
