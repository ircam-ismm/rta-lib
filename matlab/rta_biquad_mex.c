/**
 * @file   rta_biquad_mex.c
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Mon Sep  1 19:08:52 2008
 * 
 * @brief  RTA biquad filter.
 * 
 * Copyright (C) 2008 by IRCAM-Centre Georges Pompidou, Paris, France.
 * 
 */

#include "mex.h" 
#include "rta_biquad.h"

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
  double * b;
  double * a;
  double * input;
  unsigned int input_size;
  unsigned int input_m;
  unsigned int input_n;
  double * state;
  float * float_state; /* matlab float */
  int biquad_dim = 0;
  enum {df1, df2t} biquad_type = df1;
  char * biquad_name;

  /* rta inputs */
  rta_real_t * real_b;
  rta_real_t * real_a; /* [a1, a2] */
  rta_real_t * real_input;
  rta_real_t * real_state;
  unsigned int states_nb;

  /* rta outputs */
  rta_real_t * output;
  
  /* other */
  int i,j;

  /* input arguments */

  /* check proper input and output */
  if(nrhs < 3)
  {
    mexErrMsgTxt("Three inputs required.");
  }
  else if(nlhs > 2)
  {
    mexErrMsgTxt("Too many output arguments.");
  }
  
  b = mxGetData(prhs[0]);
  a = mxGetData(prhs[1]);
  input = mxGetData(prhs[2]); 
  input_size = mxGetNumberOfElements(prhs[2]);    
  input_m=mxGetM(prhs[2]);
  input_n=mxGetN(prhs[2]);

#if (RTA_REAL_TYPE == RTA_FLOAT_TYPE)
  if(mxGetClassID(prhs[0]) != mxSINGLE_CLASS)
  {
    /* input float precision conversion */
    real_b = mxMalloc( 3 * sizeof(rta_real_t)); 
    for (i=0; i<3 ;i++)
    {
      real_b[i] = (rta_real_t) b[i]; 
    }
  }
  else
  {
    /* no conversion for matlab single */
    real_b = (rta_real_t *) b;
#else
    /* no conversion for rta double */
    real_b = (rta_real_t *) b;
#endif

#if (RTA_REAL_TYPE == RTA_FLOAT_TYPE)
  }
#endif

#if (RTA_REAL_TYPE == RTA_FLOAT_TYPE)
  if(mxGetClassID(prhs[1]) != mxSINGLE_CLASS)
  {
    /* input float precision conversion */
    real_a = mxMalloc( 3 * sizeof(rta_real_t)); 
    for (i=0; i<3 ;i++)
    {
      real_a[i] = (rta_real_t) a[i]; 
    }
  }
  else
  {
    /* no conversion for matlab single */
    real_a = (rta_real_t *) a;
#else
    /* no conversion for rta double */
    real_a = (rta_real_t *) a;
#endif

#if (RTA_REAL_TYPE == RTA_FLOAT_TYPE)
  }
#endif

#if (RTA_REAL_TYPE == RTA_FLOAT_TYPE)
  if(mxGetClassID(prhs[2]) != mxSINGLE_CLASS)
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

  if(nrhs > 5)
  {
    if(mxIsChar(prhs[5]) != 1)
    {
      mexErrMsgTxt("Biquad type must be a string.");
    }

    biquad_name = mxArrayToString(prhs[5]);
    if(0 == strcmp(biquad_name, "df1"))
    {
      biquad_type = df1;
    }
    else if(0 == strcmp(biquad_name, "df2t"))
    {
      biquad_type = df2t;
    }
    else
    {
      mexErrMsgTxt("Bad biquad type.");
    }
  }
    
  switch(biquad_type)
  {
    case df1:
      states_nb = 4;
      break;

    case df2t:
      states_nb = 2;
      break;

    default:
      break;
  }

  plhs[1] = mxCreateNumericMatrix(states_nb, 1, RTA_MEX_REAL_TYPE, mxREAL);
  real_state = mxGetData(plhs[1]);

  /* copy state as it will change (and may be output) */
  if(nrhs > 3 && mxGetNumberOfElements(prhs[3]) >= states_nb)
  {
    if(mxGetClassID(prhs[3]) == mxSINGLE_CLASS)
    {
      float_state = (float *) mxGetData(prhs[3]);
      
      for(i = 0; i < states_nb; i++)
      {
        real_state[i] = (rta_real_t) float_state[i];
      }
    }
    else
    {
      state = mxGetData(prhs[3]);
      for(i = 0; i < states_nb; i++)
      {
        real_state[i] = (rta_real_t) state[i];
      }
    }
  }

  if(nrhs > 4)
  {
    biquad_dim = mxGetScalar(prhs[4]);
  }

  if(biquad_dim == 0)
  {
    if(input_m > 1)
    {
      biquad_dim = 1;
    }
    else
    {
      biquad_dim = 2;
    }    
  }
  
  plhs[0] = mxCreateNumericMatrix(input_m, input_n, RTA_MEX_REAL_TYPE, mxREAL);
  output = mxGetData(plhs[0]);

  if(biquad_dim == 1)
  {
    rta_real_t initial_state[states_nb];
    for(j = 0 ; j < states_nb ; j ++)
    {
      initial_state[j] = real_state[j];
    }

    for(i=0; i<input_n*input_m; i+=input_m)
    {
      for(j = 0 ; j < states_nb ; j ++)
      {
        real_state[j] = initial_state[j];
      }

      switch(biquad_type)
      {
        case df1:
          rta_biquad_df1_vector(output + i, real_input + i, input_m, 
                                real_b, real_a + 1, real_state);
          break;

        case df2t:
          rta_biquad_df2t_vector(output + i, real_input + i, input_m, 
                                 real_b, real_a + 1, real_state);
          break;

        default:
          break;
      }
    }
  }
  else
  {
    rta_real_t initial_state[states_nb];
    for(j = 0 ; j < states_nb ; j ++)
    {
      initial_state[j] = real_state[j];
    }

    for(i=0; i<input_m; i++)
    {
      for(j = 0 ; j < states_nb ; j ++)
      {
        real_state[j] = initial_state[j];
      }
      
      switch(biquad_type)
      {
        case df1:
          rta_biquad_df1_vector_stride(output + i, input_m, 
                                       real_input + i, input_m, input_n,
                                       real_b, 1, real_a + 1, 1,
                                       real_state, 1);
          break;

        case df2t:
          rta_biquad_df2t_vector_stride(output + i, input_m, 
                                        real_input + i, input_m, input_n,
                                        real_b, 1, real_a + 1, 1,
                                        real_state, 1);
          break;

        default:
          break;
      }
    }
  }


#if (RTA_REAL_TYPE == RTA_FLOAT_TYPE)
  if(mxGetClassID(prhs[0]) != mxSINGLE_CLASS)
  {
    mxFree(real_b);
  }

  if(mxGetClassID(prhs[1]) != mxSINGLE_CLASS)
  {
    mxFree(real_a);
  }
  
  if(mxGetClassID(prhs[2]) != mxSINGLE_CLASS)
  {
    mxFree(real_input);
  }
#endif

  return;
}

