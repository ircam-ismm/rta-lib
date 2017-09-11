/**
 * @file   rta_spectrum_to_bands_mex.c
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Thu May 15 21:12:57 2008
 * 
 * @brief  rta_spectrum_to_bands (mel) mex function
 * 
 * Copyright (C) 2008 by IRCAM-Centre Georges Pompidou, Paris, France.
 *
 */

#include "mex.h" 
#include "rta_bands.h"
#include "rta_mel.h"

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
  double * spectrum;
  unsigned int spectrum_size;
  unsigned int spectrum_m;

  double * weights;
  unsigned int bands_nb;
  unsigned int weights_size; /* spectrum_size * bands_nb */

  double * bounds;

  char * integration_name;
  
  /* rta inputs */
  rta_real_t * real_spectrum;
  rta_real_t * real_weights;
  unsigned int * uint_bounds;
  rta_spectrum_to_bands_function integration = 
    &rta_spectrum_to_bands_square_abs;
  /* rta outputs */

  /* matlab outputs */
  rta_real_t * output; /* output size is input_size */
  unsigned int output_m;
  unsigned int output_n;

  /* other */
  int i;

  /* check proper input and output */
  if(nrhs < 3)
  {
    mexErrMsgTxt("Three inputs required.");
  }
  else if(nlhs > 1)
  {
    mexErrMsgTxt("Too many output arguments.");
  }
  
  spectrum = mxGetData(prhs[0]); 
  /* take only the real spectrum + nyquist */ 
  spectrum_size = mxGetNumberOfElements(prhs[0]);    
  if(!(spectrum_size & 1))
  {
    /* if the spectrum size is even */
    /* we take only the real spectrum + nyquist */ 
    spectrum_size = spectrum_size/2 + 1;
  }

 spectrum_m=mxGetM(prhs[0]);

  weights = mxGetData(prhs[1]); 
  bands_nb = mxGetNumberOfElements(prhs[1]) / spectrum_size;
  weights_size = spectrum_size * bands_nb;

  bounds = mxGetData(prhs[2]); 

  if(spectrum_m == 1)
  {
    output_m = 1;
    output_n = bands_nb;
  }
  else
  {
    output_m = bands_nb;
    output_n = 1;
  }

  if(nrhs > 3)
  {
    if(mxIsChar(prhs[3]) != 1)
    {
      mexErrMsgTxt("Integration type must be a string.");
    }
    else
    {
      integration_name = mxArrayToString(prhs[3]);
      if(0 == strcmp(integration_name, "sqrabs"))
      {
        integration = &rta_spectrum_to_bands_square_abs;
      }
      else if(0 == strcmp(integration_name, "abs"))
      {
        integration = &rta_spectrum_to_bands_abs;
      }
      else
      {
        mexErrMsgTxt("Bad integration type.");
      }
    }
  }

#if (RTA_REAL_TYPE == RTA_FLOAT_TYPE)
  if(mxGetClassID(prhs[0]) != mxSINGLE_CLASS)
  {
    /* input float precision conversion */
    real_spectrum = mxMalloc( spectrum_size * sizeof(rta_real_t)); 
    for (i=0; i<spectrum_size ;i++)
    {
      real_spectrum[i] = (rta_real_t) spectrum[i]; 
    }
  }
  else
  {
    /* no conversion for matlab single */
    real_spectrum = (rta_real_t *) spectrum;
#else
    /* no conversion for rta double */
    real_spectrum = (rta_real_t *) spectrum;
#endif

#if (RTA_REAL_TYPE == RTA_FLOAT_TYPE)
  }
#endif


#if (RTA_REAL_TYPE == RTA_FLOAT_TYPE)
  if(mxGetClassID(prhs[1]) != mxSINGLE_CLASS)
  {
    /* input float precision conversion */
    real_weights = mxMalloc( weights_size * sizeof(rta_real_t)); 
    for (i=0; i<weights_size ;i++)
    {
      real_weights[i] = (rta_real_t) weights[i]; 
    }
  }
  else
  {
    /* no conversion for matlab single */
    real_weights = (rta_real_t *) weights;
#else
    /* no conversion for rta double */
    real_weights = (rta_real_t *) weights;
#endif

#if (RTA_REAL_TYPE == RTA_FLOAT_TYPE)
  }
#endif

  if(mxGetClassID(prhs[2]) != RTA_MEX_UNSIGNED_INTEGER_TYPE)
  {
    /* unsigned int conversion */
    uint_bounds = mxMalloc( 2 * bands_nb * sizeof(unsigned int)); 
    for (i=0; i<2*bands_nb ;i++)
    {
      uint_bounds[i] = (unsigned int) bounds[i]; 
    }
  }
  else
  {
    /* no conversion for matlab unsigned int */
    uint_bounds = (unsigned int *) bounds;
  }


  plhs[0] = mxCreateNumericMatrix(output_m, output_n, RTA_MEX_REAL_TYPE, mxREAL);
  output = mxGetData(plhs[0]);

  (*integration)(output, real_spectrum, real_weights, uint_bounds, 
                     spectrum_size, bands_nb);

#if (RTA_REAL_TYPE != RTA_DOUBLE_TYPE)
  if(mxGetClassID(prhs[0]) != mxSINGLE_CLASS)
  {
    /* free mem of tmp vec for float precision conversion */
    mxFree(real_spectrum);
  }

  if(mxGetClassID(prhs[1]) != mxSINGLE_CLASS)
  {
    /* free mem of tmp vec for float precision conversion */
    mxFree(real_weights);
  }
#endif

  if(mxGetClassID(prhs[2]) != RTA_MEX_UNSIGNED_INTEGER_TYPE)
  {
    mxFree(uint_bounds);
  }

  return;
}
