/**
 * @file   rta_bands_weights_mex.c
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Thu May 15 21:12:57 2008
 * 
 * @brief  rta_bands mex initialisation function
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
  unsigned int spectrum_size;
  unsigned int bands_number;
  char * bands_name;
  rta_real_t sample_rate = 44100.;
  rta_real_t min_freq = 0.;
  rta_real_t max_freq = 0.5 * sample_rate;

  /* rta inputs */
  rta_hz_to_mel_function hz_to_mel = &rta_hz_to_mel_slaney;
  rta_mel_to_hz_function mel_to_hz = &rta_mel_to_hz_slaney;
  rta_mel_t mel_type = rta_mel_slaney;
  

  /* rta outputs */
  int ret = 0;
  rta_real_t * bands_weights;
  unsigned int * bands_bounds;

  /* input arguments */

  /* check proper input and output */
  if(nrhs < 2)
  {
    mexErrMsgTxt("Two inputs required.");
  }
  else if(nlhs > 2)
  {
    mexErrMsgTxt("Too many output arguments.");
  }
  
  spectrum_size = mxGetScalar(prhs[0]);
  if(!(spectrum_size & 1))
  {
    /* if the spectrum size is even */
    /* we take only the real spectrum + nyquist */ 
    spectrum_size = spectrum_size/2 + 1;
  }

  bands_number = mxGetScalar(prhs[1]); 

  /* third argument is bands type*/
  if(nrhs > 2)
  {
    if(mxIsChar(prhs[2]) != 1)
    {
      mexErrMsgTxt("Bands type must be a string.");
    }
    else
    {
      bands_name = mxArrayToString(prhs[2]);
      if(0 == strcmp(bands_name, "slaney"))
      {
        mel_type = rta_mel_slaney;
        hz_to_mel = &rta_hz_to_mel_slaney;
        mel_to_hz = &rta_mel_to_hz_slaney;
      }
      else if(0 == strcmp(bands_name, "htk"))
      {
        mel_type = rta_mel_htk;
        hz_to_mel = &rta_hz_to_mel_htk;
        mel_to_hz = &rta_mel_to_hz_htk;
      }
      else
      {
        mexErrMsgTxt("Bad bands type.");
      }
    }
  }
  
  if(nrhs > 3)
  {
    sample_rate = mxGetScalar(prhs[3]);
    max_freq = 0.5 * sample_rate;
  }

  if(nrhs > 4)
  {
    min_freq = mxGetScalar(prhs[4]);
  }

  if(nrhs > 5)
  {
    max_freq = mxGetScalar(prhs[5]);
  }

  plhs[0] = mxCreateNumericMatrix(spectrum_size, bands_number, 
                                  RTA_MEX_REAL_TYPE, mxREAL);
  bands_weights = mxGetData(plhs[0]);

  plhs[1] = mxCreateNumericMatrix(2, bands_number, 
                                  RTA_MEX_UNSIGNED_INTEGER_TYPE, mxREAL);
  bands_bounds = mxGetData(plhs[1]);

  ret = rta_spectrum_to_mel_bands_weights(
    bands_weights, bands_bounds, spectrum_size,
    sample_rate, bands_number, min_freq, max_freq,
    1., hz_to_mel, mel_to_hz, mel_type);

  if(ret == 0)
  {
    mexErrMsgTxt("Error while creating the window weigths");
  }


  return;
}
