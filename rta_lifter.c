/**
 *
 *
 * @file   rta_lifter.c
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Fri Jun 15 15:29:25 2007
 * 
 * @brief  Cepstral liftering (HTK and Auditory Toolbox styles)
 * 
 * Based on Rastamat by Dan Ellis.
 * See http://www.ee.columbia.edu/~dpwe/resources/matlab/rastamat
 *
 * Copyright (C) 2007 by IRCAM-Centre Georges Pompidou, Paris, France.
 * 
 */
#include "rta_lifter.h"
#include "rta_math.h"

int rta_lifter_weights(rta_real_t * weights_vector, const unsigned int cepstrum_order,
                     const rta_real_t liftering_factor,
                     const rta_lifter_t lifter_t, const rta_lifter_mode_t lifter_m)
{
  int i;
  int ret; /* return value */

  if(lifter_t == rta_lifter_exponential)
  {
    weights_vector[0] = weights_vector[1] = 1.;
    for(i=2; i<cepstrum_order; i++)
    {
      weights_vector[i] = rta_pow(i, liftering_factor);
    }
  }
  /* HTK-style liftering */
  else if(lifter_t == rta_lifter_sinusoidal)
  {
    weights_vector[0] = 1.;
    for(i=1; i<cepstrum_order; i++)
    {
      weights_vector[i] = 1. + liftering_factor / 2. *
        rta_sin(i * M_PI / liftering_factor);
    }
  }
  else
  {
    ret = 0;
  }
  /* inversion */
  if(lifter_m == rta_lifter_mode_inverse)
  {
    for(i=1; i<cepstrum_order; i++)
    {
      weights_vector[i] = 1./weights_vector[i];
    }
  }
    
  return ret;
}

int rta_lifter_weights_stride(rta_real_t * weights_vector, const int w_stride,
                           const unsigned int cepstrum_order,
                           const rta_real_t liftering_factor,
                           const rta_lifter_t lifter_t, const rta_lifter_mode_t lifter_m)
{
  int i;
  int ret = 1; /* return value */

  if(lifter_t == rta_lifter_exponential)
  {
    weights_vector[0] = weights_vector[w_stride] = 1.;
    for(i=2; i<cepstrum_order; i++)
    {
      weights_vector[i*w_stride] = rta_pow(i, liftering_factor);
    }
  }
  /* HTK-style liftering */
  else if(lifter_t == rta_lifter_sinusoidal)
  {
    weights_vector[0] = 1.;
    for(i=1; i<cepstrum_order; i++)
    {
      weights_vector[i*w_stride] = 1. + liftering_factor / 2. *
        rta_sin(i * M_PI / liftering_factor);
    }
  }
  else
  {
    ret = 0;
  }
  /* inversion */
  if(lifter_m == rta_lifter_mode_inverse)
  {
    for(i=w_stride; i<cepstrum_order*w_stride; i+=w_stride)
    {
      weights_vector[i] = 1./weights_vector[i];
    }
  }
    
  return ret;
}


void rta_lifter_cepstrum(rta_real_t * out_cepstrum, rta_real_t * in_cepstrum,
                       const rta_real_t * weights_vector, 
                       const unsigned int cepstrum_order)
{
  int i;
  for(i=0; i<cepstrum_order; i++)
  {
    out_cepstrum[i] = in_cepstrum[i] * weights_vector[i];
  }
  return;
}

void rta_lifter_cepstrum_in_place(rta_real_t * cepstrum, const rta_real_t * weights_vector, 
                              const unsigned int cepstrum_order)
{
  int i;
  for(i=0; i<cepstrum_order; i++)
  {
    cepstrum[i] *= weights_vector[i];
  }
  return;
}

void rta_lifter_cepstrum_stride(rta_real_t * out_cepstrum, const int o_stride,
                             rta_real_t * in_cepstrum, const int i_stride,
                             const rta_real_t * weights_vector, const int w_stride,
                             const unsigned int cepstrum_order)
{
  int ii, io, iw;
  for(ii=0, io=0, iw=0;
      ii<cepstrum_order*i_stride;
      ii+=i_stride, io+=o_stride, iw+=w_stride)
  {
    out_cepstrum[io] = in_cepstrum[ii] * weights_vector[iw];
  }
  return;
}
