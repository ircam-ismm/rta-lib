/**
 *
 *
 * @file   rta_preemphasis.c
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Tue Sep  4 16:24:45 2007
 * 
 * @brief  Preemphasis filtering 
 * 
 * Simple first order difference equation
 * s(n) = s(n) - f * s(n-1) 
 *
 * Copyright (C) 2007 by IRCAM-Centre Georges Pompidou, Paris, France.
 * 
 */
#include "rta_preemphasis.h"

/* can not be in place */
/* previous_sample updated */
void rta_preemphasis_signal(rta_real_t * out_samples,
                          const rta_real_t * in_samples, const unsigned int input_size,
                          rta_real_t * previous_sample, const rta_real_t factor)
{
  int i;
  
  if(factor != 0.)
  {
    out_samples[0] = in_samples[0] - factor * (*previous_sample);
    
    for(i=1; i<input_size; i++)
    {
      out_samples[i] = in_samples[i] - factor * in_samples[i-1];
    }
    
  }
  else
  {
    for(i=0; i<input_size; i++)
    {
      out_samples[i] = in_samples[i];
    }
  }

  *previous_sample = in_samples[input_size-1];

  return;
}

/* can not be in place */
/* previous_sample updated */
void rta_preemphasis_signal_stride(rta_real_t * out_samples, const int o_stride,
                                const rta_real_t * in_samples, const int i_stride,
                                const unsigned int input_size,
                                rta_real_t * previous_sample, const rta_real_t factor)
{
  int i,o;

  if(factor != 0.)
  {
    out_samples[0] = in_samples[0] - factor * (*previous_sample);

    for(i=i_stride, o=o_stride; i<input_size*i_stride; i+=i_stride, o+=o_stride)
    {
      out_samples[i] = in_samples[i] - factor * in_samples[i-1];
    }
  }
  else
  {
    for(i=0, o=0; i<input_size*i_stride; i+=i_stride, o+=o_stride)
    {
      out_samples[i] = in_samples[i];
    } 
  }

  *previous_sample = in_samples[i-i_stride];

  return;
}
