/**
 *
 *
 * @file   rta_delta.c
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Thu Aug  2 18:39:26 2007
 * 
 * @brief  Delta (derivative for a sequence at a fixed sampling rate)
 * 
 * Simple linear slope. Each column (a scalar value during time) is
 * filtered separately.
 *
 * Copyright (C) 2007 by IRCAM-Centre Georges Pompidou, Paris, France.
 * 
 */

#include "rta_delta.h"
#include "rta_math.h"

/* <filter_size> should be odd and positive */
int rta_delta_weights(rta_real_t * weights_vector, const unsigned int filter_size)
{
  int i;
  rta_real_t filter_value;

  const rta_real_t half_filter_size = rta_floor(filter_size * 0.5);

  for(i=0, filter_value=-half_filter_size;
      i<filter_size;
      i++, filter_value+=1.)
  {
    weights_vector[i] = filter_value;
  }

  return 1;
}

/* <filter_size> should be odd and positive */
int rta_delta_weights_stride(rta_real_t * weights_vector, const int w_stride,
                          const unsigned int filter_size)
{
  int i;
  rta_real_t filter_value;

  const rta_real_t half_filter_size = rta_floor(filter_size * 0.5);

  for(i=0, filter_value=-half_filter_size;
      i<filter_size*w_stride;
      i+=w_stride, filter_value+=1.)
  {
    weights_vector[i] = filter_value;
  }

  return 1;
}


rta_real_t rta_delta_normalization_factor(const unsigned int filter_size)
{
  rta_real_t normalization = 0.;

  if(filter_size>0)
  {
    int i;
    const int half_filter_size = filter_size / 2;

    for(i=1; i<=half_filter_size; i++)
    {
      normalization += (rta_real_t) (i*i);
    }

    normalization = 0.5 / normalization;
  }
  return normalization;
}

void rta_delta(rta_real_t * delta, const rta_real_t * input_vector,
              const rta_real_t * weights_vector,
              const unsigned int filter_size)
{
  unsigned int i;
  
  *delta = 0.;

  for(i=0; i<filter_size; i++)
  {
    if(weights_vector[i] != 0.)
    {
      *delta += input_vector[i] * weights_vector[i];
    }
  }
  
  return;
}

void rta_delta_stride(rta_real_t * delta, 
                    const rta_real_t * input_vector, const int i_stride,
                    const rta_real_t * weights_vector, const int w_stride,
                    const unsigned int filter_size)
{
  unsigned int i;
  
  *delta = 0.;

  for(i=0; i<filter_size; i++)
  {
    if(weights_vector[i*w_stride] != 0.)
    {
      *delta += input_vector[i*i_stride] * weights_vector[i*w_stride];
    }
  }
  
  return;
}


void rta_delta_vector(rta_real_t * delta,
                    const rta_real_t * input_matrix, const unsigned int input_size,
                    const rta_real_t * weights_vector, const unsigned int filter_size)
{
  unsigned int i,j;
  
  for(j=0; j<input_size; j++)
  {
    delta[j] = 0.;
  }

  for(i=0; i<filter_size; i++)
  {
    if(weights_vector[i] != 0.) /* skip zeros */
    {    
      for(j=0; j<input_size; j++)
      {
        delta[j] += input_matrix[i*input_size+j] * weights_vector[i];
      }
    }
  }
  
  return;
}

void rta_delta_vector_stride(rta_real_t * delta, const int d_stride,
                       const rta_real_t * input_matrix, const int i_stride, 
                       const unsigned int input_size,
                       const rta_real_t * weights_vector, const int w_stride,
                       const unsigned int filter_size)
{
  int i,j;
  
  for(j=0; j<input_size*d_stride; j+=d_stride)
  {
    delta[j] = 0.;
  }

  for(i=0; i<filter_size; i++)
  {
    if(weights_vector[i*w_stride] != 0.) /* skip zeros */
    {    
      for(j=0; j<input_size; j++)
      {
        delta[j*d_stride] += input_matrix[(i*input_size+j)*i_stride] *
          weights_vector[i*w_stride];
      }
    }
  }
  
  return;
}
