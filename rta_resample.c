/**
 * @file   rta_resample.c
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Mon Nov 12 18:21:06 2007
 * 
 * @brief  Resample utilities
 * 
 * Copyright (C) 2007 by IRCAM-Centre Georges Pompidou, Paris, France.
 * 
 */
#include "rta_resample.h"

/* contract: factor > 0; */
/*           o_size >= i_size / factor */
void rta_downsample_int_mean(rta_real_t * output,
                             const rta_real_t * input,
                             const unsigned int i_size,
                             const unsigned int factor)
{
  const rta_real_t factor_inv = 1. / factor;
  unsigned int i,j;
  const unsigned int i_max = i_size / factor;
  
  switch(factor)
  {
    case 1:
      for(i=0; i<i_max; i++)
      {
        output[i] = input[i];
      }
      break;
      
    case 2:
      for(i=0, j=0; i<i_max; i++, j+=factor)
      {
        output[i] = factor_inv * (input[j] + input[j+1]);
      }
      break;
      
    case 3:
      for(i=0, j=0; i<i_max; i++, j+=factor)
      {
        output[i] = factor_inv * (input[j] + input[j+1] + input[j+2]);
      }
      break;

    case 4:
      for(i=0, j=0; i<i_max; i++, j+=factor)
      {
        output[i] = factor_inv * (input[j] + input[j+1] + input[j+2] + input[j+3]);
      }
      break;

    case 5:
      for(i=0, j=0; i<i_max; i++, j+=factor)
      {
        output[i] = factor_inv * (input[j] + input[j+1] + input[j+2] + input[j+3] +
          input[j+4]);
      }
      break;

    case 6:
      for(i=0, j=0; i<i_max; i++, j+=factor)
      {
        output[i] = factor_inv * (input[j] + input[j+1] + input[j+2] + input[j+3] +
          input[j+4] + input[j+5]);
      }
      break;

    case 7:
      for(i=0, j=0; i<i_max; i++, j+=factor)
      {
        output[i] = factor_inv * (input[j] + input[j+1] + input[j+2] + input[j+3] +
          input[j+4] + input[j+5] + input[j+6]);
      }
      break;

    case 8:
      for(i=0, j=0; i<i_max; i++, j+=factor)
      {
        output[i] = factor_inv * (input[j] + input[j+1] + input[j+2] + input[j+3] +
          input[j+4] + input[j+5] + input[j+6] + input[j+7]);
      }
      break;

    default:
      for(i=0, j=0; i<i_max; i++, j+=factor)
      {
        unsigned int k;
        output[i] = input[j];
        for(k=1; k<factor; k++)
        {
          output[i] += input[j+k];
        }
        output[i] *= factor_inv;
      }
  }

  return;
}

/* contract: factor > 0; */
/*           o_size >= i_size / factor */
void rta_downsample_int_mean_stride(
  rta_real_t * output, const unsigned int o_stride,
  const rta_real_t * input, const unsigned int i_stride,
  const unsigned int i_size,
  const unsigned int factor)
{
  const rta_real_t factor_inv = 1. / factor;
  unsigned int o,i;
  const unsigned int o_max = (i_size / factor) * o_stride;
  const unsigned int i_incr = factor * i_stride;

  
  switch(factor)
  {

   case 1:
      for(o=0, i=0; o<o_max; o+=o_stride, i+=i_incr)
      {
        output[o] = input[i];
      }
      break;

   case 2:
      for(o=0, i=0; o<o_max; o+=o_stride, i+=i_incr)
      {
        output[o] = factor_inv * (input[i] + input[i+i_stride]);
      }
      break;

   case 3:
      for(o=0, i=0; o<o_max; o+=o_stride, i+=i_incr)
      {
        output[o] = factor_inv *
          (input[i] + input[i+i_stride] + input[i+2*i_stride]);
      }
      break;

   case 4:
      for(o=0, i=0; o<o_max; o+=o_stride, i+=i_incr)
      {
        output[o] = factor_inv *
          (input[i] + input[i+i_stride] + input[i+2*i_stride] + 
           input[i+3*i_stride]);
      }
      break;

   case 5:
      for(o=0, i=0; o<o_max; o+=o_stride, i+=i_incr)
      {
        output[o] = factor_inv *
          (input[i] + input[i+i_stride] + input[i+2*i_stride] + 
           input[i+3*i_stride] + input[i+4*i_stride]);
      }
      break;

   case 6:
      for(o=0, i=0; o<o_max; o+=o_stride, i+=i_incr)
      {
        output[o] = factor_inv *
          (input[i] + input[i+i_stride] + input[i+2*i_stride] + 
           input[i+3*i_stride] + input[i+4*i_stride] + input[i+5*i_stride]);
      }
      break;

   case 7:
      for(o=0, i=0; o<o_max; o+=o_stride, i+=i_incr)
      {
        output[o] = factor_inv *
          (input[i] + input[i+i_stride] + input[i+2*i_stride] + 
           input[i+3*i_stride] + input[i+4*i_stride] + input[i+5*i_stride] +
           input[i+6*i_stride]);
      }
      break;

   case 8:
      for(o=0, i=0; o<o_max; o+=o_stride, i+=i_incr)
      {
        output[o] = factor_inv *
          (input[i] + input[i+i_stride] + input[i+2*i_stride] + 
           input[i+3*i_stride] + input[i+4*i_stride] + input[i+5*i_stride] +
           input[i+6*i_stride] + input[i+7*i_stride]);
      }
      break;

    default:
      for(o=0, i=0; o<o_max; o+=o_stride, i+=i_incr)
      {
        unsigned int ii;
        output[o] = input[i];
        for(ii=i_stride; ii<i_incr; ii+=i_stride)
        {
          output[o] += input[i+ii];
        }
        output[o] *= factor_inv;
      }
  }

  return;
}

/* contract: factor > 0; */
/*           o_size >= i_size / factor */
void rta_downsample_int_remove(rta_real_t * output,
                               const rta_real_t * input,
                               const unsigned int i_size,
                               const unsigned int factor)
{
  unsigned int i,j;
  const unsigned int i_max = i_size / factor;
  
  for(i=0, j=0; i<i_max; i++, j+=factor)
  {
    output[i] = input[j];
  }

  return;
}

/* contract: factor > 0; */
/*           o_size >= i_size / factor */
void rta_downsample_int_remove_stride(
  rta_real_t * output, const unsigned int o_stride,
  const rta_real_t * input, const unsigned int i_stride,
  const unsigned int i_size,
  const unsigned int factor)
{
  unsigned int o,i;
  const unsigned int o_max = (i_size / factor) * o_stride;
  const unsigned int i_incr = factor * i_stride;
  
  for(o=0, i=0; o<o_max; o+=o_stride, i+=i_incr)
  {
    output[o] = input[i];
  }

  return;
}
