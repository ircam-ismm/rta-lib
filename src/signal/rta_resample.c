/**
 * @file   rta_resample.c
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Mon Nov 12 18:21:06 2007
 *
 * @brief  Resample utilities
 *
 * @copyright
 * Copyright (C) 2007 by IRCAM-Centre Georges Pompidou, Paris, France.
 * All rights reserved.
 *
 * License (BSD 3-clause)
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include <assert.h>
#include <string.h> // for memcpy
#include <math.h> // for floor

#include "rta_resample.h"
#include "rta_util.h"	// for idefix

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
                                    rta_real_t * output, const int o_stride,
                                    const rta_real_t * input, const int i_stride,
                                    const unsigned int i_size,
                                    const unsigned int factor)
{
  const rta_real_t factor_inv = 1. / factor;
  int o,i;
  const int o_max = (i_size / factor) * o_stride;
  const int i_incr = factor * i_stride;
  
  
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
        int ii;
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
                                      rta_real_t * output, const int o_stride,
                                      const rta_real_t * input, const int i_stride,
                                      const unsigned int i_size,
                                      const unsigned int factor)
{
  int o,i;
  const int o_max = (i_size / factor) * o_stride;
  const int i_incr = factor * i_stride;
  
  for(o=0, i=0; o<o_max; o+=o_stride, i+=i_incr)
  {
    output[o] = input[i];
  }
  
  return;
}


int rta_resample_cubic (rta_real_t * out_values,
                        const rta_real_t * in_values,
                        const unsigned int i_size,
                        const unsigned int out_max_size,
                        const unsigned int i_channels,
                        const double factor)
{
  int retValue = 0;

  rta_cubic_table_init(); // conditional initialization
  
  if (factor == 1.0)
  { /* copy through */
    memcpy(out_values, in_values, i_size * i_channels * sizeof(rta_real_t));
    retValue = i_size;
  }
  else if (in_values != out_values)
  {
    int m = i_size;
    int n = i_channels;
    double inv = 1.0 / factor;
    int maxOut       = out_max_size;
    int out_m        = (int) floor((double) (m - 1) * inv) + 1;
    int out_tailm2_m = (int) floor((double) (m - 2) * inv) + 1; //WAS: floor((double) (m - 2) * inv)
    
    /* limit resampling range here? */
    if (m > 3  &&  out_m > 0)
    {
      rta_idefix_t idefix; // fractional input frame position
      rta_idefix_t incr;   // fractional input frame increment
      int i, j, onset;
      
      if(out_m > maxOut)
        out_m = maxOut;
      if(out_tailm2_m > maxOut)
        out_tailm2_m = maxOut;

      rta_idefix_set_float(&incr, factor);

      for (j = 0; j < n; j++)
      { // for all columns/channels j
        rta_idefix_set_zero(&idefix);
        
        /* copy first points with linear interpolation */
        for (i = j; (onset = rta_idefix_get_index(idefix)) < RTA_CUBIC_HEAD  &&  i < out_m * n; i += n)
        {
          float frac  = rta_idefix_get_frac(idefix);
          float left  = in_values[j + onset * n];
          float right = in_values[j + onset * n + n];
          
          //out_values[i] = rta_cubic_calc_stride_head(in_values[j + onset] * n, ft, n); // cubic interpolation without accessing element at -n
          out_values[i] = left + (right - left) * frac; // linear interpolation
          rta_idefix_incr(&idefix, incr);
        }
	assert(onset >= RTA_CUBIC_HEAD  &&  (onset <= m - RTA_CUBIC_TAIL  ||  i >= out_tailm2_m * n)); // verify that input index has advanced, since cubic interpolation accesses input sample frame at idefix.index - RTA_CUBIC_HEAD
		
        for (; i < out_tailm2_m * n; i += n)
        {
          rta_cubic_idefix_interpolate_stride(in_values + j, idefix, n, out_values + i);
          rta_idefix_incr(&idefix, incr);
        }
        
        for (; (onset = rta_idefix_get_index(idefix)) < m - 1  &&  i < out_m * n; i += n)
        {
          float frac  = rta_idefix_get_frac(idefix);
          float left  = in_values[j + onset * n];
          float right = in_values[j + onset * n + n];

          //out_values[i] = rta_cubic_calc_stride_head(in_values[j + onset] * n, ft, n);
          out_values[i] = left + (right - left) * frac;
          rta_idefix_incr(&idefix, incr);
        }
	assert(((onset == m - 2  &&  rta_idefix_get_frac(idefix) > 0)  ||
		(onset >= m - 1))  ||  factor > 2.); // we have considered all input samples (for sane factors only)

        for (; i < out_m * n; i += n) // in case fractional input position idefix had a slight overshoot
	  out_values[i] = in_values[j + onset * n]; // fill with last input sample
	
	assert(i == out_m * n + j); // we have reached the end of output (all samples filled)
      }
      retValue = out_m;
    }
  }
  
  return retValue;
}


/** EMACS **
 * Local variables:
 * mode: c
 * c-basic-offset:2
 * End:
 */
