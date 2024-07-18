/* -*- mode: c; c-basic-offset:2 -*- */
/**
 * @file   rta_histogram.c
 * @author Diemo Schwarz
 * 
 * @brief  Calculate histogram
 *
 * @copyright
 * Copyright (C) 2008 - 2021 by IRCAM-Centre Georges Pompidou, Paris, France.
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

#include <float.h>
#include "rta_histogram.h"


// init parameter struct to default values
void rta_histogram_init (rta_histogram_params_t *params)
{
  params->nhist    = 100;
  params->lo_given = 0;
  params->hi_given = 0;
  params->lo       = 0;
  params->hi       = 0;
  params->norm     = 0;
}

/* Calculate histogram */
void rta_histogram_stride (rta_histogram_params_t *params,
			   rta_real_t *input, const unsigned int i_stride, const unsigned int i_size,
			   rta_real_t *output, const unsigned int out_stride,
			   rta_real_t *bpfout, const unsigned int bpf_stride)
{
  rta_real_t one = 1;
    
  rta_histogram_weighted_stride(params,
				input, i_stride, i_size,
				&one, 0, // unweighted: all weights == 1
				output, out_stride,
				bpfout, bpf_stride); 
}

void rta_histogram_stride_multi (rta_histogram_params_t *params, int num_input,
				 rta_real_t *input[],  const unsigned int i_stride, const unsigned int i_size[],
				 rta_real_t *output,   const unsigned int out_stride,
				 rta_real_t *bpfout,   const unsigned int bpf_stride)
{
  rta_real_t one = 1;
  rta_real_t *ones[num_input]; // array of pointers to weights data

  for (int i = 0; i < num_input; i++)
    ones[i] = &one; // make then all point to 1, zero stride assures we stay there

  rta_histogram_weighted_stride_multi(params, num_input,
				      input, i_stride, i_size,
				      ones, 0, // unweighted: all weights == 1
				      output, out_stride,
				      bpfout, bpf_stride); 
}

/* Calculate weighted histogram */
void rta_histogram_weighted_stride (rta_histogram_params_t *params,
				    rta_real_t *input, const unsigned int i_stride, const unsigned int i_size,
				    rta_real_t *weights, const unsigned int w_stride,
				    rta_real_t *output, const unsigned int out_stride,
				    rta_real_t *bpfout, const unsigned int bpf_stride)
{
  rta_histogram_weighted_stride_multi(params, 1,
				      &input, i_stride, &i_size,
				      &weights, w_stride, // unweighted: all weights == 1
				      output, out_stride,
				      bpfout, bpf_stride); 
}

void rta_histogram_weighted_stride_multi (rta_histogram_params_t *params, int num_input,
					  rta_real_t *input[],   const unsigned int i_stride, const unsigned int i_size[],
					  rta_real_t *weights[], const unsigned int w_stride,
					  rta_real_t *output,    const unsigned int out_stride,
					  rta_real_t *bpfout,    const unsigned int bpf_stride)
{
  rta_real_t *ptr;
  int i, j, k;

  /* clear result matrix */
  ptr = output; 
  for (i = 0; i < params->nhist; i++, ptr += out_stride)
    *ptr = 0;

  if ((ptr = bpfout)) // clear bpf domain output if given
    for (i = 0; i < params->nhist; i++, ptr += bpf_stride)
      *ptr = 0;

  // calculate and check input data
  unsigned int numdata = 0;
  for (k = 0; k < num_input; k++)
    numdata += i_size[k];
  
  if (numdata == 0)
    return;	// no data

  /* find min/max if not given by attributes */
  if (!params->lo_given  &&  !params->hi_given)
  {
    params->lo = FLT_MAX;                /* lower histogram limit */
    params->hi = -FLT_MAX;                /* upper histogram limit */

    for (k = 0; k < num_input; k++)
      for (i = 0; i < i_size[k] * i_stride; i += i_stride)
      {
	float x = input[k][i];
	if (x < params->lo)
	  params->lo = x;
	if (x > params->hi)
	  params->hi = x;
      }
  }
  else
  {
    if (!params->lo_given)
    {
      params->lo = FLT_MAX;                /* lower histogram limit */
      
      for (k = 0; k < num_input; k++)
	for (i = 0; i < i_size[k] * i_stride; i += i_stride)
	{
	  float x = input[k][i];
	  if (x < params->lo)
	    params->lo = x;
	}
    }
    
    if (!params->hi_given)
    {
      params->hi = -FLT_MAX;                /* upper histogram limit */
      
      for (k = 0; k < num_input; k++)
	for (i = 0; i < i_size[k] * i_stride; i += i_stride)
	{
	  float x = input[k][i];
	  if (x > params->hi)
	    params->hi = x;
	}
    }
  }

  float xfact = params->nhist / (params->hi - params->lo + 1); // bin step 
    
  if (bpfout)
  { /* create bin values */
    float bfact = 1 / xfact;
    
    for (i = 0; i < params->nhist; i++)
    {
      *bpfout = i * bfact + params->lo;
      bpfout += bpf_stride;
    }
  }

  /* calculate histogram */
  for (k = 0; k < num_input; k++)
    for (i = 0, j = 0; i < i_size[k] * i_stride; i += i_stride, j += w_stride)
    {
      /* find bin index */
      int ind = (input[k][i] - params->lo) * xfact;

      /* clip in case bounds were given as attributes */
      if (ind < 0)
	ind = 0;
      if (ind >= params->nhist)
	ind = params->nhist - 1;

      output[ind * out_stride] += weights[k][j];
    }

  // normalise histogram  
  if (params->norm > 0)
  { 
    float normfact = 0;
    rta_real_t *histcol = output;
      
    if (params->norm == 1)
    { // find max
      for (i = 0; i < params->nhist; i++, histcol += out_stride)
	if (*histcol > normfact)
	  normfact = *histcol;
    }
    else if (params->norm == 2)
    { // calc sum
      for (i = 0; i < params->nhist; i++, histcol += out_stride)
	normfact += *histcol;
    }
    
    /* do normalisation */
    if (normfact != 1.0  &&  normfact != 0.0)
    {
      normfact = 1 / normfact;
      histcol  = output;
      
      for (i = 0; i < params->nhist; i++, histcol += out_stride)
	*histcol *= normfact;
    }
  }
}
