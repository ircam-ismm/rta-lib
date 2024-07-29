/* -*- mode: c; c-basic-offset:2 -*- */
/**
 * @file   rta_histogram.h
 * @author Diemo Schwarz
 * @ingroup rta_statistics
 *
 * @brief Histogram from an input vector
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

#ifndef _RTA_HISTOGRAM_H_
#define _RTA_HISTOGRAM_H_ 1

#include "rta.h"
#include <stdbool.h>	/* to get bool type in C99 */

#ifdef __cplusplus
extern "C" {
#endif


/** parameters for Histogram calculation */
typedef struct
{
  int      nhist;     /**< number of histogram bins */
  bool	   lo_given, hi_given;	/**< if true, histogram limits are given below; if false, will be calculated from data */
  float    lo, hi;    /**< histogram limits (will be calculated and returned here when not given) */
  int      norm;      /**< norm mode: 0 = off, 1 = max, 2 = sum */
} rta_histogram_params_t;


// init parameter struct to default values
void rta_histogram_init (rta_histogram_params_t *self);

     
/**
 * Calculate histogram
 *
 * @param params	pointer to histogram parameter struct
 * @param input		pointer to input data (at least i_size * i_stride elements)
 * @param i_stride	stride for input data
 * @param i_size	number of input elements
 * @param output	pointer to output data (at least params->nhist * out_stride elements)
 * @param out_stride	stride for output data
 * @param binout	NULL or pointer to bin index output data (at least params->nhist * bin_stride elements)
 * @param bin_stride	stride for bin index data
 */
void rta_histogram_stride (rta_histogram_params_t *params,
			   rta_real_t *input,  const int i_stride, const unsigned int i_size,
			   rta_real_t *output, const int out_stride,
			   rta_real_t *binout, const int bin_stride);
       
/**
 * Calculate histogram over multiple blocks of data
 *
 * @param params	pointer to histogram parameter struct
 * @param num_input	number of blocks of input data
 * @param input		array[num_input] of pointers to blocks of input data (at least i_size[i] * i_stride elements)
 * @param i_offset	offset into each block of input data
 * @param i_stride	stride for input data
 * @param i_size	array[num_input] of number of elements for each block of input data
 * @param output	pointer to output data (at least params->nhist * out_stride elements)
 * @param out_stride	stride for output data
 * @param binout	NULL or pointer to bin index output data (at least params->nhist * bpf_stride elements)
 * @param bin_stride	stride for bin index data
 */
void rta_histogram_stride_multi (rta_histogram_params_t *params, int num_input,
				 rta_real_t *input[],  const int i_offset, const int i_stride, const unsigned int i_size[],
				 rta_real_t *output,   const int out_stride,
				 rta_real_t *binout,   const int bin_stride);

/**
 * Calculate weighted histogram
 * (bins don't count occurrences, but sum weights given with each data value)
 *
 * @param params	pointer to histogram parameter struct
 * @param input		pointer to input data (at least i_size * i_stride elements)
 * @param i_stride	stride for input data
 * @param i_size	number of input elements
 * @param weights	pointer to weights data (at least i_size * w_stride elements)
 * @param w_stride	stride for weights data
 * @param output	pointer to output data (at least params->nhist * out_stride elements)
 * @param out_stride	stride for output data
 * @param binout	NULL or pointer to bin index output data (at least params->nhist * bin_stride elements)
 * @param bin_stride	stride for bin index data
 */
void rta_histogram_weighted_stride (rta_histogram_params_t *params,
				    rta_real_t *input,   const int i_stride, const unsigned int i_size,
				    rta_real_t *weights, const int w_stride,
				    rta_real_t *output,  const int out_stride,
				    rta_real_t *binout,  const int bin_stride);

  
/**
 * Calculate weighted histogram on multiple input blocks
 * (bins don't count occurrences, but sum weights given with each data value)
 *
 * @param params	pointer to histogram parameter struct
 * @param num_input	number of blocks of input data
 * @param input		array[num_input] of pointers to blocks of input data (at least i_size[i] * i_stride elements)
 * @param i_offset	offset into each block of input data
 * @param i_stride	stride for input data
 * @param i_size	array[num_input] of number of input elements for each block
 * @param weights	array[num_input] of pointers to weights data (at least i_size * w_stride elements)
 * @param w_stride	stride for weights data
 * @param output	pointer to output data (at least params->nhist * out_stride elements)
 * @param out_stride	stride for output data
 * @param binout	NULL or pointer to bin index output data (at least params->nhist * bpf_stride elements)
 * @param bin_stride	stride for bin index data
 */
void rta_histogram_weighted_stride_multi (rta_histogram_params_t *params, int num_input,
					  rta_real_t *input[],   const int i_offset, const int i_stride, const unsigned int i_size[],
					  rta_real_t *weights[], const int w_stride,
					  rta_real_t *output,    const int out_stride,
					  rta_real_t *binout,    const int bin_stride);
    
#ifdef __cplusplus
}
#endif

#endif /* _RTA_HISTOGRAM_H_ */

