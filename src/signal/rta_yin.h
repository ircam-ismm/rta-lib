/**
 * @file   rta_yin.h
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Mon Nov 12 18:21:06 2007
 * 
 * @brief  Yin algorithm for periodicity analysis
 * 
 * Based on Norbert Schnell's real-time implementation after the YIN
 * algorithm by Alain de Cheveigné and Hideki Kawahara.
 * See "YIN, a fundamental frequency estimator for speech and music"
 * published in J. Acoust. Soc. Am. 111 (4), April 2002.
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

#ifndef _RTA_YIN_H_
#define _RTA_YIN_H_ 1

#include "rta.h"

#ifdef __cplusplus
extern "C" {
#endif

/* rta_yin_setup is private */
typedef struct rta_yin_setup rta_yin_setup_t;

/** 
 * Alocate a yin setup for further yin analysis.
 *
 * \see rta_yin_setup_delete
 * 
 * @param yin_setup is an address to a pointer to a private structure
 * which may depend on actual yin implementation.
 * @param max_mins is the maximum number of minima searched by the yin
 * algorithm. 128 is a standard value for an input vector of 2048 points.
 * 
 * @return 1 on success 0 on fail. If it fails, nothing should be done
 * with 'yin_setup' (even a delete).
 */
int
rta_yin_setup_new(rta_yin_setup_t ** yin_setup, unsigned int max_mins);

/** 
 * Deallocate any (sucessfully) allocated yin setup.
 *
 * \see rta_yin_setup_new
 *
 * 
 * @param yin_setup is a pointer to the memory wich will be released.
 */
void
rta_yin_setup_delete(rta_yin_setup_t * yin_setup);


/** 
 * Yin algorithm for periodicity analysis.
 *
 * To reduce the computational cost of this algorithm, you may want to
 * downsample the signal first:
 * \see rta_downsample_int_mean
 * \see rta_downsample_int_remove
 * 
 * @param abs_min is the absolute minimum found. Its range is [0., 1.].
 * periodicity == 1. - sqrt('abs_min').
 * @param autocorrelation must be allocated before calling this
 * function. 'autocorrelation' size is 'ac_size'.
 * The energy of 'input' is
 * sqrt('autocorrelation'[0]/('input_size' - 'ac_size')). 
 * \see rta_correlation_fast
 * @param ac_size is 'autocorrelation' size. The maximum
 * searched lag is 'ac_size' - 2. 'ac_size' must be < 'input_size' 
 * and for good results it should be <= ('input_size'/2)
 * @param input size is 'input_size'
 * @param input_size is the 'input' size. It must be > 'ac_size'
 * and for good results it should be >= 2 * 'ac_size'
 * @param yin_setup must be previously allocated. \see rta_yin_setup_new
 * @param threshold under which to search the minima. It must be in
 * [0., 1.]. 0.1 is a common value. 'threshold' == (1. - confidence)^2
 * 
 * @return the lag corresponding to 'abs_min'. The corresponding
 * frequency is then sample_rate / lag.
 */
rta_real_t
rta_yin(rta_real_t * abs_min, rta_real_t * autocorrelation, 
        const unsigned int ac_size, 
        const rta_real_t * input, const unsigned int input_size,
        const rta_yin_setup_t * yin_setup,
        const rta_real_t threshold);

/** 
 * Yin algorithm for periodicity analysis.
 *
 * To reduce the computational cost of this algorithm, you may want to
 * downsample the signal first:
 * \see rta_downsample_int_mean
 * \see rta_downsample_int_remove
 * 
 * @param abs_min is the absolute minimum found. Its range is [0., 1.].
 * periodicity == 1.0 - sqrt('abs_min')
 * @param autocorrelation must be allocated before calling this
 * function. 'autocorrelation' size is 'ac_size'.
 * The energy of 'input' is
 * sqrt('autocorrelation'[0]/('input_size' - 'ac_size')).
 * \see rta_correlation_fast_stride 
 * @param ac_stride is 'autocorrelation' stride
 * @param ac_size is 'autocorrelation' size. The maximum
 * searched lag is 'ac_size' - 2. 'ac_size' must be < 'input_size' 
 * and for good results it should be <= ('input_size'/2)
 * @param input size is 'input_size'
 * @param i_stride is 'input' stride
 * @param input_size is the 'input' size. It must be > 'ac_size'
 * and for good results it should be >= 2 * 'ac_size'
 * @param yin_setup must be previously allocated. \see rta_yin_setup_new
 * @param threshold under which to search the minima. It must be in
 * [0., 1.]. 0.1 is a common value. 'threshold' == (1. - confidence)^2
 * 
 * @return the lag corresponding to 'abs_min'. The corresponding
 * frequency is then sample_rate / lag.
 */
rta_real_t
rta_yin_stride(rta_real_t * abs_min, 
               rta_real_t * autocorrelation, const int ac_stride, 
               const unsigned int ac_size, 
               const rta_real_t * input, const int i_stride, 
               const unsigned int input_size,
               const rta_yin_setup_t * yin_setup,
               const rta_real_t threshold);

#ifdef __cplusplus
}
#endif

#endif /* _RTA_YIN_H_ */

