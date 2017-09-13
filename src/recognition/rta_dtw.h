/**
 * @file rta_dtw.h
 * @author Baptiste Caramiaux
 * @date 08/10/08
 * @ingroup rta_recognition
 *
 * @brief Compute a Dynamic Time Warping
 *
 * @copyright
 * Copyright (C) 2008 - 2009 by IRCAM-Centre Georges Pompidou, Paris, France.
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

#ifndef _RTA_DTW_H_
#define _RTA_DTW_H_ 1

#include "rta.h"

#ifdef __cplusplus
extern "C" {
#endif

// #include "mnm.h"


	/*
	 * Compute temporal alignment between two signals
	 *
	 * Parameters :
	 *  @left_ptr -> pointer to the the first signal, i.e. two dimensional array
	 *  @left_m -> its first dimension, number of observations
	 *  @left_n -> its second dimension, number of variables
	 *  @right_ptr -> pointer to the the second signal, i.e. two dimensional array
	 *  @right_m -> its first dimension, number of observations
	 *  @right_n -> its second dimension, number of variables
	 *
	 *  @output_p -> first signal indices for alignment
	 *  @output_q -> second signal indices for alignment
	 *		p[i] is the sample index of the first signal corresponding
	 *		to the sample index q[i] of the second signal
	 *
	 *  @length -> new signal length
	 */
	int
	rta_dtw(float * left_ptr, int left_m, int left_n, float * right_ptr, int right_m, int right_n,
			float * output_p, float * output_q, float * output_A, float * output_B, float * output_SM, int * length);





#ifdef __cplusplus
}
#endif

#endif /* _RTA_CCA_H_ */