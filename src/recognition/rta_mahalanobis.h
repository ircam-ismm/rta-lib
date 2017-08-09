/**
 * @file rta_mahalanobis.h
 * @author Diemo Schwarz
 * @date 31.3.2008
 * @version  2.0 
 *
 * @brief  Mahalanobis distance calculation
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

#ifndef _RTA_MAHALANOBIS_H_
#define _RTA_MAHALANOBIS_H_

#include "rta.h"
#include "rta_bpf.h"

#ifdef __cplusplus
extern "C" {
#endif


/* update non-zero sigma index list */
int rta_find_nz (int n, rta_real_t *vec, int stride, int *indnz);


/** Pure Mahalanobis distance calculation
 *
 * out = sum((in - mu) .^ 2 ./ sigma)
 *
 * with
 * 
 * in    (M, N)
 * mu    (C, N)
 * sigma (C, N) or sigma(1, N)
 * out   (C, M)
 *
 * @param M	num. rows of query matrix in = num. cols of out
 * @param N	num. dimensions = num. cols of in, mu, sigma
 * @param C	num. points = num. rows of mu and maybe sigma
 * @return success
 */

int rta_mahalanobis(int M, int N, int C,
		    rta_real_t *inptr,    int instride,    int inskip, 
		    rta_real_t *muptr,    int mustride,    int muskip, 
		    rta_real_t *sigmaptr, int sigmastride, int sigmaskip, 
		    rta_real_t *outptr,   int outstride,   int outskip);


/** Mahalanbis distance calculation on non-zero dimensions
 *
 * out = sum((in - mu) .^ 2 ./ sigma)
 *
 * with
 * 
 * in    (M, N)
 * mu    (C, N)
 * sigma (C, N) or sigma(1, N)
 * out   (C, M)
 *
 * @param M	num. rows of query matrix in = num. cols of out
 * @param N	num. dimensions = num. cols of in, mu, sigma
 * @param C	num. points = num. rows of mu and maybe sigma
 * @return success
 */

int rta_mahalanobis_nz(int M, int N, int C,
		       rta_real_t *inptr,    int instride,    int inskip, 
		       rta_real_t *muptr,    int mustride,    int muskip, 
		       rta_real_t *sigmaptr, int sigmastride, int sigmaskip, 
		       rta_real_t *outptr,   int outstride,   int outskip,
		       int nnz, int *sigma_indnz, rta_bpf_t *distfuncs[]);


#ifdef __cplusplus
}
#endif

#endif /* _RTA_MAHALANOBIS_H_ */
