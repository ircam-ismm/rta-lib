/* $Id$ */
/**
@file	rta_mahalanobis.h
@author	Diemo Schwarz
@date	31.3.2008
@version	2.0 

@brief	Mahalanobis distance calculation
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
