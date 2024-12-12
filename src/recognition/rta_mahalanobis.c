/**
 * @file rta_mahalanobis.c
 * @author Norbert Schnell
 *
 * @brief  Mahalanobis distance calculation
 *
 * @copyright
 * Copyright (C) 1994, 1995, 1998, 1999, 2007 by IRCAM-Centre Georges Pompidou, Paris, France.
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

#ifdef WIN32
#include <malloc.h>
#else
#include <alloca.h>
#endif
#include <math.h>

#include "rta_mahalanobis.h"


/* update non-zero sigma index list */
int rta_find_nz (int n, rta_real_t *vec, int stride, int *indnz)
{
  int j, nnz = 0;

  for (j = 0; j < n; j++, vec += stride)
    if (*vec != 0)
      indnz[nnz++] = j;

  return nnz;
}



/** Pure Mahalanbis distance calculation
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
 * @param M     num. rows of query matrix in = num. cols of out
 * @param N     num. dimensions = num. cols of in, mu, sigma
 * @param C     num. points = num. rows of mu and maybe sigma
 * @return success
 */

int rta_mahalanobis(int M, int N, int C,
                    rta_real_t *inptr,    int instride,    int inskip,
                    rta_real_t *muptr,    int mustride,    int muskip,
                    rta_real_t *sigmaptr, int sigmastride, int sigmaskip,
                    rta_real_t *outptr,   int outstride,   int outskip)
{
  int i, j, k;

  /*  rta_post("M %d N %d C %d,  inptr %p  instride %d  inskip %d,  muptr %p  mustride %d  muskip %d,  sigmaptr %p  sigmastride %d  sigmaskip %d,  outptr %p  outstride %d  outskip %d\n",
             M, N, C, inptr, instride, inskip, muptr, mustride, muskip,
             sigmaptr, sigmastride, sigmaskip, outptr, outstride, outskip); */

  /* for each input row k / output column k */
  for (k = 0; k < M; k++, inptr += inskip, outptr += outskip)
  {
    rta_real_t *outcol = outptr;

    /* for each mu / output row i */
    for (i = 0; i < C; i++, outcol += outstride)
    {
      rta_real_t *inrow = inptr;
      rta_real_t *murow = muptr  + i * muskip;
      rta_real_t *sigmarow = sigmaptr + i * sigmaskip;
      rta_real_t  v = 0.f;

      /*
        for each input column j: calculate
         y(i, k) = sum(j=1..N) (x(k, j) - mu(i, j))^2 / sigma(i, j)^2
      */
      for (j = 0; j < N; j++, inrow += instride,
                              murow += mustride,
                              sigmarow += sigmastride)
      {
        rta_real_t x = (*inrow - *murow) / *sigmarow;
        v += x * x;

        /*
          rta_post("k %d  i %d  j %d,  inrow %d (%p)  outcol %d (%p)  murow %d (%p) sigmarow %d (%p) -- x %f  v %f\n",
          k, i, j, inrow - inptr, inrow, outcol - outptr, outcol, murow - muptr, murow, sigmarow - sigmaptr, sigmarow, x, v);
        */
      }

      *outcol = v;
    }
  }

  return 1;
}



/** Mahalanbis distance calculation on non-zero dimensions
 *
 * out = sum(distfunc(in - mu) .^ 2 ./ sigma^2)
 *
 * with
 *
 * in    (M, N)
 * mu    (C, N)
 * sigma (C, N) or sigma(1, N)
 * out   (C, M)
 *
 * @param M     num. rows of query matrix in = num. cols of out
 * @param N     num. dimensions = num. cols of in, mu, sigma
 * @param C     num. points = num. rows of mu and maybe sigma
 *
 * @param distfuncs pointer to array of N distance function bpfs, or NULL
 * @return success
 */

int rta_mahalanobis_nz(int M, int N, int C,
                       rta_real_t *inptr,    int instride,    int inskip,
                       rta_real_t *muptr,    int mustride,    int muskip,
                       rta_real_t *sigmaptr, int sigmastride, int sigmaskip,
                       rta_real_t *outptr,   int outstride,   int outskip,
                       int nnz, int *sigma_indnz, rta_bpf_t *distfuncs[])
{
  int i, j, k;

  /* for each input row k / output column k */
  for (k = 0; k < M; k++, inptr += inskip, outptr += outskip)
  {
    rta_real_t *outcol    = outptr;

    /* for each mu / output row i */
    for (i = 0; i < C; i++, outcol += outstride)
    {
      rta_real_t *inrow     = inptr;
      rta_real_t *murow     = muptr    + i * muskip;
      rta_real_t *sigmarow  = sigmaptr + i * sigmaskip;
      rta_real_t  v = 0.f;

      /*
        for each NON-zero-sigma input column jj: calculate
        y(i, k) = sum(j=1..N) (x(k, j) - mu(i, j))^2 / sigma(i, j)^2
      */
      for (j = 0; j < nnz; j++)
      {
        int jj = sigma_indnz[j];
#if RTA_USE_DISTWARP // uses rta_bpf_t, (data-compatible to FTM bpfunc_t)
        rta_real_t  d    = inrow[jj * instride] - murow[jj * mustride];
        rta_bpf_t  *dfun = distfuncs[jj];

        if (dfun)
          d = rta_bpf_get_interpolated(dfun, d);
        d /= sigmarow[jj * sigmastride];
#else
        rta_real_t d  = (inrow[jj * instride] - murow[jj * mustride]) /
                        sigmarow[jj * sigmastride];
#endif /* RTA_USE_DISTWARP */
        v += d * d;
      }

      *outcol = v;
    }
  }

  return 1;
}
