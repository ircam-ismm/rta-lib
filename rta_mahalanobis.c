/* $Id$
 *
 * FTM - RTA - Distlib
 * Copyright (C) 1994, 1995, 1998, 1999, 2007 by IRCAM-Centre Georges Pompidou, Paris, France.
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
        rta_real_t *outcol    = outptr;

        /* for each mu / output row i */
        for (i = 0; i < C; i++, outcol += outstride)
        {
          rta_real_t *inrow     = inptr;
          rta_real_t *murow     = muptr    + i * muskip;
          rta_real_t *sigmarow  = sigmaptr + i * sigmaskip;
          rta_real_t  v = 0.f;

          /* for each input column j: calculate
             y(i, k) = sum(j=1..N) (x(k, j) - mu(i, j))^2 / sigma(i, j)^2
          */
          for (j = 0; j < N; j++, inrow    += instride, 
                                  murow    += mustride, 
                                  sigmarow += sigmastride)
          {
              rta_real_t x = (*inrow - *murow) / *sigmarow;
              v += x * x;

/*            rta_post("k %d  i %d  j %d,  inrow %d (%p)  outcol %d (%p)  murow %d (%p) sigmarow %d (%p) -- x %f  v %f\n", 
              k, i, j, inrow - inptr, inrow, outcol - outptr, outcol, murow - muptr, murow, sigmarow - sigmaptr, sigmarow, x, v); */
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
 * @return success
 */

int rta_mahalanobis_nz(int M, int N, int C,
                       rta_real_t *inptr,    int instride,    int inskip, 
                       rta_real_t *muptr,    int mustride,    int muskip, 
                       rta_real_t *sigmaptr, int sigmastride, int sigmaskip, 
                       rta_real_t *outptr,   int outstride,   int outskip,
                       int nnz, int *sigma_indnz, void *distfuncs)
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

          /* for each NON-zero-sigma input column jj: calculate
             y(i, k) = sum(j=1..N) (x(k, j) - mu(i, j))^2 / sigma(i, j)^2
          */
          for (j = 0; j < nnz; j++)
          {
              int        jj    = sigma_indnz[j];
#if USE_DISTFUNC        // ARGH!!! dependency on fts_array_t and bpf_t!!!
              //TODO: replace by rta_funclib, includes bpf (data-compatible?)
              rta_real_t  d    = inrow[jj * instride] - murow[jj * mustride];
              fts_atom_t *dfun = fts_array_get_element((fts_array_t *) distfuncs, jj); //TODO: preget
              if (fts_is_a(dfun, bpf_class))
                  d = bpf_get_interpolated(fts_get_object(dfun), d);
              d /= sigmarow[jj * sigmastride];
#else
              rta_real_t d  = (inrow[jj * instride] - murow[jj * mustride]) 
                              / sigmarow[jj * sigmastride];
#endif /* USE_DISTFUNC */
              v += d * d;
          }

          *outcol = v;
        }
    }

    return 1;
}
