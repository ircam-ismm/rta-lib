/**
 * @file   rta_svd.c
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Mon Aug 18 09:58:20 2008
 * 
 * @brief  Singular Value Decomposition
 * 
 * From the TNT/Jama package jama_svd.h (Adapted from JAMA, a Java
 * Matrix Library, developed by jointly by the Mathworks and NIST; see
 * http://math.nist.gov/javanumerics/jama). 
 * 
 *  Copyright (C) 2008 by IRCAM-Centre Georges Pompidou, Paris, France.
 */

#include "rta_svd.h"
#include "rta_math.h" /* rta_abs, rta_max, rta_min, rta_hypot, rta_pow */
#include "rta_int.h" /* rta_imin, rta_imax */
#include "rta_float.h" /* RTA_REAL_EPSILON */
#include "rta_stdlib.h" /* NULL */

struct rta_svd_setup
{
  rta_svd_t svd_type;

  /* A is copied when svd_type is 'rta_svd_out_of_place'
     or n > m (transposition) */
  rta_real_t * A; /* matrix of size m x n  */
  unsigned int m;
  unsigned int n;

  /* internal workspaces */
  rta_real_t * e; /* vector of size min(m,n) */
  rta_real_t * work; /* vector of size max(m,n) */
};

int
rta_svd_setup_new(rta_svd_setup_t ** svd_setup, const rta_svd_t svd_type,
                  rta_real_t * U, rta_real_t * S, rta_real_t *  V, 
                  rta_real_t * A, const unsigned int m, const unsigned int n)
{
  int ret = 1;
  *svd_setup = (rta_svd_setup_t *) rta_malloc(sizeof(rta_svd_setup_t));
  
  if(*svd_setup == NULL)
  {
    ret = 0;
  }
  else
  {
    (*svd_setup)->svd_type = svd_type;
    (*svd_setup)->m = m;
    (*svd_setup)->n = n;
  }

  if(ret != 0)
  {
    if(svd_type == rta_svd_out_of_place || n > m)
    {
      (*svd_setup)->A = rta_malloc(m * n * sizeof(rta_real_t));
      if((*svd_setup)->A == NULL)
      {
        ret = 0;
      }
    }
    else
    {
      (*svd_setup)->A = NULL;
    }
  }

  if(ret != 0)
  {
    (*svd_setup)->e = (rta_real_t *) rta_malloc(
      rta_imin(m,n) * sizeof(rta_real_t));
    if((*svd_setup)->e == NULL)
    {
      ret = 0;
    }
  }
    
  if(ret != 0)
  {
    (*svd_setup)-> work = (rta_real_t *) rta_malloc(
      rta_imax(m,n) * sizeof(rta_real_t));
    if((*svd_setup)->work == NULL)
    {
      ret = 0;
    }
  }
  
  if(ret == 0)
  {
    rta_svd_setup_delete(*svd_setup);
  }

  return ret;
}

void
rta_svd_setup_delete(rta_svd_setup_t * svd_setup)
{
  if(svd_setup != NULL)
  {
    if((svd_setup)->A != NULL)
    {
      rta_free((svd_setup)->A);
    }

    if((svd_setup)->e != NULL)
    {
      rta_free((svd_setup)->e);
    }

    if((svd_setup)->work != NULL)
    {
      rta_free((svd_setup)->work);
    }

    rta_free(svd_setup);
  }

  return;
}

/* A = U * S * V' */

/* A is a 2D array of size m x n  */
/* U is a 2D array of size m x (m,n) */
/* S is a 1D array of size min(m,n) */
/* V is a 2D array of size n x min(n,m) */

/* 2D arrays are in row-major order */
/* A ican be modified by the computation */
/* U and V can be NULL and are not computed, then */

/* e is a 1D array of size n */
/* work is a 1D array of size m */


void
rta_svd(rta_real_t * output_U, rta_real_t * S, rta_real_t *  output_V, 
        rta_real_t * input_A, const rta_svd_setup_t * svd_setup)
{
  rta_real_t * A; /* input_A, copied or transposed into svd_setup->A */
  rta_real_t * U; /* swap with V if A is transposed */
  rta_real_t * V;
  
  rta_real_t * e = svd_setup->e; /* just to ease the reading */
  rta_real_t * work = svd_setup->work; /* just to ease the reading */

  unsigned int m = svd_setup->m;
  unsigned int n = svd_setup->n;

  int nu;
  int nct;
  int nrt;
 
  int i, j, k;
  int p, pp, iter;
  
  if(n <= m)
  {
    if(svd_setup->svd_type == rta_svd_out_of_place)
    {
      /* Use an input copy */
      A = svd_setup->A;
      j = m*n;
      for(i = 0; i<j; i++)
      {
        A[i] = input_A[i];
      }
    }
    else /* Work directly on input */
    {
      A = input_A;
    }

    U = output_U;
    V = output_V;
  }
  else
  {
    /* Use an input transposed copy */
    A = svd_setup->A;

    for(i = 0; i<m; i++)
    {
      for(j=0; j<n; j++)
      {
        A[i*n + j] = input_A[j*m + i];
      }
    }
    m = svd_setup->n;
    n = svd_setup->m;

    /* swap U and V as A is transposed */
    U = output_V;
    V = output_U;
  }

  nu = rta_imin(m,n);
  nct = rta_imin(m-1,n);
  nrt = rta_imax(0,rta_imin(n-2,m));

  /* Reduce A to bidiagonal form, storing the diagonal elements */
  /* in s and the super-diagonal elements in e. */
  for (k = 0; k < rta_imax(nct,nrt); k++) 
  {
    if (k < nct) 
    {
      /* Compute the transformation for the k-th column and */
      /* place the k-th diagonal in S[k]. */
      /* Compute 2-norm of k-th column without under/overflow. */
      S[k] = 0.0;
      for (i = k; i < m; i++) 
      {
        S[k] = rta_hypot(S[k],A[i + k*m]);
      }
      if (S[k] != 0.0) 
      {
        if (A[k + k*m] < 0.0) 
        {
          S[k] = -S[k];
        }
        for (i = k; i < m; i++) 
        {
          A[i + k*m] /= S[k];
        }
        A[k + k*m] += 1.0;
      }
      S[k] = -S[k];
    }
    for (j = k+1; j < n; j++) 
    {
      if ((k < nct) && (S[k] != 0.0))  
      {
        /* Apply the transformation. */
        rta_real_t t = 0.0;
        for (i = k; i < m; i++) 
        {
          t += A[i + k*m]*A[i + j*m];
        }
        t = -t/A[k + k*m];
        for (i = k; i < m; i++) 
        {
          A[i + j*m] += t*A[i + k*m];
        }
      }

      /* Place the k-th row of A into e for the */
      /* subsequent calculation of the row transformation. */
      e[j] = A[k + j*m];
    }
    if (U != NULL && (k < nct)) 
    {
      /* Place the transformation in U for subsequent back */
      /* multiplication. */
      for (i = k; i < m; i++) 
      {
        U[i + k*m] = A[i + k*m];
      }
    }
    if (k < nrt) 
    {
      /* Compute the k-th row transformation and place the */
      /* k-th super-diagonal in e[k]. */
      /* Compute 2-norm without under/overflow. */
      e[k] = 0.0;
      for (i = k+1; i < n; i++) 
      {
        e[k] = rta_hypot(e[k],e[i]);
      }
      if (e[k] != 0.0) 
      {
        if (e[k+1] < 0.0) 
        {
          e[k] = -e[k];
        }
        for (i = k+1; i < n; i++) 
        {
          e[i] /= e[k];
        }
        e[k+1] += 1.0;
      }
      e[k] = -e[k];
      if ((k+1 < m) && (e[k] != 0.0)) 
      {
        /* Apply the transformation. */
        for (i = k+1; i < m; i++) 
        {
          work[i] = 0.0;
        }
        for (j = k+1; j < n; j++) 
        {
          for (i = k+1; i < m; i++) 
          {
            work[i] += e[j]*A[i + j*m];
          }
        }
        for (j = k+1; j < n; j++) 
        {
          rta_real_t t = -e[j]/e[k+1];
          for (i = k+1; i < m; i++) 
          {
            A[i + j*m] += t*work[i];
          }
        }
      }
      if (V != NULL) 
      {
        /* Place the transformation in V for subsequent */
        /* back multiplication. */
        for (i = k+1; i < n; i++) 
        {
          V[i + k*n] = e[i];
        }
      }
    }
  }

  /* Set up the final bidiagonal matrix or order p. */
  p = rta_imin(n,m+1);
  if (nct < n) 
  {
    S[nct] = A[nct + nct*m];
  }
  if (m < p) 
  {
    S[p-1] = 0.0;
  }
  if (nrt+1 < p) 
  {
    e[nrt] = A[nrt + (p-1)*m];
  }
  e[p-1] = 0.0;

  /* If required, generate U. */
  if (U != NULL) 
  {
    for (j = nct; j < nu; j++) 
    {
      for (i = 0; i < m; i++) 
      {
        U[i + j*m] = 0.0;
      }
      U[j + j*m] = 1.0;
    }
    for (k = nct-1; k >= 0; k--) 
    {
      if (S[k] != 0.0) 
      {
        for (j = k+1; j < nu; j++) 
        {
          rta_real_t t = 0.0;
          for (i = k; i < m; i++) 
          {
            t += U[i + k*m]*U[i + j*m];
          }
          t = -t/U[k + k*m];
          for (i = k; i < m; i++) 
          {
            U[i + j*m] += t*U[i + k*m];
          }
        }
        for (i = k; i < m; i++ ) 
        {
          U[i + k*m] = -U[i + k*m];
        }
        U[k + k*m] = 1.0 + U[k + k*m];
        for (i = 0; i < k-1; i++) 
        {
          U[i + k*m] = 0.0;
        }
      } 
      else 
      {
        for (i = 0; i < m; i++) 
        {
          U[i + k*m] = 0.0;
        }
        U[k + k*m] = 1.0;
      }
    }
  }

  /* If required, generate V. */
  if (V != NULL) 
  {
    for (k = n-1; k >= 0; k--) 
    {
      if ((k < nrt) && (e[k] != 0.0)) 
      {
        for (j = k+1; j < nu; j++) 
        {
          rta_real_t t = 0.0;
          for (i = k+1; i < n; i++) 
          {
            t += V[i + k*n]*V[i + j*n];
          }
          t = -t/V[(k+1) + k*n];
          for (i = k+1; i < n; i++) 
          {
            V[i + j*n] += t*V[i + k*n];
          }
        }
      }
      for (i = 0; i < n; i++) 
      {
        V[i + k*n] = 0.0;
      }
      V[k + k*n] = 1.0;
    }
  }

  /* Main iteration loop for the singular values. */
  pp = p-1;
  iter = 0;

  while (p > 0) 
  {
    int k=0;
    int kase=0;

    /* Here is where a test for too many iterations would go. */

    /* This section of the program inspects for */
    /* negligible elements in the s and e arrays.  On */
    /* completion the variables kase and k are set as follows. */

    /* kase = 1     if s(p) and e[k-1] are negligible and k<p */
    /* kase = 2     if s(k) is negligible and k<p */
    /* kase = 3     if e[k-1] is negligible, k<p, and */
    /*              s(k), ..., s(p) are not negligible (qr step). */
    /* kase = 4     if e(p-1) is negligible (convergence). */

    for (k = p-2; k >= -1; k--) 
    {
      if (k == -1) 
      {
        break;
      }
      if (rta_abs(e[k]) <= RTA_REAL_EPSILON * (rta_abs(S[k]) + rta_abs(S[k+1]))) 
      {
        e[k] = 0.0;
        break;
      }
    }
    if (k == p-2) 
    {
      kase = 4;
    } 
    else 
    {
      int ks;
      rta_real_t t;
      for (ks = p-1; ks >= k; ks--) 
      {
        if (ks == k) 
        {
          break;
        }
        t = (ks != p ? rta_abs(e[ks]) : 0.) + (ks != k+1 ? rta_abs(e[ks-1]) : 0.);
        if (rta_abs(S[ks]) <= RTA_REAL_EPSILON * t)  
        {
          S[ks] = 0.0;
          break;
        }
      }
      if (ks == k) 
      {
        kase = 3;
      } 
      else if (ks == p-1) 
      {
        kase = 1;
      } 
      else 
      {
        kase = 2;
        k = ks;
      }
    }
    k++;

    /* Perform the task indicated by kase. */
    switch (kase) 
    {
      /* Deflate negligible s(p). */
      case 1: 
      {
        rta_real_t f = e[p-2];
        e[p-2] = 0.0;
        for (j = p-2; j >= k; j--) 
        {
          rta_real_t t = rta_hypot(S[j],f);
          rta_real_t cs = S[j]/t;
          rta_real_t sn = f/t;
          S[j] = t;
          if (j != k) 
          {
            f = -sn*e[j-1];
            e[j-1] = cs*e[j-1];
          }
          if (V != NULL) 
          {
            for (i = 0; i < n; i++) 
            {
              t = cs*V[i + j*n] + sn*V[i + (p-1)*n];
              V[i + (p-1)*n] = -sn*V[i + j] + cs*V[i + (p-1)*n];
              V[i + j*n] = t;
            }
          }
        }
      }
      break;

      /* Split at negligible s(k). */
      case 2: 
      {
        rta_real_t f = e[k-1];
        e[k-1] = 0.0;
        for (j = k; j < p; j++) 
        {
          rta_real_t t = rta_hypot(S[j],f);
          rta_real_t cs = S[j]/t;
          rta_real_t sn = f/t;
          S[j] = t;
          f = -sn*e[j];
          e[j] = cs*e[j];
          if (U != NULL) 
          {
            for (i = 0; i < m; i++) 
            {
              t = cs*U[i + j*m] + sn*U[i + (k-1)*m];
              U[i*n + (k-1)*m] = -sn*U[i + j*m] + cs*U[i + (k-1)*m];
              U[i + j*m] = t;
            }
          }
        }
      }
      break;

      /* Perform one qr step. */
      case 3: 
      {
        /* Calculate the shift. */
        rta_real_t scale = 
          rta_max(rta_max(rta_max(rta_max(rta_abs(S[p-1]), rta_abs(S[p-2])),
                                  rta_abs(e[p-2])), rta_abs(S[k])),rta_abs(e[k]));
        rta_real_t sp = S[p-1]/scale;
        rta_real_t spm1 = S[p-2]/scale;
        rta_real_t epm1 = e[p-2]/scale;
        rta_real_t sk = S[k]/scale;
        rta_real_t ek = e[k]/scale;
        rta_real_t b = ((spm1 + sp)*(spm1 - sp) + epm1*epm1)/2.0;
        rta_real_t c = (sp*epm1)*(sp*epm1);
        rta_real_t shift = 0.0;
        rta_real_t f;
        rta_real_t g;

        if ((b != 0.0) || (c != 0.0)) 
        {
          shift = sqrt(b*b + c);
          if (b < 0.0) 
          {
            shift = -shift;
          }
          shift = c/(b + shift);
        }
        f = (sk + sp)*(sk - sp) + shift;
        g = sk*ek;
   
        /* Chase zeros. */
        for (j = k; j < p-1; j++) 
        {
          rta_real_t t = rta_hypot(f,g);
          rta_real_t cs = f/t;
          rta_real_t sn = g/t;
          if (j != k) 
          {
            e[j-1] = t;
          }
          f = cs*S[j] + sn*e[j];
          e[j] = cs*e[j] - sn*S[j];
          g = sn*S[j+1];
          S[j+1] = cs*S[j+1];
          if (V != NULL) 
          {
            for (i = 0; i < n; i++) 
            {
              t = cs*V[i + j*n] + sn*V[i + (j+1)*n];
              V[i + (j+1)*n] = -sn*V[i + j*n] + cs*V[i + (j+1)*n];
              V[i + j*n] = t;
            }
          }
          t = rta_hypot(f,g);
          cs = f/t;
          sn = g/t;
          S[j] = t;
          f = cs*e[j] + sn*S[j+1];
          S[j+1] = -sn*e[j] + cs*S[j+1];
          g = sn*e[j+1];
          e[j+1] = cs*e[j+1];
          if (U != NULL && (j < m-1)) 
          {
            for (i = 0; i < m; i++) 
            {
              t = cs*U[i + j*m] + sn*U[i + (j+1)*m];
              U[i + (j+1)*m] = -sn*U[i + j*m] + cs*U[i + (j+1)*m];
              U[i + j*m] = t;
            }
          }
        }
        e[p-2] = f;
        iter = iter + 1;
      }
      break;

      /* Convergence. */
      case 4: 
      {
        /* Make the singular values positive. */
        if (S[k] <= 0.0) 
        {
          S[k] = (S[k] < 0.0 ? -S[k] : 0.0);
          if (V != NULL) 
          {
            for (i = 0; i <= pp; i++) 
            {
              V[i + k*n] = -V[i + k*n];
            }
          }
        }
   
        /* Order the singular values. */
        while (k < pp) 
        {
          rta_real_t t;
          if (S[k] >= S[k+1]) 
          {
            break;
          }
          t = S[k];
          S[k] = S[k+1];
          S[k+1] = t;
          if (V != NULL && (k < n-1)) 
          {
            for (i = 0; i < n; i++) 
            {
              t = V[i + (k+1)*n];
              V[i + (k+1)*n] = V[i + k*n];
              V[i + k*n] = t;
            }
          }
          if (U != NULL && (k < m-1)) 
          {
            for (i = 0; i < m; i++) 
            {
              t = U[i + (k+1)*m];
              U[i + (k+1)*m] = U[i + k*m];
              U[i + k*m] = t;
            }
          }
          k++;
        }
        iter = 0;
        p--;
      }
      break;
    }
  }
  return;
}
