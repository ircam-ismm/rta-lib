/**
 * @file rta_kdtreebuild.c
 * @author Riccardo Borghesi
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

#include "rta_kdtree.h"
#include "rta_kdtreeintern.h"
#include <stdlib.h>
#include <math.h>
#include <assert.h>


#ifndef WIN32
#include <strings.h>
#endif

#ifdef DEBUG
#define RTA_DEBUG_KDTREEBUILD 0
#else
#define RTA_DEBUG_KDTREEBUILD 0
#endif

static void compute_mean (rta_kdtree_t *t, int node, int dim)
{
  rta_kdtree_node_t *n    = &t->nodes[node];
  rta_real_t    *mean_ptr = t->mean + node * t->ndim;
  int          nstart     = n->startind;
  int          nend       = n->endind;
  int          nvector    = nend - nstart + 1; // number of vectors of processed node
  int          dstart, dend;
  int          i, j;

  if (dim < 0)
  {   /* all dimensions */
    dstart = 0;
    dend   = t->ndim;
  }
  else
  {   /* calculate mean only over given dimension */
    dstart = dim;
    dend   = dim + 1;
#if RTA_DEBUG_KDTREEBUILD
    bzero(mean_ptr, t->ndim * sizeof(rta_real_t));
#endif
  }

  for (j = dstart; j < dend; j++)
  {
    rta_real_t sum = 0;

    for (i = nstart; i <= nend; i++)
      sum += rta_kdtree_get_element(t, i, j);

    mean_ptr[j] = sum / nvector;
  }

#if RTA_DEBUG_KDTREEBUILD
  rta_post("mean vect for node %d (size %d) = ", node, nvector);
  rta_vec_post(mean_ptr, 1, t->ndim, "\n");
#endif
}


static void compute_middle (rta_kdtree_t *t, int node, int dim)
{
  rta_kdtree_node_t *n    = &t->nodes[node];
  rta_real_t    *mean_ptr = t->mean + node * t->ndim;
  int          nstart     = n->startind;
  int          nend       = n->endind;
  int          dstart, dend;
  int          i, j;

  if (dim < 0)
  {
    dstart = 0;
    dend   = t->ndim;
  }
  else
  {
    dstart = dim;
    dend   = dim + 1;
#if RTA_DEBUG_KDTREEBUILD
    bzero(mean_ptr, t->ndim * sizeof(rta_real_t));
#endif
  }

  // number of vectors from the processed node
  for (j = dstart; j < dend; j++)
  {
    rta_real_t min = MAX_FLOAT, max = -MAX_FLOAT;

    for (i = nstart; i <= nend; i++)
    {
      rta_real_t x = rta_kdtree_get_element(t, i, j);

      if (x < min)
        min = x;
      if (x > max)
        max = x;
    }

    mean_ptr[j]   = (max + min) / 2.;
  }

#if RTA_DEBUG_KDTREEBUILD
  rta_post("middle vect for node %d = ", node);
  rta_vec_post(mean_ptr, 1, t->ndim, "\n");
#endif
}

/* return 1 if node is well-behaved, 0 if node is degenerate */
static int check_node (rta_kdtree_t *t, int node, int dim)
{
  int        nstart   = t->nodes[node].startind;
  int        nend     = t->nodes[node].endind;
  rta_real_t min, max;
  int        i;

  if (nstart < nend)
    min = max = rta_kdtree_get_element(t, nstart, dim);
  else /* size == 1 */
    return 0;

  for (i = nstart + 1; i <= nend; i++)
  {
    rta_real_t x = rta_kdtree_get_element(t, i, dim);

    if (x < min)
      min = x;
    if (x > max)
      max = x;
  }

  return max != min;
}

/* compute and create node-splitting hyperplane */
static void compute_splitplane (rta_kdtree_t* t, int node, int level)
{
  rta_kdtree_node_t *n = &t->nodes[node];

 #if RTA_KDTREE_PROFILE_BUILD
  t->profile.hyperp++;
#endif
  switch (t->dmode)
  {
    case dmode_hyperplane:
    { /* compute hyperplane orthogonal to the base vector number b */
      rta_real_t *split_ptr = t->split + node * t->ndim;

#ifndef WIN32
      bzero(split_ptr, t->ndim * sizeof(rta_real_t));
#else
      memset(split_ptr, 0.0, t->ndim * sizeof(rta_real_t));
#endif
      split_ptr[n->splitdim] = 1;

#if RTA_DEBUG_KDTREEBUILD
      rta_post("Splitplane of node %i: ", node);
      rta_vec_post(split_ptr, 1, t->ndim, "\n");
#endif
    }

    /* FALLTHROUGH */
    case dmode_orthogonal:
     n->splitnorm = 1;  /* splitplane implicit orthogonal to splitdim */
  }
}


/* determine hyperplane that splits node n in two child nodes

   return: 1 if node is well-behaved, 0 if node is degenerate,
   i.e. all vectors have same distance (usually 0) to splitplane
*/
static int decompose_node (rta_kdtree_t *t, int node, int level, int use_sigma)
{
  int i, nice_node = 0, splitdim;

#if RTA_KDTREE_PROFILE_BUILD
  t->profile.mean++;
#endif
  /* determine dimension to split at
  (for the moment simply cycling through dimensions by tree level,
  skipping degenerate dimensions) */
  if (use_sigma  &&  t->sigma_nnz > 0)
  {
    splitdim = t->sigma_indnz[level % t->sigma_nnz];

    for (i = 1; i < t->sigma_nnz; i++) /* try each dim at most once */
    {
      if ((nice_node = check_node(t, node, splitdim)))
        break;

      /* step shifted by height, to use extra dimensions */
      splitdim = t->sigma_indnz[(level + i) % t->sigma_nnz];
    }
  }
  else
  {
    splitdim = level % t->ndim;

    for (i = 0; i < t->ndim; i++) /* try each dim at most once */
    {
      if ((nice_node = check_node(t, node, splitdim)))
        break;

      /* step shifted by height, to use extra dimensions */
      splitdim = (splitdim + t->height) % t->ndim;
    }
  }

  t->nodes[node].splitdim = splitdim;

#if RTA_DEBUG_KDTREEBUILD
  if (!nice_node)
    rta_post("warning: can't find non-degenerate dimension to split node %d at level %d, using dimension %d\n", node, level, splitdim);
#endif

  switch (t->mmode)
  { /* N.B.: middle and mean are only linearly  affected by sigma */
    case mmode_mean:
      if (t->dmode == dmode_orthogonal)
        compute_mean(t, node, splitdim);
      else
        compute_mean(t, node, -1);  /* all dimensions */
      break;
/*
      case mmode_median:
  compute_median(t, node);
      break;
*/
    case mmode_middle:
      if (t->dmode == dmode_orthogonal)
        compute_middle(t, node, splitdim);
      else
        compute_middle(t, node, -1);  /* all dimensions */
      break;
  }

  /* compute and create node splitting hyperplane */
  compute_splitplane(t, node, level);

  return nice_node;
}


/* vector to orthogonal plane node distance along split dimension dim */
static rta_real_t distV2orthoH (const rta_real_t* vect,
                                rta_real_t* mean, int dim,
                                rta_bpf_t  *distfunc[])
{
  return RTA_DMAP(vect[dim], mean[dim], distfunc[dim]);
}

static rta_real_t distV2orthoH_stride (const rta_real_t* vect, int stride,
                                       rta_real_t* mean, int dim,
                                       rta_bpf_t  *distfunc[])
{
  return RTA_DMAP(vect[dim * stride], mean[dim], distfunc[dim]);
}

static rta_real_t distV2orthoH_weighted (const rta_real_t* vect, int stride, rta_real_t* mean,
                                         const rta_real_t *sigma, int dim, rta_bpf_t *distfunc[])
{
#if RTA_DEBUG_KDTREEBUILD > 1
  rta_post("distV2orthoH_weighted on dim %d: (%f - %f) / %f = %f\n",
    dim, vect[dim * stride], mean[dim], sigma[dim],
    sigma[dim] > 0  ?  (vect[dim * stride] - mean[dim]) / sigma[dim]  :  0);
#endif
  return sigma[dim] > 0
       ? RTA_DMAPW(vect[dim * stride], mean[dim], sigma[dim], distfunc[dim])
       : 0;
}

/* vector to general plane node distance */
static rta_real_t distV2H (const rta_real_t* vect,
         const rta_real_t* plane,
         const rta_real_t* mean,
         int ndim, rta_real_t norm,
         rta_bpf_t *distfunc[])
{
  // standard algebra computing
  int i;
  rta_real_t dotprod = 0;

  for(i = 0; i < ndim; i++)
    dotprod += RTA_DMAP(vect[i], mean[i], distfunc[i]) * plane[i];

  return (dotprod / norm);
}

static rta_real_t distV2H_stride (const rta_real_t* vect, int stride,
                                  const rta_real_t* plane,
                                  const rta_real_t* mean,
                                  int ndim, rta_real_t norm,
                                  rta_bpf_t *distfunc[])
{
  // standard algebra computing
  int i, iv;
  rta_real_t dotprod = 0;

  for(i = 0, iv = 0; i < ndim; i++, iv += stride)
    dotprod += RTA_DMAP(vect[iv], mean[i], distfunc[i]) * plane[i];

  return (dotprod / norm);
}

static rta_real_t distV2H_weighted (const rta_real_t* vect, int stride,
                                    const rta_real_t* plane,
                                    const rta_real_t* mean,
                                    const rta_real_t *sigma,
                                    int ndim, rta_real_t norm,
                                    rta_bpf_t *distfunc[])
{
  // standard algebra computing
  int i, iv;
  rta_real_t dotprod = 0;

  for (i = 0, iv = 0; i < ndim; i++, iv += stride)
    if (sigma[i] > 0)
      dotprod += RTA_DMAPW(vect[iv], mean[i], sigma[i], distfunc[i]) * plane[i];

  return (dotprod / norm);
}


/* vector to node distance */
rta_real_t distV2N (rta_kdtree_t* t, const rta_real_t *x, const int node)
{
#if RTA_KDTREE_PROFILE_BUILD
  t->profile.v2n++;
#endif

  switch (t->dmode)
  {
    case dmode_orthogonal:
      return distV2orthoH(x, t->mean + node * t->ndim,
                          t->nodes[node].splitdim, t->dfun);
    case dmode_hyperplane:
      return distV2H(x, t->split + node * t->ndim,
                     t->mean  + node * t->ndim, t->ndim,
                     t->nodes[node].splitnorm, t->dfun);
    default:
      rta_post("error: unknown mode %d", t->dmode);
      return 0;
  }
}

rta_real_t distV2N_stride (rta_kdtree_t* t, const rta_real_t *x, int stride, const int node)
{
#if RTA_KDTREE_PROFILE_BUILD
  t->profile.v2n++;
#endif

  switch (t->dmode)
  {
    case dmode_orthogonal:
      return distV2orthoH_stride(x, stride, t->mean + node * t->ndim,
                                 t->nodes[node].splitdim, t->dfun);
    case dmode_hyperplane:
      return distV2H_stride(x, stride, t->split + node * t->ndim,
                            t->mean  + node * t->ndim, t->ndim,
                            t->nodes[node].splitnorm, t->dfun);
    default:
      rta_post("error: unknown mode %d", t->dmode);
      return 0;
  }
}

rta_real_t distV2N_weighted (rta_kdtree_t* t, const rta_real_t *x, int stride,
                             const rta_real_t *sigma, const int node)
{
  rta_kdtree_node_t *n = &t->nodes[node];
  rta_real_t *mean = t->mean + node * t->ndim;

#if RTA_KDTREE_PROFILE_BUILD
  t->profile.v2n++;
#endif
  switch (t->dmode)
  {
    case dmode_orthogonal:
      return distV2orthoH_weighted(x, stride, mean, sigma, n->splitdim, t->dfun);
    case dmode_hyperplane:
      return distV2H_weighted(x, stride, t->split + node * t->ndim, mean, sigma, t->ndim, n->splitnorm, t->dfun);
    default:
      rta_post("error: unknown mode %d", t->dmode);
      return 0;
  }
}


/* swap positions of vectors i and j: only in indirection array */
static void swap (rta_kdtree_t* t, int i, int j)
{
  rta_kdtree_object_t tmp;

  tmp = t->dataindex[i];
  t->dataindex[i] = t->dataindex[j];
  t->dataindex[j] = tmp;
}


/* (re-)insert num data vector into tree.
   If index < ndata, move to correct node,
   otherwise insert and increment ndata */
void rta_kdtree_insert (rta_kdtree_t* t, int index, int num)
{
    /* find nearest node to new vector */


    /* insert into that node: shift vector indices and node start/end
       in parent and right siblings */

    /* 1: move from insert node up to root */

    /* 2: traverse right siblings */
}


#if 0
/* regenerate nodes on already sorted data vectors */
void kdtree_rebuild (kdtree_t* t)
{
    int i = 0;
    int j = t->ndata - 1;

    /* Sort */
    while(j > i)
    {
  while((t->data[i][0] != MAX_FLOAT)&&(i < t->ndata))
  {
      i++;
  }
  while((t->data[j][0] == MAX_FLOAT)&&(j >=0 ))
  {
      j--;
  }
  if(j > i)
  {
      swap(t, i, j);
  }
    }
    t->ndata = i;
    t->nodes[0].endind = t->ndata - 1;

    kdtree_build(t, 1);
}
#endif


void rta_kdtree_build (rta_kdtree_t* t, int use_sigma)
{
  int l;    // current level number
  int n;    // current node number
  int i, j;   // loop counters

  /* Maximum length is equal to pow2(height-1) */
  if (pow2(t->height - 1) > t->ndatatot  ||  t->ndim == 0)
  {
    if (t->ndatatot == 0)
      rta_post("tree is empty!\n");
    else if (t->ndim == 0)
      rta_post("tree has 0 dimensions!  Can't build!\n");
    else
      rta_post("error: can't build this tree, try with a smaller height: %d > %d\n",
               cpow2(t->height-1), t->ndatatot);

    return;
  }

  for (l = 0; l < t->height - 1; l++)
  {   /* initialise inner nodes */
    int nstart = pow2(l)   - 1;
    int nend   = pow2(l+1) - 1;
#if RTA_DEBUG_KDTREEBUILD
    rta_post("\nLevel #%i  nodes %d..%d\n", l, nstart, nend);
#endif

    for (n = nstart; n < nend; n++)
    {   /* for all nodes at tree level l */
      int startind = t->nodes[n].startind;
      int endind   = t->nodes[n].endind;

      if (decompose_node(t, n, l, use_sigma))
      {   /* well-behaved node */
#if RTA_DEBUG_KDTREEBUILD
        rta_post("Node #%i (%i..%i): mean = ", n, startind, endind);
        rta_row_post(t->mean, n, t->ndim, "\n");
#endif
        i = startind;
        j = endind;

        while (i < j)
        { /* sort node vectors by distance to splitplane */
          while (i < j  &&  distV2N(t, kdtree_get_vector(t, i), n) <= 0)
            i++;  // if (i >= t->ndata) rta_post("n %d: i=%d\n", n, i);

          while (i < j  &&  distV2N(t, kdtree_get_vector(t, j), n) > 0)
            j--;  // if (j < 0) rta_post("n %d: j=%d\n", n, j);

          if (i < j)
            swap(t, i, j);    // rta_post("swap %i and %i\n", i ,j);
        }
      }
      else
      {
        if (startind == endind)
        {   /* singleton node: don't split, pass through to left lower (leaf) level node */
          j = startind + 1;
          i = endind + 1; /* create empty right node */
        }
        else
        { /* degenerate node: all points on splitplane -> halve */
          int middle = (startind + endind) >> 1;
          j = middle;
          i = middle + 1;
#if RTA_DEBUG_KDTREEBUILD
          rta_post("degenerate Node #%i (%i..%i): splitting at %d, %d  mean = ",
                   n, startind, endind, j, i);
          rta_row_post(t->mean, n, t->ndim, "\n");
#endif
        }
      }
#if RTA_DEBUG_KDTREEBUILD > 1
      rta_post("  --> decomposition (%i..%i), (%i..%i)\n", startind, j - 1, i, endind);
#endif

      assert(2*n+2 < t->nnodes);
      t->nodes[2*n+1].startind = startind; // start index of left child of node n
      t->nodes[2*n+1].endind   = j - 1;  // end   index of left child of node n
      t->nodes[2*n+1].size     = j - startind;

      t->nodes[2*n+2].startind = i;  // start index of right child of node n
      t->nodes[2*n+2].endind   = endind;   // end   index of right child of node n
      t->nodes[2*n+2].size     = endind - i + 1;
    }   /* end for nodes n */
  }
}

