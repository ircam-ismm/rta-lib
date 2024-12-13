/**
 * @file rta_kdtree.c
 * @author Diemo Schwarz

 * @brief  k-dimensional search tree
 *
 * Efficient multidimensional binary search tree with logarithmic time
 * complexity: From about 80-100 points, the number of comparisions
 * doesn't rise any more perceptibly.
 *
 *
 * Dimensions can be weighted while building the tree, and while
 * searching. The weight is 1 / sigma, just as in
 * mnm.mahalanobis. Sigma = 0 means: ignore this dimension.
 *
 * Call sequence for builing and using the tree:
 *
 * - 1. initialise tree structure with kdtree_init()
 *
 * - 2. set parameters like decomposition mode kdtree_t#dmode and mean
 * mode kdtree_t#mmode, tree height adaptation kdtree_t#givenheight
 *
 * - 3. set data vector with kdtree_set_data(), this returns the number
 *   of nodes the tree will build
 *
 * - 4. initialise nodes, possibly providing memory allocated according
 *   the number of nodes returned above, with kdtree_init_nodes()
 *
 * - 5. optionally set weights for building with kdtree_set_sigma()
 *
 * - 6. build the tree with kdtree_build()
 *
 * - 7. then you can query the tree with kdtree_search_knn().
 *
 * @copyright
 * Copyright (C) 1994, 1995, 1998, 1999, 2007, 2008, 2009 by IRCAM-Centre Georges Pompidou, Paris, France.
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#ifdef WIN32
#include <malloc.h>
#define snprintf sprintf_s
#else
#include <alloca.h>
#include <strings.h>
#endif

#include "rta_kdtree.h"
#include "rta_kdtreeintern.h"

const char *rta_kdtree_dmodestr[] = { "orthogonal", "hyperplane", "pca" };
const char *rta_kdtree_mmodestr[] = { "mean", "middle", "median" };


#if RTA_KDTREE_PROFILE
void rta_kdtree_profile_clear (rta_kdtree_t *t)
{
  t->profile.v2v      = 0;
  t->profile.v2n      = 0;
  t->profile.mean     = 0;
  t->profile.hyperp     = 0;
  t->profile.searches   = 0;
  t->profile.neighbours = 0;
  t->profile.maxstack   = 0;
}
#endif

#ifndef DOXYGEN_SKIP
/* debug only, not thread safe! */
void rta_vec_post (rta_real_t *v, int stride, int n, const char *suffix)
{
#define MAX_VEC_STR 2048
  static char str[MAX_VEC_STR + 1];
  int p = 0;
  int i, ii;

  for (i = 0, ii = 0; i < n; i++, ii += stride)
    p += snprintf(str + p, MAX_VEC_STR, "%s%.3f", (ii == 0  ?  "["  :  ", "), v[ii]);

  p += snprintf(str + p, MAX_VEC_STR, "]%s", suffix);
  str[p] = 0;
  rta_post(str);
}

void rta_row_post (rta_real_t *v, int row, int n, const char *suffix)
{
  rta_vec_post(v + row * n, 1, n, suffix);
}
#endif /* DOXYGEN_SKIP */

void rta_kdtree_info_display (rta_kdtree_t* t)
{
#define MB(b)  ((float) (b) / (float) (1024 * 1024))
#define FLT(b) ((b) * sizeof(rta_real_t))

  float mbdata  = MB(FLT(t->ndatatot * t->ndim));
  float mbindex = MB(t->ndatatot     * sizeof(rta_kdtree_object_t));
  float mbstack = MB(t->stack.alloc  * sizeof(rta_kdtree_stack_elem_t));
  float mbnodes = MB(t->nnodes       * sizeof(rta_kdtree_node_t));
  /* inner nodes' mean vectors and splitplanes
     (these only in hyperplane mode) */
  float mbinner = MB(t->ninner * FLT(t->ndim) *
                  (t->dmode == dmode_orthogonal ? 1 : 2));

  rta_post("\nTree Info:\n");
  rta_post("ndim        = %d\n", t->ndim);
  rta_post("ndata       = %d  (%.3f MB extern alloc size)\n", t->ndatatot, mbdata);
  rta_post("nalloc      = %d  (%.3f MB index)\n",  t->ndatatot, mbindex);
  rta_post("maxheight   = %d\n", t->maxheight);
  rta_post("givenheight = %d\n", t->givenheight);
  rta_post("height      = %d\n", t->height);
  rta_post("nnodes      = %d  (%.3f MB node struct)\n", t->nnodes, mbnodes);
  rta_post("inner nodes = %d  (%.3f MB node vectors)\n", t->ninner, mbinner);
  rta_post("stack       = %d  (%.3f MB)\n", t->stack.alloc, mbstack);
  rta_post("total size  = %.3f MB\n",
     MB(sizeof(rta_kdtree_t)) + mbnodes + mbinner + mbindex + mbstack);
  rta_post("sort mode     = %d\n", t->sort);
  rta_post("decomposition = %s\n", rta_kdtree_dmodestr[t->dmode]);
  rta_post("mean vector   = %s\n", rta_kdtree_mmodestr[t->mmode]);
}

void rta_kdtree_raw_display (rta_kdtree_t* t)
{
  int i, k;

  if (t->height == 0 || t->ndatatot == 0)
    rta_post("Empty Tree\n");

  for (k = 0; k < t->nblocks; k++)
    for (i = 0; i < t->ndata[k]; i++)
    {
      rta_post("block %d raw data vec %-3i = ", k, i);
      rta_vec_post(rta_kdtree_get_row_ptr(t, k, i), 1, t->ndim, "\n");
    }
}

void rta_kdtree_data_display (rta_kdtree_t* t, int print_data)
{
#ifdef WIN32
  rta_real_t plane[2048];
#else
  rta_real_t plane[t->ndim];
#endif
  int l, n, i;

  rta_post("\nTree Data:\n");
  if (t->height == 0 || t->ndatatot == 0)
    rta_post("Empty Tree\n");

  for (l = 0; l < t->height; l++)
  {
    rta_post("Level #%d  nodes %d..%d  splitdim %d\n",
      l, pow2(l) - 1, pow2(l+1) - 2, t->nodes[pow2(l) - 1].splitdim);

    for (n = pow2(l) - 1; n < pow2(l+1) - 1; n++)
    {
      rta_kdtree_node_t *node = &t->nodes[n];

      if (n < t->ninner)
      {
        rta_post("  inner node %d size %d <%d..%d> splitdim %d splitplane ",
           n, node->size, node->startind, node->endind, node->splitdim);

        if (t->dmode == dmode_orthogonal)
        { /* splitplane implicit orthogonal to splitdim */
#ifndef WIN32
          bzero(plane, t->ndim * sizeof(rta_real_t));
#else
          memset(plane, 0.0, t->ndim * sizeof(rta_real_t));
#endif
          plane[node->splitdim] = 1;
          rta_vec_post(plane, 1, t->ndim, "");
        }
        else
          rta_row_post(t->split, n, t->ndim, "");

        if (print_data > 1)
        {
          rta_post(" mean ");
          rta_vec_post(t->mean + n * t->ndim, 1, t->ndim, "");
        }
      }
      else
        rta_post("  leaf node %d size %d <%d..%d> ",
           n, node->size, node->startind, node->endind);

      if (print_data > 2)
      {
        rta_post(" = (");
        for (i = node->startind; i <= node->endind; i++)
        {
          rta_post("%svec (%d, %d) = ", (print_data >= 2  ?  "\n    " : ""),
            t->dataindex[i].base, t->dataindex[i].index);
          rta_vec_post(rta_kdtree_get_vector(t, i), 1, t->ndim,
            i < node->endind ? ", " : "");
        }
        rta_post(")");
      }

      rta_post("\n");
    } /* end for nodes n */
  } /* end for level l */
}




/*
 * initialisation
 */

#define rta_auto_alloc(field, in, size) do { \
  if (in == NULL) /* auto alloc */ \
    if (field) field = rta_realloc(field, size * sizeof(*field)); \
    else       field = rta_malloc(size * sizeof(*field)); \
  else field = in; /* external alloc */ } while (0)


int rta_kdtree_set_data (rta_kdtree_t *self, int nblocks, rta_real_t **data,
                         rta_kdtree_object_t *index, int *m, int n)
{
  int maxheight, givenheight, height;
  int i, j = 0, k;

  self->data      = data;
  self->nblocks   = nblocks;
  self->ndata     = m;
  self->ndim      = n;
  self->ndatatot  = 0;

  for (i = 0; i < nblocks; i++)
    self->ndatatot += self->ndata[i];

  maxheight = floor(log2(self->ndatatot));
  givenheight = self->givenheight;
  height = givenheight > 0 ? givenheight : maxheight + givenheight;

  /* clip height */
  if (height > maxheight)
    height = maxheight;

  if (height < 1)
    height = 1; /* minimum: just one node, linear search */

  self->maxheight = maxheight;
  self->height    = height;

  /* init num nodes */
  self->nnodes = pow2(height)     - 1;
  self->ninner = pow2(height - 1) - 1;

  /* init original index list */
  rta_auto_alloc(self->dataindex, index, self->ndatatot);

  if (index == NULL) /* no indices given, create them ourselves; else: use indices from outside */
    for (k = 0; k < nblocks; k++)
      for (i = 0; i < m[k]; i++, j++)
      {
        self->dataindex[j].base = k;
        self->dataindex[j].index = i;
      }

  /* init search stack size according to tree height
     (with heuristic margin of 4 times) */
  rta_kdtree_stack_grow(&self->stack, self->height * 4);

  return self->nnodes;
}


void rta_kdtree_init_nodes (rta_kdtree_t* self, rta_kdtree_node_t *nodes,
                            rta_real_t *planes, rta_real_t *means)
{
  rta_auto_alloc(self->nodes, nodes, self->nnodes);
#ifndef WIN32
  bzero(self->nodes, self->nnodes * sizeof(rta_kdtree_node_t));
#else
  memset(self->nodes, 0.0, self->nnodes * sizeof(rta_kdtree_node_t));
#endif
  rta_auto_alloc(self->mean, means, self->nnodes * self->ndim);

  if (self->dmode != dmode_orthogonal)
    rta_auto_alloc(self->split, planes, self->nnodes * self->ndim);

  if (self->nnodes > 0)
  {   /* init root node */
    self->nodes[0].startind = 0;
    self->nodes[0].endind   = self->ndatatot - 1;
    self->nodes[0].size     = self->ndatatot;
  }
}


/* update non-zero sigma index list */
int rta_kdtree_update_sigmanz (rta_kdtree_t *self)
{
  rta_real_t *sigmaptr = self->sigma;
  int j, nnz = 0;

  for (j = 0; j < self->ndim; j++)
    if (sigmaptr[j] != 0)
      self->sigma_indnz[nnz++] = j;

  self->sigma_nnz = nnz;

  return nnz;
}

void rta_kdtree_set_sigma (rta_kdtree_t *self, rta_real_t *sigma) /* todo: sigma_indnz from outside */
{
  self->sigma = sigma;
  self->sigma_indnz = rta_realloc(self->sigma_indnz, self->ndim * sizeof(*self->sigma));
  rta_kdtree_update_sigmanz(self);
}




/**********************************************************
*
*  class
*
*/

void rta_kdtree_init (rta_kdtree_t *self)
{
  /* Init */
  self->dmode       = dmode_orthogonal;
  self->mmode       = mmode_mean;
  self->sort        = 1;
  self->ndata       = NULL;
  self->ndatatot    = 0;
  self->nblocks     = 0;
  self->height      = 0;
  self->maxheight   = 0;
  self->givenheight = -1; /* -1 gives less comparisons than -2 */
  self->nnodes      = 0;
  self->ndim        = 0;
  self->dataindex   = NULL;
  self->nodes       = NULL;
  self->data        = NULL;
  self->mean        = NULL;
  self->sigma       = NULL;
  self->sigma_nnz   = 0;
  self->sigma_indnz = NULL;
  self->activecol   = -1;
  memset(self->dfun, 0, sizeof(void *) * RTA_KDTREE_MAX_DISTWARP);

  rta_kdtree_stack_init(&self->stack, 0);

#if RTA_KDTREE_PROFILE_BUILD
  rta_kdtree_profile_clear(self);
#endif
}


void rta_kdtree_free (rta_kdtree_t *self)
{
  if (self->dataindex) rta_free(self->dataindex);
  if (self->nodes) rta_free(self->nodes);
  if (self->mean) rta_free(self->mean);
  if (self->split) rta_free(self->split);
  if (self->sigma_indnz) rta_free(self->sigma_indnz);

  rta_kdtree_stack_free(&self->stack);

#if DEBUG
  self->data      = NULL;
  self->dataindex = NULL;
  self->nodes     = NULL;
#endif
}
