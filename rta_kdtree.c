/*
 * FTM 
 * Copyright (C) 1994, 1995, 1998, 1999, 2007 by IRCAM-Centre Georges Pompidou, Paris, France.
 * 
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 * 
 * See file COPYING.LIB for further informations on licensing terms.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 * 
 */

#ifdef WIN32
#include <malloc.h>
#else
#include <alloca.h>
#endif
#include <math.h>


#include "rta_kdtree.h"


const char *kdtree_dmodestr[] = { "orthogonal", "hyperplane", "pca" };
const char *kdtree_mmodestr[] = { "mean", "middle", "median" };


#if PROFILE
void profile_clear (kdtree_t *t)
{
    t->profile.v2v   	  = 0; 
    t->profile.v2n   	  = 0;
    t->profile.mean  	  = 0;
    t->profile.hyper 	  = 0;
    t->profile.searches   = 0;
    t->profile.neighbours = 0;
    t->profile.maxstack   = 0;
}
#endif


void vec_post (rta_real_t *v, int stride, int n, const char *suffix)
{
    int i, ii;

    for (i = 0, ii = 0; i < n; i++, ii += stride) 
    { 	
	rta_post("%s%.1f", (ii == 0  ?  "["  :  ", "), v[ii]);
    }
    rta_post("]%s", suffix);
}

void row_post (rta_real_t *v, int row, int n, const char *suffix)
{
    vec_post(v + row * n, 1, n, suffix);
}

void kdtree_info_display (kdtree_t* t) 
{
#   define MB(b)  ((float) (b) / (float) (1024 * 1024))
#   define FLT(b) ((b) * sizeof(rta_real_t))

    float mbdata  = MB(FLT(t->ndata));
    float mbindex = MB(t->ndata       * sizeof(int));
    float mbstack = MB(t->stack.alloc * sizeof(kdtree_stack_elem_t));
    float mbnodes = MB(t->nnodes      * sizeof(kdtree_node_t));
    /* inner nodes' mean vectors and splitplanes 
       (these only in hyperplane mode) */
    float mbinner = MB(t->ninner * FLT(t->ndim) * 
		       (t->dmode == dmode_orthogonal ? 1 : 2));

    rta_post("\nTree Info:\n");
    rta_post("ndim        = %d\n", t->ndim);
    rta_post("ndata       = %d  (%.3f MB extern alloc size)\n", t->ndata, mbdata);
    rta_post("nalloc      = %d  (%.3f MB index)\n",  t->ndata, mbindex);
    rta_post("maxheight   = %d\n", t->maxheight);
    rta_post("givenheight = %d\n", t->givenheight);
    rta_post("height      = %d\n", t->height);
    rta_post("nnodes      = %d  (%.3f MB node struct)\n", t->nnodes, mbnodes);
    rta_post("inner nodes = %d  (%.3f MB node vectors)\n", t->ninner, mbinner);
    rta_post("stack       = %d  (%.3f MB)\n", t->stack.alloc, mbstack);
    rta_post("total size  = %.3f MB\n", 
	     MB(sizeof(kdtree_t)) + mbnodes + mbinner + mbindex + mbstack);
    rta_post("sort mode     = %d\n", t->sort);
    rta_post("decomposition = %s\n", kdtree_dmodestr[t->dmode]);
    rta_post("mean vector   = %s\n", kdtree_mmodestr[t->mmode]);
}

void kdtree_raw_display (kdtree_t* t) 
{
    int i;

    if (t->height == 0 || t->ndata == 0) rta_post("Empty Tree\n");

    for (i = 0; i < t->ndata; i++) 
    {
	rta_post("raw data vec %-3i = ", i);
	vec_post(kdtree_get_row_ptr(t, i), 1, t->ndim, "\n");
    }
}

void kdtree_data_display(kdtree_t* t, int print_data) 
{
    rta_real_t plane[t->ndim];
    int l, n, i;

    rta_post("\nTree Data:\n");
    if (t->height == 0 || t->ndata == 0) rta_post("Empty Tree\n");
	
    for (l = 0; l < t->height; l++) 
    {
	rta_post("Level #%d  nodes %d..%d  splitdim %d\n", 
		 l, pow2(l) - 1, pow2(l+1) - 2, t->nodes[pow2(l) - 1].splitdim);
	for (n = pow2(l) - 1; n < pow2(l+1) - 1; n++) 
	{
	  kdtree_node_t *node = &t->nodes[n];

	  if (n < t->ninner)
	  {
	    rta_post("  inner node %d size %d <%d..%d> splitdim %d splitplane ",
		     n, node->size, node->startind, node->endind, node->splitdim);
		
	    if (t->dmode == dmode_orthogonal)
	    {	/* splitplane implicit orthogonal to splitdim */
		bzero(plane, t->ndim * sizeof(rta_real_t));
		plane[node->splitdim] = 1;
		vec_post(plane, 1, t->ndim, "");
	    }
	    else
		row_post(t->split, n, t->ndim, "");
	  }
	  else
	    rta_post("  leaf node %d size %d <%d..%d> ",
		     n, node->size, node->startind, node->endind);

	  if (print_data)
	  {
	    rta_post(" = (");
	    for (i = node->startind; i <= node->endind; i++) 
	    {
		rta_post("%svec %d = ", (print_data >= 2  ?  "\n    " : ""),
			 t->dataindex[i]);
		vec_post(kdtree_get_vector(t, i), 1, t->ndim,
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

#define auto_alloc(field, in, size) do {				\
    if (in == NULL) /* auto alloc */					\
      if (field) field = rta_realloc(field, size * sizeof(*field));	\
      else       field = rta_malloc(size * sizeof(*field));		\
    else field = in; /* external alloc */ } while (0)


int kdtree_set_data (kdtree_t *self, rta_real_t *data, int *index, int m, int n)
{
  int maxheight   = floor(log2(m));
  int givenheight = self->givenheight;
  int height      = givenheight > 0  ?  givenheight :  maxheight + givenheight;
  int i;

  /* clip height */
  if (height > maxheight)
      height = maxheight;
  if (height < 1)	
      height = 1;	/* minimum: just one node, linear search */

  self->maxheight = maxheight;
  self->height    = height;

  self->data      = data;
  self->ndata     = m;
  self->ndim      = n;

  /* init num nodes */
  self->nnodes = pow2(height)     - 1;
  self->ninner = pow2(height - 1) - 1;

  /* init index list */
  auto_alloc(self->dataindex, index, m);
  for (i = 0; i < m; i++)
      self->dataindex[i] = i;

  /* init search stack size according to tree height 
     (with heuristic margin of 4 times) */
  kdtree_stack_grow(&self->stack, self->height * 4);
  
  return self->nnodes;
}


void kdtree_init_nodes (kdtree_t* self, kdtree_node_t *nodes, rta_real_t *planes, rta_real_t *means) 
{
  auto_alloc(self->nodes, nodes, self->nnodes);
  bzero(self->nodes, self->nnodes * sizeof(kdtree_node_t));

  auto_alloc(self->mean,  means,  self->nnodes * self->ndim);

  if (self->dmode != dmode_orthogonal)
      auto_alloc(self->split, planes, self->nnodes * self->ndim);

  if (self->nnodes > 0)
  {   /* init root node */
      self->nodes[0].startind = 0;		   
      self->nodes[0].endind   = self->ndata - 1;
      self->nodes[0].size     = self->ndata;
  }
}


/* update non-zero sigma index list */
int kdtree_update_sigmanz (kdtree_t *self)
{
    rta_real_t *sigmaptr = self->sigma;
    int j, nnz = 0;

    for (j = 0; j < self->ndim; j++)
	if (sigmaptr[j] != 0)
	    self->sigma_indnz[nnz++] = j;
    self->sigma_nnz = nnz;
    return nnz;
}

void kdtree_set_sigma (kdtree_t *self, rta_real_t *sigma) /* todo: sigma_indnz from outside */
{
    self->sigma = sigma;
    self->sigma_indnz = rta_realloc(self->sigma_indnz, 
				self->ndim * sizeof(*self->sigma));
    kdtree_update_sigmanz(self);
}




/**********************************************************
*
*  class
*
*/

void kdtree_init (kdtree_t *self)
{
  /* Init */
  self->dmode       = dmode_orthogonal;
  self->mmode       = mmode_mean;
  self->sort        = 1;
  self->ndata  	    = 0;
  self->height 	    = 0;
  self->maxheight   = 0;
  self->givenheight = -1;	/* -1 gives less comparisons than -2 */
  self->nnodes 	    = 0;
  self->ndim   	    = 0;
  self->dataindex   = NULL;
  self->nodes  	    = NULL;
  self->data   	    = NULL;
  self->sigma  	    = NULL;
  self->sigma_nnz   = 0;
  self->sigma_indnz = NULL;

  kdtree_stack_init(&self->stack, 0);

#if PROFILE_BUILD
  profile_clear(self);
#endif
}


void kdtree_free (kdtree_t *self)
{
  if (self->dataindex)	 rta_free(self->dataindex);
  if (self->nodes)	 rta_free(self->nodes);
  if (self->mean)	 rta_free(self->mean);
  if (self->split)	 rta_free(self->split);
  if (self->sigma_indnz) rta_free(self->sigma_indnz);

  kdtree_stack_free(&self->stack);

#if DEBUG
  self->data      = NULL;
  self->dataindex = NULL;
  self->nodes     = NULL;
#endif
}
