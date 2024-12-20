/* -*-mode:c; c-basic-offset: 2-*- */
/** 
 * @file rta_kdtreesearch.c
 * @author Diemo Schwarz
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

#include <stdlib.h>
#include <math.h>
#include "rta_kdtree.h"
#include "rta_kdtreeintern.h"


#ifdef DEBUG
#define RTA_DEBUG_KDTREESEARCH 0 // 1
#else
#define RTA_DEBUG_KDTREESEARCH 0
#endif

/*
 *  optimised search stack
 */

void rta_kdtree_stack_init (rta_kdtree_stack_t *s, int size)
{
  s->alloc  = size;
  s->buffer = (rta_kdtree_stack_elem_t *) rta_malloc(s->alloc * sizeof(rta_kdtree_stack_elem_t));
  s->size = 0;
}

#define stack_clear(s) ((s)->size = 0)

void rta_kdtree_stack_free (rta_kdtree_stack_t *s)
{
  if (s->buffer)
    rta_free(s->buffer);
}

static void stack_realloc (rta_kdtree_stack_t *stack, int alloc)
{
#if RTA_DEBUG_KDTREESEARCH
    rta_post("kdtree: grow stack from d to %d\n", stack->alloc, alloc);
#endif
    stack->buffer = (rta_kdtree_stack_elem_t *) rta_realloc( stack->buffer, alloc * sizeof(rta_kdtree_stack_elem_t));
    stack->alloc = alloc;
}

void rta_kdtree_stack_grow (rta_kdtree_stack_t *stack, int alloc)
{
  if (alloc > stack->alloc)
    stack_realloc (stack, alloc);
}

#define stack_push(s, n, d) ((((s)->size >= (s)->alloc) ? stack_realloc((s), 2 * (s)->alloc) : 0), \
                            ((s)->buffer)[(s)->size].node = (n), \
                            ((s)->buffer)[(s)->size++].dist = (d))

#define stack_pop_safe(s, elem) (((s)->size > 0) ? (elem = s->buffer[--((s)->size)]) : 0)

#define stack_pop(s, elem) (*elem = s->buffer[--((s)->size)])

#define stack_empty(s) ((s)->size == 0)

static void kdtree_stack_display (rta_kdtree_stack_t *s)
{
  int i;

  if (s->size == 0)
    rta_post("stack empty\n");
  else
    for (i = s->size - 1; i >= 0; i--)
      rta_post("    stack pos %d:  node %3d, square dist %f\n",
               i, s->buffer[i].node, s->buffer[i].dist);
}



/*
 * support routines
 */

#define MAX(a, b) ((a) > (b) ? (a) : (b))

static int maxArr (rta_real_t* array, int size)
{
  int index = 0;
  rta_real_t max = array[0];
  int i;

  for (i = 1; i < size; i++)
  {
    if (array[i] > max)
    {
      index = i;
      max = array[i];
    }
  }

  return index;
}

rta_real_t rta_euclidean_distance (rta_real_t* v1, int stride1,
                                   rta_real_t* v2, int dim,
                                   rta_bpf_t  *distwarp[])
{
  int i, i1;
  rta_real_t sum = 0;

  for (i = 0, i1 = 0; i < dim; i++, i1 += stride1)
  {
    rta_real_t diff = v2[i] - v1[i1];
#if RTA_USE_DISTWARP // uses rta_bpf_t, (data-compatible to FTM bpfunc_t)
    rta_bpf_t *dfun = distwarp[i];

    if (dfun)
      diff = rta_bpf_get_interpolated(dfun, diff);
#endif /* RTA_USE_DISTWARP */
    sum += diff * diff;
  }

  return sum; // returns square distance
}

rta_real_t rta_weighted_euclidean_distance_stride (rta_real_t* v1, int stride1,
                                                   rta_real_t* v2,
                                                   rta_real_t *sigma, int ndim,
                                                   rta_bpf_t *distwarp[])
{
  int i, i1;
  rta_real_t sum = 0;

  for (i = 0, i1 = 0; i < ndim; i++, i1 += stride1)
    if (sigma[i] > 0)
    {
#if RTA_USE_DISTWARP // uses rta_bpf_t, (data-compatible to FTM bpfunc_t)
      rta_real_t diff = v2[i] - v1[i1];
      rta_bpf_t *dfun = distwarp[i];

      if (dfun)
        diff = rta_bpf_get_interpolated(dfun, diff);

      diff /= sigma[i];
#else
      rta_real_t diff = (v2[i] - v1[i1]) / sigma[i];
#endif /* RTA_USE_DISTWARP */
      sum += diff * diff;
    }

  return sum; // returns square distance
}


rta_real_t rta_euclidean_distance_Linf (rta_real_t* v1, int stride1,
					rta_real_t* v2, int dim,
					rta_bpf_t  *distwarp[])
{
  int i, i1;
  rta_real_t max = 0;

  for (i = 0, i1 = 0; i < dim; i++, i1 += stride1)
  {
    rta_real_t diff = v2[i] - v1[i1];
#if RTA_USE_DISTWARP // uses rta_bpf_t, (data-compatible to FTM bpfunc_t)
    rta_bpf_t *dfun = distwarp[i];

    if (dfun)
      diff = rta_bpf_get_interpolated(dfun, diff);
#endif /* RTA_USE_DISTWARP */
    diff = fabs(diff);
    if (diff > max)
      max = diff;
  }

  return max; // returns max distance
}

rta_real_t rta_weighted_euclidean_distance_stride_Linf (rta_real_t* v1, int stride1,
							rta_real_t* v2,
							rta_real_t *sigma, int ndim,
							rta_bpf_t *distwarp[])
{
  int i, i1;
  rta_real_t max = 0;

  for (i = 0, i1 = 0; i < ndim; i++, i1 += stride1)
    if (sigma[i] > 0)
    {
#if RTA_USE_DISTWARP // uses rta_bpf_t, (data-compatible to FTM bpfunc_t)
      rta_real_t diff = v2[i] - v1[i1];
      rta_bpf_t *dfun = distwarp[i];

      if (dfun)
        diff = rta_bpf_get_interpolated(dfun, diff);

      diff /= sigma[i];
#else
      rta_real_t diff = (v2[i] - v1[i1]) / sigma[i];
#endif /* RTA_USE_DISTWARP */
    diff = fabs(diff);
    if (diff > max)
      max = diff;
    }

  return max; // returns max distance
}


/** set column index to use for filtering (in)active units, or -1 to disable filtering */
void rta_kdtree_set_activecolumn (rta_kdtree_t *t, int col)
{
  if (col >= t->ndim)
    col = t->ndim - 1;  // clip to num. columns, < 0 means don't use active column
  
  t->activecol = col;
}


#if RTA_KDTREE_USE_ACTIVE         
static bool get_active (rta_kdtree_t *t, int index)
{
  return t->activecol >= 0  ?  (rta_kdtree_get_element(t, index, t->activecol) != 0)  :  true;
}
#else
#  define get_active(t, i)  (true)
#endif


/* Perform search in kd-tree structure t
   params:
     vector of ndim elements to search nearest neighbours of
     stride in input vector
     k max number of neighbours to find (actual number can be lower)
     r max squared distance of neighbours to find (r = 0 means no limit)
     use_sigma flag to use weights
   out:
     indx[K] = (base, element) index of the Kth nearest neighbour
     dist[K] = squared distance of the Kth nearest neighbour
     return: actual number of found neighbours 
*/
int rta_kdtree_search_knn (rta_kdtree_t *t, rta_real_t* vector, int stride,
                           int k, const rta_real_t r, int use_sigma,
                 /* out */ rta_kdtree_object_t *indx, rta_real_t *dist)
{
  int kmax = 0; /* index of current kth neighbour */
  int leaves_start = t->ninner; /* first leaf node */
  rta_real_t sentinel = (r == 0 ? MAX_FLOAT : r);
  rta_real_t *sigmaptr = t->sigma;
  rta_real_t dxx; /* distance between 2 vectors */
  int i; /* current processed vector */

  rta_kdtree_stack_t *s = &t->stack;
  rta_kdtree_stack_elem_t cur; /* current (node, dist) couple */

  if (t->ndatatot == 0)
    return 0;

  if (k < 1)
    k = 1;

  // Init distances
  for (i = 0; i < k; i++)
    dist[i] = sentinel;

  // Init Search Stack
  stack_clear(s);
  stack_push(s, 0, 0);

  while (!stack_empty(s))
  {
#if RTA_DEBUG_KDTREESEARCH >= 2
    kdtree_stack_display(s);
#endif
#if RTA_KDTREE_PROFILE_SEARCH
    if (s->size > t->profile.maxstack)
      t->profile.maxstack = s->size;
#endif
    stack_pop(s, &cur);

    if (cur.dist <= dist[kmax])  // elimination rule
    {
      if (cur.node >= leaves_start)
      {   /* leaf node: search through vectors linearly */
        int istart = t->nodes[cur.node].startind;
        int iend = t->nodes[cur.node].endind;
        int i;

#if RTA_DEBUG_KDTREESEARCH
        rta_post("Leaf node p = %d  cur.dist %f\n", cur.node, cur.dist);
#endif
        for (i = istart; i <= iend; i++)
        {
          const bool isactive = get_active(t, i);

          if (isactive)
          {
            switch (use_sigma)
	    {
	    case 1:
              dxx = rta_weighted_euclidean_distance_stride(
                      vector, stride, rta_kdtree_get_vector(t, i), sigmaptr, t->ndim, t->dfun);
	    break;
	      
            case 0:
              dxx = rta_euclidean_distance(
                      vector, stride, rta_kdtree_get_vector(t, i), t->ndim, t->dfun);
	    break;

	    case 3:
              dxx = rta_weighted_euclidean_distance_stride_Linf(
                      vector, stride, rta_kdtree_get_vector(t, i), sigmaptr, t->ndim, t->dfun);
	    break;
	      
            case 2:
              dxx = rta_euclidean_distance_Linf(
                      vector, stride, rta_kdtree_get_vector(t, i), t->ndim, t->dfun);
	    break;
	    }
	    
#if RTA_KDTREE_PROFILE_SEARCH
	    t->profile.v2v++;
#endif
#if RTA_DEBUG_KDTREESEARCH
	    rta_post("  distance = %f between vector %d (elem %d, %d) ", dxx, i, t->dataindex[i].base, t->dataindex[i].index);
	    rta_vec_post(rta_kdtree_get_vector(t, i), 1, t->ndim, " and x ");
	    rta_vec_post(vector, stride,             t->ndim, "\n");
#endif
	    if (dxx <= dist[kmax])
	    {   /* return original index in data and distance */
	      if (k == 1)
	      {
		indx[kmax] = t->dataindex[i];
		dist[kmax] = dxx;
	      }
	      else if (t->sort)
	      {
		int pos = kmax; /* where to insert */
		
		if (kmax < k - 1)
		{   /* first move or override */
		  dist[kmax + 1] = dist[kmax];
		  indx[kmax + 1] = indx[kmax];
		  kmax++;
		}
		
		/* insert into sorted list of distance */
		while (pos > 0  &&  dxx < dist[pos - 1])
		{   /* move up */
		  dist[pos] = dist[pos - 1];
		  indx[pos] = indx[pos - 1];
		  pos--;
		}
		
		indx[pos] = t->dataindex[i];
		dist[pos] = dxx;
	      }
	      else
	      {
		indx[kmax] = t->dataindex[i];
		dist[kmax] = dxx;
		kmax = maxArr(dist, k);
	      }
	    } // end if dxx <= dist[kmax]
	  } // end if vector isactive
	} // end for all vectors in node
      } // end if leaf node
      else
      { // branched node
        rta_real_t d; // signed distance of target vector to node

        if (use_sigma)
          d = distV2N_weighted(t, vector, stride, sigmaptr, cur.node);
        else
          d = distV2N_stride(t, vector, stride, cur.node);

#if RTA_DEBUG_KDTREESEARCH
        rta_post("Inner node %d  d %f  cur.dist %f --> push max %f\n",
                 cur.node, d, cur.dist, MAX(cur.dist, d*d));
#endif
        if (d < 0)
        {
          stack_push(s, 2*cur.node+2, MAX(cur.dist, d*d));
          stack_push(s, 2*cur.node+1, cur.dist);
        }
        else
        {
          stack_push(s, 2*cur.node+1, MAX(cur.dist, d*d));
          stack_push(s, 2*cur.node+2, cur.dist);
        }
      }
    } // end if (cur.dist <= dist[kmax])  // elimination rule
#if RTA_DEBUG_KDTREESEARCH
    else /* node distance to target > than max distance found: can be eliminated from search */
    {
      rta_post("eliminate node %d (size %d): cur.dist %f > dist[kmax=%d] = %f\n",
               cur.node, t->nodes[cur.node].size, cur.dist, kmax, dist[kmax]);
    }
#endif
  } // end while stack is not empty
#if RTA_KDTREE_PROFILE_SEARCH
  t->profile.searches++;
  t->profile.neighbours += kmax + 1;
#endif
#if RTA_DEBUG_KDTREESEARCH
  rta_post("kdtree_search found %d vectors < radius %f\n",
           kmax + (dist[kmax] < sentinel), r);
#endif

  /* return actual number of found neighbours, can be less than k,
     then kmax is the index of the next one to find */
  return kmax + (dist[kmax] < sentinel);
}
