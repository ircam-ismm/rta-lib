
#include <stdlib.h>
#include <math.h>
#include "rta_kdtree.h"
#include "rta_kdtreeintern.h"


#ifdef DEBUG
#define DEBUG_KDTREESEARCH 0
#else
#define DEBUG_KDTREESEARCH 0
#endif

#define MAX_FLOAT 0x7FFFFFFF


/*
 *  optimised search stack
 */

void kdtree_stack_init (kdtree_stack_t *s, int size)
{
  s->alloc  = size;
  s->buffer = (kdtree_stack_elem_t *) rta_malloc(s->alloc * sizeof(kdtree_stack_elem_t));
  s->size = 0;
}

#define stack_clear(s) ((s)->size = 0)

void kdtree_stack_free (kdtree_stack_t *s)
{
    if (s->buffer)
	rta_free(s->buffer);
}

static void stack_realloc (kdtree_stack_t *stack, int alloc)
{
#if DEBUG_KDTREESEARCH      
    rta_post("kdtree: grow stack from d to %d\n", stack->alloc, alloc);
#endif
    stack->buffer = (kdtree_stack_elem_t *) rta_realloc( stack->buffer, alloc * sizeof(kdtree_stack_elem_t));
    stack->alloc = alloc;
}

void kdtree_stack_grow (kdtree_stack_t *stack, int alloc)
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

static void kdtree_stack_display (kdtree_stack_t *s)
{
    int i;

    if (s->size == 0)
	rta_post("stack empty\n");
    else
	for (i = s->size - 1; i >= 0; i--)
	{
	    rta_post("    stack pos %d:  node %3d, square dist %f\n", 
		     i, s->buffer[i].node, s->buffer[i].dist);
	}
}



/*
 * support routines
 */

#define MAX(a, b) ((a) > (b) ? (a) : (b))

static int maxArr (rta_real_t* array, int size) 
{
    int        index = 0;
    rta_real_t max   = array[0];
    int i;
	
    for (i = 1; i < size; i++) 
    {
	if (array[i] > max)
	{
	    index = i;
	    max   = array[i];
	}
    }
	
    return index;
}

rta_real_t rta_euclidean_distance (rta_real_t* v1, int stride1, 
				   rta_real_t* v2, int dim,
				   rta_bpf_t  *distfunc[]) 
{
    int i, i1;
    rta_real_t sum = 0;

    for (i = 0, i1 = 0; i < dim; i++, i1 += stride1) 
    {
#if USE_DISTFUNC // uses rta_bpf_t, (data-compatible to FTM bpfunc_t)
      rta_real_t diff = v2[i] - v1[i1];
      rta_bpf_t *dfun = distfunc[i];
      if (dfun)
	diff = rta_bpf_get_interpolated(dfun, diff);
#else
      rta_real_t diff = v2[i] - v1[i1];
#endif /* USE_DISTFUNC */
      sum += diff * diff;
    }

    return sum;
}


rta_real_t rta_weighted_euclidean_distance (rta_real_t* v1, 
					    rta_real_t* v2, 
					    rta_real_t *sigma, int ndim,
					    rta_bpf_t  *distfunc[]) 
{
    int i;
    rta_real_t sum = 0, sqrtsum = 0;

  for (i = 0; i < ndim; i++) 
    if (sigma[i] > 0)
    {
#if USE_DISTFUNC // uses rta_bpf_t, (data-compatible to FTM bpfunc_t)
      rta_real_t diff = v2[i] - v1[i];
      rta_bpf_t *dfun = distfunc[i];
      if (dfun)
	diff = rta_bpf_get_interpolated(dfun, diff);
      diff /= sigma[i];
#else
      rta_real_t diff = (v2[i] - v1[i]) / sigma[i];
#endif /* USE_DISTFUNC */
      sum += diff * diff;
      
      rta_post("rta_weighted_euclidean_distance %d (%f - %f)  ->  %f sum d^2 %f sum %f\n",
               i, v2[i], v1[i], diff, sum, sqrt(sum));
    }
  sqrtsum = sqrt(sum); 
  return sqrtsum;
}

rta_real_t rta_weighted_euclidean_distance_stride (rta_real_t* v1, int stride1,
						   rta_real_t* v2, 
						   rta_real_t *sigma, int ndim,
						   rta_bpf_t  *distfunc[]) 
{
    int i, i1;
    rta_real_t sum = 0;

    for (i = 0, i1 = 0; i < ndim; i++, i1 += stride1) 
	if (sigma[i] > 0)
	{
#if USE_DISTFUNC // uses rta_bpf_t, (data-compatible to FTM bpfunc_t)
	  rta_real_t diff = v2[i] - v1[i1];
	  rta_bpf_t *dfun = distfunc[i];
	  if (dfun)
	    diff = rta_bpf_get_interpolated(dfun, diff);
	  diff /= sigma[i];
#else
	  rta_real_t diff = (v2[i] - v1[i1]) / sigma[i];
#endif /* USE_DISTFUNC */
	  sum += diff * diff;
	}

    return sum;
}


/* out:    indx[K] = index of the Kth nearest neighbour
	   dist[K] = distance of the Kth nearest neighbour
   return: actual number of found neighbours */
int kdtree_search_knn (kdtree_t *t, rta_real_t* vector, int stride, 
		       int k, const rta_real_t r, int use_sigma, 
		       /* out */ kdtree_object_t *indx, rta_real_t *dist) 
{
    int		kmax         = 0; 	  /* index of current kth neighbour */
    int		leaves_start = t->ninner; /* first leaf node */
    rta_real_t  sentinel     = (r == 0 ? MAX_FLOAT : r);
    rta_real_t *sigmaptr     = t->sigma;
    rta_real_t  dxx;	  	     	  /* distance between 2 vectors */
    int		i;  	     		  /* current processed vector */

    kdtree_stack_t     *s    = &t->stack;
    kdtree_stack_elem_t cur;		  /* current (node, dist) couple */

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
#if DEBUG_KDTREESEARCH >= 2
	kdtree_stack_display(s);
#endif
#if KDTREE_PROFILE_SEARCH
	if (s->size > t->profile.maxstack)
	    t->profile.maxstack = s->size;
#endif
	stack_pop(s, &cur); 

	if (cur.dist <= dist[kmax])  // elimination rule
	{  
	    if (cur.node >= leaves_start) 
	    {   /* leaf node: search through vectors linearly */
		int istart = t->nodes[cur.node].startind;
		int iend   = t->nodes[cur.node].endind;
		int i;

#if DEBUG_KDTREESEARCH
		rta_post("Leaf node p = %d  cur.dist %f\n", cur.node, cur.dist);
#endif
		for (i = istart; i <= iend; i++)
		{
		    if (use_sigma)
			dxx = rta_weighted_euclidean_distance_stride(vector, stride, 
				kdtree_get_vector(t, i), sigmaptr, t->ndim, t->dfun);
		    else
			dxx = rta_euclidean_distance(vector, stride, 
				kdtree_get_vector(t, i), t->ndim, t->dfun);
#if KDTREE_PROFILE_SEARCH
		    t->profile.v2v++;
#endif
#if DEBUG_KDTREESEARCH
		    rta_post("  distance = %f between vector %d ", dxx, i);
		    vec_post(kdtree_get_vector(t, i), 1, t->ndim, " and x ");
		    vec_post(vector, stride,             t->ndim, "\n");
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
			    int pos = kmax;	/* where to insert */
			
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
		    }
		}
	    }
	    else 
	    {	// branched node
		rta_real_t d;

		if (use_sigma)
		    d = distV2N_weighted(t, vector, stride, sigmaptr, cur.node);
		else
		    d = distV2N_stride(t, vector, stride, cur.node);

#if DEBUG_KDTREESEARCH
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
	}
#if DEBUG_KDTREESEARCH
	else /* node can be eliminated from search */
	{
	    rta_post("eliminate node %d (size %d): cur.dist %f > dist[kmax=%d] = %f\n",
		     cur.node, t->nodes[cur.node].size, cur.dist, kmax, dist[kmax]);
	}
#endif

    }
#if KDTREE_PROFILE_SEARCH
    t->profile.searches++;
    t->profile.neighbours += kmax + 1;
#endif
#if DEBUG_KDTREESEARCH
    rta_post("kdtree_search found %d vectors < radius %f\n", 
	     kmax + (dist[kmax] < sentinel), r);
#endif

    /* return actual number of found neighbours, can be less than k,
       then kmax is the index of the next one to find */
    return kmax + (dist[kmax] < sentinel);	
}
