
#include "rta_kdtree.h"
#include <stdlib.h>
#include <math.h>


#ifdef DEBUG
#define DEBUG_KDTREESEARCH 0
#else
#define DEBUG_KDTREESEARCH 0
#endif


/*
 *  optimised search stack
 */

void kdtree_stack_init (kdtree_stack_t *s, int size)
{
  s->alloc  = size;
  s->buffer = (kdtree_stack_elem_t *) fts_malloc(s->alloc * sizeof(kdtree_stack_elem_t));
  s->size = 0;
}

#define stack_clear(s) ((s)->size = 0)

void kdtree_stack_free (kdtree_stack_t *s)
{
    if (s->buffer)
	fts_free(s->buffer);
}

static void stack_realloc (kdtree_stack_t *stack, int alloc)
{
#if DEBUG_KDTREESEARCH      
    fts_post("kdtree: grow stack from d to %d\n", stack->alloc, alloc);
#endif
    stack->buffer = (kdtree_stack_elem_t *) fts_realloc( stack->buffer, alloc * sizeof(kdtree_stack_elem_t));
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
	fts_post("stack empty\n");
    else
	for (i = s->size - 1; i >= 0; i--)
	{
	    fts_post("    stack pos %d:  node %3d, square dist %f\n", 
		     i, s->buffer[i].node, s->buffer[i].dist);
	}
}



/*
 * support routines
 */

#define MAX(a, b) ((a) > (b) ? (a) : (b))

static int maxArr (float* array, int size) 
{
    int   index = 0;
    float max   = array[0];
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

static float euclidean_distance (float* v1, int stride1, float* v2, int dim) 
{
    int i, i1;
    float sum = 0;

    for (i = 0, i1 = 0; i < dim; i++, i1 += stride1) 
    {
	float diff = v2[i] - v1[i1];
	sum += diff * diff;
    }

    return sum;
}

static float weighted_euclidean_distance (float* v1, int stride1, float* v2, float *sigma, int ndim) 
{
    int i, i1;
    float sum = 0;

    for (i = 0, i1 = 0; i < ndim; i++, i1 += stride1) 
	if (sigma[i] > 0)
	{
	    float diff = (v2[i] - v1[i1]) / sigma[i];
	    sum += diff * diff;
	}

    return sum;
}


/* out:    y[K] = index of the Kth nearest neighbour (in float for interfacing reasons)
	   d[K] = distance of the Kth nearest neighbour
   return: actual number of found neighbours */
int kdtree_search_knn (kdtree_t *t, float* vector, int stride, int k, const float r, 
		       int use_sigma, /* out */ float *indx, float *dist) 
{
    int    kmax         = 0; 		/* index of current kth neighbour */
    int    leaves_start = t->ninner;	/* first leaf node */
    float  sentinel     = (r == 0 ? MAX_FLOAT : r);
    float *sigmaptr     = t->sigma;
    float  dxx;	  	     		/* distance between 2 vectors */
    int    i;	  	     		/* current processed vector */

    kdtree_stack_t     *s = &t->stack;
    kdtree_stack_elem_t cur;		/* current (node, dist) couple */

    if (t->ndata == 0) 
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
#if PROFILE_SEARCH
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
		fts_post("Leaf node p = %d  cur.dist %f\n", cur.node, cur.dist);
#endif
		for (i = istart; i <= iend; i++)
		{
		    if (use_sigma)
			dxx = weighted_euclidean_distance(vector, stride, 
				kdtree_get_vector(t, i), sigmaptr, t->ndim);
		    else
			dxx = euclidean_distance(vector, stride, 
				kdtree_get_vector(t, i), t->ndim);
#if PROFILE_SEARCH
		    t->profile.v2v++;
#endif
#if DEBUG_KDTREESEARCH
		    fts_post("  distance = %f between vector %d ", dxx, i);
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
		float d;

		if (use_sigma)
		    d = distV2N_weighted(t, vector, stride, sigmaptr, cur.node);
		else
		    d = distV2N_stride(t, vector, stride, cur.node);

#if DEBUG_KDTREESEARCH
		fts_post("Inner node %d  d %f  cur.dist %f --> push max %f\n", 
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
	    fts_post("eliminate node %d (size %d): cur.dist %f > dist[kmax=%d] = %f\n",
		     cur.node, t->nodes[cur.node].size, cur.dist, kmax, dist[kmax]);
	}
#endif

    }
#if PROFILE_SEARCH
    t->profile.searches++;
    t->profile.neighbours += kmax + 1;
#endif
#if DEBUG_KDTREESEARCH
    fts_post("kdtree_search found %d vectors < radius %f\n", 
	     kmax + (dist[kmax] < sentinel), r);
#endif

    /* return actual number of found neighbours, can be less than k,
       then kmax is the index of the next one to find */
    return kmax + (dist[kmax] < sentinel);	
}
