
#include "kdtree.h"
#include <stdlib.h>
#include <math.h>


#ifdef DEBUG
#define DEBUG_KDTREEBUILD 0
#else
#define DEBUG_KDTREEBUILD 0
#endif


static void compute_mean (kdtree_t *t, int node, int dim) 
{
    node_t *n        = &t->nodes[node];
    float  *mean_ptr = fmat_get_ptr(n->mean);
//  float  *data_ptr = fmat_get_ptr(t->data); 
    int     nstart   = n->startind;
    int     nend     = n->endind;
    int     nvector  = nend - nstart + 1; // number of vectors of processed node
    int     dstart, dend;
    int     i, j;

    if (dim < 0)
    {   /* all dimensions */
	dstart = 0;
	dend   = t->ndim;
    }
    else
    {   /* calculate mean only over given dimension */
	dstart = dim;
	dend   = dim + 1;
#if DEBUG_KDTREEBUILD
	bzero(mean_ptr, t->ndim * sizeof(float));
#endif
    }

    for (j = dstart; j < dend; j++) 
    {
	float sum = 0;
	
	for (i = nstart; i <= nend; i++) 
	{
	    sum += kdtree_get_element(t, i, j);
	}
	mean_ptr[j] = sum / nvector;
    }

#if DEBUG_KDTREEBUILD
    fts_post("mean vect for node %d (size %d) = ", node, nvector);
    vec_post(mean_ptr, 1, t->ndim, "\n");
#endif
}


static void compute_middle (kdtree_t *t, int node, int dim) 
{
    node_t *n        = &t->nodes[node];
    float  *mean_ptr = fmat_get_ptr(n->mean);
    int     nstart   = n->startind;
    int     nend     = n->endind;
    int     dstart, dend;
    int     i, j;

    if (dim < 0)
    {
	dstart = 0;
	dend   = t->ndim;
    }
    else
    {
	dstart = dim;
	dend   = dim + 1;
#if DEBUG_KDTREEBUILD
	bzero(mean_ptr, t->ndim * sizeof(float));
#endif
    }

    // number of vectors from the processed node
    for (j = dstart; j < dend; j++) 
    {
	float min = MAX_FLOAT, max = -MAX_FLOAT;

	for (i = nstart; i <= nend; i++) 
	{
	    float x = kdtree_get_element(t, i, j);

	    if (x < min)
		min = x;
	    if (x > max)
		max = x;
	}
	mean_ptr[j]   = (max + min) / 2.;
    }

#if DEBUG_KDTREEBUILD
    fts_post("middle vect for node %d = ", node);
    vec_post(mean_ptr, 1, t->ndim, "\n");
#endif
}

/* return 1 if node is well-behaved, 0 if node is degenerate */
static int check_node (kdtree_t *t, int node, int dim) 
{
    int     nstart   = t->nodes[node].startind;
    int     nend     = t->nodes[node].endind;
    float   min, max;
    int     i;

    min = max = kdtree_get_element(t, nstart, dim);

    for (i = nstart + 1; i <= nend; i++) 
    {
	float x = kdtree_get_element(t, i, dim);

	if (x < min)
	    min = x;
	if (x > max)
	    max = x;
    }

    return max != min;
}


static int compare(float* currentVect, float* mean, int component) 
{
// fts_post("compare %f and %f\n", mean[component], currentVect[component]);
    return (currentVect[component] - mean[component]); 
}

/* compute and create node-splitting hyperplane */
static void compute_splitplane (kdtree_t* t, int node, int level) 
{
    node_t *n = &t->nodes[node];

 #if PROFILE_BUILD
    t->profile.hyper++;
#endif
    switch (t->dmode)
    {
    case dmode_hyperplane:
	/* compute hyperplane orthogonal to the base vector number b */
	n->split = fmat_create(t->ndim, 1);
	fts_object_refer((fts_object_t *) n->split);
	
	fmat_set_const  (n->split, 0);
	fmat_set_element(n->split, n->splitdim, 0, 1);

#if DEBUG_KDTREEBUILD
    fts_post("Splitplane of node %i: ", node);
    vec_post(fmat_get_ptr(n->split), 1, t->ndim, "\n");
#endif
    /* FALLTHROUGH */
    case dmode_orthogonal:
	n->splitnorm = 1;	/* splitplane implicit orthogonal to splitdim */
    }
}


/* determine hyperplane that splits node n in two child nodes

   return: 1 if node is well-behaved, 0 if node is degenerate,
   i.e. all vectors have same distance (usually 0) to splitplane 
*/
static int decompose_node (kdtree_t *t, int node, int level, int use_sigma) 
{
    int i, nice_node = 0, splitdim;
    float *sigmaptr  = use_sigma  ?  fmat_get_ptr(t->sigma)  :  NULL;

#if PROFILE_BUILD
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
    if (!nice_node) fts_post("warning: can't find non-degenerate dimension to split node %d at level %d, using dimension %d\n", node, level, splitdim);

    /* always create node mean vector */
    t->nodes[node].mean = fmat_create(t->ndim, 1);
    fts_object_refer((fts_object_t *) t->nodes[node].mean);

    switch (t->mmode)
    { /* N.B.: middle and mean are only linearly  affected by sigma */
      case mmode_mean:
	  if (t->dmode == dmode_orthogonal)
	      compute_mean(t, node, splitdim);
	  else
	      compute_mean(t, node, -1);	/* all dimensions */
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
	      compute_middle(t, node, -1);	/* all dimensions */
      break;
    }

    /* compute and create node splitting hyperplane */
    compute_splitplane(t, node, level);

    return nice_node;
}


/* vector to orthogonal plane node distance along split dimension dim */
static float distV2orthoH (const float* vect, float* mean, int dim) 
{
    return vect[dim] - mean[dim];
}
static float distV2orthoH_stride (const float* vect, int stride, float* mean, int dim) 
{
    return vect[dim * stride] - mean[dim];
}
static float distV2orthoH_weighted (const float* vect, int stride, 
				    float* mean, const float *sigma, int dim) 
{
#if DEBUG_KDTREEBUILD > 1
    fts_post("distV2orthoH_weighted on dim %d: (%f - %f) / %f = %f\n",
	dim, vect[dim * stride], mean[dim], sigma[dim],
	sigma[dim] > 0  ?  (vect[dim * stride] - mean[dim]) / sigma[dim]  :  0);
#endif
    return sigma[dim] > 0  ?  (vect[dim * stride] - mean[dim]) / sigma[dim]  
	                   :  0;
}


/* vector to general plane node distance */
static float distV2H (const float* vect, float* plane, float* mean, 
		      int ndim, float norm) 
{
    // standard algebra computing
    int i;
    float dotprod = 0;

    for(i = 0; i < ndim; i++) 
    {
	dotprod += (vect[i] - mean[i]) * plane[i];
    }
    return (dotprod / norm);
}
static float distV2H_stride (const float* vect, int stride, float* plane, 
			     float* mean, int ndim, float norm) 
{
    // standard algebra computing
    int i, iv;
    float dotprod = 0;

    for(i = 0, iv = 0; i < ndim; i++, iv += stride) 
    {
	dotprod += (vect[iv] - mean[i]) * plane[i];
    }
    return (dotprod / norm);
}
static float distV2H_weighted (const float* vect, int stride, float* plane, 
			       float* mean, const float *sigma, int ndim, float norm) 
{
    // standard algebra computing
    int i, iv;
    float dotprod = 0;

    for (i = 0, iv = 0; i < ndim; i++, iv += stride) 
	if (sigma[i] > 0)
	    dotprod += (vect[iv] - mean[i]) / sigma[i] * plane[i];

    return (dotprod / norm);
}


/* vector to node distance */
float distV2N (kdtree_t* t, const float *x, const int node)
{
    node_t *n = &t->nodes[node];

#if PROFILE_BUILD
    t->profile.v2n++;
#endif

    switch (t->dmode)
    {
    case dmode_orthogonal:
	return distV2orthoH(x, fmat_get_ptr(n->mean), n->splitdim);
    case dmode_hyperplane:
	return distV2H(x, fmat_get_ptr(n->split), 
		          fmat_get_ptr(n->mean), t->ndim, n->splitnorm);
    default:
	fts_post("error: unknown mode %d", t->dmode);
	return 0;
    }
}

float distV2N_weighted (kdtree_t* t, const float *x, int stride, 
			const float *sigma, const int node)
{
    node_t *n    = &t->nodes[node];
    float  *mean = fmat_get_ptr(n->mean);

#if PROFILE_BUILD
    t->profile.v2n++;
#endif
    switch (t->dmode)
    {
    case dmode_orthogonal:
	return distV2orthoH_weighted(x, stride, mean, sigma, n->splitdim);
    case dmode_hyperplane:
	return distV2H_weighted(x, stride, fmat_get_ptr(n->split), mean, sigma, t->ndim, n->splitnorm);
    default:
	fts_post("error: unknown mode %d", t->dmode);
	return 0;
    }
}


/* swap positions of vectors i and j: only in indirection array */
static void swap (kdtree_t* t, int i, int j) 
{
    int tmp;

    tmp = t->dataindex[i];
    t->dataindex[i] = t->dataindex[j];
    t->dataindex[j] = tmp;
}


/* (re-)insert num data vector into tree.  
   If index < ndata, move to correct node, 
   otherwise insert and increment ndata */
void kdtree_insert (kdtree_t* t, int index, int num) 
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


void kdtree_build (kdtree_t* t, int use_sigma) 
{			
    int l;		// current level number
    int n;		// current node number
    int	i, j;		// loop counters
	
    /* Maximum length is equal to pow2(height-1) */ 
    if (pow2(t->height - 1) > t->ndata  ||  t->ndim == 0) 
    {
	if (t->ndata == 0) 
	    fts_post("tree is empty!\n");
	else if (t->ndim == 0) 
	    fts_post("tree has 0 dimensions!  Can't build!\n");
	else
	    fts_post("error: can't build this tree, try with a smaller height: %d > %d\n", pow2(t->height-1), t->ndata);

	return;
    }

    for (l = 0; l < t->height - 1; l++) 
    {   /* initialise inner nodes */
	int nstart = pow2(l)   - 1;
	int nend   = pow2(l+1) - 1;
#if DEBUG_KDTREEBUILD
	fts_post("\nLevel #%i  nodes %d..%d\n", l, nstart, nend);
#endif
	for (n = nstart; n < nend; n++) 
	{   /* for all nodes at tree level l */
	    if (decompose_node(t, n, l, use_sigma))
	    {   /* well-behaved node */
#if DEBUG_KDTREEBUILD
		fts_post("Node #%i (%i..%i): mean = ", 
			 n, t->nodes[n].startind, t->nodes[n].endind); 
		vec_post(fmat_get_ptr(t->nodes[n].mean), 1, t->ndim, "\n");
#endif
		i = t->nodes[n].startind; 
		j = t->nodes[n].endind;
		
		while (i < j) 
		{ /* sort node vectors by distance to splitplane */
		    while (distV2N(t, kdtree_get_vector(t, i), n) <= 0) 
		    {
			i++;	if (i >= t->ndata) fts_post("n %d: i=%d\n", n, i);
		    }
		    while (distV2N(t, kdtree_get_vector(t, j), n) > 0) 
		    {
			j--;	if (j < 0) fts_post("n %d: j=%d\n", n, j);
		    }
		    if (i < j) 
		    { 
			swap(t, i, j);    // fts_post("swap %i and %i\n", i ,j);
		    }
		}
	    }
	    else
	    {   /* degenerate node: all points on splitplane -> halve */
		int middle = (t->nodes[n].startind + t->nodes[n].endind) >> 1;
		j = middle;
		i = middle + 1;
#if DEBUG_KDTREEBUILD
		fts_post("degenerate Node #%i (%i..%i): splitting at %d, %d  mean = ", 
			 n, t->nodes[n].startind, t->nodes[n].endind, j, i); 
		vec_post(fmat_get_ptr(t->nodes[n].mean), 1, t->ndim, "\n");
#endif		
	    }

	    assert(2*n+2 < t->nnodes);
	    t->nodes[2*n+1].startind = t->nodes[n].startind; // start index of left child of node n
	    t->nodes[2*n+1].endind   = j;	             // end   index of left child of node n
	    t->nodes[2*n+1].size     = j - t->nodes[n].startind + 1;

	    t->nodes[2*n+2].startind = i;	             // start index of right child of node n
	    t->nodes[2*n+2].endind   = t->nodes[n].endind;   // end   index of right child of node n
	    t->nodes[2*n+2].size     = t->nodes[n].endind - i + 1;
	}
    }
}


void kdtree_init_data (kdtree_t* t, int h, int vect_num, int dim) 
{
    int i;

    t->ndim   = dim;
    t->ndata  = vect_num;
    t->height = h;

    /* init nodes */
    t->nnodes = pow2(h)     - 1;
    t->ninner = pow2(h - 1) - 1;
    t->nodes  = fts_realloc(t->nodes, t->nnodes * sizeof(node_t));
    bzero(t->nodes, t->nnodes * sizeof(node_t));

    if (t->nnodes > 0)
    {   /* init root node */
	t->nodes[0].startind = 0;		   
	t->nodes[0].endind   = t->ndata - 1;
	t->nodes[0].size     = t->ndata;
    }

    /* init index list */
    t->dataalloc = t->data->alloc / dim; /* index alloc size = data alloc rows */
    t->dataindex = fts_realloc(t->dataindex, t->dataalloc * sizeof(int));

    for (i = 0; i < t->ndata; i++)
	t->dataindex[i] = i;

    /* init search stack size according to tree height 
       (with heuristic margin of 4 times) */
    kdtree_stack_grow(&t->stack, t->height * 4);
}
