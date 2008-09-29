/*
 * RTA
 * Copyright (C) 2008 by IRCAM-Centre Georges Pompidou, Paris, France.
 */

#ifndef _DATA_KDTREE_H_
#define _DATA_KDTREE_H_

#include "fts.h"


#define PROFILE_BUILD  1
#define PROFILE_SEARCH 1
#define PROFILE	       (PROFILE_BUILD || PROFILE_SEARCH)

#define MAX_FLOAT 0x7FFFFFFF  

#define MAX_NODES 100

/* stack elements for search algorithm */
typedef struct _elem 
{
    int   node;
    float dist;
} kdtree_stack_elem_t;

typedef struct _stack
{
    int size;
    int alloc;
    kdtree_stack_elem_t *buffer;
} kdtree_stack_t;


typedef struct node
{
    int startind;	/* index of first vector in node in dataindex array */
    int endind;		/* index of last vector in node in dataindex array */
    int size;		/* number of vectors in node */
    int splitdim;	/* for dmode othogonal, dimension along which node is split*/

    fmat_t *mean;	/* mean vector (todo: median), always present */
    fmat_t *split;	/* hyperplane A1*X1 + A2*X2 +...+ An*Xn + An+1 = 0, 
			   or NULL in dmode_orthogonal */
    float   splitnorm;	/* length of split vector */
} node_t;


/* decomposition mode */
typedef enum 
{ 
    dmode_orthogonal, 	/* optimised othogonal to axes */
    dmode_hyperplane, 	/* hyperplane orthogonal to axes */
    dmode_pca	      	/* hyperplane along principal components */
} kdtree_dmode_t;

/* pivot (mean vector to split at) calculation mode */
typedef enum 
{
    mmode_mean,		/* mean of values / distances to splitplane */
    mmode_middle,	/* middle between min/max */
    mmode_median	/* true median (guarantees equal number of points left
			   and right and thus a well-balanced and optimal tree) */
} kdtree_mmode_t;

typedef struct _kdtree
{
    kdtree_dmode_t dmode; /* decomposition mode */
    kdtree_mmode_t mmode; /* pivot (mean vector to split at) calculation mode */

    int     ndim;	  /* Dimension of vectors */
    int     ndata;	  /* Number of vectors */
    fmat_t *data;	  /* data matrix (ndata, ndim) */
    int	    dataalloc;	  /* allocated size of dataindex */
    int    *dataindex;	  /* data vector indirection array (ndata) */

    fmat_t *sigma;	  /* weight, 0 == inf */
    int     sigma_nnz;    /* number of non-zero sigma */
    int    *sigma_indnz;  /* non-zero sigma lines */

    int     height;	  /* Height of the kdtree */
    int     maxheight;	  /* Maximal height of the kdtree */
    int     givenheight;  /* Height given by user, subtract from max if <= 0 */
    int     nnodes; 	  /* Number of nodes (must be a power of 2) */
    int     ninner; 	  /* Number of inner nodes (=index of first leaf node)*/
    node_t *nodes;	  /* nodes (nnodes) */

    int	    sort;	  /* sort search result by distance */
    kdtree_stack_t stack;

    struct {
	int v2v;	  /* vector to vector distances */
	int v2n;	  /* vector to node distances */
	int mean;	  /* mean vector calculations */
	int hyper;	  /* split plane calculations */
	int searches;	  /* searches performed */
	int neighbours;	  /* neighbours found */
	int maxstack;	  /* highest stack size */
    } profile;
} kdtree_t;


const char *kdtree_dmodestr[];
const char *kdtree_mmodestr[];


/* get data element via indirection order array */
#define kdtree_get_element(t, i, j)  fmat_get_element(t->data, t->dataindex[i], j)

/* get data vector via indirection order array */
#define kdtree_get_vector(t, i)     (fmat_get_ptr(t->data) + t->dataindex[i] * t->ndim)

#define fmat_get_row_ptr(f, i)    (fmat_get_ptr(f) + (i) * fmat_get_n(f))

#define pow2(x)  (1 << (x))


void vec_post (float *v, int stride, int n, const char *suffix);
void kdtree_info_display (kdtree_t* t);
void kdtree_raw_display  (kdtree_t* t);
void kdtree_data_display (kdtree_t* t, int print_data);
void profile_clear (kdtree_t *t);

void kdtree_stack_init (kdtree_stack_t *s, int size);
void kdtree_stack_free (kdtree_stack_t *s);
void kdtree_stack_grow (kdtree_stack_t *stack, int alloc);

/* vector to node distance */
float distV2N (kdtree_t* t, const float *x, const int node);
float distV2N_stride (kdtree_t* t, const float *x, int stride, const int node);
float distV2N_weighted (kdtree_t* t, const float *x, int stride, const float *sigma, const int node);

void kdtree_set_decomposition(kdtree_t *t, kdtree_dmode_t mode, void *param);

void kdtree_clear_nodes (kdtree_t *self);
void kdtree_set   (kdtree_t *self, fmat_t *data, fmat_t *sigma, int use_sigma);
void kdtree_set_sigma (kdtree_t *self, fmat_t *sigma);
int  kdtree_update_sigmanz (kdtree_t *self);
void kdtree_build (kdtree_t* t, int use_sigma);
void kdtree_init_data (kdtree_t* t, int h, int vect_num, int dim);

int  kdtree_search_knn (kdtree_t *t, float* x, int stride, int k, float r, float *y, float *d);

#endif
