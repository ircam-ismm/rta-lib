/**
 * @file	rta_kdtree.h
 * @author	Diemo Schwarz
 * @date	29.9.2008
 * @version	1.0 
 *
 * @brief	k-dimensional search tree
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
 * Copyright (C) 2008 by IRCAM-Centre Georges Pompidou, Paris, France.
 */

#ifndef _RTA_KDTREE_H_
#define _RTA_KDTREE_H_

#include "rta.h"


#define PROFILE_BUILD  1
#define PROFILE_SEARCH 1
#define PROFILE	       (PROFILE_BUILD || PROFILE_SEARCH)

#define MAX_FLOAT 0x7FFFFFFF  

#define MAX_NODES 100

/** stack element for search algorithm */
typedef struct _elem_struct
{
    int   node;
    float dist;
} kdtree_stack_elem_t;

/** optimised stack for search algorithm */
typedef struct _stack_struct
{
    int size;
    int alloc;
    kdtree_stack_elem_t *buffer;
} kdtree_stack_t;

/** one node of the kd-tree */
typedef struct _node_struct
{
    int startind;	/**< index of first vector in node in dataindex array */
    int endind;		/**< index of last vector in node in dataindex array */
    int size;		/**< number of vectors in node */
    int splitdim;	/**< for dmode othogonal, dimension along which node is split*/

    float *mean;	/**< mean vector (todo: median), always present */
    float *split;	/**< hyperplane A1*X1 + A2*X2 +...+ An*Xn + An+1 = 0, 
			   or NULL in dmode_orthogonal */
    float  splitnorm;	/**< spatial length of split vector */
} kdtree_node_t;


/** decomposition mode
 *
 * Different ways to split data space at each tree level: 
 * \li orthogonal to dimension (fastest tree-building speed), 
 * \li by arbitrary hyperplane (for testing only), or 
 * \li by PCA (optimal decomposition)
 */
typedef enum 
{ 
    dmode_orthogonal, 	/**< optimised othogonal to axes */
    dmode_hyperplane, 	/**< hyperplane orthogonal to axes */
    dmode_pca	      	/**< hyperplane along principal components */
} kdtree_dmode_t;

/** pivot calculation mode 
 *
 * the \em pivot is the mean vector to split each tree node at 
 */
typedef enum 
{
    mmode_mean,		/**< mean of values / distances to splitplane */
    mmode_middle,	/**< middle between min/max */
    mmode_median	/**< true median (guarantees equal number of points left
			   and right and thus a well-balanced and optimal tree) */
} kdtree_mmode_t;


/** k-dimensional search tree data structure */
typedef struct _kdtree_struct
{
    kdtree_dmode_t dmode; /**< decomposition mode */
    kdtree_mmode_t mmode; /**< pivot calculation mode */

    int     ndim;	  /**< Dimension of vectors */
    int     ndata;	  /**< Number of vectors */
    float  *data;	  /**< data matrix (ndata, ndim) */
    int	    dataalloc;	  /**< allocated size of dataindex */
    int    *dataindex;	  /**< data vector indirection array (ndata) */

    float  *sigma;	  /**< 1/weight, 0 == inf */
    int     sigma_nnz;    /**< number of non-zero sigma */
    int    *sigma_indnz;  /**< non-zero sigma lines */

    int     height;	  /**< Height of the kdtree */
    int     maxheight;	  /**< Maximal height of the kdtree */
    int     givenheight;  /**< Height given by user, gives tree height
			     if positive, subtract from maxheight if negative */
    int     nnodes; 	  /**< Number of nodes (must be a power of 2) */
    int     ninner; 	  /**< Number of inner nodes (=index of first leaf node)*/
    kdtree_node_t *nodes; /**< nodes (nnodes) */

    int	    sort;	  /**< sort search result by distance */
    kdtree_stack_t stack;

    /** profiling data: count internal operations */
    struct {
	int v2v;	  /**< vector to vector distances */
	int v2n;	  /**< vector to node distances */
	int mean;	  /**< mean vector calculations */
	int hyper;	  /**< split plane calculations */
	int searches;	  /**< searches performed */
	int neighbours;	  /**< neighbours found */
	int maxstack;	  /**< highest stack size */
    } profile;
} kdtree_t;


const char *kdtree_dmodestr[];
const char *kdtree_mmodestr[];


/** get data element via indirection order array 
 *
 * This macro returns the element of the kdtree_t#data space by ordered row
 * index \p i, that is mapped to the original row ordering by the
 * indirection array kdtree_t#dataindex.
 *
 * @param t	kd-tree structure
 * @param i	ordered row index (0 <= i < kdtree_t#ndata)
 * @param j	column index      (0 <= j < kdtree_t#ndim)
 * @return	data element at row \p i, column \p j
 */
#if DOXYGEN_FUNCTIONS
float   kdtree_get_element(kdtree_t *t, int i, int j);
#else
#define kdtree_get_element(t, i, j)  ((t)->data[(t)->dataindex[i] * (t)->ndim + (j)])
#endif

/** get data vector via indirection order array
 *
 * This macro returns a pointer to the \p i'th data vector of the
 * kdtree_t#data space by ordered row index \p i, that is mapped to
 * the original row ordering by the indirection array
 * kdtree_t#dataindex.
 *
 * @param t	kd-tree structure
 * @param i	ordered row index (0 <= i < kdtree_t#ndata)
 * @return	pointer to data row \p i
 */
#if DOXYGEN_FUNCTIONS
float  *kdtree_get_vector(kdtree_t *t, int i);
#else
#define kdtree_get_vector(t, i)     ((t)->data + (t)->dataindex[i] * (t)->ndim)
#endif

#define kdtree_get_row_ptr(t, i)    ((t)->data + (i) * (t)->ndim)

#define pow2(x)  (1 << (x))

/** helper function to print a vector \p v of length \p n with stride \p stride to the console */
void vec_post (float *v, int stride, int n, const char *suffix);

/** print only tree info to console */
void kdtree_info_display (kdtree_t* t);
/** print tree raw unsorted data to console */
void kdtree_raw_display  (kdtree_t* t);
/** print tree data to console
    \p print_data = 1 or 2 controls verbosity */
void kdtree_data_display (kdtree_t* t, int print_data);

/** set all counters in kdtree_t#profile to zero */
void profile_clear (kdtree_t *t);

void kdtree_stack_init (kdtree_stack_t *s, int size);
void kdtree_stack_free (kdtree_stack_t *s);
void kdtree_stack_grow (kdtree_stack_t *stack, int alloc);

/** vector to node distance */
float distV2N (kdtree_t* t, const float *x, const int node);
/** vector to node distance with stride */
float distV2N_stride (kdtree_t* t, const float *x, int stride, const int node);
/** vector to node distance with stride and weights 1/sigma */
float distV2N_weighted (kdtree_t* t, const float *x, int stride, const float *sigma, const int node);

void kdtree_set_decomposition(kdtree_t *t, kdtree_dmode_t mode, void *param);

/** free nodes */
void kdtree_clear_nodes (kdtree_t *self);

/** set new data and build tree */
void kdtree_set   (kdtree_t *self, float *data, int m, int n, float *sigma, int use_sigma);

/** set new pointer to weight vector sigma of length kdtree_t#ndim */
void kdtree_set_sigma (kdtree_t *self, float *sigma);

/** check for changes in weights 

Dimensions can be weighted while building the tree, and while
searching. The weight is 1 / sigma, just as in mnm.mahalanobis. Sigma
= 0 means: ignore this dimension.
*/
int  kdtree_update_sigmanz (kdtree_t *self);

/** init tree according to data */
void kdtree_init_data (kdtree_t* t, int h, int vect_num, int dim);

/** build tree according to data 
 *
 * @param t		kd-tree structure
 * @param use_sigma	use weights for distance calculations while building tree
*/
void kdtree_build (kdtree_t* t, int use_sigma);

/** rebuild search tree from changed data or weights 
 *
 * @param t		kd-tree structure
 * @param use_sigma	use weights for distance calculations while rebuilding tree*/
void kdtree_rebuild (kdtree_t* t, int use_sigma);

/** (re-)insert data vectors into tree.  
 *
 * New or changed vectors are already within or appended to
 * kdtree_t::data.  If index < ndata, move to correct node, otherwise
 * insert and increment ndata.
 *
 * @param t	kd-tree structure
 * @param index	start index of vectors to insert
 * @param num	number of vectors to insert
 */
void kdtree_insert (kdtree_t* t, int index, int num);

/** remove data vectors from tree.  
 *
 * Signal removal of rows in kdtree_t::data from search tree.  
 *
 * @param t	kd-tree structure
 * @param index	start index of vectors to remove
 * @param num	number of vectors to remove
 */
void kdtree_delete (kdtree_t* t, int index, int num);


/** Perform search in kd-tree structure \p t.
 *
 * Return the \p k nearest neighbours in multi-dimensional data using
 * the search tree \p t.  The output vector \p d contains the squared
 * distances to the \p r closest data vectors in kdtree_t::data
 * according to the formula
 \f[ 
     d_j = \sum_{j=0}^{ndim} \left( \frac{data_j - x_j}{\sigma} \right)^2
 \f]
 *
 * @param t	kd-tree structure
 * @param x	vector of kdtree_t#ndim elements to search nearest neighbours of
 * @param stride stride in vector \p x
 * @param k	max number of neighbours to find (actual number can be lower)
 * @param r	max squared distance of neighbours to find (\p r = 0 means no limit)
 * @param y	output vector (size == \p r <= \p k) of indices into original data kdtree_t#data 
 * @param d	output vector (size == \p r <= \p k) of squared distances to data vectors 
 * @return \p r = the number of actual neighbours found, 0 <= \p r <= \p k
 */
int  kdtree_search_knn (kdtree_t *t, float* x, int stride, int k, float r, float *y, float *d);

#endif
