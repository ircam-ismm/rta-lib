/**
 * @file rta_kdtree.h
 * @author Diemo Schwarz
 * @date 29.9.2008
 * @version  1.0
 * @ingroup rta_recognition
 *
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

// build doc with:
// cd distlib; make doc

#ifndef _RTA_KDTREE_H_
#define _RTA_KDTREE_H_

#include "rta.h"
#include "rta_bpf.h"

#ifdef __cplusplus
extern "C" {
#endif


#define RTA_KDTREE_PROFILE_BUILD  1
#define RTA_KDTREE_PROFILE_SEARCH 1
#define RTA_KDTREE_PROFILE (RTA_KDTREE_PROFILE_BUILD || RTA_KDTREE_PROFILE_SEARCH)


#define RTA_USE_DISTFUNC 1
#define RTA_KDTREE_MAX_DISTFUNC 256


#if RTA_USE_DISTFUNC // uses rta_bpf_t, (data-compatible to FTM bpfunc_t)
#  define RTA_DMAPW(x, y, s, dfun) ((dfun) ? rta_bpf_get_interpolated(dfun, ((x) - (y))) / (s) : ((x) - (y)) / (s))
#else
#  define RTA_DMAPW(x, y, s, dfun) ((x) - (y)) / (s)
#endif /* RTA_USE_DISTFUNC */

#if RTA_USE_DISTFUNC // uses rta_bpf_t, (data-compatible to FTM bpfunc_t)
#  define RTA_DMAP(x, y, dfun) ((dfun) ? rta_bpf_get_interpolated(dfun, ((x) - (y))) : ((x) - (y)))
#else
#  define RTA_DMAP(x, y, dfun) ((x) - (y))
#endif /* RTA_USE_DISTFUNC */



/** stack element for search algorithm */
typedef struct _kdtree_stack_elem_struct
{
  int node;
  rta_real_t dist;
} rta_kdtree_stack_elem_t;

/** optimised stack for search algorithm */
typedef struct _kdtree_stack_struct
{
  int size;
  int alloc;
  rta_kdtree_stack_elem_t *buffer;
} rta_kdtree_stack_t;


/** decomposition mode
 *
 * Different ways to split data space at each tree level:
 * \li orthogonal to dimension (fastest tree-building speed),
 * \li by arbitrary hyperplane (for testing only), or
 * \li by PCA (optimal decomposition)
 */
typedef enum
{
  dmode_orthogonal,   /**< optimised orthogonal to axes */
  dmode_hyperplane,   /**< hyperplane orthogonal to axes */
  dmode_pca         /**< hyperplane along principal components */
} rta_kdtree_dmode_t;

/** pivot calculation mode
 *
 * the \em pivot is the mean vector to split each tree node at
 */
typedef enum
{
  mmode_mean,   /**< mean of values / distances to splitplane */
  mmode_middle, /**< middle between min/max */
  mmode_median  /**< true median (guarantees equal number of points left
                   and right and thus a well-balanced and optimal tree) */
} rta_kdtree_mmode_t;


/** one node of the kd-tree */
typedef struct _kdtree_node_struct
{
  int startind; /**< index of first vector in node in dataindex array */
  int endind;   /**< index of last vector in node in dataindex array */
  int size;   /**< number of vectors in node */
  int splitdim; /**< for dmode orthogonal, dimension along which node is split*/

  rta_real_t  splitnorm;  /**< spatial length of split vector */
} rta_kdtree_node_t;

typedef struct _kdtree_profile_struct
{
  int v2v;    /**< vector to vector distances */
  int v2n;    /**< vector to node distances */
  int mean;   /**< mean vector calculations */
  int hyperp;    /**< split plane calculations */
  int searches;   /**< searches performed */
  int neighbours;   /**< neighbours found */
  int maxstack;   /**< highest stack size */
} rta_kdtree_profile_t;

/** struct holding a block index and an element index (matrix row) */
typedef struct _kdtree_object_struct
{
  int base;   /**< number of block (in data array) */
  int index;    /**< row vector index inside block */
} rta_kdtree_object_t;

/** k-dimensional search tree data structure */
typedef struct _kdtree_struct
{
  rta_kdtree_dmode_t dmode; /**< decomposition mode */
  rta_kdtree_mmode_t mmode; /**< pivot calculation mode */

  int     ndim;   /**< Dimension of vectors */
  int     ndatatot;   /**< Number of total vectors */
  int     nblocks;    /**< Number of blocks of data */
  int     *ndata;   /**< Number of vectors per block*/
  rta_real_t **data;    /**< nblocks pointers to data matrices (ndata, ndim) */
  rta_kdtree_object_t *dataindex;   /**< data vector indirection array (ndata):
           original index of data vector at tree array position  */

  rta_real_t *sigma;    /**< 1/weight, 0 == inf */
  int     sigma_nnz;    /**< number of non-zero sigma */
  int    *sigma_indnz;  /**< non-zero sigma lines */
  rta_bpf_t *dfun[RTA_KDTREE_MAX_DISTFUNC]; /* distance transfer functions */

  int     height;   /**< Height of the kdtree */
  int     maxheight;    /**< Maximal height of the kdtree */
  int     givenheight;  /**< Height given by user, gives tree height
         if positive, subtract from maxheight if negative */
  int     nnodes;     /**< Number of nodes (must be a power of 2) */
  int     ninner;     /**< Number of inner nodes (=index of first leaf node)*/
  rta_kdtree_node_t *nodes; /**< nodes (nnodes) */
  rta_real_t *mean;   /**< mean vectors in nnodes rows (todo: median), always present */
  rta_real_t *split;    /**< hyperplanes A1*X1 + A2*X2 +...+ An*Xn + An+1 = 0,
         in nnodes rows or NULL in dmode_orthogonal */

  int     sort;   /**< sort search result by distance */
  rta_kdtree_stack_t stack;

    /** profiling data: count internal operations */
  rta_kdtree_profile_t profile;

} rta_kdtree_t;


extern const char *rta_kdtree_dmodestr[];
extern const char *rta_kdtree_mmodestr[];


/** get data element via indirection order array
 *
 * This macro returns the element of the kdtree_t#data space by ordered row
 * index \p i, that is mapped to the original row ordering by the
 * indirection array kdtree_t#dataindex.
 *
 * @param t kd-tree structure
 * @param i ordered row index (0 <= i < kdtree_t#ndata)
 * @param j column index      (0 <= j < kdtree_t#ndim)
 * @return  data element at row \p i, column \p j
 */
#if DOXYGEN_FUNCTIONS
rta_real_t rta_kdtree_get_element(rta_kdtree_t *t, int i, int j);
#else
#define rta_kdtree_get_element(t, i, j) ((t)->data[(t)->dataindex[i].base][(t)->dataindex[i].index * (t)->ndim + (j)])
#endif

/** get data vector via indirection order array
 *
 * This macro returns a pointer to the \p i'th data vector of the
 * kdtree_t#data space by ordered row index \p i, that is mapped to
 * the original row ordering by the indirection array
 * kdtree_t#dataindex.
 *
 * @param t kd-tree structure
 * @param i ordered row index (0 <= i < kdtree_t#ndata)
 * @return  pointer to data row \p i
 */
#if DOXYGEN_FUNCTIONS
rta_real_t *rta_kdtree_get_vector(rta_kdtree_t *t, int i);
#else
#define rta_kdtree_get_vector(t, i) ((t)->data[(t)->dataindex[i].base] + (t)->dataindex[i].index * (t)->ndim)
#endif

/* @param k number of data block */
#define rta_kdtree_get_row_ptr(t, k, i) ((t)->data[k] + (i) * (t)->ndim)
// old: #define kdtree_get_row_ptr(t, i)        kdtree_get_row_ptr2(t, 0, i)

/** print only tree info to console */
void rta_kdtree_info_display (rta_kdtree_t* t);

/** print tree raw unsorted data to console */
void rta_kdtree_raw_display  (rta_kdtree_t* t);
  
/** print tree data to console
    \p print_data = 1 or 2 controls verbosity */
void rta_kdtree_data_display (rta_kdtree_t* t, int print_data);

/** set all counters in kdtree_t#profile to zero */
void rta_kdtree_profile_clear (rta_kdtree_t *t);

/** set decomposition mode */
void rta_kdtree_set_decomposition (rta_kdtree_t *t, rta_kdtree_dmode_t mode, void *param);

/** set pivot mode */
void rta_kdtree_set_pivot (rta_kdtree_t *t, rta_kdtree_mmode_t mode);

/** initialise tree structure */
void rta_kdtree_init (rta_kdtree_t *self);

/** free auto-allocated tree memory

    Call this ONLY when you gave the required node memory pointers as NULL.
*/
void rta_kdtree_free (rta_kdtree_t *self);

/** set new data vector and size

    @param self   kd-tree structure
    @param nblocks  number of data matrices to index
    @param data   pointer to an array(\p nblocks) of pointers to \p m * \p n real values of data to be searched
    @param index  if not NULL, must point to space for \p m int values with row indices to use (1..mdata), otherwise the library will auto-allocate
    @param m    pointer to an array(\p nblocks) of the numbers of data vectors (rows) in each data matrix
    @param n    dimension of data vectors (columns)

    @return the number of nodes the tree will build

    \p data and \p index must point to arrays in memory that will stay
    valid throughout the existence of the tree!
*/
int rta_kdtree_set_data (rta_kdtree_t *self, int nblocks, rta_real_t **data, rta_kdtree_object_t *index, int *m, int n);

/** build tree only on m lines of data listed in ind */
// int kdtree_set_data_ind (kdtree_t *self, rta_real_t *data, int *index, int m, int n, int *ind);

/** set new pointer to weight vector sigma of length kdtree_t#ndim

Dimensions can be weighted while building the tree, and while
searching. The weight is 1 / sigma, just as in mnm.mahalanobis. Sigma
= 0 means: ignore this dimension.
*/
void rta_kdtree_set_sigma (rta_kdtree_t *self, rta_real_t *sigma);

/** check for changes in weights

    This functions searches the non-zero elements in kdtree_t#sigma.
    Only these dimensions will be taken into account for searching.
*/
int rta_kdtree_update_sigmanz (rta_kdtree_t *self);


/** initialise tree nodes

    This function must be called after kdtree_set_data() and before
    kdtree_build().  It initialises the tree nodes and the vectors
    decomposing the search space.

    @param self   kd-tree structure
    @param nodes  space for tree nodes or NULL for automatic allocation
    @param means  space for mean vectors or NULL for automatic allocation
    @param planes space for hyperplanes or NULL for automatic allocation

    \em Prerequisites:
    - tree data(m, n) must have been set with rta_kdtree_set_data(), this returned nnodes
    - if \p nodes is not NULL, it must point to space for the tree nodes of size nnodes * sizeof(kdtree_node_t), which must have been allocated outside of the library
    - if \p means is not NULL, it must point to space for the mean vectors of size nnodes / 2 * ndim, which must have been allocated outside of the library
    - if \p planes is not NULL and decomposition mode is not dmode_orthogonal, it must point to space for split hyperplane base vectors of size nnodes / 2 * ndim, which must have been allocated outside of the library
*/
void rta_kdtree_init_nodes (rta_kdtree_t *self, rta_kdtree_node_t *nodes, rta_real_t *means, rta_real_t *planes);

/** build tree

    @param self   kd-tree structure
    @param use_sigma  use weights for distance calculations while building tre
    \em Prerequisites:
    - tree data(m, n) must have been set with kdtree_set_data(), this returned nnodes
    - nodes must have been initialised with kdtree_init_nodes()
    - if \p use_sigma is on, sigma must have been set with kdtree_set_sigma().
*/
void rta_kdtree_build (rta_kdtree_t *self, int use_sigma);

/** TBI: rebuild search tree from changed data or weights
 *
 * @param t   kd-tree structure
 * @param use_sigma use weights for distance calculations while rebuilding tree*/
void rta_kdtree_rebuild (rta_kdtree_t* t, int use_sigma);

/** TBI: (re-)insert data vectors into tree.
 *
 * New or changed vectors are already within or appended to
 * kdtree_t::data.  If index < ndata, move to correct node, otherwise
 * insert and increment ndata.
 *
 * @param t kd-tree structure
 * @param index start index of vectors to insert
 * @param num number of vectors to insert
 */
void rta_kdtree_insert (rta_kdtree_t* t, int index, int num);

/** TBI: remove data vectors from tree.
 *
 * Signal removal of rows in kdtree_t::data from search tree.
 *
 * @param t kd-tree structure
 * @param index start index of vectors to remove
 * @param num number of vectors to remove
 */
void rta_kdtree_delete (rta_kdtree_t* t, int index, int num);

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
 * @param t kd-tree structure
 * @param x vector of kdtree_t#ndim elements to search nearest neighbours of
 * @param stride stride in vector \p x
 * @param k max number of neighbours to find (actual number can be lower)
 * @param r max squared distance of neighbours to find (\p r = 0 means no limit)
 * @param use_sigma use weights set by #rta_kdtree_set_sigma
 * @param y output vector (size == \p r <= \p k) of (base, element) indices into original data kdtree_t#data
 * @param d output vector (size == \p r <= \p k) of squared distances to data vectors
 * @return \p n = the number of actual neighbours found, 0 <= \p n <= \p k
 */
int rta_kdtree_search_knn (rta_kdtree_t *t, rta_real_t* x, int stride, int k, const rta_real_t r, int use_sigma,
                           /*out*/ rta_kdtree_object_t *y, rta_real_t *d);

/**
 * Weighted squared vector distance (v1 - v2)^2
 */
rta_real_t rta_euclidean_distance (rta_real_t* v1, int stride1,
                                   rta_real_t* v2, int dim,
                                   rta_bpf_t  *distfunc[]);

rta_real_t rta_weighted_euclidean_distance (rta_real_t* v1, rta_real_t* v2,
                                            rta_real_t *sigma, int ndim,
                                            rta_bpf_t  *distfunc[]);

rta_real_t rta_weighted_euclidean_distance_stride (rta_real_t* v1, int stride1,
                                                   rta_real_t* v2,
                                                   rta_real_t *sigma, int ndim,
                                                   rta_bpf_t  *distfunc[]);


#ifdef __cplusplus
}
#endif

#endif /* _RTA_KDTREE_H_ */
