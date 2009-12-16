/**
@file	mif.h
@author	Diemo Schwarz
@date	21.11.2008
@version 0.1 

@brief	Metric Inverted File Index Structure



Copyright (C) 2008 - 2009 by IRCAM-Centre Georges Pompidou, Paris, France.
 */


/* build doc with:
   cd distlib; make doc
*/


#ifndef _RTA_MIF_H_
#define _RTA_MIF_H_

#include "rta.h"

#ifdef __cplusplus
extern "C" {
#endif


#define MIF_PROFILE_BUILD  1
#define MIF_PROFILE_SEARCH 1
#define MIF_PROFILE	       (MIF_PROFILE_BUILD || MIF_PROFILE_SEARCH)


/** an "object" stored in the mif index */
typedef struct _mif_object
{
    void *base;		/**< base pointer of the object */
    int   index;	/**< index of the object relative to base */
} mif_object_t;

/** linked list of objects */
typedef struct _mif_pl_entry
{
    struct _mif_pl_entry *next;	/**< next entry */
    mif_object_t obj;		/**< object */
} mif_pl_entry_t;


/** one posting list of ki objects
    indexed (sorted) by start of runs of ref. object order */
typedef struct _mif_postinglist
{
    int 	    size;	/**< number of objects stored */
    mif_pl_entry_t **entries;    /**< array[ki] entries indexed by sort order */
} mif_postinglist_t;


typedef struct _mif_profile_struct
{
	int v2v;	  /**< vector to vector distances */
	int searches;	  /**< searches performed */
} mif_profile_t;


/** distance function pointer */
typedef rta_real_t (*mif_distance_function_t) (mif_object_t *a, mif_object_t *b);


/** k-dimensional search tree data structure */
typedef struct _mif_index
{
    /* parameters */
    mif_distance_function_t distance;  /**< domain distance function */
    int		numobj;	/**< number of data objects */
    int		numref;	/**< number of reference objects */
    int		ki;	/**< number of reference objects used for indexing */
    int		ks;	/**< number of reference objects used for searching */
    int		mpd;	/**< maximum allowed position difference */    

    mif_object_t *refobj;  /**< array[numref] of reference objects */
    mif_postinglist_t *pl; /**< array[numref] of postinglists */

    mif_profile_t profile; /**< profiling data: count internal operations */
} mif_index_t;



/** print only tree info to console */
 void mif_print (mif_index_t* t, int verbosity);


/** set all counters in mif_index_t#profile to zero */
void mif_profile_clear (mif_profile_t *t);


/** initialise index structure */
void mif_init (mif_index_t *self, mif_distance_function_t distfunc, int nr, int ki);

/** free allocated memory */
void mif_free (mif_index_t *self);

/** bulk load new data and index it
    
    @param self		mif structure
    @param base		base pointer to data
    @param numobj	number of data objects on this base pointer

    (@return the number of nodes the tree will build)
*/
int mif_add_data  (mif_index_t *self, int numbase, void **base, int *numbaseobj);


/** build tree 
    
    @param self		kd-tree structure
    @param use_sigma	use weights for distance calculations while building tre    
    \em Prerequisites: 
    - tree data(m, n) must have been set with mif_set_data(), this returned nnodes
    - nodes must have been initialised with mif_init_nodes()
    - if \p use_sigma is on, sigma must have been set with mif_set_sigma().
*/
void mif_build (mif_index_t *self, int use_sigma);


/** TBI: rebuild search tree from changed data or weights 
 *
 * @param t		kd-tree structure
 * @param use_sigma	use weights for distance calculations while rebuilding tree*/
void mif_rebuild (mif_index_t* t, int use_sigma);

/** TBI: (re-)insert data vectors into tree.  
 *
 * New or changed vectors are already within or appended to
 * mif_index_t::data.  If index < ndata, move to correct node, otherwise
 * insert and increment ndata.
 *
 * @param t	kd-tree structure
 * @param index	start index of vectors to insert
 * @param num	number of vectors to insert
 */
void mif_insert (mif_index_t* t, int index, int num);

/** TBI: remove data vectors from tree.  
 *
 * Signal removal of rows in mif_index_t::data from search tree.  
 *
 * @param t	kd-tree structure
 * @param index	start index of vectors to remove
 * @param num	number of vectors to remove
 */
void mif_delete (mif_index_t* t, int index, int num);


/** Perform search in mif index structure
 *
 * Return the \p k nearest neighbours in multi-dimensional data using
 * the search tree \p t.  The output vector \p d contains the squared
 * distances to the \p r closest data vectors in mif_index_t::data
 * according to the reference objects.
 *
 * @param mif	mif structure
 * @param x	vector of mif_index_t#ndim elements to search nearest neighbours of
 * @param ks	number of nearest reference objects to search
 * @param mpd	max position distance ref. objects to take into account
 * @param k	max number of neighbours to find (actual number can be lower)
 * @param r	max squared distance of neighbours to find (\p r = 0 means no limit)
 * @param y	output vector (size == \p r <= \p k) of indices into original data mif_index_t#data 
 * @param d	output vector (size == \p r <= \p k) of squared distances to data vectors 
 * @return \p r = the number of actual neighbours found, 0 <= \p r <= \p k
 */
int mif_search (mif_index_t *mif, mif_object_t* x, int ks, int mpd, int k, const rta_real_t r, 
		/*out*/ rta_real_t *y, rta_real_t *d);


#ifdef __cplusplus
}
#endif

#endif /* _RTA_MIF_H_ */
