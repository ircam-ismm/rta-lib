/**
@file	mif.h
@author	Diemo Schwarz
@date	21.11.2009
@brief	Metric Inverted File Index Structure

Copyright (C) 2008 - 2009 by IRCAM-Centre Georges Pompidou, Paris, France.
 */

/** \mainpage

@author	 Diemo Schwarz
@date	 21.11.2009
@version 0.5

@brief	Library for an index on metric spaces using inverted files

The ::mif_index_t and the related functions in mif.h implement efficient
indexing and similarity search algorithms on objects in a metric
space, where only distances are known.  These distances are calculated
by a caller-defined distance function with prototype
::mif_distance_function_t, that operates on anonymous objects represented as ::mif_object_t.

It implements the algorithm described in Amato G., Savino P., <i>Approximate similarity search in metric spaces using inverted files</i>, 
Third Interational ICST Conference on Scalable Information Systems, 2008.

Copyright (C) 2009 by IRCAM-Centre Georges Pompidou, Paris, France.



\section parcompl Parameters and Complexities

This section explains the connection of the algorithm's parameters to
the computational complexity.

There are 4 parameters, the first two used for building the index, the
last two used for searching it.  According to Amato and Savino (2008),
the bounds should be (depending on the number of data vectors \f$ n_{data} \f$):

- mif_index_t::numref number of reference objects 			 
	\f$ n_{ref} \ge  2 \sqrt{n_{data}} \f$ 
- mif_index_t::ki     number of reference objects used for indexing
 	 \f$ k_i    \le  n_{ref} \f$
- mif_index_t::ks     number of reference objects used for searching
 	 \f$ k_s    \le k_i \f$
- mif_index_t::mpd    maximum allowed position difference
	 \f$ mpd    \le k_i \f$

The complexities are then 
(using a balancing factor \f$\alpha \ge 2\f$ for \f$ n_{ref} = \alpha \sqrt{n_{data}} \f$):

- Number of distance calculations for building: 
	\f$  n_{data} \cdot n_{ref}  =   \alpha \, n_{data}^{\frac{3}{2}} \f$
- Number of posting list entries after building: 
	\f$ O\left( k_i \cdot n_{data} \right) \f$

- Number of distance calculations for searching: 
	\f$ n_{ref} =  \alpha \sqrt{n_{data}} \f$
- Number of posting list accesses for searching: 
	\f$  O\left( \frac{k_s \cdot (2 mpd + 1) \cdot n_{data}}{n_{ref}} \right) 
	  =  O\left( k_s \cdot (2 mpd + 1) \cdot \frac{1}{\alpha} \sqrt{n_{data}} \right) \f$

*/


/* build doc with:
   cd distlib; make doc
*/


#ifndef _RTA_MIF_H_
#define _RTA_MIF_H_

#ifdef __cplusplus
extern "C" {
#endif

#include "rta.h"


#define MIF_PROFILE_BUILD  1
#define MIF_PROFILE_SEARCH 1
#define MIF_PROFILE	       (MIF_PROFILE_BUILD || MIF_PROFILE_SEARCH)


/** structure to collect profiling data */
typedef struct _mif_profile_struct
{
    int o2o;		/**< number of object to object distance calculations */
    int searches;	/**< number of searches performed */
    int placcess;	/**< number of posting lists accessed */
    int plbinaccess;	/**< number of posting list bins accessed */
    int indexaccess;	/**< number of accesses to the index entries in the posting lists */
    int numhashobj;	/**< number of hashed objects used for searching */
    int numhashalloc;	/**< number of hashed objects allocated for searching */
    int numhashbin;	/**< number of hashed table bins allocated for searching */
} mif_profile_t;


/** representation of an object to index 

    A data "object" stored in the mif index is represented as an anonymous base
    pointer and an index relative to it.  The pointer can for instance
    refer to a matrix and the index to the rows, or the pointer points
    to a sound file structure, and the index refers to analysis
    windows or frames.  The index code never looks at these values,
    but only passes them to the distance function, defined by the
    caller.
*/
typedef struct _mif_object
{
    int base;	/**< base ID of the object */
    int index;	/**< index of the object relative to base */
} mif_object_t;


/** array of objects */
typedef struct _mif_pl_bin
{
    int		  num;	 /**< number of objects in bin */
    int		  alloc; /**< number of objects allocated while building index, byte size of compressed blob while searching */
    mif_object_t *obj;	 /**< pointer to array(num) objects */
} mif_pl_bin_t;


/** one posting list of mif_index_t::ki objects
    indexed (sorted) by start of runs of ref. object order */
typedef struct _mif_postinglist
{
    int 	  size;   /**< number of objects stored */
    mif_pl_bin_t *bin;   /**< array[ki] bins indexed by sort order */
} mif_postinglist_t;

typedef struct _mif_files
{
    int		nbase;			    /**< number of base files */
    int		ndim;  			    /**< copy of number of dims (to avoid access to memory-mapped file header) */
    int 	descrid;		    /**< descriptor ID */
    char      **filename;		    /**< array(nbase) of file names */
    void      **base;			    /**< array(nbase) of file base pointer */
    int	       *numbaseobj;		    /**< array(nbase) of num obj. per base */
} mif_files_t;


/** distance function 

    This defines a pointer to a function prototype that must be
    defined by the user of the MIF index.  The function is passed two
    pointers to external mif_object_t objects and returns a distance value. */
typedef rta_real_t (*mif_distance_function_t) (void *private, mif_object_t *a, mif_object_t *b);

/** init/deinit of temporary data for distance function */
typedef void (*mif_manage_function_t) (void *private, mif_files_t *database);


/** Metric inverted file index data structure 

    This is the main data structure holding the MIF index.  Note that
    the objects themselves are not stored in the index, only mif_index_t::numref reference objects.
*/
typedef struct _mif_index
{
    /* parameters */
    mif_distance_function_t distance;	    /**< domain distance function */
    mif_manage_function_t   distance_init;  /**< domain init     function */
    mif_manage_function_t   distance_free;  /**< domain cleanup  function */
    void       *distance_private;	    /**< private data for distance function */
    mif_files_t *files;		    /**< link to database */
/* todo: void *database */  

    int		numobj;	/**< total number of data objects */
    int		numref;	/**< number of reference objects */
    int		ki;	/**< number of reference objects used for indexing */
    int		ks;	/**< number of reference objects used for searching */
    int		mpd;	/**< maximum allowed position difference */		

    mif_object_t *refobj;  /**< array[numref] of reference objects */
    mif_postinglist_t *pl; /**< array[numref] of postinglists */

    mif_profile_t profile; /**< profiling data: count internal operations */
} mif_index_t;



/** print index info to console */
void mif_print (mif_index_t* t, int verbosity, char *msg);


/** set all counters in mif_index_t#profile to zero */
void mif_profile_clear (mif_profile_t *t);

/** print profile info */
void mif_profile_print (mif_profile_t *t);


/** initialise index structure */
void mif_init (mif_index_t *self, mif_distance_function_t distfunc, 
	       mif_manage_function_t distinit, mif_manage_function_t distfree, 
	       void *distprivate, int nr, int ki);

/** free allocated memory */
void mif_free (mif_index_t *self);

/* init distance function for query */ 
void mif_init_index (mif_index_t *self, mif_files_t *db);

/** bulk load new data and index it
    
    @param self		mif structure
    @param base		base pointer to data
    @param numobj	number of data objects on this base pointer

    @return the number of objects in the data
*/
int mif_add_data  (mif_index_t *self, mif_files_t *db);


/** TBI: rebuild search index from changed data or weights 
 *
 * @param t		index structure
 * @param use_sigma	use weights for distance calculations while rebuilding index*/
void mif_rebuild (mif_index_t* t, int use_sigma);

/** TBI: (re-)insert data vectors into index.  
 *
 * New or changed vectors are already within or appended to
 * mif_index_t::data.  If index < ndata, move to correct node, otherwise
 * insert and increment ndata.
 *
 * @param t	index structure
 * @param index	start index of vectors to insert
 * @param num	number of vectors to insert
 */
void mif_insert (mif_index_t* t, int index, int num);

/** TBI: remove data vectors from index.  
 *
 * Signal removal of rows in mif_index_t::data from search index.  
 *
 * @param t	index structure
 * @param index	start index of vectors to remove
 * @param num	number of vectors to remove
 */
void mif_delete (mif_index_t* t, int index, int num);


/** Perform search in mif index structure
 *
 * Return the \p k nearest neighbours in MIF index.
 * The output vector \p dist contains the transformed
 * distances to the \p r closest data objects in \p obj
 * according to the reference objects.
 *
 * The search uses these parameters: the mif_index_t::ks number of nearest reference objects to search and the 
 * mif_index_t::mpd max position distance ref. objects to take into account.
 *
 * @param mif	mif structure
 * @param query	object to search nearest neighbours of
 * @param k	max number of neighbours to find (actual number can be lower)
 * @param obj	output vector (size == \p r <= \p k) of mif_object_t object structures
 * @param dist	output vector (size == \p r <= \p k) of transformed distances to data objects
 * @return \p r = the number of actual neighbours found, 0 <= \p r <= \p k
 */
int mif_search_knn (mif_index_t *self, mif_object_t *query, int k, 
		    /*out*/ mif_object_t *obj, int *dist) ;

/* @param ks	number of nearest reference objects to search
 * @param mpd	max position distance ref. objects to take into account
 * @param r	max squared distance of neighbours to find (\p r = 0 means no limit) */
//int mif_search (mif_index_t *mif, mif_object_t* x, int ks, int mpd, int k, const rta_real_t r, 
//		/*out*/ rta_real_t *y, rta_real_t *d);


#ifdef __cplusplus
}
#endif

#endif /* _RTA_MIF_H_ */
