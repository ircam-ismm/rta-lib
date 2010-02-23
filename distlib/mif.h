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
@version 0.8

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


\section versions Version History

\subsection v0_8 v0.8 23.02.2010
	- Use zlib for compression of postinglists in sqlite database, reduces db size to 1/3
	- DB is read entirely into memory by the mifquery program (uses up to 1.2 kB/obj, i.e. 1.2 GB of memory for 1M objects).

\subsection v0_7 v0.7 22.02.2010
	- persistent storage of MIF index in sqlite database

\todo
	- don't load db into memory for query, but access postinglist bins from db
	- don't open all data disco files at once but only when accessed (for indexing query object by refobjs)
	- even better: store refobj data vectors in a separate data file, so NO data file needs to be accessed for querying!


\section performance Performance Mesurements v0.7 vs. v0.8 (...z.log)

\subsection build Time for Building Index

<pre>
testindex10k.log:    time for building of index = 3.420000 s, 0.000342 s / obj
testindex10kz.log:   time for building of index = 3.410000 s, 0.000341 s / obj
testindex100k.log:   time for building of index = 157.059998 s, 0.001571 s / obj
testindex100kz.log:  time for building of index = 171.610001 s, 0.001716 s / obj
testindex500k.log:   time for building of index = 1832.589966 s, 0.003665 s / obj
testindex500kz.log:  time for building of index = 1876.540039 s, 0.003753 s / obj
testindex1000k.log:  time for building of index = 6316.720215 s, 0.006317 s / obj
testindex1000kz.log: time for building of index = 6290.030273 s, 0.006290 s / obj
</pre>


\subsection dump Time for Dumping Index to Database

<pre>
testindex10k.log:    time for dumping index = 0.180000 s, 0.000018 s / obj
testindex10kz.log:   time for dumping index = 1.630000 s, 0.000163 s / obj
testindex100k.log:   time for dumping index = 3.380000 s, 0.000034 s / obj
testindex100kz.log:  time for dumping index = 31.250000 s, 0.000312 s / obj
testindex500k.log:   time for dumping index = 7.570048 s, 0.000015 s / obj
testindex500kz.log:  time for dumping index = 67.499901 s, 0.000135 s / obj
testindex1000k.log:  time for dumping index = 21.319679 s, 0.000021 s / obj
testindex1000kz.log: time for dumping index = 178.099716 s, 0.000178 s / obj
</pre>


\subsection load Time for Loading Index from Database

<pre>
testquery10k.log:   time for loading index = 0.030000 s, 0.000003 s / obj
testquery10z.log:   time for loading index = 0.120000 s, 0.000012 s / obj
testquery100k.log:  time for loading index = 0.970000 s, 0.000010 s / obj
testquery100z.log:  time for loading index = 2.820000 s, 0.000028 s / obj
testquery500k.log:  time for loading index = 2.730000 s, 0.000005 s / obj
testquery500z.log:  time for loading index = 6.590000 s, 0.000013 s / obj
testquery1000k.log: time for loading index = 7.210000 s, 0.000007 s / obj
testquery1000z.log: time for loading index = 11.100000 s, 0.000011 s / obj
</pre>


\subsection query Time for Querying 20-NN

<pre>
testquery10k.log:   time for 100 queries = 0.340000 s, 0.003400 s / queryobj
testquery10z.log:   time for 100 queries = 0.300000 s, 0.003000 s / queryobj
testquery100k.log:  time for 100 queries = 7.570000 s, 0.075700 s / queryobj
testquery100z.log:  time for 100 queries = 7.490000 s, 0.074900 s / queryobj
testquery500k.log:  time for 100 queries = 17.570000 s, 0.175700 s / queryobj
testquery500z.log:  time for 100 queries = 15.470000 s, 0.154700 s / queryobj
testquery1000k.log: time for 100 queries = 68.769997 s, 0.687700 s / queryobj
testquery1000z.log: time for 100 queries = 60.410000 s, 0.604100 s / queryobj
</pre>


\subsection acc Simulated Number of Pages Accessed for Querying

<pre>
testquery10z.log:   #bytes accessed in zipped index:    12.159480 MB = 2968 blocks of 4096
testquery10z.log:   #bytes accessed in optimal index:   20.146384 MB = 4918 blocks of 4096
testquery10z.log:   #bytes accessed in index:           40.212768 MB = 9817 blocks of 4096
testquery100z.log:  #bytes accessed in zipped index:   116.033352 MB = 28328 blocks of 4096
testquery100z.log:  #bytes accessed in optimal index:  186.884948 MB = 45626 blocks of 4096
testquery100z.log:  #bytes accessed in index:          373.517096 MB = 91190 blocks of 4096
testquery500z.log:  #bytes accessed in zipped index:   187.737304 MB = 45834 blocks of 4096
testquery500z.log:  #bytes accessed in optimal index:  314.145432 MB = 76695 blocks of 4096
testquery500z.log:  #bytes accessed in index:          628.130864 MB = 153352 blocks of 4096
testquery1000z.log: #bytes accessed in zipped index:   789.029648 MB = 192634 blocks of 4096
testquery1000z.log: #bytes accessed in optimal index: 2022.288396 MB = 493722 blocks of 4096
testquery1000z.log: #bytes accessed in index:         4044.336792 MB = 987386 blocks of 4096
</pre>


\subsection space Size of Index Structure in Memory

<pre>
testindex10kz.log:   TOTAL    4.004000 MB =    4004000 B (20020 B/refobj, 400.399994 B/obj)
testindex100kz.log:  TOTAL  125.830408 MB =  125830408 B (199098 B/refobj, 1258.304077 B/obj)
testindex500kz.log:  TOTAL  398.893288 MB =  398893288 B (282102 B/refobj, 797.786560 B/obj)
testindex1000kz.log: TOTAL 1198.917552 MB = 1198917552 B (599458 B/refobj, 1198.917480 B/obj)
</pre>


\subsection dbsize Size of Index in Database File

<pre>
   5.2M data/testindex10k.db
   1.8M data/testindex10k.dbz
 146.1M data/testindex100k.db
  58.3M data/testindex100k.dbz
 429.1M data/testindex500k.db
 153.5M data/testindex500k.dbz
1224.7M data/testindex1000k.db
 228.5M data/testindex1000k.dbz
</pre>

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
    int indexaccess;	/**< number of accesses to the index entries in the posting list bins */
    int indexaccessbytes; /**< number of (compressed) bytes accessed above */
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
    int		  alloc; /**< number of objects allocated while building index, 
			      byte size of compressed blob rounded to object size while searching */
    mif_object_t *obj;	 /**< pointer to array(num) objects */
} mif_pl_bin_t;


/** one posting list of mif_index_t::ki objects
    indexed (sorted) by start of runs of ref. object order */
typedef struct _mif_postinglist
{
    int 	  size;  /**< number of objects stored */
    mif_pl_bin_t *bin;   /**< array[ki] bins indexed by sort order */
} mif_postinglist_t;

typedef struct _mif_files
{
    int		 nbase;		/**< number of base files */
    int		 ndim;  	/**< copy of number of dims (to avoid access to memory-mapped file header) */
    int 	 descrid;	/**< descriptor ID */
    const char **filename;	/**< array(nbase) of file names */
    void       **base;		/**< array(nbase) of file base pointer */
    int	        *numbaseobj;	/**< array(nbase) of num obj. per base */
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

/** init distance function for query */ 
void mif_init_index (mif_index_t *self, mif_files_t *db);

/** init and allocate posting lists */
void mif_pl_init(mif_postinglist_t *pl, int ki, int kpl);

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
