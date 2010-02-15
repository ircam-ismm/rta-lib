/**
@file	mif.c
@author	Diemo Schwarz
@date	30.11.2008
@version 0.1 

@brief	Metric Inverted File Index Structure

Copyright (C) 2008 - 2009 by IRCAM-Centre Georges Pompidou, Paris, France.
*/



#ifdef WIN32
#include <malloc.h>
static double log2(double x){ return log(x)/log(2);}
#else
#include <alloca.h>
#endif

#include <assert.h>
#include <string.h>

#include "mif.h"
#include "mifhash.h"
#include "rta_math.h"
#include "rta_util.h"

#define MAX_FLOAT 0x7FFFFFFF


/*
 * posting list handling 
 */

#define PLBLOCKSIZE 4
#define PLBLOCKS(n) ((int) (n / PLBLOCKSIZE))

void mif_pl_init(mif_postinglist_t *pl, int ki, int kpl)
{
    int i;

    pl->size = 0;
    pl->bin  = rta_malloc(ki * sizeof(mif_pl_bin_t));

    for (i = 0; i < ki; i++)
    {
	kpl = PLBLOCKS(kpl) * PLBLOCKSIZE;
	pl->bin[i].num   = 0;
	pl->bin[i].alloc = kpl;
	pl->bin[i].obj   = rta_malloc(kpl * sizeof(mif_object_t));
    }
}

void mif_pl_free(mif_postinglist_t *pl, int ki)
{
    int i;

    for (i = 0; i < ki; i++)
	rta_free(pl->bin[i].obj);

    rta_free(pl->bin);
}

/* insert object into correct bin of posting list pl of a reference object with order k */
void mif_pl_insert (mif_postinglist_t *pl, mif_object_t *newobj, int k)
{
    mif_pl_bin_t *bin = &pl->bin[k];

    /* grow list if necessary */
    if (bin->num + 1 >= bin->alloc)
    {
	bin->alloc += PLBLOCKSIZE;
	bin->obj    = realloc(bin->obj, bin->alloc * sizeof(mif_object_t));
	assert(bin->obj);
    }

    bin->obj[bin->num].base  = newobj->base;
    bin->obj[bin->num].index = newobj->index;

    bin->num++;
    pl->size++;
}



/*
 *  index structure handling
 */

/** initialise index structure */
void mif_init (mif_index_t *self, mif_distance_function_t distfunc, 	    
	       mif_manage_function_t distinit, mif_manage_function_t distfree, 
	       void *distprivate,  int nr, int ki)
{
    self->distance_private = distprivate;
    self->distance	   = distfunc;
    self->distance_init    = distinit;
    self->distance_free    = distfree;
    self->numobj = 0;
    self->numref = nr;
    self->ki     = ki;
    self->ks     = 0;
    self->refobj = rta_malloc(nr * sizeof(mif_object_t));
    self->pl     = rta_zalloc(nr * sizeof(mif_postinglist_t));

    mif_profile_clear(&self->profile);
}


/** free allocated memory */ 
void mif_free (mif_index_t *self)
{
    int i;

    if (self->distance_free)
	(*self->distance_free)(self->distance_private, NULL);

    /* free posting lists */
    for (i = 0; i < self->numref; i++)
	mif_pl_free(&self->pl[i], self->ki);

    rta_free(self->refobj); 
    rta_free(self->pl);
}


void mif_print_pl(mif_index_t *self, int i)
{
    mif_postinglist_t *pl = &self->pl[i];
    int j, k;

    rta_post("pl %d  refobj %d  size %d:\n", i, self->refobj[i].index, pl->size);

    for (k = 0; k < self->ki; k++)
    {
	mif_pl_bin_t *bin = &pl->bin[k];

	rta_post("  <%d: ", k);
	for (j = 0; j < bin->num; j++)
	{
	    rta_post("%d ", bin->obj[j].index);
	    //rta_post("%p.%d ", entry->obj.base, entry->obj.index);
	}
	rta_post(">\n");
    }
    rta_post("\n");
}

void mif_print (mif_index_t *self, int verb)
{
    rta_post("\nMIF index info:\n");
    rta_post("n_obj    = %d\n", self->numobj);
    rta_post("n_ref    = %d\n", self->numref);
    rta_post("k_i      = %d\n", self->ki);
    rta_post("k_s      = %d\n", self->ks);
    rta_post("mpd      = %d\n", self->mpd);
 
    if (verb >= 1)
    {	/* space used */
	int i, npe = 0;
	unsigned long spl     = self->numref * sizeof(*self->pl);
	unsigned long srefobj = self->numref * sizeof(mif_object_t);
	unsigned long stotal;

	for (i = 0; i < self->numref; i++)
	    npe += self->pl[i].size;

	spl += npe * sizeof(mif_pl_bin_t);
	stotal = srefobj + spl;

	rta_post("\nspace for struct              = %14lu B\n", sizeof(mif_index_t));
	rta_post("%9d ref.obj.            = %14lu B\n", self->numref, srefobj);
	rta_post("%9d postinglist entries = %14lu B", npe, spl);
	rta_post(" (%d #/refobj, %f #/bin)\n", 
		 npe / self->numref, (float) npe / self->numref / self->ki);
	rta_post("TOTAL %20f MB = %14lu B (%lu B/refobj, %f B/obj)\n", 
		 stotal / 1e6, stotal, stotal / self->numref, (float) stotal / self->numobj);
    }

    if (verb >= 2)
    {
	int i;

	rta_post("\npostinglists: ro -> <order: object...>  i.e.  for object, ro is order-closest ref.obj\n");

	for (i = 0; i < self->numref; i++)
	    mif_print_pl(self, i);
    }
}

/** set all counters in mif_index_t#profile to zero */ 
void mif_profile_clear (mif_profile_t *t)
{
    t->o2o = 0;
    t->searches = 0;
    t->placcess = 0;
    t->plbinaccess = 0;
    t->indexaccess = 0;
    t->numhashobj  = 0;
    t->numhashalloc  = 0;
    t->numhashbin  = 0;
}

/** print profile info */ 
void mif_profile_print (mif_profile_t *t)
{
    int n = t->searches ? t->searches : 0x7fffffff;
    rta_post("\nProfile: (total / per search)\n"
	     "#searches:    %9d\n" 
	     "#dist:        %9d %6d\n"
	     "#placcess:    %9d %6d\t(%ld bytes each)\n"
	     "#binaccess:   %9d %6d\t(%ld bytes each)\n"
	     "#entryaccess: %9d %6d\t(%ld bytes each)\n"
	     "#hashobj:     %9d %6d\t(%ld bytes each)\n"
	     "#hashalloc:   %9d %6d\t(%ld bytes each)\n"
	     "#hashbin:     %9d %6d\t(%ld bytes each)\n",
	     t->searches, 
	     t->o2o,	     t->o2o / n, 	     
	     t->placcess,    t->placcess / n,    sizeof(mif_postinglist_t), 
	     t->plbinaccess, t->plbinaccess / n, sizeof(mif_pl_bin_t *), 
	     t->indexaccess, t->indexaccess / n, sizeof(mif_pl_bin_t), 
	     t->numhashobj,  t->numhashobj  / n, HASHOBJSIZE,
	     t->numhashalloc, t->numhashalloc / n, HASHOBJSIZE,
	     t->numhashbin,  t->numhashbin  / n, sizeof(unsigned int));
}


/*
 *  build index
 */

/*choose reference objects and fill refobj array in mif struct */ 
static int mif_choose_refobj (mif_index_t *self, mif_files_t *files)
{	
    int numbase = files->nbase;
    int numobj = 0;
    int *cumobj = alloca((numbase + 1) * sizeof(int));
    int *sample = alloca(self->numref * sizeof(int));
    int i;

    /* how many elements are given */
    for (i = 0, cumobj[0] = 0; i < numbase; i++)
    {
	cumobj[i + 1] = numobj += files->numbaseobj[i];
    }
	
    /* choose reference objects: draw numref random indices */
    rta_choose_k_from_n(self->numref, numobj, sample);
	
    /* lookup base and relative index */
    for (i = 0; i < self->numref; i++)
    {
	int objind  = sample[i];
	int baseind = rta_find_int(objind, numbase, cumobj + 1);
	int relind  = objind - cumobj[baseind];
	    
	self->refobj[i].base  = baseind;
	self->refobj[i].index = relind;
    }

    return numobj;
}

/** Find k closest reference objects for object <base[b], index> and their ordering 

  newobj	object to index
  indx[k]	out: ref. obj. index is list of ref.obj. index ordered by distance dist
  dist[k]	out: distance between obj and refobj indx[i]

  @return	number k of reference objects found (k <= ki)
*/ 
static int mif_index_object (mif_index_t *self, mif_object_t *newobj, int k,
			     /*out*/ int *indx, rta_real_t *dist)
{
    int r, kmax = 0; /* numfound */

    /* init distances */ 
    for (r = 0; r < k; r++) 
	dist[r] = MAX_FLOAT;

    for (r = 0; r < self->numref; r++)
    {
	rta_real_t d = (*self->distance)(self->distance_private, &self->refobj[r], newobj);
	    
	if (d <= dist[kmax]) 
	{   /* return original index in data and distance */
	    int pos = kmax;	/* where to insert */

	    if (kmax < k - 1)
	    {   /* first move or override */
		dist[kmax + 1] = dist[kmax];
		indx[kmax + 1] = indx[kmax];
		kmax++;
	    }

	    /* insert into sorted list of distance */
	    while (pos > 0  &&  d < dist[pos - 1])
	    {   /* move up */
		dist[pos] = dist[pos - 1];
		indx[pos] = indx[pos - 1];
		pos--;
	    }

	    indx[pos] = r;
	    dist[pos] = d;
	}
    }

    /* return actual number of found objects, can be less than k,
       then kmax is the index of the next one to find */
    return kmax + (dist[kmax] < MAX_FLOAT);
}


/* index data objects by reference objects, build postinglists */ 
static void mif_build_index (mif_index_t *self, mif_files_t *files)
{  
    int          *indx = alloca(self->ki * sizeof(*indx));	/* ref. obj. index */
    rta_real_t   *dist = alloca(self->ki * sizeof(*dist));	/* distance to obj */
    mif_object_t  newobj;	/* object to index */
    int		  numbase = files->nbase;
    int		  b, i, k, kfound;

    /* init distance function with file list */
    (*self->distance_init)(self->distance_private, files);

    /* for each object */
    for (b = 0; b < numbase; b++)
    {
	newobj.base = b;

	for (i = 0; i < files->numbaseobj[b]; i++)
	{
	    newobj.index = i;
	    	    	
	    /* find ki closest reference objects for object <base[b], i> and their ordering */
	    kfound = mif_index_object(self, &newobj, self->ki, indx, dist);
	    /* now indx is list of ref.obj. ordered by distance dist */

	    /* insert object into posting lists of ki closest reference objects with order k */
	    for (k = 0; k < kfound; k++)
	    {
		mif_pl_insert(&self->pl[indx[k]], &newobj, k);
	    }

#if MIF_PROFILE_BUILD
	    self->profile.o2o += self->numref;
#endif
	}
    }
}


/** bulk load new data and index it
    
    @param self		mif structure
    @param numbase	number of base data blocks
    @param base		array[numbase] of base pointer to data
    @param numobj	array[numbase] of number of data objects on this base pointer

    (@return the number of nodes the tree will build)
*/ 
int mif_add_data (mif_index_t *self, mif_files_t *db)
{
    int i;

    self->files = db;

    /* init posting lists */
    for (i = 0; i < self->numref; i++)
	mif_pl_init(&self->pl[i], self->ki, self->numobj / self->numref);

    /* choose numref reference objects and fill refobj array in mif struct */
    self->numobj = mif_choose_refobj(self, db);

    /* index data */
    mif_build_index(self, db);

    return self->numobj;
}



/******************************************************************************
 *
 *	Searching
 *
 */

#define MIF_DEBUG_SEARCH 0

/* out:    indx[K] = index of the Kth nearest neighbour
   return: actual number of found neighbours */
int mif_search_knn (mif_index_t *self, mif_object_t *query, int k, 
		    /* out */ mif_object_t *obj, int *dist) 
{
    rta_real_t *qdist = alloca(self->ks * sizeof(*qdist));	/* distance of query to refobj */
    int        *qind  = alloca(self->ks * sizeof(*qind));	/* index of closest refobj */
    int r, kq, kmax = 0;
    int hashpredict = self->ks * (2 * self->mpd + 1) * self->numobj / self->numref;

    /* preliminary array to allow iterating over hash.  WARNING: not reentrant */
    static int		   hash_allocated = 0;

    /* hash of objects seen in ks closest ref. objects posting lists */
    static fts_hashtable_t hash;

    if (hash_allocated == 0)
    {
	fts_hashtable_init(&hash, hashpredict, self->ks * self->mpd);
	hash_allocated = 1;
    }
    else
    {
	fts_hashtable_clear(&hash);
    }

    /* index query object by ks-closest ref. objects, sorted by
        distance to query object. (cost: nref distance calculations) */
    kq = mif_index_object(self, query, self->ks, qind, qdist);

#if MIF_DEBUG_SEARCH >= 2
    rta_post("query %d indexed %d: (indx, dist) ", query->index, kq);
    for (r = 0; r < kq; r++)
	rta_post("(ro %d = obj %d, %f) ", 
		 qind[r], self->refobj[qind[r]].index, qdist[r]);
    rta_post("\n");
#endif

    /* induced limited spearman footrule distance:
       for all ks-closest ref. obj. (sort order position r) */
    for (r = 0; r < kq; r++)
    {   /* */
#undef rta_max
#undef rta_min
#define rta_max(a,b) (((a)>(b))? (a) : (b))
#define rta_min(a,b) (((a)<(b))? (a) : (b))

	int minp = rta_max(0,            r - self->mpd);
	int maxp = rta_min(self->ki - 1, r + self->mpd);
	mif_postinglist_t *pl = &self->pl[qind[r]];
	int p;

//	rta_post("  r %d  refobj #%d=obj %d:  MPD range %d..%d\n", 
//		 r, qind[r], self->refobj[qind[r]].index, minp, maxp);
#if MIF_PROFILE_SEARCH
	    self->profile.placcess++;
#endif

	/* go through posting list order range between minp and maxp */
	for (p = minp; p <= maxp; p++)
	{
	    mif_pl_bin_t *bin = &pl->bin[p];
	    int i;

#if MIF_PROFILE_SEARCH
	    self->profile.plbinaccess++;
#endif
#if MIF_DEBUG_SEARCH
	    rta_post("  accessing bin %d of refobj %d (pl size %d).\n", p, r, pl->size);
#endif
	    for (i = 0; i < bin->num; i++)
	    {
		int *accumulator;

		if (fts_hashtable_put(&hash, &bin[i].obj, &accumulator) == 0)
		{ /* new occurence of obj: init entry inserted into accumulator hash */
		    *accumulator = self->ki * self->ks;	/* largest possible distance */
		}

		/* correct accumulated distance by actual transformed distance */
		*accumulator += abs(r - p) - self->ki;
		//rta_post("  dist %d  accu %s = %d\n", p, keybuf, *accu);

#if MIF_PROFILE_SEARCH
		self->profile.indexaccess++;
#endif
	    }   /* end while posting list entry bin not empty */
	}
    }


    /* init distances */ 
    for (r = 0; r < k; r++) 
	dist[r] = MAX_FLOAT;

#if MIF_DEBUG_SEARCH >= 2   
    rta_post("  dist of %d hashed obj: ", numhashobj);
#endif

    /* iterate through hash and pick k lowest distances */
    for (r = 1; r <= hash.count; r++)
    {
	int  dtrans;

	dtrans = hash.cell_heap[r].accumulator;

#if MIF_DEBUG_SEARCH >= 2
	rta_post("(%d, %d)  ", hash.cell_heap[r].obj->index, dtrans);
#endif
	    
	if (dtrans <= dist[kmax]) 
	{   /* return original index in data and distance */
	    int pos = kmax;	/* where to insert */

	    if (kmax < k - 1)
	    {   /* first move or override */
		dist[kmax + 1] = dist[kmax];
		obj [kmax + 1] = obj [kmax];
		kmax++;
	    }

	    /* insert into sorted list of distance */
	    while (pos > 0  &&  dtrans < dist[pos - 1])
	    {   /* move up */
		dist[pos] = dist[pos - 1];
		obj [pos] = obj [pos - 1];
		pos--;
	    }

	    obj [pos] = *hash.cell_heap[r].obj;  /* de-hash */
	    dist[pos] = dtrans;
	}
    }
#if MIF_DEBUG_SEARCH >= 2   
    rta_post("\n");
#endif

#if MIF_PROFILE_SEARCH
    self->profile.numhashobj   += hash.count;
    self->profile.numhashalloc += hash.alloc;
    self->profile.numhashbin   += hash.length;
    self->profile.o2o	       += self->numref;
    self->profile.searches++;
#endif

    /* return actual number of found objects, can be less than k,
       then kmax is the index of the next one to find */
    return kmax + (dist[kmax] < MAX_FLOAT);
}


#if MIF_BUILD_TEST

#include "rta_kdtree.h"

#define NDIM 2

static rta_real_t mif_euclidean_distance (mif_object_t *a, mif_object_t *b)
{
/*    return rta_euclidean_distance(a->base + a->index * NDIM, 1,
				  b->base + b->index * NDIM, NDIM);
*/
    rta_real_t* v1 = (rta_real_t *) a->base + a->index * NDIM;
    rta_real_t* v2 = (rta_real_t *) b->base + b->index * NDIM;
    rta_real_t sum = 0;
    int i;

    for (i = 0; i < NDIM; i++) 
    {
	rta_real_t diff = v2[i] - v1[i];
	sum += diff * diff;
    }

    return sum;
}


int main (int argc, char *argv[])
{
#   define      nrow 3
#   define	ki   3
#   define	K    3

    mif_index_t mif;
    int		nref, nobj = 10; 	/* (nrow * 2 + (nrow - 2))*/
    rta_real_t  *data = rta_malloc(sizeof(rta_real_t) * nobj * NDIM);
    int i, j, kfound;

    /* create some test data */
    for (i = 0; i < nobj; i++)
    {
	data[i * NDIM] = i;
	data[i * NDIM + 1] = 0;
    }

    /* num. ref. obj. >= 2 * sqrt(nobj) */
    nref = 2 * sqrt(nobj);

    mif_init(&mif, mif_euclidean_distance, nref, ki);
    mif_add_data(&mif, 1, (void **) &data, &nobj);
    mif.ks  = 3;
    mif.mpd = 3;
    mif_print(&mif, 1);
    mif_profile_print(&mif.profile);
    mif_profile_clear(&mif.profile);

    { /* test searching */
	mif_object_t *obj  = alloca(K * sizeof(*obj));	/* objects */
	int          *dist = alloca(K * sizeof(*dist));	/* transformed distance to obj */
	mif_object_t  query;

	for (i = 0; i < nobj; i++)
	{
	    query.base = data;
	    query.index = i;
	    
	    kfound = mif_search_knn(&mif, &query, K, obj, dist);

	    rta_post("--> %d-NN of query obj %d (found %d):  ", K, i, kfound);
	    for (j = 0; j < kfound; j++)
		rta_post("%d ", obj[j].index);
	    rta_post("\n\n");
	}
    }

    mif_profile_print(&mif.profile);
    mif_free(&mif);

    rta_free(data);

    return 0;
}
#endif
