/**
@file	mif.c
@author	Diemo Schwarz
@date	30.11.2008
@version 0.1 

@brief	Metric Inverted File Index Structure

Copyright (C) 2008 - 2009 by IRCAM-Centre Georges Pompidou, Paris, France.
*/


#include "rta.h"
#include "mif.h"

/*
 * posting list handling 
 */

void mif_pl_init(mif_postinglist_t *pl, int alloc, int ki)
{
    pl->size     = 0;
    pl->capacity = alloc;

    pl->order    = rta_alloc(alloc * sizeof(*pl->order));
    pl->entryind = rta_alloc(alloc * sizeof(*pl->entryind));
    pl->entries  = rta_alloc(ki    * sizeof(*pl->entries));
}

void mif_pl_free(mif_postinglist_t *pl)
{
    pl->size     = 0;
    pl->capacity = 0;

    rta_free(pl->order);
    rta_free(pl->entryind);
    rta_free(pl->entries);
}



/** initialise index structure */
void mif_init (mif_t *self, mif_distance_function_t distfunc, int nr, int ki)
{
    int i;

    self->dist	 = distfunc;
    self->numobj = 0;
    self->numref = nr;
    self->ki     = ki;
    self->ks     = 0;
    self->refobj       = rta_alloc(numref * sizeof(mif_object_t));
    self->postinglists = rta_alloc(numref * sizeof(mif_postinglist_t));

    /* init posting lists */
    for (i = 0; i < nr; i++)
	mif_pl_init(&self->postinglists [i], ki, ki);

    mif_profile_clear(&self->profile);
}


/** free allocated memory */
void mif_free (mif_t *self)
{
    int i;

    /* free posting lists */
    for (i = 0; i < nr; i++)
	mif_pl_free(&self->postinglists[i]);

    rta_free(self->refobj); 
    rta_free(self->postinglists);
}


/** bulk load new data and index it
    
    @param self		mif structure
    @param numbase	number of base data blocks
    @param base		array[numbase] of base pointer to data
    @param numobj	array[numbase] of number of data objects on this base pointer

    (@return the number of nodes the tree will build)
*/
int mif_add_data (mif_t *self, int numbase, void **base, int *numbaseobj)
{
    int i, k, b;

    {	/*choose reference objects and fill refobj array in mif struct */

	int numobj = 0;
	int *cumobj = alloca(numbase * sizeof(int));
	int *sample = alloca(self->numref);
	
	/* how many elements are given */
	for (i = 0; i < numbase; i++)
	{
	    cumobj[i] = numobj += numbaseobj[i];
	}
	
	/* choose reference objects: draw numref random indices */
	rta_choose_k_from_n(self->numref, numobj, sample);
	
	/* lookup base and relative index */
	for (i = 0; i < self->numref; i++)
	{
	    int objind  = sample[i];
	    int baseind = rta_find_int(objind, numbase, cumobj);
	    int relind  = objind - cumobj[baseind];
	    
	    self->refobj[i].base  = base[baseind];
	    self->refobj[i].index = relind;
	}
    }
     
    {   /* index data */
	int          *indx = alloca(self->ki * sizeof(*indx));	/* ref. obj. index */
	rta_real_t   *dist = alloca(self->ki * sizeof(*dist));	/* distance to obj */
	mif_object_t  newobj;	/* object to index */

	/* for each object */
    for (b = 0; b < numbase; b++)
    {
	newobj.base = base[b];

	for (i = 0; i < numbaseobj[b]; i++)
	{
	    newobj.index = i;
	    	
	    /* find ki closest reference objects for object <base[b], i> and their ordering */
	    for (r = 0; r < self->numref; r++)
	    {
		int kmax = 0;	/* numfound */
		rta_real_t d = self->distance(&self->refobj[r], &newobj);
	    
		if (d <= dist[kmax]) 
		{   /* return original index in data and distance */
		    int pos = kmax;	/* where to insert */
			
		    if (kmax < self->ki - 1)
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
	    } /* now indx is list of ref.obj. ordered by distance dist */

	    /* insert object unsorted into posting lists of ki closest reference objects */
	    for (k = 0; k < self->ki; k++)
	    {
		mif_pl_append(self->pl[indx[k]], &newobj, k);
	    }
	}
    }

    /* sort posting lists and create order-lookup index */
    for (i = 0; i < self->numref; i++)
    {
	qsort(self->pl[i].entries, self->ki, sizeof(*self->pl->entries), comppl);
    }
}

static int compobj (const void *a, const void *b)
{
    return *(int *) a - *(int *) b;
}
