
#ifndef _MIF_HASH_H_
#define _MIF_HASH_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <assert.h>

#include "mif.h"
#include "rta.h"

typedef struct _mif_hash_entry_t 
{
    const mif_object_t *obj;		/* key */
    int		  accumulator;	/* value */
    unsigned int  next;
    // struct	 _mif_hash_entry_t *next;
} mif_hash_entry_t;


#if MIF_USE_SYSTEM_HASH

#include <search.h>

/* hash */
#define MAXLEN   (sizeof(void *) * 2 + 4 + 1) /* 32 bit pointer + 16 bit index in hex */

/* for profile output */
#define HASHOBJSIZE sizeof(mif_hash_entry_t) + MAXLEN


static void mif_object_hash_create(mif_index_t *self)
{
    hcreate(self->numobj);	/* todo: ki * ks */
}

static void mif_object_hash_destroy(mif_index_t *self)
{
    hdestroy();
}

static void mif_object_hash_makekey(mif_object_t *obj, char *key)
{
    int len = sprintf(key, "%lx%x", (long unsigned int) obj->base, obj->index);
    assert(len < MAXLEN - 1);
}

static int mif_object_hash_new(char *key, int entry)
{
    ENTRY item, *ret;

    item.key  = strdup(key);
    item.data = (void *) entry;

    ret = hsearch(item, ENTER);
    return (ret  ?  (int) ret->data  :  -1);
}

static int mif_object_hash_get(char *key)
{
    ENTRY item, *ret;

    item.key  = key;
    ret = hsearch(item, FIND);
    return (ret  ?  (int) ret->data  :  -1);
}

#else

/* from ftmlib/hashtable.c */

/**
 * The FTS hashtable.
 *
 * An FTS hashtable maps keys to values, where both keys and values are FTS atoms.
 *
 * FTS hashtable performs automatic rehashing when the ration between number of inserted keys and the
 * capacity is greater than a "load factor" (typically 0.75).
 *
 * The initial capacity of the hashtable can be given by an argument at the initialisation
 * (FTS_HASHTABLE_SMALL, FTS_HASHTABLE_MEDIUM or FTS_HASHTABLE_BIG).
 * To work with a big hashtable, it is more efficient to create the hashtable with a big capacity 
 * to avoid intermediate automatic rehashing.
 *
 * @defgroup fts_hashtable FTS hashtable
 * @ingroup fts_struct
 */

/**
 * @brief 
 * FTS hashtable data type.
 *
 * @ingroup fts_hashtable
 */
typedef struct fts_hashtable fts_hashtable_t;

typedef struct _mif_hash_entry_t fts_hashtable_cell_t;

struct fts_hashtable
{
    unsigned int length;	/* table size */
    unsigned int alloc;		/* cell_heap size */
    unsigned int incr;		/* cell_heap alloc increment */
    int count;
    int rehash_count;
    unsigned int *table;	/* length indices into cell_heap */
    //fts_hashtable_cell_t **table;
    fts_hashtable_cell_t *cell_heap;
};

#define FTS_HASHTABLE_STANDARD_LOAD_FACTOR  0.75

/* for profile output */
#define HASHOBJSIZE sizeof(mif_hash_entry_t)


static const unsigned int primes_suite[] = {
  7,
  17,
  31,
  67,
  127,
  257,
  521,
  1031,
  2053,
  4099,
  8191,
  16411,
  32771,
  65537,
  131071,
  262147,
  524287,
  1048583,
  2097169,
  4194319,
  8388617,
  16777259,
  33554467,
};

static unsigned int hash (const mif_object_t *p)
{
    /* fine with pointers */
    return (unsigned int) p->base + p->index;
}


static int equals (const mif_object_t *p1, const mif_object_t *p2)
{
    return p1->base == p2->base  &&  p1->index == p2->index;
}


static unsigned int new_length (unsigned int length)
{
  unsigned int i;
  
  for ( i = 0; i < sizeof (primes_suite) / sizeof (unsigned int); i++)
    if (length < primes_suite[i])
      return primes_suite[i];
  
  return primes_suite[i-1];
}


void fts_hashtable_init (fts_hashtable_t *h, int initial_capacity, int alloc_increment)
{
    h->length = new_length(initial_capacity);
    h->alloc  = initial_capacity;
    h->incr   = alloc_increment;
    h->count  = 0;
    h->rehash_count = (int)(h->length * FTS_HASHTABLE_STANDARD_LOAD_FACTOR);

    h->table     = (unsigned int *) rta_zalloc((h->length + 1) * sizeof(unsigned int));
    //h->table     = (fts_hashtable_cell_t **) rta_calloc(h->length * sizeof(fts_hashtable_cell_t *));
    h->cell_heap = (fts_hashtable_cell_t *)  rta_malloc(h->alloc  * sizeof(fts_hashtable_cell_t));
}

void fts_hashtable_clear (fts_hashtable_t *h)
{
  unsigned int i;
  
  for (i = 0; i <= h->length; i++)
      h->table[i] = 0;
  
  h->count = 0;
}


void fts_hashtable_free (fts_hashtable_t *h)
{
  fts_hashtable_clear( h);
  rta_free(h->cell_heap);
}


// static fts_hashtable_cell_t **
static unsigned int *lookup_cell (const fts_hashtable_t *h, const mif_object_t *key)
{
    unsigned int *c = &h->table[hash(key) % h->length + 1];
    // fts_hashtable_cell_t **c = &h->table[hash( key) % h->length]; 
  
    while (*c  && !equals( h->cell_heap[*c].obj, key))
	c = &h->cell_heap[*c].next;

    // while (*c && !equals( &(*c)->obj, key))
        // c = &(*c)->next;
  
  return c;
}


int fts_hashtable_get (const fts_hashtable_t *h, const mif_object_t *key, /*out*/ int **value)
{
    unsigned int *c = lookup_cell(h, key);
    // fts_hashtable_cell_t **c = lookup_cell( h, key);
  
    if (*c > 0)
       //if(*c != NULL)
    {
	*value = &h->cell_heap[*c].accumulator;
	// *value = &(*c)->accumulator;
	return 1;
    }
  
    return 0;
}


int fts_hashtable_put (fts_hashtable_t *h, const mif_object_t *key, /*inout*/ int **value)
{
    unsigned int *c = lookup_cell(h, key);
    fts_hashtable_cell_t *cc;
    //fts_hashtable_cell_t **c = lookup_cell( h, key);
  
    if (*c > 0)
       //if(*c != NULL)
    {
	*value = &h->cell_heap[*c].accumulator;
	// *value = &(*c)->accumulator;
	return 1;
    } /* else: new object */
  
    if (h->count + 1 >= h->alloc)
    {
	h->alloc += h->incr;
	rta_post("reallocate hash object store from %d to %d obj (size %lu)\n", 
		 h->alloc - h->incr, h->alloc, h->alloc * sizeof(fts_hashtable_cell_t));
	h->cell_heap = rta_realloc(h->cell_heap, h->alloc * sizeof(fts_hashtable_cell_t));
	assert(h->cell_heap);
    }

    /* append after existing cell */
    *c = h->count + 1;
    cc = &h->cell_heap[*c];
    // *c = &h->cell_heap[h->count];
    (cc)->obj = key;
    // (*c)->accumulator = *value;
    *value = &(cc)->accumulator;
    (cc)->next = c; //NULL;
    h->count++;
    
    if ( h->count == h->rehash_count) /* usually: >= */
    {
	rta_post("hashtable (count %d) wants to rehash table of length %d\n", h->count, h->length);
	// rehash( h);
    }  
    return 0;
}


#if 0
static int hashtable_iterator_has_more( fts_iterator_t *iter)
{
  fts_hashtable_iterator_t *i = (fts_hashtable_iterator_t *) iter->data;
  
  if (i->cell)
    return 1;
  
  while ( i->index-- )
  {
    if ( i->table[i->index] )
    {
      i->cell = i->table[i->index];
      return 1;
    }
  }
  
  fts_heap_free( iter->data, iterator_heap);
  
  return 0;
}

static void hashtable_iterator_next( fts_iterator_t *iter, fts_atom_t *a)
{
  fts_hashtable_iterator_t *i = (fts_hashtable_iterator_t *) iter->data;
  
  if ( !i->cell)
  {
    while( i->index--)
    {
      if (i->table[i->index])
      {
	i->cell = i->table[i->index];
	break;
      }
    }
  }
  
  if (!i->cell)
    return;
  
  *a = (i->keys) ? i->cell->key : i->cell->value;
  
  i->cell = i->cell->next;
}

static void hashtable_iterator_get( const fts_hashtable_t *h, fts_iterator_t *i, int keys)
{
  fts_hashtable_iterator_t *hiter;

  hiter = (fts_hashtable_iterator_t *)fts_heap_alloc( iterator_heap);

  hiter->table = h->table;
  hiter->keys = keys;
  hiter->cell = NULL;
  hiter->index = h->length;

  i->has_more = hashtable_iterator_has_more;
  i->next = hashtable_iterator_next;
  i->data = hiter;
}
#endif

#endif

#ifdef __cplusplus
}
#endif

#endif /* _MIF_HASH_H_ */
