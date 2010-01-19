

typedef struct _mif_hash_entry_t 
{
    mif_object_t *obj;		/* key */
    int		  accumulator;	/* value */
    struct	 _mif_hash_entry_t *next;
} mif_hash_entry_t;


#if 1 || MIF_USE_SYSTEM_HASH

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

static unsigned int 
hash( const mif_object_t *p)
{
    return (unsigned int) p->base + p->index;
}


static int 
equals( const mif_object_t *p1, const mif_object_t *p2)
{
//    return p1->base == p2->base  &&  p1->index == p2->index;
    return *p1 == *p2;
}


static unsigned int 
new_length( unsigned int length)
{
  unsigned int i;
  
  for ( i = 0; i < sizeof (primes_suite) / sizeof (unsigned int); i++)
    if (length < primes_suite[i])
      return primes_suite[i];
  
  return primes_suite[i-1];
}


static fts_hashtable_cell_t **
lookup_cell( const fts_hashtable_t *h, const mif_object_t *key)
{
  fts_hashtable_cell_t **c = &h->table[hash( key) % h->length];
  
  while(*c && !equals( &(*c)->key, key))
    c = &(*c)->next;
  
  return c;
}


void 
fts_hashtable_init( fts_hashtable_t *h, int initial_capacity)
{
    h->length = new_length(initial_capacity);
    h->count = 0;
    h->rehash_count = (int)(h->length * FTS_HASHTABLE_STANDARD_LOAD_FACTOR);

    h->table = (fts_hashtable_cell_t **) rta_calloc( h->length * sizeof( fts_hashtable_cell_t *));
}


void 
fts_hashtable_clear( fts_hashtable_t *h)
{
  unsigned int i;
  
  for ( i = 0; i < h->length; i++)
  {
	fts_hashtable_cell_t *c, *next;
	
	for( c = h->table[i]; c ; c = next)
	{
	  next = c->next;
	  fts_heap_free( c, cell_heap);
	}
	
	h->table[i] = 0;
  }
  
  h->count = 0;
}



#endif
