
#include "rta.h"
#include "mif.h"

typedef struct _disco_KLS_private
{
    int		N;	/**< number of dimensions */
    int		r1, r2;	/**< temporary vars */
    rta_real_t *D;	/**< space for vector distance */
} disco_KLS_private_t;


/** precalc and preallocate everything */
void disco_KLS_init (void *private, mif_object_t *obj);

/** precalc and preallocate everything */
void disco_KLS_free (void *private, mif_object_t *base);

rta_real_t disco_KLS (void *private, mif_object_t *a, mif_object_t *b);
