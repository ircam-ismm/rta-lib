
#include <math.h>
#include "mif.h"
#include "discofile.h"
#include "discodist.h"


/* rta_real_t *DISCO_VECTOR(mif_files_t *db, mif_object_t *o) 
   calculate pointer to start of frame vector (skipping time) */
#define DISCO_VECTOR(db, o) (rta_real_t *) (disco_file_data((disco_file_header_t *) db->base[o->base]) \
					    + 1 /* skip first element in feature vector (contains time) */ \
					    + o->index * db->ndim)

/** precalc and preallocate everything */
void disco_KLS_init (void *private, mif_files_t *db)
{
    disco_KLS_private_t *kls = private;

    kls->db = db;
    kls->N  = (int) ((-1 + sqrt(1 + 8 * (db->ndim - 1))) / 4);
    kls->r1 = kls->N;
    kls->r2 = kls->N * kls->N + kls->N;

    kls->D =  (rta_real_t *) rta_malloc(sizeof(rta_real_t) * kls->N);
}

/** precalc and preallocate everything */
void disco_KLS_free (void *private, mif_files_t *db)
{
    free(((disco_KLS_private_t *) private)->D);
}

rta_real_t disco_KLS (void *private, mif_object_t *a, mif_object_t *b)
{
    disco_KLS_private_t *kls = private;
    rta_real_t* v1 = DISCO_VECTOR(kls->db, a);
    rta_real_t* v2 = DISCO_VECTOR(kls->db, b);
    rta_real_t dist = 0;
    int i,j;

/*  for (i=0; i < 10; i++)
       rta_post("%f ", v1[i]);
    rta_post("\n");
*/
    double T1=0;
    double T2=0;
    double S1=0;
    double S2=0;
    double tmp1=0;
    double tmp2=0;

    for (i=0 ; i < kls->N; i++)
    {
	kls->D[i]=v1[i]-v2[i];
    }

    for (i=0 ; i < kls->N; i++)
    {
	tmp1=0;
	tmp2=0;

	for (j=0; j < kls->N; j++)
	{
	    T1 += v2[j + i*kls->N + kls->r2]*v1[i + j*kls->N + kls->r1];
	    T2 += v2[j + i*kls->N + kls->r1]*v1[i + j*kls->N + kls->r2];

	    tmp1 += v2[i + j*kls->N + kls->r2]*kls->D[j];
	    tmp2 += v1[i + j*kls->N + kls->r2]*kls->D[j];
	}

	S1 += tmp1*kls->D[i]-1;
	S2 += tmp2*kls->D[i]-1;
    }

    dist=(S1 + S2 + T1 + T2)/4 ;

    if (dist <0)
	dist=0;

    return dist;
}

