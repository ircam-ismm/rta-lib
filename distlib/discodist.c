
#include <math.h>
#include "mif.h"
#include "discofile.h"
#include "discodist.h"

/* calculate pointer to start of frame vector (skipping time) */
#define DISCO_NDIM(o) ((int *) o->base)[1]
#define DISCO_VECTOR(o) (rta_real_t *) ((float *) (o->base + 3 * sizeof(int)) + 1 + o->index * DISCO_NDIM(o))


/** precalc and preallocate everything */
void disco_KLS_init (void *private, mif_object_t *obj)
{
    disco_KLS_private_t *kls = private;

    kls->N  = (int) ((-1 + sqrt(1 + 8 * (DISCO_NDIM(obj) - 1))) / 4);
    kls->r1 = kls->N;
    kls->r2 = kls->N * kls->N + kls->N;

    kls->D =  (rta_real_t *) rta_malloc(sizeof(rta_real_t) * kls->N);
}

/** precalc and preallocate everything */
void disco_KLS_free (void *private, mif_object_t *base)
{
    free(((disco_KLS_private_t *) private)->D);
}

rta_real_t disco_KLS (void *private, mif_object_t *a, mif_object_t *b)
{
    disco_KLS_private_t *kls = private;
    rta_real_t* v1 = DISCO_VECTOR(a);
    rta_real_t* v2 = DISCO_VECTOR(b);
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

