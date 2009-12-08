
/* Hybrid multi-dimensional scaling (Chalmers, Ross, Morrison algorithm) 
   
   Based on a mass-string-damper model (rta_msdr_t),
   pre-place a random sample of points,
   then run iterations until stress < threshold or maxiter reached:
   - for each point, find ns random points,
     - insert into neighbour list if close, into random list otherwise
   - run one iteration of msdr model
 */

#include "rta_mds_hybrid.h"

/* local abbreviation for number of dimensions */
#define NDIM	RTA_MDS_NDIM

/* index to i'th row of NDIM matrix */
#define ROWIND(i) ((i) * NDIM * sizeof(float))

/* minimum distance to avoid numeric overflow */
#define EPS 1e-9	


/* place a sample of ns points randomly 
float rta_mds_hybrid_preplace_random (rta_mds_hybrid_t *sys, int ns)
{
}
*/

//    rta_mds_hybrid_preplace_random(sys, nsamp, sampind, ndata, ndim, data, weights);
    //- place data[sampind] into layout space
    //- add full links with length = distance_weighted(p1, p2, weights)
    //- run system until delta(stress) < thresh or i > maxiter:
    //  - stress = rta_msdr_update(sys);

static float *vectors_alloc (int n)
{
    int    size = n * NDIM * sizeof(float);
    float *vect = (float *) rta_malloc(size);
    bzero(vect, size);
    return vect;
}

#define vect_alloc(n)    ((float *) rta_malloc((n) * sizeof(float)))
