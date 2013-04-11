
/* Hybrid multi-dimensional scaling 
   (Chalmers, Ross, Morrison algorithm) 
*/

#ifndef _RTA_MDS_HYBRID_H_
#define _RTA_MDS_HYBRID_H_

#include "rta.h"
#include "rta_msdr.h"

#ifdef __cplusplus
extern "C" {
#endif


/** number of dimensions D, redefine before including*/
#ifndef RTA_MDS_NDIM
#define RTA_MDS_NDIM      2
#define RTA_MDS_NDIM_STR "2"
#endif


rta_real_t weighted_euclidean_distance (rta_real_t* v1, rta_real_t* v2, rta_real_t *sigma, int ndim);


#ifdef __cplusplus
}
#endif

#endif /* _RTA_MDS_HYBRID_H_ */
