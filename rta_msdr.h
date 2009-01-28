
/* Mass Spring Damper Repulsion model */


#ifndef _RTA_MSDR_H_
#define _RTA_MSDR_H_

#include "rta.h"

#ifdef __cplusplus
extern "C" {
#endif

/** number of dimensions D, redefine before including*/
#ifndef RTA_MSDR_NDIM
#define RTA_MSDR_NDIM	   2
#define RTA_MSDR_NDIM_STR "2"
#endif


/** max. number of link categories */
#define RTA_MSDR_MAXCAT	2


/** struct for bounding values */
typedef struct _rta_msdr_limits {
    int	  unlimited;		/** flag to test for limits at all */
    float min [RTA_MSDR_NDIM];  /** lower limit for each dimension */
    float max [RTA_MSDR_NDIM];  /** upper limit for each dimension */
} rta_msdr_limits_t;


/** struct representing one mass
    N.B.: position, force, velocity stored as vectors in msdr struct */
typedef struct _rta_msdr_mass {
    int   nlinks [RTA_MSDR_MAXCAT]; /** number of links per link category */
    int  *links  [RTA_MSDR_MAXCAT];	/** indices of links */
    float maxdist[RTA_MSDR_MAXCAT];  /** distance of farthest neighbour */
} rta_msdr_mass_t;


/** struct representing one link between two masses
    N.B.: allmost all parameters are stored as vectors in msdr struct */
typedef struct _rta_msdr_link {
    int m1; /** index of mass 1 */
    int m2; /** index of mass 2 */
} rta_msdr_link_t;


typedef struct _rta_msdr {
    /** masses */
    int nmasses;		/** number of masses in system */
    int massalloc;		/** allocated number of masses in system */
    int linklistalloc;
    rta_msdr_mass_t *masses;	/** non-vector mass data */

    /** vectors of mass data */
    float *pos;			/** (nmass, D) position */
    float *invmass;		/** (nmass, 1) inverse mass, 0 = fixed */
    float *force;		/** (nmass, D) force vectors */
    float *speed;		/** (nmass, D) speed vectors */
    float *pos2;		/** (nmass, D) old position */
    float *outforce;		/** (nmass, D) force output */

    rta_msdr_limits_t poslim;	/** position limits for masses */
    rta_msdr_limits_t speedlim;	/** speed limits for masses */

    /** links */
    int nlinks;			/** number of links in system */
    int linksalloc;		/** allocated number of links in system */
    rta_msdr_link_t *links;

    float	*K1;		/** rigidity(nlink) */
    float	*D1;		/** damping1(nlink) */
    float	*D2;		/** damping2(nlink) friction? */
    float	*Rt;		/** repulsion length(nlink) */
    float	*Rf;		/** repulsion force(nlink) */

    float	*l0;		/** nominal length(nlink) */
    float	*lcurr;		/** current length(nlink) */
    float	*lprev;		/** previous length(nlink) */
    float       *stress;     	/** link stress(nlink, 1) */
    float       *forceabs;     	/** link force(nlink, 1) */
} rta_msdr_t;


int rta_msdr_add_link (rta_msdr_t *sys, int m1, int m2, float len, int cat,
		       float K1, float D1, float D2, float Rt, float Rf);

void rta_msdr_clear_links (rta_msdr_t *sys);

/* update links and mass positions
   return total stress */
float rta_msdr_update (rta_msdr_t *sys);

/* update links and mass positions
   return total stress */
float rta_msdr_update_ind (rta_msdr_t *sys, int nind, int *ind);

/* update links and mass positions, no inertia
   return total stress */
float rta_msdr_update_limp_ind (rta_msdr_t *sys, int nind, int *ind);

/** update nind masses in index list ind */
void  rta_msdr_update_masses_ind (rta_msdr_t *sys, int nind, int *ind);

/** compute damping and friction forces */
float rta_msdr_update_links_damping (rta_msdr_t *sys);

/** compute link elasticity forces */
float rta_msdr_update_links_elasticity (rta_msdr_t *sys);

void rta_msdr_init (rta_msdr_t *sys, int maxmass, int maxlinkstotal, int maxlinkscat);

void rta_msdr_free (rta_msdr_t *sys);

/* set mass data (pos and masses) and reset to initial state */
void rta_msdr_set (rta_msdr_t *sys, int nmasses, float *pos, float *invmass);

/* copy masses pos only, n rows from index i  */
void rta_msdr_set_pos (rta_msdr_t *sys, int i, int n, float *pos);

/* copy inv. mass only, n rows from index i */
void rta_msdr_set_mass (rta_msdr_t *sys, int i, int n, float *invmass);


/* set vector to recieve current link forces */
void rta_msdr_set_outforce (rta_msdr_t *sys, float *outforce);

/* set rigidity parameter for all links */
void rta_msdr_set_K1 (rta_msdr_t *sys, float k1);

/* set damping parameter for all links */
void rta_msdr_set_D1 (rta_msdr_t *sys, float d1);

/* set friction parameter for all links */
void rta_msdr_set_D2 (rta_msdr_t *sys, float d2);

/* set repulsion threshold parameter for all links */
void rta_msdr_set_Rt (rta_msdr_t *sys, float rt);

/* set repulsion force parameter for all links */
void rta_msdr_set_Rf (rta_msdr_t *sys, float rf);

/* get number of masses in system */
int rta_msdr_get_num_masses (rta_msdr_t *sys);

/* get number of links in system */
int rta_msdr_get_num_links (rta_msdr_t *sys);

/* copy ncol columns if link data to out(nlinks, 8):
   masses id(2), masses pos (4), stress, force */
int rta_msdr_get_links (rta_msdr_t *sys, float *out);

/* get total movement after last masses update */
float rta_msdr_get_movement (rta_msdr_t *sys);

/* copy force vector (after links update) */
float rta_msdr_get_force (rta_msdr_t *sys);


/** set masses position limits if id == 0, speed limits if id == 1 */
void rta_msdr_set_limits (rta_msdr_t *sys, int id, float *min, float *max);

/** set to no limits */
void rta_msdr_set_unlimited (rta_msdr_t *sys, int id);


#ifdef __cplusplus
}
#endif

#endif /* _RTA_MSDR_H_ */
