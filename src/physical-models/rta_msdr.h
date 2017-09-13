/**
 * @file rta_msdr.h
 * @author Norbert Schnell
 * @ingroup rta_physical_models
 *
 * @brief Mass Spring Damper Repulsion model
 *
 * @copyright
 * Copyright (C) 2008 - 2009 by IRCAM-Centre Georges Pompidou, Paris, France.
 * All rights reserved.
 *
 * License (BSD 3-clause)
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

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
#define RTA_MSDR_MAXCAT	3


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
    int ind1; /** index in link list of mass 1 UNUSED??? */
    int ind2; /** index in link list of mass 2 UNUSED??? */
    int cat;  /** category of link list of masses */
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
    int nlinks [RTA_MSDR_MAXCAT];		/** number of links in system */
    int nactive [RTA_MSDR_MAXCAT];		/** number of links active */
    int linksalloc;				/** allocated number of links in system */
    rta_msdr_link_t *links [RTA_MSDR_MAXCAT];

    float	*K1 [RTA_MSDR_MAXCAT];		/** rigidity(nlink) */
    float	*D1 [RTA_MSDR_MAXCAT];		/** damping1(nlink) */
    float	*D2 [RTA_MSDR_MAXCAT];		/** damping2(nlink) friction? */
    float	*Rt [RTA_MSDR_MAXCAT];		/** repulsion length(nlink) */
    float	*Rf [RTA_MSDR_MAXCAT];		/** repulsion force(nlink) */

    float	*l0 [RTA_MSDR_MAXCAT];		/** nominal length(nlink) */
    float	*lcurr [RTA_MSDR_MAXCAT];	/** current length(nlink) */
    float	*lprev [RTA_MSDR_MAXCAT];	/** previous length(nlink) */
    float       *stress [RTA_MSDR_MAXCAT];     	/** link stress(nlink, 1) */
    float       *forceabs [RTA_MSDR_MAXCAT];   	/** link force(nlink, 1) */
} rta_msdr_t;


/*
 *  links
 */

int rta_msdr_add_link (rta_msdr_t *sys, int m1, int m2, float len, int cat,
		       float K1, float D1, float D2, float Rt, float Rf);

int rta_msdr_insert_link (rta_msdr_t *sys, int m1, int m2, float len, int cat,
			  float K1, float D1, float D2, float Rt, float Rf, int max);

void rta_msdr_clear_links (rta_msdr_t *sys);


/*
 *  masses
 */

/* set mass data (pos and masses) and reset to initial state */
void rta_msdr_set (rta_msdr_t *sys, int nmasses, float *pos, float *invmass);

/* copy masses pos only, n rows from index i  */
void rta_msdr_set_pos (rta_msdr_t *sys, int i, int n, float *pos);

/* copy inv. mass only, n rows from index i */
void rta_msdr_set_mass (rta_msdr_t *sys, int i, int n, float *invmass);


/*
 *  system update
 */

/* update links and mass positions
   return total stress */
float rta_msdr_update (rta_msdr_t *sys);

float rta_msdr_update_limp (rta_msdr_t *sys);

/* update links and mass positions
   return total stress */
float rta_msdr_update_ind (rta_msdr_t *sys, int nind, int *ind);

/* update links and mass positions, no inertia
   return total stress */
float rta_msdr_update_limp_ind (rta_msdr_t *sys, int nind, int *ind);



void rta_msdr_init (rta_msdr_t *sys, int maxmass, int maxlinkstotal, int maxlinkscat);

void rta_msdr_free (rta_msdr_t *sys);


/* set vector to recieve current link forces */
void rta_msdr_set_outforce (rta_msdr_t *sys, float *outforce);

/* set rigidity parameter for all links or in category */
void rta_msdr_set_K1 (rta_msdr_t *sys, int cat, float k1);

/* set damping parameter for all links */
void rta_msdr_set_D1 (rta_msdr_t *sys, int cat, float d1);

/* set friction parameter for all links */
void rta_msdr_set_D2 (rta_msdr_t *sys, int cat, float d2);

/* set repulsion threshold parameter for all links */
void rta_msdr_set_Rt (rta_msdr_t *sys, int cat, float rt);

/* set repulsion force parameter for all links */
void rta_msdr_set_Rf (rta_msdr_t *sys, int cat, float rf);

/* set length by index */
void rta_msdr_set_link_length (rta_msdr_t *sys, int i, int c, float L);
void rta_msdr_set_link_K1 (rta_msdr_t *sys, int i, int c, float K1);
void rta_msdr_set_link_D1 (rta_msdr_t *sys, int i, int c, float D1);
void rta_msdr_set_link_D2 (rta_msdr_t *sys, int i, int c, float D2);
void rta_msdr_set_link_Rt (rta_msdr_t *sys, int i, int c, float rt);
void rta_msdr_set_link_Rf (rta_msdr_t *sys, int i, int c, float rf);


/* get number of masses in system */
int rta_msdr_get_num_masses (rta_msdr_t *sys);

/* get number of links in system */
int rta_msdr_get_num_links (rta_msdr_t *sys);
int rta_msdr_get_num_active_links (rta_msdr_t *sys);

int rta_msdr_get_num_links_cat (rta_msdr_t *sys, int c);

int rta_msdr_get_mass_num_links(rta_msdr_t *sys, int massi, int cat);

void rta_msdr_clear_cat_links (rta_msdr_t *sys, int cat);


/* copy ncol columns if link data to out(nlinks, 8):
   masses id(2), masses pos (4), stress, force */
int rta_msdr_get_links (rta_msdr_t *sys, float *out);

/* get max link distance of category */
float rta_msdr_get_mass_maxdist(rta_msdr_t *sys, int massi, int cat);

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
