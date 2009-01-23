

#include "rta_msdr.h"

/* local abbreviation for number of dimensions */
#define NDIM	RTA_MSDR_NDIM

/* index to i'th row of NDIM matrix */
#define ROWIND(i) ((i) * NDIM * sizeof(float))

/* minimum distance to avoid numeric overflow */
#define EPS 1e-9	


/*
 * physics engine
 */


/* return value x clipped to bounds for dimension d */
static float limit (rta_msdr_limits_t *lim, int d, float x)
{
    if (x < lim->min[d])
	return lim->min[d];
    else if (x > lim->max[d])
	return lim->max[d];
    else
	return x;
}


/** euclidean distance between two vectors */
static float distance (float *x1, float *x2)
{
    double dist = 0;
    int   d;

    for (d = 0; d < NDIM; d++)
    {
	double distorth = x1[d] - x2[d];
	dist += distorth * distorth;
    }

    return (float) sqrt(dist);
}


// #define SIGN(x)	((x) >= 0 ? 1 : -1)
#define SIGN(x)	((x) >= 0 ? 1 : -1)

/* compute link forces */
static float update_links (rta_msdr_t *sys)
{
    float totalstress = 0;
    int i;

    for (i = 0; i < sys->nlinks; i++)
    {
	int    m1       = sys->links[i].m1;
	int    m2       = sys->links[i].m2;
	float *m1ptr    = sys->pos + m1 * NDIM;
	float *m2ptr    = sys->pos + m2 * NDIM;

	/* calculate absolute force */
	float dist     = distance(m1ptr, m2ptr);
	int d;

	sys->stress[i] = dist - sys->l0[i];

	/* if using pow(stress): if (stress > 0) */
	if (dist > 0)
	{
	    float viscodamping = sys->D1[i] * (dist - sys->lprev[i]);
	    sys->forceabs[i] = (sys->K1[i] * sys->stress[i] + viscodamping) / dist;
	}
	else
	    sys->forceabs[i] = 0;

	/* add repulsion term if close enough */
	if (dist < sys->Rt[i])
	    sys->forceabs[i] -= (1 - dist / sys->Rt[i]) * sys->Rf[i];

	/* apply link force to mass force vector */
	for (d = 0; d < NDIM; d++)
	{
	    float forcenorm = sys->forceabs[i] * (m1ptr[d] - m2ptr[d]);
	
	    /* apply force normal vector to both masses, add friction */
	    sys->force[m1 * NDIM + d] 
		-= forcenorm + sys->D2[i] * sys->speed[m1 * NDIM + d];
	    sys->force[m2 * NDIM + d] 
		+= forcenorm - sys->D2[i] * sys->speed[m2 * NDIM + d];
	}

/*	fts_post("link %2d: dist %f len %f -> stress %f forceabs %f\n", 
		 i, dist, sys->l0[i], sys->stress[i], sys->forceabs[i]);
*/
	sys->lprev[i] = dist;
	totalstress += fabs(sys->stress[i]);
    }

    return totalstress;
}


/* compute damping and friction forces */
static float update_links_damping (rta_msdr_t *sys)
{
    int i;

    for (i = 0; i < sys->nlinks; i++)
    {
	int    m1       = sys->links[i].m1;
	int    m2       = sys->links[i].m2;
	float *m1ptr    = sys->pos + m1 * NDIM;
	float *m2ptr    = sys->pos + m2 * NDIM;

	/* calculate absolute force */
	float dist     = distance(m1ptr, m2ptr);
	int d;

	if (dist > 0)  /* viscosity damping */
	    sys->forceabs[i] = sys->D1[i] * (dist - sys->lprev[i]) / dist;
	else
	    sys->forceabs[i] = 0;

	/* apply link force to mass force vector */
	for (d = 0; d < NDIM; d++)
	{
	    float forcenorm = sys->forceabs[i] * (m1ptr[d] - m2ptr[d]);
	
	    /* apply force normal vector to both masses, add friction */
	    sys->force[m1 * NDIM + d] 
		-= forcenorm + sys->D2[i] * sys->speed[m1 * NDIM + d];
	    sys->force[m2 * NDIM + d] 
		+= forcenorm - sys->D2[i] * sys->speed[m2 * NDIM + d];
	}
    }
}

/* compute link elasticity forces */
static float update_links_elasticity (rta_msdr_t *sys)
{
    float totalstress = 0;
    int i;

    for (i = 0; i < sys->nlinks; i++)
    {
	int    m1       = sys->links[i].m1;
	int    m2       = sys->links[i].m2;
	float *m1ptr    = sys->pos + m1 * NDIM;
	float *m2ptr    = sys->pos + m2 * NDIM;

	/* calculate absolute force */
	float dist     = distance(m1ptr, m2ptr);
	int d;

	sys->stress[i] = dist - sys->l0[i];

	if (dist > 0)
	    sys->forceabs[i] = sys->K1[i] * sys->stress[i] / dist;
	else
	    sys->forceabs[i] = 0;

	/* add repulsion term if close enough */
	if (dist < sys->Rt[i])
	    sys->forceabs[i] -= (1 - dist / sys->Rt[i]) * sys->Rf[i];

	/* apply link force to mass force vector */
	for (d = 0; d < NDIM; d++)
	{
	    float forcenorm = sys->forceabs[i] * (m1ptr[d] - m2ptr[d]);
	
	    /* apply force normal vector to both masses */
	    sys->force[m1 * NDIM + d] -= forcenorm;
	    sys->force[m2 * NDIM + d] += forcenorm;
	}

	sys->lprev[i] = dist;
	totalstress += fabs(sys->stress[i]);
    }

    return totalstress;
}


/* update masses from link forces */
static float update_masses (rta_msdr_t *sys, float *totalmovement)
{
    int i, d;

    for (i = 0; i < sys->nmasses; i++)
    {
	int    massidx  = i * NDIM;

	for (d = 0; d < NDIM; d++)
	{
	    int   dd   = massidx + d;
	    float pold = sys->pos[dd];

	    /* apply force and current speed */
	    float pnew = sys->force[dd] * sys->invmass[i] + 2 * pold - sys->pos2[dd];
	    /* set and clip speed, recalc new position */
	    sys->speed[dd] = limit(&sys->speedlim, d, pnew - pold);
	    pnew           = limit(&sys->poslim,   d, pold + sys->speed[dd]);

	    /* cycle state: clear applied force */
	    sys->pos   [dd] = pnew;
	    sys->pos2  [dd] = pold;
	    if (sys->outforce)
		sys->outforce[dd] = sys->force[dd];
	    sys->force [dd] = 0;
	}
    }

    if (totalmovement != NULL)
    {	/* if we want to know total movement, run expensive sqr/sqrt
	   over all lines again */
	int dd = 0;
	*totalmovement = 0;

	for (i = 0; i < sys->nmasses; i++)
	{
	    float  massmove = 0;

	    for (d = 0; d < NDIM; d++, dd++)
		massmove += sys->speed[dd] * sys->speed[dd];

	    *totalmovement += sqrt(massmove);
	}

	return *totalmovement;
    }
    else
	return 0;
}

/* update links and mass positions
   return total stress, on demand return total movement magnitude */
float rta_msdr_update (rta_msdr_t *sys, float *move)
{
    float stress, move1;

/* single-step:
    stress = update_links(sys);
    update_masses(sys, move);
*/
    /* two-step: */
    update_links_damping(sys);
    update_masses(sys, &move1);
    stress = update_links_elasticity(sys);
    update_masses(sys, move);

    fts_post("update: damp move %f  elastic move %f  stress %f\n",
	     move1, move ? *move : 0, stress);

    return stress;
}


/*
 *  querying the system
 */

int rta_msdr_get_num_masses (rta_msdr_t *sys)
{
    return sys->nmasses;
}

int rta_msdr_get_num_links (rta_msdr_t *sys)
{
    return sys->nlinks;
}

/* set vector to recieve current link forces */
void rta_msdr_set_outforce (rta_msdr_t *sys, float *outforce)
{
    sys->outforce = outforce;
}

/* copy ncol columns if link data to out(nlinks, 8):
   masses id(2), masses pos (4), stress, force */
int rta_msdr_get_links (rta_msdr_t *sys, float *out)
{
    int i;

    for (i = 0; i < sys->nlinks; i++)
    {
	rta_msdr_link_t *link = &sys->links[i];
	out[0] = link->m1;
	out[1] = link->m2;
	out[2] = sys->pos[link->m1 * NDIM];
	out[3] = sys->pos[link->m1 * NDIM + 1];
	out[4] = sys->pos[link->m2 * NDIM];
	out[5] = sys->pos[link->m2 * NDIM + 1];
	out[6] = sys->stress[i];
	out[7] = sys->forceabs[i];
	out += 8;
    }
}


/*
 *  constructing the system
 */

/* set masses position limits */
void rta_msdr_set_limits (rta_msdr_t *sys, int id, float *min, float *max)
{
    rta_msdr_limits_t *lim = (id == 1  ?  &sys->speedlim  :  &sys->poslim);
    int i;

    for (i = 0; i < NDIM; i++)
    {
	lim->min[i] = min[i];
	lim->max[i] = max[i];
    }
}


/** add link between masses m1, m2
    m1, m2 must be valid mass ids, cat must be valid link category */
int rta_msdr_add_link (rta_msdr_t *sys, int m1, int m2, float len, int cat,
		   float K1, float D1, float D2, float Rt, float Rf)
{
    if (sys->nlinks < sys->linksalloc  &&  m1 != m2)
    {
	rta_msdr_mass_t *m1ptr, *m2ptr;
	int i = sys->nlinks++;

	sys->links[i].m1 = m1;
	sys->links[i].m2 = m2;
	sys->l0[i] = len;	/* nominal length */
	sys->lprev[i] = len;

	sys->K1[i] = K1;
	sys->D1[i] = D1;
	sys->D2[i] = D2;
	sys->Rt[i] = Rt;
	sys->Rf[i] = Rf;

	/* record linkage in both masses' links lists */
	m1ptr = &sys->masses[m1];
	m2ptr = &sys->masses[m2];

	/* if (nlinks[cat] < maxlinks[cat]) */
	/* todo: insert into sorted list; if already present: overwrite link */
	m1ptr->links[cat][m1ptr->nlinks[cat]++] = i;
	m2ptr->links[cat][m2ptr->nlinks[cat]++] = i;

	return i;
    }
    else
	return -1;
}

/* set rigidity parameter for all links */
void rta_msdr_set_K1 (rta_msdr_t *sys, float k1)
{
    int i;

    for (i = 0; i < sys->nlinks; i++)
    {
	sys->K1[i] = k1;
    }
}

/* set damping parameter for all links */
void rta_msdr_set_D1 (rta_msdr_t *sys, float d1)
{
    int i;

    for (i = 0; i < sys->nlinks; i++)
    {
	sys->D1[i] = d1;
    }
}

/* set friction parameter for all links */
void rta_msdr_set_D2 (rta_msdr_t *sys, float d2)
{
    int i;

    for (i = 0; i < sys->nlinks; i++)
    {
	sys->D2[i] = d2;
    }
}

/* set repulsion threshold parameter for all links */
void rta_msdr_set_Rt (rta_msdr_t *sys, float rt)
{
    int i;

    for (i = 0; i < sys->nlinks; i++)
    {
	sys->Rt[i] = rt;
    }
}

/* set repulsion force parameter for all links */
void rta_msdr_set_Rf (rta_msdr_t *sys, float rf)
{
    int i;

    for (i = 0; i < sys->nlinks; i++)
    {
	sys->Rf[i] = rf;
    }
}


/* copy masses positions */
void rta_msdr_set_pos (rta_msdr_t *sys, int ind, int n, float *pos)
{   /* reset also old pos and speed, leave forces! */
    memcpy(sys->pos   + ind * NDIM, pos, n * NDIM * sizeof(float));
    memcpy(sys->pos2  + ind * NDIM, pos, n * NDIM * sizeof(float));
    bzero (sys->speed + ind * NDIM,      n * NDIM * sizeof(float));
}

/* set inv. mass only */
void rta_msdr_set_mass (rta_msdr_t *sys, int i, int n, float *invmass)
{
    memcpy(sys->invmass + i, invmass, n * sizeof(float));
}



static float *vectors_alloc (int n)
{
    int    size = n * NDIM * sizeof(float);
    float *vect = (float *) rta_malloc(size);
    bzero(vect, size);
    return vect;
}

#define vect_alloc(n)    ((float *) rta_malloc((n) * sizeof(float)))


void rta_msdr_init (rta_msdr_t *sys, int maxmass, int maxlinkstotal, int maxlinkscat)
{
    int i, c;

    /* init masses */
    sys->nmasses       = 0;
    sys->massalloc     = maxmass;
    sys->linklistalloc = maxlinkscat;
    
    sys->masses = rta_malloc(maxmass * sizeof(rta_msdr_mass_t));
    sys->force  = vectors_alloc(maxmass);
    sys->speed  = vectors_alloc(maxmass);
    sys->pos2   = vectors_alloc(maxmass);
    sys->outforce = NULL;

    /* init masses' links lists */
    for (i = 0; i < maxmass; i++)
    {
	rta_msdr_mass_t *m = sys->masses + i;

	for (c = 0; c < RTA_MSDR_MAXCAT; c++)
	{
	    m->nlinks[c]  = 0;
	    m->maxdist[c] = 0;
	    m->links[c]   = rta_malloc(maxlinkscat * sizeof(int));
	}
    }

    /* init links */
    sys->nlinks     = 0;
    sys->linksalloc = maxlinkstotal;

    sys->links = rta_malloc(maxlinkstotal * sizeof(rta_msdr_link_t));
    sys->K1    = vect_alloc(maxlinkstotal);
    sys->D1    = vect_alloc(maxlinkstotal);
    sys->D2    = vect_alloc(maxlinkstotal);
    sys->Rt    = vect_alloc(maxlinkstotal);
    sys->Rf    = vect_alloc(maxlinkstotal);
    sys->l0    = vect_alloc(maxlinkstotal);
    sys->lcurr = vect_alloc(maxlinkstotal);
    sys->lprev = vect_alloc(maxlinkstotal);
    sys->stress = vect_alloc(maxlinkstotal);
    sys->forceabs = vect_alloc(maxlinkstotal);
}


/* set mass data (pos and masses) and reset to initial state */
void rta_msdr_set (rta_msdr_t *sys, int nmasses, float *pos, float *invmass)
{
    int i, c;
    
    sys->nmasses = nmasses;
    sys->nlinks  = 0;
    sys->pos     = pos;
    sys->invmass = invmass;

    memcpy(sys->pos2, sys->pos, sys->nmasses * NDIM * sizeof(float));
    bzero(sys->force, sys->nmasses * NDIM * sizeof(float));
    bzero(sys->speed, sys->nmasses * NDIM * sizeof(float));

    /* init masses' links lists */
    for (i = 0; i < sys->nmasses; i++)
    {
	rta_msdr_mass_t *m = sys->masses + i;

	for (c = 0; c < RTA_MSDR_MAXCAT; c++)
	{
	    m->nlinks[c]  = 0;
	    m->maxdist[c] = 0;
	}
    }
}


void rta_msdr_free (rta_msdr_t *sys)
{
    int i, c;

    /* free masses' links lists */
    for (i = 0; i < sys->nmasses; i++)
    {
	rta_msdr_mass_t *m = sys->masses + i;

	for (c = 0; c < RTA_MSDR_MAXCAT; c++)
	{
	    rta_free(m->links[c]);
	}
    }
    
    /* free masses' vectors (pos and invmass given from outside) */
    rta_free(sys->masses);
    rta_free(sys->force);
    rta_free(sys->speed);
    rta_free(sys->pos2);

    rta_free(sys->K1);	 
    rta_free(sys->D1);	 
    rta_free(sys->D2);	 
    rta_free(sys->Rt);	 
    rta_free(sys->Rf);	 
    rta_free(sys->l0);   
    rta_free(sys->lcurr);
    rta_free(sys->lprev);
    rta_free(sys->stress);
    rta_free(sys->forceabs);
 }
