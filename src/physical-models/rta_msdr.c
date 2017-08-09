
#include "math.h"
#include "rta_msdr.h"
#include <float.h>

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
    if (lim->unlimited)
	return x;
    else if (x < lim->min[d])
	return lim->min[d];
    else if (x > lim->max[d])
	return lim->max[d];
    else
	return x;
}


/** euclidean distance between two vectors */
static void distance (float *x1, float *x2, float *ret)
{
    float dist = 0;
    int   d;

    for (d = 0; d < NDIM; d++)
    {
	float distorth = x1[d] - x2[d];
	dist += distorth * distorth;
    }

    //fts_post("distance x1 %p: %f %f - x2 %p: %f %f -> d^2 %f  d %f\n", x1, x1[0], x1[1], x2, x2[0], x2[1], dist, (double) sqrtf(dist));

    *ret = sqrtf(dist);
}


/* compute link forces */
static float update_links (rta_msdr_t *sys)
{
    float totalstress = 0;
    int i, c;

    for (c = 0; c < RTA_MSDR_MAXCAT; c++)
    {
	for (i = 0; i < sys->nlinks[c]; i++)
	{
	    int    m1       = sys->links[c][i].m1;
	    int    m2       = sys->links[c][i].m2;

	    if (m1 >= 0)
	    {
		float *m1ptr    = sys->pos + m1 * NDIM;
		float *m2ptr    = sys->pos + m2 * NDIM;

		/* calculate absolute force */
		float dist     = 0;
		int d;

		distance(m1ptr, m2ptr, &dist);
		sys->stress[c][i] = dist - sys->l0[c][i];

		/* if using pow(stress): if (stress > 0) */
		if (dist > 0)
		{
		    float viscodamping = sys->D1[c][i] * (dist - sys->lprev[c][i]);
		    sys->forceabs[c][i] = (sys->K1[c][i] * sys->stress[c][i] + viscodamping) / dist;
		}
		else
		    sys->forceabs[c][i] = 0;

		/* add repulsion term if close enough */
		if (dist < sys->Rt[c][i])
		    sys->forceabs[c][i] -= (1 - dist / sys->Rt[c][i]) * sys->Rf[c][i];

		/* apply link force to mass force vector */
		for (d = 0; d < NDIM; d++)
		{
		    float forcenorm = sys->forceabs[c][i] * (m1ptr[d] - m2ptr[d]);
	
		    /* apply force normal vector to both masses, add friction */
		    sys->force[m1 * NDIM + d] 
			-= forcenorm + sys->D2[c][i] * sys->speed[m1 * NDIM + d];
		    sys->force[m2 * NDIM + d] 
			+= forcenorm - sys->D2[c][i] * sys->speed[m2 * NDIM + d];
		}

//	fts_post("update link %2d: dist %f len %f -> stress %f forceabs %f\n",  		 i, dist, sys->l0[c][i], sys->stress[c][i], sys->forceabs[c][i]);

		sys->lprev[c][i] = dist;
		totalstress += fabs(sys->stress[c][i]);
	    }
	}
    }

    return totalstress;
}


/* compute damping and friction forces */
static float rta_msdr_update_links_damping (rta_msdr_t *sys)
{
    int i, c, d;

    for (c = 0; c < RTA_MSDR_MAXCAT; c++)
    {
	for (i = 0; i < sys->nlinks[c]; i++)
	{
	    int    m1       = sys->links[c][i].m1;
	    int    m2       = sys->links[c][i].m2;

	    if (m1 >= 0)
	    {
		float *m1ptr    = sys->pos + m1 * NDIM;
		float *m2ptr    = sys->pos + m2 * NDIM;
		float  dist     = 0;

		distance(m1ptr, m2ptr, &dist);

		if (dist > 0)  /* viscosity damping */
		    sys->forceabs[c][i] = sys->D1[c][i] * (dist - sys->lprev[c][i]) / dist;
		else
		    sys->forceabs[c][i] = 0;

		/* apply link force to mass force vector */
		for (d = 0; d < NDIM; d++)
		{
		    float forceorth = sys->forceabs[c][i] * (m1ptr[d] - m2ptr[d]);
	
		    /* apply force normal vector to both masses, add friction */
		    sys->force[m1 * NDIM + d] 
			-= forceorth + sys->D2[c][i] * sys->speed[m1 * NDIM + d];
		    sys->force[m2 * NDIM + d] 
			+= forceorth - sys->D2[c][i] * sys->speed[m2 * NDIM + d];
		}
//	fts_post("update link damping %2d: len0 %f  lprev %f  lcurr %f  -> forceabs %f\n", i, sys->l0[c][i], sys->lprev[c][i], dist, sys->forceabs[c][i]);
	    }
	}
    }
}

/* compute link elasticity forces */
static float rta_msdr_update_links_elasticity (rta_msdr_t *sys)
{
    float totalstress = 0;
    int i, c, d;

    for (c = 0; c < RTA_MSDR_MAXCAT; c++)
    {
	for (i = 0; i < sys->nlinks[c]; i++)
	{
	    int    m1       = sys->links[c][i].m1;
	    int    m2       = sys->links[c][i].m2;

	    if (m1 >= 0)
	    { 
		float *m1ptr    = sys->pos + m1 * NDIM;
		float *m2ptr    = sys->pos + m2 * NDIM;
		float  dist     = 0;

		distance(m1ptr, m2ptr, &dist);
		sys->stress[c][i] = dist - sys->l0[c][i];

		if (dist > 0)
		    sys->forceabs[c][i] = sys->K1[c][i] * sys->stress[c][i] / dist;
		else
		    sys->forceabs[c][i] = 0;

		/* add repulsion term if close enough */
		if (dist < sys->Rt[c][i])
		    sys->forceabs[c][i] -= (1 - dist / sys->Rt[c][i]) * sys->Rf[c][i] / dist;

		/* apply link force to mass force vector */
		for (d = 0; d < NDIM; d++)
		{
		    float forceorth = sys->forceabs[c][i] * (m1ptr[d] - m2ptr[d]);
	
		    /* apply force normal vector to both masses */
		    sys->force[m1 * NDIM + d] -= forceorth;
		    sys->force[m2 * NDIM + d] += forceorth;
		}

//	fts_post("update link elasticity %2d: len0 %f  lprev %f  lcurr %f  -> forceabs %f  stress %f\n", i, sys->l0[c][i], sys->lprev[c][i], dist, sys->forceabs[c][i], sys->stress[c][i]);

		sys->lprev[c][i] = dist;
		totalstress += fabs(sys->stress[c][i]);
	    }
	}
    }

    return totalstress;
}


/* update masses from link forces */
static void update_masses (rta_msdr_t *sys)
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
	    sys->force [dd] = 0;
	}
    }
}


/* update masses from link forces */
static void rta_msdr_update_masses_ind (rta_msdr_t *sys, int nind, int *ind)
{
    int i, d;

    for (i = 0; i < nind; i++)
    {
	int    dd   = ind[i] * NDIM;	/* element index of row ind[i] */
	float  mass = sys->invmass[ind[i]];

	for (d = 0; d < NDIM; d++, dd++)
	{
	    /* apply force and current speed */
	    float pold = sys->pos[dd];
	    float pnew = sys->force[dd] * mass + 2 * pold - sys->pos2[dd];

	    /* set and clip speed, recalc new position */
	    sys->speed[dd] = limit(&sys->speedlim, d, pnew - pold);
	    pnew           = limit(&sys->poslim,   d, pold + sys->speed[dd]);

//	    fts_post("update_masses_ind i %d -> ind %d  mass %f  dim %d:  force %f -> speed %f  pold %f -> pnew %f\n", i, ind[i], 1. / mass, d, sys->force [dd], sys->speed[dd], pold, pnew);

	    /* cycle state: clear applied force */
	    sys->pos   [dd] = pnew;
	    sys->pos2  [dd] = pold;
	    sys->force [dd] = 0;
	}
    }
}


/* update masses from link forces without inertia */
static void rta_msdr_update_masses_limp (rta_msdr_t *sys)
{
    int i, d;

    for (i = 0; i < sys->nmasses; i++)
    {
	int    dd   = i * NDIM;	/* element index of row i */
	float  mass = sys->invmass[i];

	for (d = 0; d < NDIM; d++, dd++)
	{
	    /* apply force, DON'T apply last speed, because no inertia! */
	    float pold = sys->pos[dd];
	    float pnew = sys->force[dd] * mass + pold;

	    /* calc and clip speed, recalc new position */
	    sys->speed[dd] = limit(&sys->speedlim, d, pnew - pold);
	    pnew           = limit(&sys->poslim,   d, pold + sys->speed[dd]);

	    /* cycle state: clear applied force */
	    sys->pos   [dd] = pnew;
	    sys->force [dd] = 0;
	}
    }
}

/* update indexed masses from link forces without inertia */
static void rta_msdr_update_masses_limp_ind (rta_msdr_t *sys, int nind, int *ind)
{
    int i, d;

    for (i = 0; i < nind; i++)
    {
	int    dd   = ind[i] * NDIM;	/* element index of row ind[i] */
	float  mass = sys->invmass[ind[i]];

	for (d = 0; d < NDIM; d++, dd++)
	{
	    /* apply force, DON'T apply last speed, because no inertia! */
	    float pold = sys->pos[dd];
	    float pnew = sys->force[dd] * mass + pold;

	    /* calc and clip speed, recalc new position */
	    sys->speed[dd] = limit(&sys->speedlim, d, pnew - pold);
	    pnew           = limit(&sys->poslim,   d, pold + sys->speed[dd]);

//	    fts_post("update_masses_ind i %d -> ind %d  mass %f  dim %d:  force %f -> speed %f  pold %f -> pnew %f\n", i, ind[i], 1. / mass, d, sys->force [dd], sys->speed[dd], pold, pnew);

	    /* cycle state: clear applied force */
	    sys->pos   [dd] = pnew;
	    sys->force [dd] = 0;
	}
    }
}


/* update links and mass positions
   return total stress */
float rta_msdr_update (rta_msdr_t *sys)
{
    float stress;

    /* two-step: */
    rta_msdr_update_links_damping(sys);
    update_masses(sys);
    stress = rta_msdr_update_links_elasticity(sys);
    if (sys->outforce)
	rta_msdr_get_force(sys);
    update_masses(sys);

    return stress;
}


/* update links and mass positions
   return total stress */
float rta_msdr_update_limp (rta_msdr_t *sys)
{
    float stress = update_links(sys);

    if (sys->outforce)
	rta_msdr_get_force(sys);

    rta_msdr_update_masses_limp(sys);

    return stress;
}


/* update links and mass positions
   return total stress */
float rta_msdr_update_limp_ind (rta_msdr_t *sys, int nind, int *ind)
{
    float stress = update_links(sys);

    if (sys->outforce)
	rta_msdr_get_force(sys);

    rta_msdr_update_masses_limp_ind(sys, nind, ind);

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
    int n = 0, c;

    for (c = 0; c < RTA_MSDR_MAXCAT; c++)
    {
	n += sys->nlinks[c];
    }

    return n;
}

int rta_msdr_get_num_active_links (rta_msdr_t *sys)
{
    int n = 0, c;

    for (c = 0; c < RTA_MSDR_MAXCAT; c++)
    {
	n += sys->nactive[c];
    }

    return n;
}

int rta_msdr_get_num_links_cat (rta_msdr_t *sys, int c)
{
    return sys->nlinks[c];
}

/* set vector to recieve current link forces */
void rta_msdr_set_outforce (rta_msdr_t *sys, float *outforce)
{
    sys->outforce = outforce;
}

/* copy ncol columns of link data to out(nlinks, 11):
   masses id(2), masses pos (4), l0, lprev, stress, force, category */
int rta_msdr_get_links (rta_msdr_t *sys, float *out)
{
    int i, c, k = 0;

    for (c = 0; c < RTA_MSDR_MAXCAT; c++)
    {
	for (i = 0; i < sys->nlinks[c]; i++)
	{
	    rta_msdr_link_t *link = &sys->links[c][i];

	    if (link->m1 >= 0)
	    {
		out[0] = link->m1;
		out[1] = link->m2;
		out[2] = sys->pos[link->m1 * NDIM];
		out[3] = sys->pos[link->m1 * NDIM + 1];
		out[4] = sys->pos[link->m2 * NDIM];
		out[5] = sys->pos[link->m2 * NDIM + 1];
		out[6] = sys->l0[c][i];
		out[7] = sys->lprev[c][i];
		out[8] = sys->stress[c][i];
		out[9] = sys->forceabs[c][i];
		out[10]= c;
		out += 11;
		k++;
	    }
	}
    }

    return k;
}

/* return distance to farthest neighbour, inf if no neighbours */
float rta_msdr_get_mass_maxdist(rta_msdr_t *sys, int massi, int cat)
{
    //fts_post("mass_maxdist %d cat %d -> %f\n", massi, cat, sys->masses[massi].maxdist[cat]);	    
    return sys->masses[massi].maxdist[cat];
}

int rta_msdr_get_mass_num_links(rta_msdr_t *sys, int massi, int cat)
{
    return sys->masses[massi].nlinks[cat];
}


/* if we want to know total masses movement, run expensive sqr/sqrt
   over all lines again */
float rta_msdr_get_movement (rta_msdr_t *sys)
{
    float totalmovement = 0;
    float *speedy = sys->speed;
    int i, d;

    for (i = 0; i < sys->nmasses; i++)
    {
	float  massmove = 0;

	for (d = 0; d < NDIM; d++, speedy++)
	    massmove += *speedy * *speedy;

	totalmovement += sqrtf(massmove);
    }

    return totalmovement;
}


/* copy force vector (after links update) */
float rta_msdr_get_force (rta_msdr_t *sys)
{
    float totalforce = 0;
    float *forcey = sys->force;
    float *forceo = sys->outforce;
    int i, d;

    for (i = 0; i < sys->nmasses; i++)
    {
	float mag2 = 0;

	for (d = 0; d < NDIM; d++, forcey++, forceo++)
	{
	    float orth = *forceo = *forcey;
	    mag2      += orth * orth;
	}

	totalforce += sqrtf(mag2);
    }

    return totalforce;
}



/*
 *  constructing the system
 */

/* set masses position limits */
void rta_msdr_set_limits (rta_msdr_t *sys, int id, float *min, float *max)
{
    rta_msdr_limits_t *lim = (id == 1  ?  &sys->speedlim  :  &sys->poslim);
    int i;

    lim->unlimited = 0;
    for (i = 0; i < NDIM; i++)
    {
	lim->min[i] = min[i];
	lim->max[i] = max[i];
    }
}


/* set to no limits */
void rta_msdr_set_unlimited (rta_msdr_t *sys, int id)
{
#   define             NO_LIMIT FLT_MAX
    rta_msdr_limits_t *lim = (id == 1  ?  &sys->speedlim  :  &sys->poslim);
    int                i;
    
    /* clear limits */
    lim->unlimited = 1;
    for (i = 0; i < NDIM; i++)
    {
	lim->min[i] = -NO_LIMIT;
	lim->max[i] =  NO_LIMIT;
    }
}


void rta_msdr_set_link_length (rta_msdr_t *sys, int i, int c, float L)
{
    sys->l0[c][i] 	     = L; /* set length by index */
}
void rta_msdr_set_link_K1 (rta_msdr_t *sys, int i, int c, float K1)
{
    sys->K1[c][i] 	     = K1; /* set stiffness */
}
void rta_msdr_set_link_D1 (rta_msdr_t *sys, int i, int c, float D1)
{
    sys->D1[c][i] 	     = D1; /* set viscosity */
}
void rta_msdr_set_link_D2 (rta_msdr_t *sys, int i, int c, float D2)
{
    sys->D2[c][i] 	     = D2; /* set friction */
}
void rta_msdr_set_link_Rf (rta_msdr_t *sys, int i, int c, float rf)
{
    sys->Rf[c][i] 	     = rf; /* set repulsion force */
}
void rta_msdr_set_link_Rt (rta_msdr_t *sys, int i, int c, float rt)
{
    sys->Rt[c][i] 	     = rt; /* set repulsion distance */
}


/* set link data only */
static void set_link (rta_msdr_t *sys, int i, int c, int m1, int m2, float len,
		      float K1, float D1, float D2, float Rt, float Rf, float *ret)
{
    float *m1pos = sys->pos + m1 * NDIM;
    float *m2pos = sys->pos + m2 * NDIM;
    float  dist = 0;

    distance(m1pos, m2pos, &dist); /* previous len = current dist */
    
    sys->links[c][i].m1 = m1;
    sys->links[c][i].m2 = m2;
    sys->l0[c][i]       = len;		       /* nominal length */
    sys->lprev[c][i]    = dist;
    sys->K1[c][i] 	     = K1; 
    sys->D1[c][i] 	     = D1;
    sys->D2[c][i] 	     = D2;
    sys->Rt[c][i] 	     = Rt;
    sys->Rf[c][i] 	     = Rf;

    *ret = dist;
}

/* set backpointers to linked masses in link */
static void set_link_masses (rta_msdr_t *sys, int i, int m1, int m2, int cat, float dist)
{
	rta_msdr_mass_t *m1ptr, *m2ptr;

	/* record linkage in both masses' links lists */
	m1ptr = &sys->masses[m1];
	m2ptr = &sys->masses[m2];

	sys->links[cat][i].ind1 = m1ptr->nlinks[cat];
	sys->links[cat][i].ind2 = m2ptr->nlinks[cat];
	sys->links[cat][i].cat  = cat;

	/* nlinks[cat] < linklistalloc tested in caller */
	/* todo: insert into sorted list; if already present: overwrite link */
	m1ptr->links[cat][m1ptr->nlinks[cat]++] = i;
	m2ptr->links[cat][m2ptr->nlinks[cat]++] = i;
 
	/* update max distance of link list */
	if (dist > m1ptr->maxdist[cat])
	    m1ptr->maxdist[cat] = dist;
	if (dist > m2ptr->maxdist[cat])
	    m2ptr->maxdist[cat] = dist;
}

/** add link between masses m1, m2
    m1, m2 must be valid mass ids, cat must be valid link category */
int rta_msdr_add_link (rta_msdr_t *sys, int m1, int m2, float len, int cat,
		       float K1, float D1, float D2, float Rt, float Rf)
{
    rta_msdr_mass_t *m1ptr, *m2ptr;
    
    m1ptr = &sys->masses[m1];
    m2ptr = &sys->masses[m2];

    if (m1 != m2  &&
	sys->nlinks[cat]   < sys->linksalloc  &&
	m1ptr->nlinks[cat] < sys->linklistalloc  &&
	m2ptr->nlinks[cat] < sys->linklistalloc) 
    {
	int i = sys->nlinks[cat]++;
	float dist = 0;
	set_link(sys, i, cat, m1, m2, len, K1, D1, D2, Rt, Rf, &dist);
	set_link_masses(sys, i, m1, m2, cat, dist);
	sys->nactive[cat]++;

	//rta_post("rta_msdr_add_link %d %d  len %f  dist %f  cat %d -> %d nlinks %d %d\n", 
	//	 m1, m2, len, sys->lprev[cat][i], cat, i, m1ptr->nlinks[cat], m2ptr->nlinks[cat]);

	return i;
    }
    else
    {
	rta_post("REFUSED!!!  rta_msdr_add_link %d %d len %f  cat %d -> nlinks %d %d\n", 
		 m1, m2, len, cat, m1ptr->nlinks[cat], m2ptr->nlinks[cat]);
	return -1;
    }
}


void rta_msdr_remove_link (rta_msdr_t *sys, int i, int cat)
{
    sys->links[cat][i].m1 = -1;
    sys->links[cat][i].m2 = -1;
    sys->nactive[cat]--;
}


static void insert_mass_link(rta_msdr_t *sys, rta_msdr_mass_t *m1ptr, int i, int cat, 
			     float dist, int max)
{
    /* insert from behind into sorted list of distance */
    int    pos   = m1ptr->nlinks[cat]; /* where to insert in masses links list */
    int   *links = m1ptr->links[cat];  /* indices of links in mass */
    float *ldist = sys->l0[cat];       /* link distances */

    while (pos > 0  &&  dist < ldist[links[pos - 1]])
    {   /* move up */
	if (pos < max)
	    links[pos] = links[pos - 1];
	pos--;
    }
    
    if (pos < max)
    {
	if (m1ptr->nlinks[cat] < max)
	    m1ptr->nlinks[cat]++;
	else
	{   /* links list was full, remove link to be thrown out */
	    /*fts_post("insert_mass_link: pos %d  remove link %d  id %d  cat %d\n", 
	      pos, max-1, links[max-1], cat);*/
	    rta_msdr_remove_link(sys, links[max-1], cat);
	}

	links[pos] = i;

	/* update max distance of link list: take last element */
	m1ptr->maxdist[cat] = ldist[links[m1ptr->nlinks[cat] - 1]];
    }
    /* else: nothing to insert, dist was higher than maxdist */
}


/** insert link between masses m1, m2 into list sorted by distance, 
    throw out farthest link if not enough space, update maxdist
    m1, m2 must be valid mass ids, cat must be valid link category */
int rta_msdr_insert_link (rta_msdr_t *sys, int m1, int m2, float len, int cat,
			  float K1, float D1, float D2, float Rt, float Rf, int max)
{
    rta_msdr_mass_t *m1ptr, *m2ptr;
    
    /* record linkage in both masses' links lists */
    m1ptr = &sys->masses[m1];
    m2ptr = &sys->masses[m2];

    if (m1 != m2  &&  
	sys->nlinks[cat]   < sys->linksalloc)
    {
	int   i    = sys->nlinks[cat]++;
	float dist = 0;
	set_link(sys, i, cat, m1, m2, len, K1, D1, D2, Rt, Rf, &dist);
	sys->nactive[cat]++;

	/* if (nlinks[cat] < maxlinks[cat]) */
	/* todo: insert into sorted list; if already present: overwrite link */
	insert_mass_link(sys, m1ptr, i, cat, dist, max);
	insert_mass_link(sys, m2ptr, i, cat, dist, max);

	/*fts_post("rta_msdr_insert_link %d %d len %f  prev %f  dist %f  cat %d -> %d nlinks %d %d\n", 
	  m1, m2, len, sys->lprev[cat][i], dist, cat, i, m1ptr->nlinks[cat], m2ptr->nlinks[cat]);*/

	return i;
    }
    else
    {
	fts_post("REFUSED!!!  rta_msdr_insert_link %d %d len %f  cat %d -> nlinks %d %d\n", 
		 m1, m2, len, cat, m1ptr->nlinks[cat], m2ptr->nlinks[cat]);
	return -1;
    }
}



/* set rigidity parameter for all links */
void rta_msdr_set_K1 (rta_msdr_t *sys, int cat, float k1)
{
    int i, c, c0 = 0, c1 = RTA_MSDR_MAXCAT - 1;

    if (cat != -1)
	c0 = c1 = cat;

    for (c = c0; c <= c1; c++)
	for (i = 0; i < sys->nlinks[c]; i++)
	{
	    sys->K1[c][i] = k1;
	}
}

/* set damping parameter for all links */
void rta_msdr_set_D1 (rta_msdr_t *sys, int cat, float d1)
{
    int i, c, c0 = 0, c1 = RTA_MSDR_MAXCAT - 1;

    if (cat != -1)
	c0 = c1 = cat;

    for (c = c0; c <= c1; c++)
	for (i = 0; i < sys->nlinks[c]; i++)
	{
	    sys->D1[c][i] = d1;
	}
}

/* set friction parameter for all links */
void rta_msdr_set_D2 (rta_msdr_t *sys, int cat, float d2)
{
    int i, c, c0 = 0, c1 = RTA_MSDR_MAXCAT - 1;

    if (cat != -1)
	c0 = c1 = cat;

    for (c = c0; c <= c1; c++)
	for (i = 0; i < sys->nlinks[c]; i++)
	{
	    sys->D2[c][i] = d2;
	}
}

/* set repulsion threshold parameter for all links */
void rta_msdr_set_Rt (rta_msdr_t *sys, int cat, float rt)
{
    int i, c, c0 = 0, c1 = RTA_MSDR_MAXCAT - 1;

    if (cat != -1)
	c0 = c1 = cat;

    for (c = c0; c <= c1; c++)
	for (i = 0; i < sys->nlinks[c]; i++)
	{
	    sys->Rt[c][i] = rt;
	}
}

/* set repulsion force parameter for all links */
void rta_msdr_set_Rf (rta_msdr_t *sys, int cat, float rf)
{
    int i, c, c0 = 0, c1 = RTA_MSDR_MAXCAT - 1;

    if (cat != -1)
	c0 = c1 = cat;

    for (c = c0; c <= c1; c++)
	for (i = 0; i < sys->nlinks[c]; i++)
	{
	    sys->Rf[c][i] = rf;
	}
}


/* copy masses positions */
void rta_msdr_set_pos (rta_msdr_t *sys, int ind, int n, float *pos)
{   /* reset also old pos and speed, leave forces! */
    memcpy(sys->pos   + ind * NDIM, pos, n * NDIM * sizeof(float));
    memcpy(sys->pos2  + ind * NDIM, pos, n * NDIM * sizeof(float));
#ifndef WIN32
    bzero (sys->speed + ind * NDIM,      n * NDIM * sizeof(float));
#else
	memset(sys->speed + ind * NDIM, 0.0, n * NDIM * sizeof(float));
#endif
}

/* set inv. mass only */
void rta_msdr_set_mass (rta_msdr_t *sys, int i, int n, float *invmass)
{
    memcpy(sys->invmass + i, invmass, n * sizeof(float));
}


/* clear all links */
void rta_msdr_clear_links (rta_msdr_t *sys)
{
    int i, c;

    for (c = 0; c < RTA_MSDR_MAXCAT; c++)
    {
        sys->nlinks[c] = 0;
        sys->nactive[c] = 0;      
    }
    
    for (i = 0; i < sys->massalloc; i++)
    {
	rta_msdr_mass_t *m = sys->masses + i;

	for (c = 0; c < RTA_MSDR_MAXCAT; c++)
	{
	    m->nlinks[c]  = 0;
	    m->maxdist[c] = FLT_MAX;
	}
    }
}

/* clear all links of category cat in all masses*/
void rta_msdr_clear_cat_links (rta_msdr_t *sys, int cat)
{
    int i;

    sys->nlinks[cat] = 0;
    sys->nactive[cat] = 0;

    for (i = 0; i < sys->massalloc; i++)
    {
	rta_msdr_mass_t *m = sys->masses + i;

	m->nlinks[cat]  = 0;
	m->maxdist[cat] = FLT_MAX;
    }
}


static float *vectors_alloc (int n)
{
    int    size = n * NDIM * sizeof(float);
    float *vect = (float *) rta_malloc(size);
#ifndef WIN32
    bzero(vect, size);
#else
	memset(vect, 0.0, size);
#endif
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

    /* allocate masses' links lists */
    for (i = 0; i < sys->massalloc; i++)
    {
	for (c = 0; c < RTA_MSDR_MAXCAT; c++)
	     sys->masses[i].links[c] = rta_malloc(maxlinkscat * sizeof(int));
    }

    /* allocate links */
    sys->linksalloc = maxlinkstotal;
  
    for (c = 0; c < RTA_MSDR_MAXCAT; c++)
    {
	sys->links[c] = rta_malloc(maxlinkstotal * sizeof(rta_msdr_link_t));
	sys->K1[c]    = vect_alloc(maxlinkstotal);
	sys->D1[c]    = vect_alloc(maxlinkstotal);
	sys->D2[c]    = vect_alloc(maxlinkstotal);
	sys->Rt[c]    = vect_alloc(maxlinkstotal);
	sys->Rf[c]    = vect_alloc(maxlinkstotal);
	sys->l0[c]    = vect_alloc(maxlinkstotal);
	sys->lcurr[c] = vect_alloc(maxlinkstotal);
	sys->lprev[c] = vect_alloc(maxlinkstotal);
	sys->stress[c] = vect_alloc(maxlinkstotal);
	sys->forceabs[c] = vect_alloc(maxlinkstotal);
    }

    /* init masses' links lists and links */
    rta_msdr_clear_links (sys);

    rta_msdr_set_unlimited(sys, 0);
    rta_msdr_set_unlimited(sys, 1);
}


/* set mass data (pos and masses) and reset to initial state */
void rta_msdr_set (rta_msdr_t *sys, int nmasses, float *pos, float *invmass)
{
    int i, c;
    
    sys->nmasses = nmasses;
    sys->pos     = pos;
    sys->invmass = invmass;

//    memcpy(sys->pos1, sys->pos, sys->nmasses * NDIM * sizeof(float));
    memcpy(sys->pos2, sys->pos, sys->nmasses * NDIM * sizeof(float));
#ifndef WIN32
    bzero(sys->force, sys->nmasses * NDIM * sizeof(float));
    bzero(sys->speed, sys->nmasses * NDIM * sizeof(float));
#else
	memset(sys->force, 0.0, sys->nmasses * NDIM * sizeof(float));
    memset(sys->speed, 0.0, sys->nmasses * NDIM * sizeof(float));
#endif

    /* init masses' links lists */
    rta_msdr_clear_links (sys);
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

    for (c = 0; c < RTA_MSDR_MAXCAT; c++)
    {
	rta_free(sys->links[c]);
	rta_free(sys->K1[c]);	 
	rta_free(sys->D1[c]);	 
	rta_free(sys->D2[c]);	 
	rta_free(sys->Rt[c]);	 
	rta_free(sys->Rf[c]);	 
	rta_free(sys->l0[c]);   
	rta_free(sys->lcurr[c]);
	rta_free(sys->lprev[c]);
	rta_free(sys->stress[c]);
	rta_free(sys->forceabs[c]);
    }
}
