#ifndef _RTA_CONFIGURATION_H_
#define _RTA_CONFIGURATION_H_ 1

//#define rta_post post
#define rta_post printf

/** fts memory allocation */
#undef rta_malloc
#define rta_malloc malloc

/** fts memory reallocation */
#undef rta_realloc
#define rta_realloc realloc

/** fts memory deallocation */
#undef rta_free
#define rta_free free

/** simple floating point precision */
#undef RTA_REAL_TYPE
#define RTA_REAL_TYPE RTA_FLOAT_TYPE

/* Apple VecLib for float and double */
#if defined(__APPLE__) && defined(__MACH__) && \
(RTA_REAL_TYPE == RTA_FLOAT_TYPE || RTA_REAL_TYPE == RTA_DOUBLE_TYPE)
#define RTA_USE_VECLIB 1
#endif

// rta_complex.h is entirely determined by a #ifdef WIN32 directive
// ----------------------------------------------------------------

// FOUND IN rta_biquad:
// --------------------

#ifdef WIN32
#define inline
#endif

// FOUND IN rta_kdtree.c:
// ----------------------

// beginning:
#ifdef WIN32
#include <malloc.h>
#define snprintf sprintf_s
#else
#include <alloca.h>
#include <strings.h>
#endif

// kdtree_data_display:
#ifdef WIN32
  rta_real_t plane[2048];
#else
  rta_real_t plane[t->ndim];
#endif

#ifndef WIN32
    bzero(plane, t->ndim * sizeof(rta_real_t));
#else
    memset(plane, 0.0, t->ndim * sizeof(rta_real_t));
#endif

// kdtree init nodes:
#ifndef WIN32
  bzero(self->nodes, self->nnodes * sizeof(kdtree_node_t));
#else
  memset(self->nodes, 0.0, self->nnodes * sizeof(kdtree_node_t));
#endif

// FOUND IN rta_kdtreebuild.c:
// ---------------------------

// beginning:
#ifndef WIN32
#include <strings.h>
#endif

// compute_splitplane:
#ifndef WIN32
  bzero(split_ptr, t->ndim * sizeof(rta_real_t));
#else
  memset(split_ptr, 0.0, t->ndim * sizeof(rta_real_t));
#endif

// FOUND IN rta_mahalanobis.h:
// ---------------------------

#ifdef WIN32
#include <malloc.h>
#else
#include <alloca.h>
#endif

// rta_math.h is half determined by a #ifdef WIN32 directive
// ---------------------------------------------------------

// FOUND IN rta_mds_hybrid.h:
// --------------------------
#ifndef WIN32
  bzero(vect, size);
#else
  memset(vect, 0.0, size);
#endif

// FOUND IN rta_msdr.h:
// --------------------

// rta_msdr_set_pos:
#ifndef WIN32
  bzero (sys->speed + ind * NDIM,      n * NDIM * sizeof(float));
#else
  memset(sys->speed + ind * NDIM, 0.0, n * NDIM * sizeof(float));
#endif

// vectors_alloc:
#ifndef WIN32
    bzero(vect, size);
#else
  memset(vect, 0.0, size);
#endif

// rta_msdr_set:
#ifndef WIN32
  bzero(sys->force, sys->nmasses * NDIM * sizeof(float));
  bzero(sys->speed, sys->nmasses * NDIM * sizeof(float));
#else
  memset(sys->force, 0.0, sys->nmasses * NDIM * sizeof(float));
  memset(sys->speed, 0.0, sys->nmasses * NDIM * sizeof(float));
#endif

// FOUND in rta_onepole.h:
// -----------------------
#ifdef WIN32
#define inline
#endif

// FOUND in rta_selection.h:
// -------------------------
#ifdef WIN32
#define inline
#endif

// FOUND IN rta_util.h:
// --------------------
#ifdef WIN32
static long random(){return rand();}
#endif


#endif /* _RTA_CONFIGURATION_H_ */
