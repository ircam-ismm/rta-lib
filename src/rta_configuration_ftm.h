/**
 * @file   rta_configuration_ftm.h
 * @author Jean-Philippe Lambert
 * @date   Mon Sep 10 11:05:09 2007
 * 
 * @brief  Custom configuration for the rta library
 *
 * 
 * Copyright (C) 2007 by IRCAM-Centre Georges Pompidou, Paris, France.
 * 
 */


#ifndef _RTA_CONFIGURATION_H_
#define _RTA_CONFIGURATION_H_ 1

/** FTM configuration */
#include "ftmlib.h"

/** console output */
#undef rta_post
#define rta_post fts_post

/** fts memory allocation */
#undef rta_malloc
#define rta_malloc fts_malloc

/** fts memory reallocation */
#undef rta_realloc
#define rta_realloc fts_realloc

/** fts memory deallocation */
#undef rta_free
#define rta_free fts_free

/** simple floating point precision */
#undef RTA_REAL_TYPE
#define RTA_REAL_TYPE RTA_FLOAT_TYPE
/* end of FTM configuration */

/* Apple VecLib for float and double */
#if defined(__APPLE__) && defined(__MACH__) && \
  (RTA_REAL_TYPE == RTA_FLOAT_TYPE || RTA_REAL_TYPE == RTA_DOUBLE_TYPE)
#define RTA_USE_VECLIB 1
#endif

#endif /* _RTA_CONFIGURATION_H_ */
