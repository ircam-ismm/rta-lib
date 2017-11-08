/*
 *  rta_configuration.h
 *  @brief  Custom configuration for the rta library
 *  Copyright 2008-2014 by IRCAM – Centre Pompidou. All rights reserved.
 *
 */

#ifndef _RTA_CONFIGURATIONx_H_
#define _RTA_CONFIGURATIONx_H_ 1

#define rta_post post

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

#ifdef WIN32
#define snprintf sprintf_s
#endif

#endif
