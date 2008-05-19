/**
 * @file   rta_configuration.h
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Thu May 15 21:58:21 2008
 * 
 * @brief  Matlab configuration (simple precision)
 * 
 * 
 */

#ifndef _RTA_CONFIGURATION_H_
#define _RTA_CONFIGURATION_H_ 1

/* Matlab Mex */
#include "mex.h"

/* Matlab memory management */
#undef rta_malloc
#define rta_malloc mxMalloc

#undef rta_realloc
#define rta_realloc mxRealloc

#undef rta_free
#define rta_free mxFree

/* simple floating point precision */
#undef RTA_REAL_TYPE
#define RTA_REAL_TYPE RTA_FLOAT_TYPE



#endif /* _RTA_CONFIGURATION_H_ */
