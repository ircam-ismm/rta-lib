/**
 * @file   rta_float.h
 * @author Jean-Philippe Lambert
 * @date   Mon Sep 10 11:05:09 2007
 * 
 * @brief  rta_real_t type wrapper for <float.h>
 *
 * Default is RTA_REAL_FLOAT. Define your RTA_REAL_TYPE (into
 * rtaConfiguration.h) to override these.
 * 
 * Copyright (C) 2007 by IRCAM-Centre Georges Pompidou, Paris, France.
 * 
 */
#ifndef _RTA_FLOAT_H_
#define _RTA_FLOAT_H_ 1

#include "rta.h"

/** @name Real type constants */
/** @{ */
#if (RTA_REAL_TYPE == RTA_FLOAT_TYPE)

#define RTA_REAL_DIG FLT_DIG
#define RTA_REAL_EPSILON FLT_EPSILON
#define RTA_REAL_MANT_DIG FLT_MANT_DIG
#define RTA_REAL_MAX FLT_MAX
#define RTA_REAL_MAX_10_EXP FLT_MAX_10_EXP
#define RTA_REAL_MAX_EXP FLT_MAX_EXP
#define RTA_REAL_MIN FLT_MIN
#define RTA_REAL_MIN_10_EXP FLT_MIN_10_EXP
#define RTA_REAL_MIN_EXP FLT_MIN_EXP

#endif /* float */

#if (RTA_REAL_TYPE == RTA_DOUBLE_TYPE)

#define RTA_REAL_DIG DBL_DIG
#define RTA_REAL_EPSILON DBL_EPSILON
#define RTA_REAL_MANT_DIG DBL_MANT_DIG
#define RTA_REAL_MAX DBL_MAX
#define RTA_REAL_MAX_10_EXP DBL_MAX_10_EXP
#define RTA_REAL_MAX_EXP DBL_MAX_EXP
#define RTA_REAL_MIN DBL_MIN
#define RTA_REAL_MIN_10_EXP DBL_MIN_10_EXP
#define RTA_REAL_MIN_EXP DBL_MIN_EXP

#endif /* double */

#if (RTA_REAL_TYPE == RTA_LONG_DOUBLE_TYPE)

#define RTA_REAL_DIG LDBL_DIG
#define RTA_REAL_EPSILON LDBL_EPSILON
#define RTA_REAL_MANT_DIG LDBL_MANT_DIG
#define RTA_REAL_MAX LDBL_MAX
#define RTA_REAL_MAX_10_EXP LDBL_MAX_10_EXP
#define RTA_REAL_MAX_EXP LDBL_MAX_EXP
#define RTA_REAL_MIN LDBL_MIN
#define RTA_REAL_MIN_10_EXP LDBL_MIN_10_EXP
#define RTA_REAL_MIN_EXP LDBL_MIN_EXP

#endif /* long double */
/** @} */
/* Real type constants */

#include <float.h>

#endif /* _RTA_FLOAT_H_ */
