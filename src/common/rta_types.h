/**
 * @file   rta_types.h
 * @author Jean-Philippe Lambert
 * @date   Mon Sep 10 11:05:09 2007
 * 
 * @brief  rta_real_t type warper for float, double or long double
 *
 * Default is RTA_REAL_FLOAT. Define your RTA_REAL_TYPE to override these.
 * 
 * Copyright (C) 2007 by IRCAM-Centre Georges Pompidou, Paris, France.
 * 
 */


#ifndef _RTA_TYPES_H_
#define _RTA_TYPES_H_ 1

/** default floating point precision */
#ifndef RTA_REAL_TYPE
#define RTA_REAL_TYPE RTA_FLOAT_TYPE
#endif

#undef RTA_FLOAT_TYPE
#define RTA_FLOAT_TYPE 1
#undef RTA_DOUBLE_TYPE
#define RTA_DOUBLE_TYPE 2
#undef RTA_LONG_DOUBLE_TYPE
#define RTA_LONG_DOUBLE_TYPE 3

#if (RTA_REAL_TYPE == RTA_FLOAT_TYPE)
#undef rta_real_t
#define rta_real_t float
#endif

#if (RTA_REAL_TYPE == RTA_DOUBLE_TYPE)
#undef rta_real_t
#define rta_real_t double
#endif

#if (RTA_REAL_TYPE == RTA_LONG_DOUBLE_TYPE)
#undef rta_real_t
#define rta_real_t long double
#endif

#endif /* _RTA_TYPES_H_ */
