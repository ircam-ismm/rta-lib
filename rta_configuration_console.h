/**
 * @file   rta_configuration.h
 * @author Diemo Schwarz
 * @date   8.12.2009
 * 
 * @brief  console configuration for the rta library
 *
 * 
 * Copyright (C) 2007 by IRCAM-Centre Georges Pompidou, Paris, France.
 * 
 */


#ifndef _RTA_CONFIGURATION_H_
#define _RTA_CONFIGURATION_H_ 1


#include "rta_stdio.h"
#include "rta_stdlib.h"


/** redefine console printing */
#undef rta_post
#define rta_post(args...) fprintf(stderr, ##args)
#include <stdio.h>

/** simple floating point precision */
#undef RTA_REAL_TYPE
#define RTA_REAL_TYPE RTA_FLOAT_TYPE

/* Apple VecLib for float and double */
#if defined(__APPLE__) && defined(__MACH__) && \
  (RTA_REAL_TYPE == RTA_FLOAT_TYPE || RTA_REAL_TYPE == RTA_DOUBLE_TYPE)
#define RTA_USE_VECLIB 1
#endif

#endif /* _RTA_CONFIGURATION_H_ */
