/**
 * @file   rta_stdlib.h
 * @author Jean-Philippe Lambert
 * @date   Mon Sep 10 11:05:09 2007
 * 
 * @brief  Default defines for the rta library
 *
 * Define your own to override these.
 * \see rta_configuration.h
 * 
 * Copyright (C) 2007 by IRCAM-Centre Georges Pompidou, Paris, France.
 * 
 */


#ifndef _RTA_STDLIB_H_
#define _RTA_STDLIB_H_ 1

/** default define for console printing */
#ifndef rta_post
#define rta_post printf
#include <stdlib.h>
#endif

/** default define for memory allocation */
#ifndef rta_malloc
#define rta_malloc malloc
#include <stdlib.h>
#endif

/** default define for memory reallocation */
#ifndef rta_realloc
#define rta_realloc realloc
#include <stdlib.h>
#endif

/** default define for memory deallocation */
#ifndef rta_free
#define rta_free free
#include <stdlib.h>
#endif

#ifndef NULL
#define NULL 0
#endif


#endif /* _RTA_STDLIB_H_ */
