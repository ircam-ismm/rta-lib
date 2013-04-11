/**
 * @file   rta_stdio.h
 * @author Diemo Schwarz
 * @date   30.10.2008
 * 
 * @brief  Default defines for the rta library
 *
 * Define your own to override these.
 * \see rta_configuration.h
 * 
 * Copyright (C) 2008 by IRCAM-Centre Georges Pompidou, Paris, France.
 * 
 */


#ifndef _RTA_STDIO_H_
#define _RTA_STDIO_H_ 1

/** default define for console printing */
#ifndef rta_post
#define rta_post printf
#include <stdio.h>
#endif

#endif /* _RTA_STDIO_H_ */
