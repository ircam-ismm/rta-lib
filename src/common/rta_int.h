/** 
 * @file   rta_int.h
 * @author Jean-Philippe Lambert
 * @date   Thu Sep 12 18:10:41 2007
 * 
 * @brief  Integer mathematical functions
 * 
 *
 * Copyright (C) 2007 - 2009 by IRCAM-Centre Georges Pompidou, Paris, France.
 * 
 */
#ifndef _RTA_INT_H_
#define _RTA_INT_H_ 1

#include "rta.h"

#ifdef __cplusplus
extern "C" {
#endif

/** 
 * Integer version of the log2 function (rounded down)
 * 
 * @param n 
 * 
 * @return log2('n') or 0 if 'n' == 0
 */
unsigned int rta_ilog2(unsigned int n);

/** 
 * Integer version of the maximum
 * 
 * @param m 
 * @param n
 * 
 * @return the smallest integer between 'm' and 'n'
 */
extern inline int rta_imax(int m, int n);

/** 
 * Integer version of the minimum
 * 
 * @param m 
 * @param n 
 * 
 * @return the smallest integer between 'm' and 'n'
 */
extern inline int rta_imin(int m, int n);


/** 
 * Next (or equal) 'n' power of 2
 * 
 * @param n 
 * 
 * @return minimum p such as 2^p >= 'n'
 */
unsigned int rta_inextpow2(unsigned int n);

#ifdef __cplusplus
}
#endif

#endif /* _RTA_INT_H_ */
