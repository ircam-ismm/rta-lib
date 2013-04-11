/**
 * @file   rta_util.h
 * @author Diemo Schwarz
 * @date   1.12.2009
 * 
 * @brief  file with common support functions
 *
 * 
 * Copyright (C) 2007 by IRCAM-Centre Georges Pompidou, Paris, France.
 * 
 */

#ifndef _RTA_UTIL_H_
#define _RTA_UTIL_H_

#include "rta.h"

#ifdef __cplusplus
extern "C" {
#endif


/** generate k random indices out of 0..n-1 in vector sample(k) */
void rta_choose_k_from_n (int k, int n, int *sample);


/** return index i such that i is the highest i for which x < arr[i]
   num !> 0 */
int rta_find_int (int x, int num, int *arr);

#ifdef __cplusplus
}
#endif

#endif /* _RTA_UTIL_H_ */
