/*
 *  rta_cca.h
 *
 *  Created by Caramiaux Baptiste on 08/10/08.
 *  Copyright (C) 2008 by IRCAM-Centre Georges Pompidou, Paris, France.
 *
 *	Using GSL (GNU Scientific Library) under GPL Licence
 */
#ifndef _RTA_CCA_H_
#define _RTA_CCA_H_ 1


#ifdef __cplusplus
extern "C" {
#endif
	
#include "mnm.h"
	
	
	
	/*
	 * Compute Canonical Correlation Analysis between X and Y
	 *		i.e. X.A and Y.B are maximally correlated with
	 *		correlation coefficients in C
	 *
	 * Parameters :
	 *  @left_ptr -> pointer to the the first signal, i.e. two dimensional array
	 *  @left_m -> its first dimension, number of observations
	 *  @left_n -> its second dimension, number of variables
	 *  @right_ptr -> pointer to the the second signal, i.e. two dimensional array
	 *  @right_m -> its first dimension, number of observations
	 *  @right_n -> its second dimension, number of variables
	 *
	 *  @output_A -> projection matrix applied to the first signal
	 *					in order to create first canonical components
	 *  @output_q ->  projection matrix applied to the second signal
	 *					in order to create second canonical components
	 *
	 *		column (X.A)_:,j is correlated to (Y.B)_:,j with a maximal correlation
	 *		coefficient C_j
	 *		
	 *  @self -> cca structure
	 */
	int
	rta_cca(float * left_ptr, int left_m, int left_n,
			float * right_ptr, int right_m, int right_n,
			float * output_A, float * output_B, float * output_C, ftmext_t * self);
	
	
	
	
#ifdef __cplusplus
}
#endif

#endif /* _RTA_CCA_H_ */