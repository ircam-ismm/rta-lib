/*
 *  rta_dtw.h
 *  
 *	Compute a Dynamic Time Warping
 *
 *  Created by Caramiaux Baptiste on 08/10/08.
 *  Copyright (C) 2008 - 2009 by IRCAM-Centre Georges Pompidou, Paris, France.
 *
 */



#ifndef _RTA_DTW_H_
#define _RTA_DTW_H_ 1

#include "rta.h"

#ifdef __cplusplus
extern "C" {
#endif
	
#include "mnm.h"
	
	
	/*
	 * Compute temporal alignment between two signals
	 *
	 * Parameters :
	 *  @left_ptr -> pointer to the the first signal, i.e. two dimensional array
	 *  @left_m -> its first dimension, number of observations
	 *  @left_n -> its second dimension, number of variables
	 *  @right_ptr -> pointer to the the second signal, i.e. two dimensional array
	 *  @right_m -> its first dimension, number of observations
	 *  @right_n -> its second dimension, number of variables
	 *
	 *  @output_p -> first signal indices for alignment
	 *  @output_q -> second signal indices for alignment
	 *		p[i] is the sample index of the first signal corresponding 
	 *		to the sample index q[i] of the second signal
	 *		
	 *  @length -> new signal length
	 */
	int
	rta_dtw(float * left_ptr, int left_m, int left_n, float * right_ptr, int right_m, int right_n,
			float * output_p, float * output_q, float * output_A, float * output_B, float * output_SM, int * length);
	
	
	
	
	
#ifdef __cplusplus
}
#endif

#endif /* _RTA_CCA_H_ */