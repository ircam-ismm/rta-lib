/**
 * @file rta_cca.h
 * @author Baptiste Caramiaux
 * @date 08/10/08
 *
 * @brief Canonical Correlation Analysis
 *
 * @copyright
 * Copyright (C) 2008 - 2009 by IRCAM-Centre Georges Pompidou, Paris, France.
 * All rights reserved.
 *
 *  Using GSL (GNU Scientific Library) under GPL License
 *
 * License : GPL-v3
 *
 * This file is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _RTA_CCA_H_
#define _RTA_CCA_H_ 1


#ifdef __cplusplus
extern "C" {
#endif

// #include "mnm.h"

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