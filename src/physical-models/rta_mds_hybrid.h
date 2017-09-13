/**
 * @file rta_mds_hybrid.h
 * @author Norbert Schnell
 * @ingroup rta_physical_models
 *
 * @brief Hybrid multi-dimensional scaling (Chalmers, Ross, Morrison algorithm)
 *
 * Based on a mass-string-damper model (rta_msdr_t),
 * pre-place a random sample of points,
 * then run iterations until stress < threshold or maxiter reached:
 * - for each point, find ns random points,
 * - insert into neighbour list if close, into random list otherwise
 * - run one iteration of msdr model
 *
 * @copyright
 * Copyright (C) 2008 - 2009 by IRCAM-Centre Georges Pompidou, Paris, France.
 * All rights reserved.
 *
 * License (BSD 3-clause)
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef _RTA_MDS_HYBRID_H_
#define _RTA_MDS_HYBRID_H_

#include "rta.h"
#include "rta_msdr.h"

#ifdef __cplusplus
extern "C" {
#endif


/** number of dimensions D, redefine before including*/
#ifndef RTA_MDS_NDIM
#define RTA_MDS_NDIM      2
#define RTA_MDS_NDIM_STR "2"
#endif


rta_real_t weighted_euclidean_distance (rta_real_t* v1, rta_real_t* v2, rta_real_t *sigma, int ndim);


#ifdef __cplusplus
}
#endif

#endif /* _RTA_MDS_HYBRID_H_ */
