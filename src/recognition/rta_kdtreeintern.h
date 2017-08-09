/**
 * @file rta_kdtreeintern.h
 * @author Diemo Schwarz
 * @date 30.10.2008
 * @version 1.0 
 *
 * @brief  private definitions and declarations for k-dimensional search tree
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

#ifndef _RTA_KDTREEINTERN_H_
#define _RTA_KDTREEINTERN_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <float.h>
#define MAX_FLOAT FLT_MAX
  
#define pow2(x)  (1 << (x))


/** helper function to print a vector \p v of length \p n with stride \p stride to the console */
void vec_post (rta_real_t *v, int stride, int n, const char *suffix);

/** helper function to print row \p i of matrix \p m of length \p n to the console */
void row_post (rta_real_t *m, int i, int n, const char *suffix);


void kdtree_stack_init (kdtree_stack_t *s, int size);
void kdtree_stack_free (kdtree_stack_t *s);
void kdtree_stack_grow (kdtree_stack_t *stack, int alloc);


/** vector to node distance */
rta_real_t distV2N (kdtree_t* t, const rta_real_t *x, const int node);
/** vector to node distance with stride */
rta_real_t distV2N_stride (kdtree_t* t, const rta_real_t *x, int stride, const int node);
/** vector to node distance with stride and weights 1/sigma */
rta_real_t distV2N_weighted (kdtree_t* t, const rta_real_t *x, int stride, const rta_real_t *sigma, const int node);




#ifdef __cplusplus
}
#endif

#endif /* _RTA_KDTREEINTERN_H_ */
