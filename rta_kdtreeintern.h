/**
@file	rta_kdtreeintern.h
@author	Diemo Schwarz
@date	30.10.2008
@version 1.0 

@brief	private definitions and declarations for k-dimensional search tree
*/


#ifndef _RTA_KDTREEINTERN_H_
#define _RTA_KDTREEINTERN_H_

#ifdef __cplusplus
extern "C" {
#endif


#define MAX_FLOAT 0x7FFFFFFF  

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
