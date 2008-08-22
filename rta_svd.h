/**
 * @file   rta_svd.h
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Mon Aug 18 09:58:20 2008
 * 
 * @brief  Singular Value Decomposition
 * 
 * Copyright (C) 2008 by IRCAM-Centre Georges Pompidou, Paris, France.
 */
#ifndef _RTA_SVD_H_
#define _RTA_SVD_H_ 1

#include "rta.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum
{
  rta_svd_out_of_place = 1,  /**< the input matrix A is first copied */
  rta_svd_in_place = 2  /**< the input matrix A can be modified */
} rta_svd_t;


/* rta_svd_setup is private (depends on implementation) */
typedef struct rta_svd_setup rta_svd_setup_t;


/** 
 * Allocate an svd setup for further svd decomposition
 *
 * Any 2D array is in row-major order. To use this function in a
 * column-major order, swap U and V, and swap m an n (also for
 * rta_svd function).
 * <pre>
 * rta_svd_setup_new(&svd_setup, svd_type, V, S, U, A, n, m);
 * rta_svd(V, S, U, A, svd_setup);
 * </pre>
 *
 * \see rta_svd_setup_delete
 * \see rta_svd
 * 
 * @param svd_setup is an address to a pointer to a private structure
 * which may depend on actual svd implementation.
 * @param svd_type can be
 *  'rta_svd_out_of_place': the matrix 'A' is first copied;
 *  'rta_svd_in_place': 'A' may be modified. This is the most
 *     efficient if m >= n. (If n > m, A is copied and transposed.)
 * @param U is a 2D array of size 'm' x 'n', or a 'NULL' pointer (it
 * is not calculated, then).
 * @param S is a 1D array of size 'n', representing the diagonal matrix
 * of the singular values of 'A'.
 * @param V is a 2D array f size 'n' x 'n', or a 'NULL' pointer (it is
 * not calculated, then).  
 * @param A is a 2D array of size 'm' x 'n'. 'A' can be modified by the
 * computation if 'svd_type' is 'rta_svd_in_place'.
 * @param m is the first dimension of 'A'.
 * @param n is the second dimension of 'A'.
 * 
 * @return 1 on success 0 on fail. If it fails, nothing should be done
 * with 'svd_setup' (even a delete).
 */
int
rta_svd_setup_new(rta_svd_setup_t ** svd_setup, const rta_svd_t svd_type,
                  rta_real_t * U, rta_real_t * S, rta_real_t *  V, 
                  rta_real_t * A, const unsigned int m, const unsigned int n);

/** 
 * Deallocate any (sucessfully) allocated svd setup.
 *
 * \see rta_svd_setup_new
 * 
 * @param svd_setup is a pointer to the memory wich will be released.
 */
void
rta_svd_setup_delete(rta_svd_setup_t * svd_setup);

/** 
 *
 * For an m-by-n matrix A with m >= n, the singular value
 * decomposition is an m-by-n orthogonal matrix U, an n-by-n diagonal
 * matrix S, and an n-by-n orthogonal matrix V so that A = U*S*V'. 
 * (V' is the conjugate transpose of V and only the diagonal vector of
 * S is represented here.) 
 *
 * The singular values, S[k] are ordered so that
 * S[0] >= S[1] >= ... >= S[n-1].
 *
 * The singular value decompostion always exists, so the constructor will
 * never fail. The matrix condition number and the effective numerical
 * rank can be computed from the singlar values.
 *
 * Any array must be allocated before calling this function.
 *
 * Any 2D array is in row-major order. To use this function in a
 * column-major order, swap U and V (and swap m an n during the setup,
 * too). 
 * <pre>
 * rta_svd_setup_new(&svd_setup, svd_type, V, S, U, A, n, m);
 * rta_svd(V, S, U, A, svd_setup);
 * </pre>
 * \see rta_svd_setup_new
 * 
 * @param U is a 2D array of size m x n, or a NULL pointer (it
 * is not calculated, then).
 * @param S is a 1D array of size n, representing the diagonal matrix
 * of the singular values of 'A'.
 * @param V is a 2D array f size n x n, or a NULL pointer (it is not
 * calculated, then).  
 * @param A is a 2D array of size m x n. A is modified by the
 * computation
 * @param svd_setup is a previously allocated setup (by rta_setup_new).
 * 
 */
void
rta_svd(rta_real_t * U, rta_real_t * S, rta_real_t *  V, rta_real_t * A,
        const rta_svd_setup_t * svd_setup);

#ifdef __cplusplus
}
#endif

#endif /* _RTA_SVD_H_ */
