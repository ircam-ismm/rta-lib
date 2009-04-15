/**
 *  rta_cca.c
 *
 *  Created by Caramiaux Baptiste on 08/10/08.
 *  Copyright (C) 2008 - 2009 by IRCAM-Centre Georges Pompidou, Paris, France.
 *
 *	Using GSL (GNU Scientific Library) under GPL Licence
 */




#include<gsl/gsl_matrix.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_blas.h>
#include<gsl/gsl_statistics_float.h>
#include<gsl/gsl_math.h>

#include<math.h>
#include "rta_cca.h"



/*
 * Center Variables
 */
static int
center_variables(gsl_matrix * X, gsl_matrix * Y)
{
	int i, j, k;
	
	int left_m = X->size1;
	int left_n = X->size2;
	
	int right_m = Y->size1;
	int right_n = Y->size2;
	
	
	gsl_vector * left_means = gsl_vector_alloc( left_n);
	gsl_vector_set_zero( left_means);
	gsl_vector * right_means = gsl_vector_alloc( right_n);
	gsl_vector_set_zero(right_means);
	
	
	if (left_m != 0 && right_m != 0)
	{
		for (i = 0; i < left_m; i++)
		{
			for (j = 0; j < left_n; j++)
				(left_means->data)[j]	+= (X->data)[i * left_n + j];
			
			for (k = 0; k < right_n; k++)
				(right_means->data)[k]	+= (Y->data)[i * right_n + k];
		}
		for (i = 0; i < left_m; i++)
		{
			for (j = 0; j < left_n; j++)
				(X->data)[i * left_n + j] -= (left_means->data)[j] / left_m;
			
			for (k = 0; k < right_n; k++)
				(Y->data)[i * right_n + k] -= (right_means->data)[k] / right_m;
		}
	}
	else
	{
		return -1;
	}
	
	
	free(left_means);
	free(right_means);
}


/* 
 * Compute the product A * B = C 
 *
 * BLAS matrix product, available in GSL
 * has produced numerical errors
 */
static int
matrix_product( gsl_matrix * A, gsl_matrix * B, gsl_matrix * C)
{
	
	int i, j, k;
	double tmp;
	int mA = A->size1; int nA = A->size2;
	int mB = B->size1; int nB = B->size2;
	
	if ( nA != mB )
		return -1;
	else
	{
		for (i = 0 ; i < mA; i++)
			for (j = 0 ; j < nB; j++)
			{
				tmp = 0.;
				
				for (k = 0; k < nA; k++)
					tmp += gsl_matrix_get(A, i, k) * gsl_matrix_get(B, k, j);
				
				gsl_matrix_set( C, i, j, tmp);
			}
	}
	return 1;
}




/* 
 * QR Decomposition
 *
 * M matrix is decomposed in Q and R. P is a permutation matrix which
 * sort diagonal elements in R with decreasing absolute values
 *
 */
static void
sorted_QR_decomposition( gsl_matrix * M, gsl_matrix * Q, gsl_matrix * R, gsl_permutation * P)
{
	int m		= M->size1;
	int n		= M->size2;
	int min_mn	= GSL_MIN_INT( m, n);
	
	gsl_vector * Tau	= gsl_vector_alloc( min_mn);
	gsl_vector * Norms	= gsl_vector_alloc( n);
	int * Signs			= (int *)malloc( n * sizeof(int));
	
	
	//QR decomposition with permutations in P
	gsl_linalg_QRPT_decomp( M, Tau, P, Signs, Norms);
	
	//Extract from QR Q and R matrices
	gsl_linalg_QR_unpack( M, Tau, Q, R);
	
	
	gsl_vector_free(Tau);
	gsl_vector_free(Norms);
	free(Signs);
	
}


/* 
 * Compute Rank
 *
 * Return matrix rank in the case of a triangular matrix M
 */
static int 
compute_rank(gsl_matrix * M)
{
	
	//M is a triangular matrix
	int n = M->size2;
	int i;
	int rank = 0;
	float epsilon	= 2e-12;
	
	
	for (i = 0; i < n; i++)
	{
		if ( fabs((M->data)[i * n + i]) > epsilon )
			rank += 1;
	}
	
	return rank;
}


/* 
 * Compute Pseudo-Inverse
 *
 * Return M pseudo-inverse in invM using Singular
 * Value Decomposition available in GSL
 */
static void
compute_pseudo_inverse( gsl_matrix * M, gsl_matrix * invM)
{
	
	int i,j,k;
	int m = M->size1;
	int n = M->size2;
	double tmp;
	double epsilon = 2e-10;
	
	gsl_matrix * V = gsl_matrix_alloc( n, n);
	gsl_vector * S = gsl_vector_alloc( n);
	gsl_vector * work = gsl_vector_alloc( n);
	
	//Singular Value Decomposition
	// M = U.S.Vt
	gsl_linalg_SV_decomp( M, V, S, work);
	
	
	/*
	 * Pseudo-Inverse is invM = V.inv(S).Ut
	 */
	
	//Compute matrix product V.inv(S)
	for (i = 0; i < n; i++)
		for (j = 0; j < n; j++)
		{
			if ( fabs((S->data)[j]) < epsilon)
				gsl_matrix_set( V, i, j, 0.);
			else
			{
				tmp = gsl_matrix_get( V, i, j);
				tmp /= gsl_vector_get( S, j);
				gsl_matrix_set( V, i, j, tmp);
			}	
		}
	
	//Compute matrix product V.inv(S) with Ut	
	for (i = 0; i < n; i++)
		for (j = 0; j < m; j++)
		{
			tmp = 0.;
			for (k = 0; k < n; k++)
				tmp += (V->data)[i * n + k] * (M->data)[j * n + k];
			
			gsl_matrix_set( invM, i, j, tmp);
		}
	
	
	gsl_matrix_free(V);
	gsl_vector_free(S);
	gsl_vector_free(work);
	
}




/*
 * Canonical Correlation Analysis routine
 *	see rta_cca.h for specifications
 */
int
rta_cca(float * left_ptr, int left_m, int left_n,
		float * right_ptr, int right_m, int right_n,
		float * output_A, float * output_B, float * output_C, ftmext_t * self)
{
	
	
	
	int i, j, k;
	
	//Condition : same number of observations
	if ( left_m != right_m )
		return 1;
	
	
	
	//X : first multivariate table
	//Y : second multivariate table
	gsl_matrix * X = gsl_matrix_alloc( left_m, left_n);
	gsl_matrix * Y = gsl_matrix_alloc( right_m, right_n);
	
	//Fill matrices
	for (i = 0; i < left_m; i++)
	{
		for (j = 0; j < left_n; j++)
			(X->data)[i * left_n + j] = (double) left_ptr[i * left_n + j];
		
		for (k = 0; k < right_n; k++)
			(Y->data)[i * right_n + k] = (double) right_ptr[i * right_n + k];
	}
	
	//Center Variables
	if ( center_variables( X, Y) == -1 )
		return 2;
	
	
	

	/* 
	 * QR Decomposition
	 *
	 * X -> QX, RX
	 * Y -> QY, RY
	 */
	
	gsl_matrix * QX			= gsl_matrix_alloc( left_m, left_m);
	gsl_matrix * RX			= gsl_matrix_alloc( left_m, left_n);
	gsl_permutation * PX	= gsl_permutation_alloc( left_n);
	
	sorted_QR_decomposition(X, QX, RX, PX);

	
	gsl_matrix * QY			= gsl_matrix_alloc( right_m, right_m);
	gsl_matrix * RY			= gsl_matrix_alloc( right_m, right_n);
	gsl_permutation * PY	= gsl_permutation_alloc( right_n);
	
	sorted_QR_decomposition(Y, QY, RY, PY);

	
	
	/*
	 * Following RX and RY ranks, algorithm
	 * keeps the first "rank" variables
	 */
	
	// Computing ranks
	int rankX		= compute_rank( RX);
	int rankY		= compute_rank( RY);
	int min_ranks	= GSL_MIN_INT( rankX, rankY);
	
	// Modifying matrices
	gsl_matrix * subQX = gsl_matrix_alloc( left_m, rankX);
	gsl_matrix * subRX = gsl_matrix_alloc( rankX, rankX);
	gsl_matrix * subQY = gsl_matrix_alloc( right_m, rankY);
	gsl_matrix * subRY = gsl_matrix_alloc( rankY, rankY);
	
	if (rankX == 0)
		return 3;
	else
	{
		for(i = 0; i < rankX; i++)
			for(j = 0; j < rankX; j++)
			{
				gsl_matrix_set(subQX, i, j, gsl_matrix_get(QX, i, j));
				gsl_matrix_set(subRX, i, j, gsl_matrix_get(RX, i, j));
			}
		for(i = rankX; i < left_m; i++)
			for(j = 0; j < rankX; j++)
			{
				gsl_matrix_set(subQX, i, j, gsl_matrix_get(QX, i, j));
			}
	}
	
	if (rankY == 0)
		return 3;
	else
	{
		for(i = 0; i < rankY; i++)
			for(j = 0; j < rankY; j++)
			{
				gsl_matrix_set(subQY, i, j, gsl_matrix_get(QY, i, j));
				gsl_matrix_set(subRY, i, j, gsl_matrix_get(RY, i, j));
			}
		for(i = rankY; i < right_m; i++)
			for(j = 0; j < rankY; j++)
			{
				gsl_matrix_set(subQY, i, j, gsl_matrix_get(QY, i, j));
			}
		
	}
	
	gsl_matrix_free(RX);
	gsl_matrix_free(QX);
	gsl_matrix_free(RY);
	gsl_matrix_free(QY);
	
	
	
	
	/* 
	 * Singular Value Decomposition of QXt.QY
	 * 
	 * Two cases :
	 *	- if QXt number of rows > QY number of columns
	 *	- if QXt number of rows < QY number of columns
	 */
	int svdm;
	int svdn;
	double tmp_sum = 0.;
	
	gsl_matrix * M_svd;
	gsl_vector * S;
	gsl_matrix * V;
	gsl_vector * work;
	gsl_matrix * L2;
	gsl_matrix * L1;
	
	
	//Compute QXtQY or QYtQX depending to the case
	if (rankX <= rankY)
	{
		svdm = rankY;
		svdn = rankX;
		M_svd = gsl_matrix_alloc( svdm, svdn);
		
		for (i = 0 ; i < subQY->size2; i++)
			for (j = 0 ; j < subQX->size2; j++)
			{
				tmp_sum = 0.;
				
				for (k = 0; k < subQY->size1; k++)
					tmp_sum += gsl_matrix_get(subQY, k, i) * gsl_matrix_get(subQX, k, j);
				
				gsl_matrix_set(M_svd, i, j, tmp_sum);
			}

	}
	else
	{
		svdm = rankX;
		svdn = rankY;
		M_svd = gsl_matrix_alloc( svdm, svdn);
		
		for (i = 0 ; i < subQX->size2; i++)
			for (j = 0 ; j < subQY->size2; j++)
			{
				tmp_sum = 0.;
				
				for (k = 0; k < subQX->size1; k++)
					tmp_sum += gsl_matrix_get(subQX, k, i) * gsl_matrix_get(subQY, k, j);
				
				gsl_matrix_set(M_svd, i, j, tmp_sum);
			}
		
	}
	
	
	//Compute Singular Value Decomposition
	V		= gsl_matrix_alloc( svdn, svdn);
	S		= gsl_vector_alloc( svdn);
	work	= gsl_vector_alloc( svdn);

	gsl_linalg_SV_decomp( M_svd, V, S, work);
	
	
	
	
	/*
	 * Keep first "rank" variables in eigenvector matrices
	 * produced by SVD
	 *
	 * U -> L1
	 * V -> L2
	 */
	if (rankX <= rankY)
	{
		L1	= gsl_matrix_alloc( svdn, min_ranks);
		L2	= gsl_matrix_alloc( svdm, min_ranks);
		
		for(j = 0; j < min_ranks; j++)
		{
			for(i = 0; i < svdn; i++)
				gsl_matrix_set(L1, i, j, gsl_matrix_get(V, i, j));
			
			for(k = 0; k < svdm; k++)
				gsl_matrix_set(L2, k, j, gsl_matrix_get(M_svd, k, j));
		}
	}
	else
	{
		L1	= gsl_matrix_alloc( svdm, min_ranks);
		L2	= gsl_matrix_alloc( svdn, min_ranks);
		
		for(j = 0; j < min_ranks; j++)
		{
			for(i = 0; i < svdn; i++)
				gsl_matrix_set(L2, i, j, gsl_matrix_get(V, i, j));
			
			for(k = 0; k < svdm; k++)
				gsl_matrix_set(L1, k, j, gsl_matrix_get(M_svd, k, j));
		}
	}
	
	
	/*
	 * Compute full rank projection matrices : ResX and ResY
	 *
	 * pseudo-inverses of subRX and subRY : invsubRX, invsubRY
	 * ResX = invsubRX.L1
	 * ResY = invsubRY.L2
	 *
	 */

	//Pseudo-inverse of subRX
	gsl_matrix * ResX = gsl_matrix_alloc( rankX, min_ranks);
	gsl_matrix * invsubRX = gsl_matrix_alloc( rankX, rankX);
	compute_pseudo_inverse( subRX, invsubRX);
	
	//Compute ResX
	matrix_product( invsubRX, L1, ResX);
	gsl_matrix_scale( ResX, sqrt((float)(left_m - 1)));
	

	//Pseudo-inverse of subRY
	gsl_matrix * ResY = gsl_matrix_alloc( rankY, min_ranks);
	gsl_matrix * invsubRY = gsl_matrix_alloc( rankY, rankY);
	compute_pseudo_inverse( subRY, invsubRY);
	
	//Compute ResY
	matrix_product( invsubRY, L2, ResY);
	gsl_matrix_scale( ResY, sqrt((float)(right_m - 1)));	
	
	
	
	
	/* 
	 * Computing A and B : projection matrices
	 * 
	 * Computation from ResX and ResY with real row order
	 * thanks to permutation matrices
	 */
	
	gsl_matrix * preA = gsl_matrix_alloc( left_n, min_ranks);
	gsl_matrix_set_zero(preA);
	
	for(i = 0; i < rankX; i++)
		for(j = 0; j < min_ranks; j++)
			gsl_matrix_set( preA, i, j, gsl_matrix_get( ResX, i, j));
	
	for(i = 0; i < left_n; i++)
		for(j = 0; j < min_ranks; j++)
			output_A[(int) gsl_permutation_get( PX, i) * min_ranks + j] = (float) gsl_matrix_get( preA, i, j);
	
	
	
	gsl_matrix * preB = gsl_matrix_alloc(right_n, min_ranks);
	gsl_matrix_set_zero(preB);
	
	for(i = 0; i < rankY; i++)
		for(j = 0; j < min_ranks; j++)
			gsl_matrix_set( preB, i, j, gsl_matrix_get( ResY, i, j));
	
	
	
	//Outputs ....
	for(i = 0; i < right_n; i++)
		for(j = 0; j < min_ranks; j++)
			output_B[(int) gsl_permutation_get( PY, i) * min_ranks + j] = (float) gsl_matrix_get( preB, i, j);

	for(i = 0; i < S->size; i++)
		output_C[i] = (float) gsl_vector_get( S, i);
	
	
	
	gsl_matrix_free(preA);
	gsl_matrix_free(preB);
	
}