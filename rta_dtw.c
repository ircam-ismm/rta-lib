/**
 *  rta_dtw.c
 *
 *  Created by Caramiaux Baptiste on 08/10/08.
 *  Copyright (C) 2008 - 2009 by IRCAM-Centre Georges Pompidou, Paris, France.
 *
 */


#include<math.h>

#include "rta_math.h"
#include "rta_dtw.h"

#define INF HUGE_VAL



/*
 * Construct Score Matrix
 *
 * Score matrix is based on euclidean distance
 * Signals must have the same number of observations
 */
static int
prepare_score_matrix(float * A, int m_A, float * B, int m_B, int n, float * out_R)
{
	int		i, j, k;
	float	max_d = 0.;
	float	tmp = 0.;
	
	//Compute euclidean distance between i-th first signal observation
	//and i-th second signal observation for all i
	for ( i = 0; i < m_A; i++ )
		for ( j = 0; j < m_B; j++ )
		{
			out_R[i * m_B + j] = 0.0;
			
			for ( k = 0; k < n; k++ )
				out_R[i * m_B + j] += ( A[i * n + k] - B[j * n + k] ) * ( A[i * n + k] - B[j * n + k] );
			
			if ( fabs(out_R[i * m_B + j]) > max_d)
				max_d = fabs(out_R[i * m_B + j]);
		}
	
	//Normalization
	for ( i = 0; i < m_A; i++ )
		for ( j = 0; j < m_B; j++ )
		{
			tmp = (out_R[i * m_B + j] / max_d);
			out_R[i * m_B + j] = tmp;
		}	
	
	return 0;
}



/*
 * Dynamic Programming Routine
 * 
 * Copyright (c) 2003-05 Dan Ellis <dpwe@ee.columbia.edu>
 * released under GPL - see file COPYRIGHT
 */
static int
dpcore(float * pM, int rows, int cols, int * pC, int crows, int ccols, float * pD, int * pP)
{
	
	
	//int rows, cols, 
	int i, j, k, tb;
	
	float d1, d2, v;
	float *costs;
	int *steps;
	int ncosts;
	int ii;
	
	ncosts = crows;
	costs = (float *)malloc( ncosts * sizeof(float));
	steps = (int *)malloc( ncosts * 2 * sizeof(int));
	
	for (ii = 0; ii < ncosts; ++ii) 
	{
		steps[2*ii] = (int)(pC[ii * ccols]);
		steps[2*ii+1] = (int)(pC[ii * ccols + 1]);
		costs[ii] = (float)(pC[ii * ccols + 2]);
	}
	
	
	
	/* do dp */
	v = 0;	
	tb = 1;	/* value to use for 0,0 */
	
	for (j = 0; j < cols; ++j) 
	{
		for (i = 0; i < rows; ++i) 
		{
			d1 = pM[i * cols + j];
			
			for (k = 0; k < ncosts; ++k) 
			{
				if ( i >= steps[2*k] && j >= steps[2*k+1] ) 
				{
					d2 = costs[k]*d1 + pD[(i-steps[2*k]) * cols + (j-steps[2*k+1])];
					
					if (d2 < v) {
						v = d2;
						tb = k;
					}
				}
			}
			
			pD[i * cols + j] = v;
			pP[i * cols + j] = tb;
			v = INF;
		}
	}
	free((void *)costs);
	free((void *)steps);
	
}





/*
 *
 * dpfast fill index tables p and q using the dpcore function
 * dpcore function is copyrighted by Dan Ellis, under GPL licence
 *
 */
static int
dpfast(float * SM, int m, int n, int T1, int T2, int * p, int * q, int * length)
{
	
	
	int i	= m-1;
	int	j	= n-1;
	int cpt = 0;
	int k, tb;
	
	
	/*
	 * C is the step and weight matrix
	 * two first columns : first and second dimension steps
	 * thirs columns : the corresponding weight at (i,j)
	 */
	int * C = malloc( 9 * sizeof(int));
	C[0] = 1; C[1] = 1; C[2] = 1;
	C[3] = 1; C[4] = 0; C[5] = 1;
	C[6] = 0; C[7] = 1; C[8] = 1;
	
	//cost matrix
	float * pD	= malloc( m * n * sizeof(float));
	
	//index matrix
	int * pP	= malloc( m * n * sizeof(int));
	
	//initialization
	for( k = 0; k < m * n; k++)
	{
	    pD[k] = 0.0; pP[k] = 0;
	}
	
	
	//computation of the cost matrix
	dpcore( SM, m, n, C, 3, 3, pD, pP);

	
	//build index tables corresponding to the
	//first and the second signal
	p[0] = i;
	q[0] = j;
	
	while ( i!= 0 || j!= 0)
	{
		tb = pP[i * n + j];
			
		i -= C[tb * 3];
		j -= C[tb * 3 + 1];
		
		cpt+=1;
		
		p[cpt] = i;
		q[cpt] = j;
	}
	
	*length = cpt + 1;
}



/*
 * rta_dtw : main function
 *	see rta_dtw.h
 */
int
rta_dtw(float * left_ptr, int left_m, int left_n, float * right_ptr, int right_m, int right_n,
		float * output_p, float * output_q, int * length)
{
	
	int i;
	
	
	//Fill cost (or score) matrix : SM
	float * SM	= malloc( left_m * right_m * sizeof(float));
	prepare_score_matrix( left_ptr, left_m, right_ptr, right_m, left_n, SM);
	
	
	//Compute index tables p,q using dynamic programming
	int * p		= malloc( (left_m + right_m + 1) * sizeof(int));
	int * q		= malloc( (left_m + right_m + 1) * sizeof(int));
	dpfast(SM, left_m, right_m, 1, 1, p, q, length);

	
	//Build outputs
	for (i = 0; i < *length; i++)
	{
			output_p[i] = p[(*length - 1) - i];
			output_q[i] = q[(*length - 1) - i];
	}
	
}