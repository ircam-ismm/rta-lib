/**
@file	discodist.h
@author	Diemo Schwarz
@date	19.1.2010
@brief	Distance function for Metric Inverted File Index

This file defines the KLS distance function for the MIF index.

Copyright (C) 2008 - 2009 by IRCAM-Centre Georges Pompidou, Paris, France.
*/

#include "rta.h"
#include "mif.h"

/** private precalculated data for KLS distance function */
typedef struct _disco_KLS_private
{
    mif_files_t *db;	/**< copy of files list, contains number of vector dimensions */
    int		 N;	/**< size of KL matrix */
    int		 r1, r2;/**< temporary vars */
    rta_real_t  *D;	/**< space for vector distance */
} disco_KLS_private_t;


/** precalc and preallocate everything */
void disco_KLS_init (void *private, mif_files_t *db);

/** free preallocated memory */
void disco_KLS_free (void *private, mif_files_t *db);

/** DISCO KLS distance function by Christophe Charbuillet */
rta_real_t disco_KLS (void *private, mif_object_t *a, mif_object_t *b);
