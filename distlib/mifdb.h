/**
@file	mifdb.h
@author	Diemo Schwarz
@date	11.2.2010
@brief	Metric Inverted File Index Database

This file defines the interface for persistent storage of a MIF index.

Copyright (C) 2008 - 2009 by IRCAM-Centre Georges Pompidou, Paris, France.
*/

#ifndef _RTA_DISCODB_H_
#define _RTA_DISCODB_H_

#ifdef __cplusplus
extern "C" {
#endif

#include <sqlite3.h>
#include "mif.h"	// just for mif_object_t

/* sqlite-specific db struct */
typedef struct _mifdbsqlite
{
    sqlite3	 *sqlite3db;
    sqlite3_stmt *insfile, *insrefobj, *inspl;
    sqlite3_stmt *getfile, *getrefobj, *getpl;
} mifdbsqlite_t;


/* superclass: general db struct */
typedef struct _mifdb
{
    mifdbsqlite_t mifdbsqlite;
    mif_files_t   files;
} mifdb_t;


int mifdb_open (mifdb_t *database, const char *dbname);
int mifdb_close (mifdb_t *database);
int mifdb_begin_read (mifdb_t *database);
int mifdb_end_read (mifdb_t *database);

int mifdb_create(mifdb_t *database, const char *dbname, int nref, int ki, int ndim, int descrid);

/* get parameters from db */
int mifdb_get_params  (mifdb_t *database, int *nref, int *ki, int *ndim, int *descrid, int *nfiles, int *ndata); 

int mifdb_begin_transaction (mifdb_t *database);
int mifdb_commit_transaction (mifdb_t *database);

int mifdb_add_file (mifdb_t *database, int index, const char *filename, int nobj);
int mifdb_get_file (mifdb_t *database, int *index, const char **filename, int *nobj);

int mifdb_add_refobj (mifdb_t *database, int index, const mif_object_t *obj);
int mifdb_get_refobj (mifdb_t *database, int *index, mif_object_t *obj);

int mifdb_add_postinglist (mifdb_t *database, int index, int binindex, int size, const mif_object_t *obj);
int mifdb_get_postinglist (mifdb_t *database, int *index, int *binindex, int *size, int *bytes, mif_object_t **entries);


#ifdef __cplusplus
}
#endif

#endif /* _RTA_DISCODB_H_ */
