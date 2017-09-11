/**
@file	mifdbsqlite3.c
@brief	Metric Inverted File Index Database
@author	Diemo Schwarz
@date	11.2.2010

Implementation of the interface for persistent storage of a MIF index on sqlite3.

Creates an sqlite database.  

Note that postinglists are stored as a (compressed) blob of
mif_object_t in field \a entries.
This eliminates the memory disadvantage of the design decision to
identify objects by a base (sound file) index and an object (audio
frame vector) index in mif_object_t.

\section queries Useful queries for sqlite3 testindex.db:

<code>
.headers on
.mode columns
select *, count(numobj), sum(numobj) from indexparams, discofile;
select * from discofile;
select sum(size), min(size), max(size), avg(size) from postinglist;

select refobjid, binindex, size, quote(entries) from postinglist where refobjid = 0;
select refobjid, binindex, size, length(entries), round(length(entries) / (size * 8.0) * 100) as compression, quote(entries) from postinglist where refobjid <= 1;
</code>

\todo move mifdb_dump, mifdb_read into new file mifdb.c (mif index-struct dependent code), keep mifdbsqlite3.c independent from mif.h

Copyright (C) 2008 - 2009 by IRCAM-Centre Georges Pompidou, Paris, France.
*/


#include <assert.h>
#include <string.h>
#include <sqlite3.h>
#include "mifdb.h"

#if USE_ZLIB
#include <zlib.h>
#endif


#define DEBUG_MIFDB	0


/** sql commands to create schema */
static const char *create = "\n\
drop table if exists indexparams;					\n\
create table indexparams ( /* parameters of index building */		\n\
  name		text,							\n\
  date		time,							\n\
  mifversion	decimal,						\n\
  nrefobj	integer, 						\n\
  ki		integer, 						\n\
  ndim		integer, 						\n\
  descrid	integer 						\n\
);									\n\
drop table if exists discofile;						\n\
create table discofile ( /* data file names and number of vectors */	\n\
  fileid	integer primary key,					\n\
  filename	text,    						\n\
  numobj	integer 						\n\
);									\n\
drop table if exists refobj;						\n\
create table refobj ( /* reference objects 0..nrefobj-1 */		\n\
  refobjid	integer primary key,					\n\
  fileid	integer references discofile (fileid),			\n\
  objid		integer 						\n\
);									\n\
drop table if exists postinglist;					\n\
create table postinglist ( /* posting list bins per refobj and order */ \n\
  refobjid	integer references refobj (refobjid),			\n\
  binindex	integer, 						\n\
  size		integer, /* number of entries */			\n\
  entries	blob,    /* (compressed) array of (fileid, index) */	\n\
  unique (refobjid, binindex)						\n\
);									\n\
";

/* single sql commands to insert data */
static const char *insparam = "\
insert into indexparams values (?, datetime('now'), ?, ?, ?, ?, ?);";

static const char *insfile = "\
insert into discofile values (?, ?, ?);";

static const char *insrefobj = "\
insert into refobj values (?, ?, ?);";

static const char *inspl = "\
insert into postinglist values (?, ?, ?, ?);";

/* single sql commands to query data */
static const char *getfile = "\
select fileid, filename, numobj from discofile;";

static const char *getrefobj = "\
select refobjid, fileid, objid from refobj;";

static const char *getpl = "\
select refobjid, binindex, size, entries from postinglist;";


/* convenience wrappers */
static int mifdbsqlite3_prepare (mifdbsqlite_t *db, const char *query, 
				 /*out*/ sqlite3_stmt **stmt)
{
    const char *tail   = NULL;

    if (sqlite3_prepare(db->sqlite3db, query, strlen(query) + 1, stmt, &tail) != SQLITE_OK)
    {
	rta_post("***sqlite query preparation error %d: '%s'\n  query '%s'\n",
		 sqlite3_errcode(db->sqlite3db), sqlite3_errmsg(db->sqlite3db), query);
	return 0;
    }
#if DEBUG_MIFDB
    else
	rta_post("mifdbsqlite3_prepare: stmt %p query '%s'\n  tail '%s'\n", 
		 *stmt, query, tail);
#endif

    return 1;
}

/* execute and then reset prepared query stmt that doesn't return a result */
static int mifdbsqlite3_update (mifdbsqlite_t *db, sqlite3_stmt *stmt)
{
    if (sqlite3_step(stmt) != SQLITE_DONE)
    {
	rta_post("***sqlite update execution error %d: '%s'\n",
		 sqlite3_errcode(db->sqlite3db), 
		 sqlite3_errmsg(db->sqlite3db));
	return 0;
    }
#if DEBUG_MIFDB
    else
	rta_post("mifdbsqlite3_update: step\n");
#endif

    if (sqlite3_reset(stmt) != SQLITE_OK)
    {
	rta_post("***sqlite query reset error %d: '%s'\n",
		 sqlite3_errcode(db->sqlite3db), 
		 sqlite3_errmsg(db->sqlite3db));
	return 0;
    }
#if DEBUG_MIFDB
    else
	rta_post("mifdbsqlite3_update: reset\n");
#endif

    return 1;
}


/* execute prepared query stmt for one row */
static int mifdbsqlite3_step (mifdbsqlite_t *db, sqlite3_stmt *stmt)
{
    int status = sqlite3_step(stmt);

    switch (status)
    {
    case SQLITE_ROW:
    case SQLITE_DONE:
    case SQLITE_OK:
#if DEBUG_MIFDB
	rta_post("mifdbsqlite3_step: stmt %p status %d\n", stmt, status);
#endif
	return 1;

    default:
	rta_post("***sqlite query stmt %p execution step error %d (status %d): '%s'\n",
		 stmt, sqlite3_errcode(db->sqlite3db), status,
		 sqlite3_errmsg(db->sqlite3db));
	return 0;
    }
}


/* open db */
int mifdb_open (mifdb_t *database, const char *dbname)
{
    mifdbsqlite_t *db = (mifdbsqlite_t *) database;

    if (sqlite3_open(dbname, &db->sqlite3db) != SQLITE_OK)
	/* todo: check errors */
	return 0;

    return 1;
}

/* close db */
int mifdb_close (mifdb_t *database)
{
    mifdbsqlite_t *db = (mifdbsqlite_t *) database;

    /* todo: clean up queries with sqlite3_finalize() */

    return sqlite3_close(db->sqlite3db);
    /* check errors */
}


int mifdb_begin_read (mifdb_t *database)
{
    mifdbsqlite_t *db = (mifdbsqlite_t *) database;

    /* compile statements for later */
    if (!mifdbsqlite3_prepare(db, getfile,   &db->getfile))	return 0;
    if (!mifdbsqlite3_prepare(db, getrefobj, &db->getrefobj))	return 0;
    if (!mifdbsqlite3_prepare(db, getpl,     &db->getpl))	return 0;

    return mifdb_begin_transaction(database);
}

int mifdb_end_read (mifdb_t *database)
{
    mifdbsqlite_t *db = (mifdbsqlite_t *) database;

    if (sqlite3_finalize(db->getfile)   != SQLITE_OK)	return 0;
    if (sqlite3_finalize(db->getrefobj) != SQLITE_OK)	return 0;
    if (sqlite3_finalize(db->getpl)     != SQLITE_OK)	return 0;

    return mifdb_commit_transaction(database);
}

int mifdb_begin_write (mifdb_t *database, int maxpl)
{
    mifdbsqlite_t *db = (mifdbsqlite_t *) database;

    db->complen = maxpl * sizeof(mif_object_t) + 8; /* 8 B buffer for zlib markage */
    db->compbuf = malloc(db->complen);

    /* compile statements for later */
    if (!mifdbsqlite3_prepare(db, insfile,   &db->insfile))	return 0;
    if (!mifdbsqlite3_prepare(db, insrefobj, &db->insrefobj))	return 0;
    if (!mifdbsqlite3_prepare(db, inspl,     &db->inspl))	return 0;

    return mifdb_begin_transaction(database);
}

int mifdb_end_write (mifdb_t *database)
{
    mifdbsqlite_t *db = (mifdbsqlite_t *) database;

    if (sqlite3_finalize(db->insfile)   != SQLITE_OK)	return 0;
    if (sqlite3_finalize(db->insrefobj) != SQLITE_OK)	return 0;
    if (sqlite3_finalize(db->inspl)     != SQLITE_OK)	return 0;

    free (db->compbuf);

    return mifdb_commit_transaction(database);
}

int mifdb_begin_transaction (mifdb_t *database)
{
    mifdbsqlite_t *db = (mifdbsqlite_t *) database;
    const char *query = "begin;";
    char *errmsg;

    if (sqlite3_exec(db->sqlite3db, query, NULL, NULL, &errmsg) != SQLITE_OK)
    {
	rta_post("***sqlite query execution error %d in begin: '%s'\n  query '%s'\n",
		 sqlite3_errcode(db->sqlite3db), errmsg, query);
	return 0;
    }   
    return 1;
}

int mifdb_commit_transaction(mifdb_t *database)
{
    mifdbsqlite_t *db = (mifdbsqlite_t *) database;
    const char *query = "commit;";
    char *errmsg;
    int status = sqlite3_exec(db->sqlite3db, query, NULL, NULL, &errmsg);

    if (status != SQLITE_OK  &&  status != SQLITE_DONE)
    {
	rta_post("***sqlite query execution error %d in commit: '%s'\n  query '%s'\n",
		 sqlite3_errcode(db->sqlite3db), errmsg, query);
	return 0;
    }   
    return 1;
}


/* create schema on open db, store params */
int mifdb_create (mifdb_t *database, const char *dbname, int nref, int ki, int ndim, int descrid)
{
    mifdbsqlite_t *db = (mifdbsqlite_t *) database;
    sqlite3_stmt  *stmt;
    char *errmsg;

    /* call multiple sql commands to create base */
    if (sqlite3_exec(db->sqlite3db, create, NULL, NULL, &errmsg) != SQLITE_OK)
    {
	rta_post("***sqlite query execution error %d: '%s'\n  query '%s'\n",
		 sqlite3_errcode(db->sqlite3db), errmsg, create);
	return 0;
    }

    if (!mifdbsqlite3_prepare(db, insparam, &stmt))		return 0;

    sqlite3_bind_text(stmt, 1, dbname, strlen(dbname), SQLITE_STATIC);
    sqlite3_bind_double(stmt, 2, MIFVERSION);
    sqlite3_bind_int (stmt, 3, nref);
    sqlite3_bind_int (stmt, 4, ki);
    sqlite3_bind_int (stmt, 5, ndim);
    sqlite3_bind_int (stmt, 6, descrid);

    if (!mifdbsqlite3_update(db, stmt))				return 0;
    return sqlite3_finalize(stmt) == SQLITE_OK;
}

/* get parameters from db */
int mifdb_get_params  (mifdb_t *database, int *nref, int *ki, int *ndim, int *descrid, int *nfiles, int *ndata)
{
    mifdbsqlite_t *db = (mifdbsqlite_t *) database;
    sqlite3_stmt  *stmt;
    char	  *getparam = "\
	select nrefobj, ki, ndim, descrid, count(numobj), sum(numobj) \
	from indexparams, discofile;";

    if (!mifdbsqlite3_prepare(db, getparam, &stmt))	return 0;
    if (!mifdbsqlite3_step(db, stmt))			return 0;
    
    *nref    = sqlite3_column_int(stmt, 0);
    *ki      = sqlite3_column_int(stmt, 1);
    *ndim    = sqlite3_column_int(stmt, 2);
    *descrid = sqlite3_column_int(stmt, 3);
    *nfiles  = sqlite3_column_int(stmt, 4);
    *ndata   = sqlite3_column_int(stmt, 5);

    /* free prepared statement and do last error check */
    return sqlite3_finalize(stmt) == SQLITE_OK;
}


/* add file to open db */
int mifdb_add_file (mifdb_t *database, int index, const char *filename, int nobj)
{
    mifdbsqlite_t *db = (mifdbsqlite_t *) database;
    
    sqlite3_bind_int (db->insfile, 1, index);
    sqlite3_bind_text(db->insfile, 2, filename, strlen(filename), SQLITE_STATIC);
    sqlite3_bind_int (db->insfile, 3, nobj);
 
    if (!mifdbsqlite3_update(db, db->insfile))		return 0;
   
    return 1;
}
/* get file from open db */
int mifdb_get_file (mifdb_t *database, int *index, const char **filename, int *nobj)
{
    mifdbsqlite_t *db = (mifdbsqlite_t *) database;
    sqlite3_stmt  *stmt = db->getfile;
    
    if (mifdbsqlite3_step(db, stmt))
    {
	*index    = sqlite3_column_int (stmt, 0);
	*filename = (const char *) sqlite3_column_text(stmt, 1);
	*nobj     = sqlite3_column_int (stmt, 2);
	return 1;
    }
    else
	return 0;
}

/* add refobj to open db */
int mifdb_add_refobj (mifdb_t *database, int index, const mif_object_t *obj) //int fileid, int objid)
{
    mifdbsqlite_t *db = (mifdbsqlite_t *) database;
    
    sqlite3_bind_int (db->insrefobj, 1, index);
    sqlite3_bind_int (db->insrefobj, 2, obj->base); //fileid);
    sqlite3_bind_int (db->insrefobj, 3, obj->index); //objid);
 
    if (!mifdbsqlite3_update(db, db->insrefobj))	return 0;

    return 1;
}
/* get refobj from open db */
int mifdb_get_refobj (mifdb_t *database, int *index, mif_object_t *obj)
{
    mifdbsqlite_t *db = (mifdbsqlite_t *) database;
    sqlite3_stmt  *stmt = db->getrefobj;

    if (mifdbsqlite3_step(db, stmt))
    {
	*index     = sqlite3_column_int(stmt, 0);
	obj->base  = sqlite3_column_int(stmt, 1);
	obj->index = sqlite3_column_int(stmt, 2);
	return 1;
    }
    else
	return 0;
}


/* add postinglist to open db */
int mifdb_add_postinglist (mifdb_t *database, int index, int binindex, int size, const mif_object_t *obj)
{
    mifdbsqlite_t *db = (mifdbsqlite_t *) database;
    
    sqlite3_bind_int (db->inspl, 1, index);
    sqlite3_bind_int (db->inspl, 2, binindex);
    sqlite3_bind_int (db->inspl, 3, size);
#if USE_ZLIB
    {
	uLongf complen = db->complen;
	if (compress(db->compbuf, &complen, (const Bytef *) obj, size * sizeof(mif_object_t)) != Z_OK)
	{
	    rta_post("compression error!\n");
	    return 0;
	}
	sqlite3_bind_blob(db->inspl, 4, db->compbuf, complen, SQLITE_STATIC);
    }
#else
    sqlite3_bind_blob(db->inspl, 4, obj, size * sizeof(mif_object_t), SQLITE_STATIC);
#endif 
    if (!mifdbsqlite3_update(db, db->inspl))		return 0;

    return 1;
}

/* get postinglist from open db */
int mifdb_get_postinglist (mifdb_t *database, int *index, int *binindex, int *num, int *bytes, mif_object_t **entries)
{
    mifdbsqlite_t *db = (mifdbsqlite_t *) database;
    sqlite3_stmt  *stmt = db->getpl;
    const char	  *blob;

    if (mifdbsqlite3_step(db, stmt))
    {
	*index     = sqlite3_column_int(stmt, 0);
	*binindex  = sqlite3_column_int(stmt, 1);
	*num       = sqlite3_column_int(stmt, 2);
	blob       = sqlite3_column_blob(stmt, 3);
	*bytes     = sqlite3_column_bytes(stmt, 3);
	*entries   = rta_malloc(*num * sizeof(mif_object_t));
#if USE_ZLIB
	{   /* decompress into allocated pl bin */
	    uLongf complen = *num * sizeof(mif_object_t);
	    if (uncompress((Bytef *) *entries, &complen, (const Bytef *) blob, *bytes) != Z_OK)
		return 0;
	    assert(complen == *num * sizeof(mif_object_t));
	}
#else
	/* copy given num. objects */
	memcpy(*entries, blob, *num * sizeof(mif_object_t));
#endif 	
	return 1;
    }
    else
	return 0;
}

