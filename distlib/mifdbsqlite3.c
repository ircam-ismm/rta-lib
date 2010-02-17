/**
@file	mifdbsqlite3.c
@author	Diemo Schwarz
@date	11.2.2010
@brief	Metric Inverted File Index Database

Implementation of the interface for persistent storage of a MIF index on sqlite3.

Creates an sqlite database.  

Useful queries in sqlite3 testindex.db:

.headers on
.mode tabs
select * from indexparams;
select sum(size), min(size), max(size), avg(size) from postinglist;

select * from postinglist where refobjid = 0 and binindex = 0;

select nrefobj, ki, ndim, descrid, count(numobj), sum(numobj) from indexparams, discofile;


Copyright (C) 2008 - 2009 by IRCAM-Centre Georges Pompidou, Paris, France.
*/


#include <string.h>
#include <sqlite3.h>
#include "mifdb.h"


#define DEBUG_MIFDB	0


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
	rta_post("mifdbsqlite3_prepare: query '%s'\n  tail '%s'\n", query, tail);
#endif

    return 1;
}

/* execute and then reset prepared query stmt that doesn't return a result */
static int mifdbsqlite3_update (mifdbsqlite_t *db, sqlite3_stmt *stmt)
{
    if (sqlite3_step(stmt) != SQLITE_DONE)
    {
	rta_post("***sqlite query execution error %d: '%s'\n",
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

    if (status == SQLITE_ROW)
	return 1;
    else if (status != SQLITE_DONE)
	rta_post("***sqlite query execution error %d: '%s'\n",
		 sqlite3_errcode(db->sqlite3db), 
		 sqlite3_errmsg(db->sqlite3db));
    
    return 0;
}


/* open db */
int mifdb_open (mifdb_t *database, const char *dbname)
{
    mifdbsqlite_t *db = (mifdbsqlite_t *) database;

    return sqlite3_open(dbname, &db->sqlite3db) == SQLITE_OK;
    /* check errors */
}

/* close db */
int mifdb_close (mifdb_t *database)
{
    mifdbsqlite_t *db = (mifdbsqlite_t *) database;

    /* todo: clean up queries with sqlite3_finalize() */

    return sqlite3_close(db->sqlite3db);
    /* check errors */
}

int mifdb_begin_transaction (mifdb_t *database)
{
    mifdbsqlite_t *db = (mifdbsqlite_t *) database;
    const char *query = "begin;";
    char *errmsg;

    if (sqlite3_exec(db->sqlite3db, query, NULL, NULL, &errmsg) != SQLITE_OK)
    {
	rta_post("***sqlite query execution error %d: '%s'\n  query '%s'\n",
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

    if (sqlite3_exec(db->sqlite3db, query, NULL, NULL, &errmsg) != SQLITE_OK)
    {
	rta_post("***sqlite query execution error %d: '%s'\n  query '%s'\n",
		 sqlite3_errcode(db->sqlite3db), errmsg, query);
	return 0;
    }   
    return 1;
}


/* sql commands to create schema */
static const char *create = "\n\
drop table if exists discofile;						\n\
create table discofile (						\n\
  fileid	integer primary key,					\n\
  filename	text,    						\n\
  numobj	integer 						\n\
);									\n\
drop table if exists refobj;						\n\
create table refobj (							\n\
  refobjid	integer primary key,					\n\
  fileid	integer references discofile (fileid),			\n\
  objid		integer 						\n\
);									\n\
drop table if exists postinglist;					\n\
create table postinglist (						\n\
  refobjid	integer references refobj (refobjid),			\n\
  binindex	integer, 						\n\
  size		integer, /* number of entries */			\n\
  entries	blob,    /* (compressed) array of (fileid, index) */	\n\
  unique (refobjid, binindex)						\n\
);									\n\
drop table if exists indexparams;					\n\
create table indexparams (						\n\
  name		text,							\n\
  nrefobj	integer, 						\n\
  ki		integer, 						\n\
  ndim		integer, 						\n\
  descrid	integer 						\n\
);									\n\
";

/* single sql commands to insert data */
static const char *insparam = "\
insert into indexparams values (?, ?, ?, ?, ?);";

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
    sqlite3_bind_int (stmt, 2, nref);
    sqlite3_bind_int (stmt, 3, ki);
    sqlite3_bind_int (stmt, 4, ndim);
    sqlite3_bind_int (stmt, 5, descrid);

    if (!mifdbsqlite3_update(db, stmt))				return 0;
    sqlite3_finalize(stmt);

    /* compile statements for later */
    if (!mifdbsqlite3_prepare(db, insfile,   &db->insfile))	return 0;
    if (!mifdbsqlite3_prepare(db, insrefobj, &db->insrefobj))	return 0;
    if (!mifdbsqlite3_prepare(db, inspl,     &db->inspl))	return 0;
    if (!mifdbsqlite3_prepare(db, getfile,   &db->getfile))	return 0;
    if (!mifdbsqlite3_prepare(db, getrefobj, &db->getrefobj))	return 0;
    if (!mifdbsqlite3_prepare(db, getpl,     &db->getpl))	return 0;

    return 1;
}

/* get parameters from db */
int mifdb_get_params  (mifdb_t *database, int *nref, int *ki, int *ndim, int *descrid, int *nfiles, int *ndata)
{
    mifdbsqlite_t *db = (mifdbsqlite_t *) database;
    sqlite3_stmt  *stmt;

    if (!mifdbsqlite3_prepare(db, "select nrefobj, ki, ndim, descrid, count(numobj), sum(numobj) from indexparams, discofile;", &stmt))	
	return 0;

    if (!mifdbsqlite3_step(db, stmt))
	return 0;
    
    *nref    = sqlite3_column_int(stmt, 1);
    *ki      = sqlite3_column_int(stmt, 2);
    *ndim    = sqlite3_column_int(stmt, 3);
    *descrid = sqlite3_column_int(stmt, 4);
    *nfiles  = sqlite3_column_int(stmt, 5);
    *ndata   = sqlite3_column_int(stmt, 6);

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
	*index    = sqlite3_column_int(stmt, 1);
	*filename = sqlite3_column_text(stmt, 2);
	*nobj     = sqlite3_column_int(stmt, 3);
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
	*index     = sqlite3_column_int(stmt, 1);
	obj->base  = sqlite3_column_int(stmt, 2);
	obj->index = sqlite3_column_int(stmt, 3);
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
    sqlite3_bind_blob(db->inspl, 4, obj, size * sizeof(mif_object_t), SQLITE_STATIC);
 
    if (!mifdbsqlite3_update(db, db->inspl))		return 0;

    return 1;
}

/* get postinglist from open db */
int mifdb_get_postinglist (mifdb_t *database, int *index, int *binindex, int *size, int *bytes, mif_object_t **entries)
{
    mifdbsqlite_t *db = (mifdbsqlite_t *) database;
    sqlite3_stmt  *stmt = db->getpl;
    
    if (mifdbsqlite3_step(db, stmt))
    {
	*index     = sqlite3_column_int(stmt, 1);
	*binindex  = sqlite3_column_int(stmt, 2);
	*size      = sqlite3_column_int(stmt, 3);
	*entries   = (mif_object_t *) sqlite3_column_blob(stmt, 4);
	*bytes     = sqlite3_column_bytes(stmt, 4);

	return 1;
    }
    else
	return 0;
}

