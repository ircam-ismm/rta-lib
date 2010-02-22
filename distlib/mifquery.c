
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "mif.h"
#include "mifdb.h"
#include "discofile.h"
#include "discodist.h"


#define BLOCKSIZE 4096



static void usage()
{
    fprintf(stderr, "\
Query nquery first vectors from query-file against index in index-file, write results to result-file or stdout.\n\
Query-file is in DISCO format, index is an sqlite database generated with mifindex.\n\
\n\
usage: mifquery ks mpd  index-file-name query-file-name [nquery] [K to query]  [result-file-name]\n\
(-1 for defaults)\n\
\n\
result-file-name will be created as a text file which contains nquery lines of results: \n\
[audioID.segID of the query1] [K] [audioID1.segID1] [dist] ... [audioIDk.segIDk] [dist] \n\
[audioID.segID of the query2] [K] [audioID1.segID1] [dist] ... [audioIDk.segIDk] [dist] \n\
...\n\
[audioID.segID of the queryN] [K] [audioID1.segID1] [dist] ... [audioIDk.segIDk] [dist] \n\
\n");

    exit(1);
}


static int mifdb_read (mifdb_t *mifdb, mif_index_t *mif)
{
    int i, j, ok = 1;
    
    ok = mifdb_begin_read(mifdb);
    mif->numobj = 0;

    /* read file names */
    for (i = 0; ok  &&  i < mif->files->nbase; i++)
    {
	int ind, nobj;
	const char *name;

	if (!(ok &= mifdb_get_file(mifdb, &ind, &name, &nobj)))  break;
	mif->files->filename[ind]   = strdup(name);
	mif->files->numbaseobj[ind] = nobj;
	mif->numobj += nobj;
    }

    /* read refobj */
    for (i = 0; ok  &&  i < mif->numref; i++)
    {
	int ind;
	mif_object_t obj;

	if (!(ok &= mifdb_get_refobj(mifdb, &ind, &obj)))  break;
	mif->refobj[ind].base  = obj.base;
	mif->refobj[ind].index = obj.index;
	mif->pl[ind].size = 0;
    }

    /* init posting lists */
    for (i = 0; i < mif->numref; i++)
	mif_pl_init(&mif->pl[i], mif->ki, -1);

    /* read postinglists */
    for (i = 0; ok  &&  i < mif->numref; i++)
    {
	for (j = 0; ok  &&  j < mif->ki; j++)
	{
	    int indref, indbin, num, bytes;
	    mif_object_t *entries;

	    if ((ok &= mifdb_get_postinglist(mifdb, &indref, &indbin, &num, &bytes, &entries)))
	    {
		mif_pl_bin_t *bin = &mif->pl[indref].bin[indbin];
		    
		mif->pl[indref].size += num;
		bin->num   = num;
		bin->alloc = bytes / sizeof(mif_object_t);  /* store used space for blob (rounded to object size) */
		bin->obj   = rta_malloc(bytes);
		memcpy(bin->obj, entries, bytes);
	    }
	}
    }

    ok &= mifdb_end_read(mifdb);

    return ok;
}


int main (int argc, char *argv[])
{
    const char  *outname = NULL, *dbname = NULL, *queryname = NULL;
    int		 ks = 0, mpd = 0, nquery = -1, K = 5;
    mif_index_t  mif;
    mifdb_t	 mifdb;
    FILE	*outfile;
    disco_file_t qfile, *dbfile;
    disco_KLS_private_t kls;
    int		 i, nref, ki, nfiles, ndata, ndim, descrid;
    clock_t	 startload, startquery;
    float	 loadtime, querytime;

    /* args: nref ki ks mpd  db-disco-file-name query-disco-file-name nquery k [result-file-name] */
    switch (argc)
    {	/* fallthrough */
    case 8:	outname	  = argv[7];	  
    case 7:	K	  = atoi(argv[6]);
    case 6:	nquery	  = atoi(argv[5]);
    case 5:	queryname = argv[4];	  
    case 4:	dbname	  = argv[3];	  
    case 3:	mpd	  = atoi(argv[2]);
    case 2:	ks	  = atoi(argv[1]);
    break;

    default:	/* wrong number or no args */
    case 1:
	usage();
    break;
    }

    /* open DISCO db file, get parameters */
    if (dbname)
    {
	if (mifdb_open(&mifdb, dbname)  &&
	    mifdb_get_params(&mifdb, &nref, &ki, &ndim, &descrid, &nfiles, &ndata))
	{
	    fprintf(stderr, "opened index database file '%s'\n", dbname);
	}
	else
	{
	    fprintf(stderr, "error: can't open index database file '%s'\n", dbname);
	    return -3;
	}
    }
    else
    {
	fprintf(stderr, "error: no input index database file name given\n");
	usage();
    }

    /* open DISCO query file */
    if (queryname)
    {
	disco_file_map(&qfile, queryname, nquery);

	if (qfile.base == NULL)
	    return -2;
	else
	    fprintf(stderr, "mapped query file '%s' length %ld to pointer %p,\n ndata %d  ndim %d  descr %d\n",
		    queryname, qfile.len, qfile.base, qfile.base->ndata, qfile.base->ndim, qfile.base->descrid);
    }
    else
    {
	fprintf(stderr, "no query file name given\n");
	usage();
    }

    /* check compatibility*/
    if (ndim != qfile.base->ndim  ||  descrid != qfile.base->descrid)
    {
	fprintf(stderr, "error: index database (ndim %d, descrid %d) and query file (ndim %d, descrid %d) are incompatible!\n",
		ndim, descrid, qfile.base->ndim, qfile.base->descrid);
	return -4;
    }

    /* init index with parameters */
    mif_init(&mif, disco_KLS, disco_KLS_init, disco_KLS_free, &kls, nref, ki);

    /* read index from database */
    startload = clock();

    /* set up files list (allocate one more space for query file) */
    mifdb.files.nbase      = nfiles;
    mifdb.files.ndim       = ndim;
    mifdb.files.descrid    = descrid;
    mifdb.files.base       =           alloca((nfiles+1) * sizeof(void *));
    mifdb.files.filename   = (const char **) alloca((nfiles+1) * sizeof(char *));
    mifdb.files.numbaseobj =   (int *) alloca((nfiles+1) * sizeof(int));
    dbfile =          (disco_file_t *) alloca((nfiles+1) * sizeof(disco_file_t));
    mif.files = &mifdb.files;

    if (!mifdb_read(&mifdb, &mif))
    {
	fprintf(stderr, "error: can't read from index database file '%s'\n", dbname);
	return -5;
    }
    else
    {
	loadtime = ((float) clock() - startload) / (float) CLOCKS_PER_SEC;
	fprintf(stderr, "time for loading index = %f s, %f s / obj\n\n", 
		loadtime, loadtime / mif.numobj);
	fflush(stderr);
    }

    /* open DISCO database files */
    for (i = 0; i < nfiles; i++)
    {	
	disco_file_map(&dbfile[i], mifdb.files.filename[i], 0);
	mifdb.files.base[i] = dbfile[i].base;	/* start of mapped file */
    }

    /* add query file to files list */
    mifdb.files.base[nfiles]	   = qfile.base;	/* start of mapped file */
    mifdb.files.filename[nfiles]   = qfile.filename;	/* start of mapped file */
    mifdb.files.numbaseobj[nfiles] = qfile.base->ndata;	/* start of mapped file */

    /* open result file */
    if (outname)
	outfile = fopen(outname, "w");
    else
	outfile = stdout;

    /* check args and set for search */
    if (ks <= 0)
	ks   = nref / 4;
    if (mpd <= 0)
	mpd   = 5;

    /* init distance function for query */ 
    mif_init_index(&mif, &mifdb.files);

    mif.ks  = ks;
    mif.mpd = mpd;
    mif_print(&mif, 1, "after loading");
    mif_profile_print(&mif.profile);
    mif_profile_clear(&mif.profile);

    {   /* query with first nquery frames from query file */
	int i, j, kfound;
	mif_object_t *obj  = alloca(K * sizeof(*obj));	/* query objects */
	int          *dist = alloca(K * sizeof(*dist));	/* transformed distance to obj */
	mif_object_t  query;
	size_t	      bytesaccpl, bytesaccessed, bytesaccopt, bytesacczip;

	startquery  = clock();
	query.base = nfiles;

	if (nquery < 0)
	    nquery = qfile.base->ndata;

	for (i = 0; i < nquery; i++)
	{
	    query.index = i;
	    
	    fprintf(stderr, "querying %d-NN of query obj %d ", K, i);
	    kfound = mif_search_knn(&mif, &query, K, obj, dist);
	    fprintf(stderr, "--> found %d NN, best is db obj %d with dist %d\n", 
		    kfound, obj[0].index, dist[0]);

	    /* write result file line (query file as no. nfiles, db file with their index) */
	    fprintf(outfile, "%d.%d %d ", query.base, query.index, K);
	    for (j = 0; j < kfound; j++)
		fprintf(outfile, "%d.%d %d ", obj[j].base, obj[j].index, dist[j]);
	    fprintf(outfile, "\n");
	    fflush(outfile);
	}

	querytime = ((float) clock() - startquery) / (float) CLOCKS_PER_SEC;
	mif_print(&mif, 0, "after query");
	mif_profile_print(&mif.profile);
	rta_post("#accesses predicted:    %6d\n",
		 mif.ks * (mif.mpd * 2 + 1) * mif.numobj / mif.numref); 
	
	bytesaccpl    = mif.profile.placcess * sizeof(mif_postinglist_t);
	bytesaccessed = bytesaccpl + mif.profile.indexaccess * sizeof(mif_object_t);
	bytesaccopt   = bytesaccpl + mif.profile.indexaccess * 4;
	bytesacczip   = bytesaccpl + mif.profile.indexaccessbytes;
	fprintf(stderr, 
		"#bytes accessed in index:         %f MB = %d blocks of %d\n"
		"#bytes accessed in optimal index: %f MB = %d blocks of %d\n"
		"#bytes accessed in zipped index:  %f MB = %d blocks of %d\n", 
		bytesaccessed / 1e6, (int) ceil(bytesaccessed / BLOCKSIZE), BLOCKSIZE,
		bytesaccopt   / 1e6, (int) ceil(bytesaccopt   / BLOCKSIZE), BLOCKSIZE,
		bytesacczip   / 1e6, (int) ceil(bytesacczip   / BLOCKSIZE), BLOCKSIZE);
	fprintf(stderr, "time for %d queries = %f s, %f s / queryobj\n\n", 
		nquery, querytime, querytime / nquery);
    }

    /* cleanup and close files */
    mif_free(&mif);
    mifdb_close(&mifdb);
    disco_file_unmap(&qfile);

    for (i = 0; i < nfiles; i++)
    {
	disco_file_unmap(&dbfile[i]);
	free((char *) mifdb.files.filename[i]);
    }

    if (outfile != stdout)
	fclose(outfile);

    return 0;
}
