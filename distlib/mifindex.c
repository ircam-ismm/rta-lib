
#include <stdio.h>
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
Index input files, store index in sqlite database.\n\
Data files are in DISCO format.\n\
\n\
usage: mifindex nref ki database-name input-file-names...\n\
(-1 for defaults)\n\
\n");

    exit(1);
}

static int mifdb_dump (mifdb_t *mifdb, mif_index_t *mif, disco_file_t *infile)
{
    int i, j, ok;
    
    ok = mifdb_begin_transaction(mifdb);

    for (i = 0; ok  &&  i < mif->files->nbase; i++)
	ok &= mifdb_add_file(mifdb, i, infile[i].filename, mif->files->numbaseobj[i]);

    for (i = 0; ok  &&  i < mif->numref; i++)
	ok &= mifdb_add_refobj(mifdb, i, &mif->refobj[i]);

    for (i = 0; ok  &&  i < mif->numref; i++)
	for (j = 0; ok  &&  j < mif->ki; j++)
	    ok &= mifdb_add_postinglist(mifdb, i, j, 
					mif->pl[i].bin[j].num, 
					mif->pl[i].bin[j].obj);

    ok &= mifdb_commit_transaction(mifdb);

    return ok;
}


int main (int argc, char *argv[])
{
    const char	       	 *inname = NULL, *dbname = NULL;
    int		       	  nref = 0, ki = 0;
    disco_KLS_private_t	  kls;
    disco_file_t       	 *infile;
    disco_file_header_t	**base;
    int  	       	 *ndata;
    mif_index_t	       	  mif;
    mifdb_t	       	  mifdb;
    int  	       	  ntotal = 0;	/* total number of objects */
    int  	       	  nfiles, ndim = 0, descrid = 0;

    clock_t startbuild, startdump;
    float   buildtime, dumptime;
    int i;

    /* args: nref ki dbname disco-file-names... */
    if (argc < 4)
	usage();	/* wrong number or no args */
    else
    {
	nref   = atoi(argv[1]);
	ki     = atoi(argv[2]);
	dbname = argv[3];
	argc  -= 4; argv += 4;	/* eat args and executable name */
    }

    /* open database file */
    mifdb_open(&mifdb, dbname);

    nfiles = argc;
    infile = (disco_file_t *)	     alloca(nfiles * sizeof(disco_file_t));
    base   = (disco_file_header_t *) alloca(nfiles * sizeof(void *));
    ndata  = (int *)		     alloca(nfiles * sizeof(int));

    /* open DISCO input files */
    for (i = 0; i < nfiles; i++)
    {	
	inname = argv[i];
	disco_file_map(&infile[i], inname, 0);

	if (i == 0)
	{   /* first file determines dims */
	    ndim    = infile[i].base->ndim;
	    descrid = infile[i].base->descrid;
	}

	/* store base pointer and num. obj. per base in array */
	base[i]  = infile[i].base;	/* start of mapped file */
	ndata[i] = infile[i].base->ndata;
	ntotal  += infile[i].base->ndata;

	if (infile[i].base == NULL)
	{   /* todo: rollback on error */ 
	    return -3;
	}
	else
	    fprintf(stderr, "mapped database file '%s' length %ld to pointer %p,\n ndata %d  ndim %d  descr %d\n",
		    inname, infile[i].len, infile[i].base, infile[i].base->ndata, infile[i].base->ndim, infile[i].base->descrid);
    
	/* check compatibility*/
	if (ndim != infile[i].base->ndim  ||  descrid != infile[i].base->descrid)
	{
	    fprintf(stderr, "input file '%s' uses different descriptor or dimensions!\n",
		    inname);
	    /* todo: rollback on error */ 
	    return -4;
	}
    }

    /* determine parameters if not given */
    if (nref <= 0)  /* num. ref. obj. >= 2 * sqrt(nobj) */
	nref = 2 * sqrt(ntotal);
    if (ki <= 0)
	ki   = nref / 4;

    /* set up files list */
    mifdb.files.nbase    = nfiles;
    mifdb.files.ndim     = ndim;
    mifdb.files.descrid  = descrid;
    mifdb.files.base     = (void **) base;
    mifdb.files.numbaseobj = ndata;

    /* init index */
    mif_init(&mif, disco_KLS, disco_KLS_init, disco_KLS_free, &kls, nref, ki);

    /* index all input objects (= models) */
    startbuild = clock();
    mif_add_data(&mif, &mifdb.files);
    buildtime = ((float) clock() - startbuild) / (float) CLOCKS_PER_SEC;
    mif_print(&mif, 0);
    mif_profile_print(&mif.profile);
    mif_profile_clear(&mif.profile);

    fprintf(stderr, "time for building of index = %f s, %f s / obj\n\n", 
	    buildtime, buildtime / mif.numobj);
    fflush(stderr);
    
    /* dump whole index to db */
    startdump = clock();
    if (mifdb_create(&mifdb, dbname, nref, ki, ndim, descrid))
	mifdb_dump(&mifdb, &mif, infile);
    mifdb_close(&mifdb);

    dumptime = ((float) clock() - startdump) / (float) CLOCKS_PER_SEC;
    fprintf(stderr, "time for dumping index = %f s, %f s / obj\n\n", 
	    dumptime, dumptime / mif.numobj);
    fflush(stderr);

    /* cleanup and close files */
    mif_free(&mif);
    for (i = 0; i < nfiles; i++)
	disco_file_unmap(&infile[i]);

    return 0;
}
