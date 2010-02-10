
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "mif.h"
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


int main (int argc, char *argv[])
{
    const char *inname = NULL, *dbname = NULL;
    int   nref = 0, ki = 0;
    disco_file_t dbfile, qfile;
    mif_index_t mif;
    int  *base = NULL;
    int   ndata, ndim, descrid;
    clock_t startbuild, startquery;
    float   buildtime, querytime;

    /* args: nref ki ks mpd  db-disco-file-name query-disco-file-name nquery k [result-file-name] */
    if (argc < 4)
	usage();	/* wrong number or no args */
    else
    {
	nref = atoi(argv[1]);
	ki   = atoi(argv[2]);
	dbname	  = argv[3];
	argc -= 4; argv += 4;
    }

    /* open database file */



    for (i = 0; i < argc; i++)
    {	/* open DISCO input file */
	inname = argv[i];
	disco_file_map(&infile, inname, 0);

	if (i == 0)
	{
	    ndim  = infile.base->ndim;
	    ndata = infile.base->ndata;
	    descr = infile.base->descrid;
	}
	base  = infile.base;

	if (infile.base == NULL)
	    return -3;
	else
	    fprintf(stderr, "mapped database file '%s' length %ld to pointer %p,\n ndata %d  ndim %d  descr %d\n",
		    inname, infile.len, infile.base, infile.base->ndata, infile.base->ndim, infile.base->descrid);
    

	/* check compatibility*/
	if (ndim != infile.base->ndim  ||  descrid != infile.base->descrid)
	{
	    fprintf(stderr, "input file '%s' uses different descriptor or dimensions!\n",
		    inname);
	    return -4;
	}



	{   /* index every frame (= model) */
	    /* check args */
	    if (nref <= 0)  /* num. ref. obj. >= 2 * sqrt(nobj) */
		nref = 2 * sqrt(ndata);
	    if (ki <= 0)
		ki   = nref / 4;

	    mif_init(&mif, disco_KLS, nref, ki);
	    startbuild = clock();
	    mif_add_data(&mif, 1, (void **) &base, &ndata);
	    buildtime = ((float) clock() - startbuild) / (float) CLOCKS_PER_SEC;
	    mif.ks  = ks;
	    mif.mpd = mpd;
	    mif_print(&mif, 0);
	    mif_profile_print(&mif.profile);
	    mif_profile_clear(&mif.profile);

	    fprintf(stderr, "time for building of index = %f s, %f s / obj\n\n", 
		    buildtime, buildtime / mif.numobj);
	    fflush(stderr);
	}
    
	disco_file_unmap(&infile);
    }

    /* cleanup and close database */
    mif_free(&mif);

    return 0;
}
