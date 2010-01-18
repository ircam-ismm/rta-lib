
#include <stdio.h>
#include <math.h>
#include <time.h>

#include "mif.h"
#include "discofile.h"


#define BLOCKSIZE 4096


/* calculate pointer to start of frame vector (skipping time) */
#define DISCO_NDIM(o) ((int *) o->base)[1]
#define DISCO_VECTOR(o) (rta_real_t *) ((float *) (o->base + 3 * sizeof(int)) + 1 + o->index * DISCO_NDIM(o))


static rta_real_t disco_euclidean_distance (mif_object_t *a, mif_object_t *b)
{
/*    return rta_euclidean_distance(a->base + a->index * NDIM, 1,
				  b->base + b->index * NDIM, NDIM);
*/
    rta_real_t* v1 = DISCO_VECTOR(a);
    rta_real_t* v2 = DISCO_VECTOR(b);
    rta_real_t sum = 0;
    int i, ndim = DISCO_NDIM(a);

    for (i = 0; i < ndim; i++) 
    {
	rta_real_t diff = v2[i] - v1[i];
	sum += diff * diff;
    }

    return sum;
}


static rta_real_t disco_KLS (mif_object_t *a, mif_object_t *b)
{
   rta_real_t* v1 = DISCO_VECTOR(a);
   rta_real_t* v2 = DISCO_VECTOR(b);
   rta_real_t* D;
   rta_real_t dist = 0;
   int i,j, ndim = DISCO_NDIM(a);

/*   for (i=0; i < 10; i++)
       rta_post("%f ", v1[i]);
   rta_post("\n");
*/
   int r1,r2,N;
   double T1=0;
   double T2=0;
   double S1=0;
   double S2=0;
   double tmp1=0;
   double tmp2=0;

   N=  (int) ((-1 + sqrt(1 + 8 *(ndim - 1))) / 4);
   r1=N;
   r2=N*N+N;

    D =  (rta_real_t*)alloca(sizeof(rta_real_t)*N);

   for (i=0 ; i < N; i++){
       D[i]=v1[i]-v2[i];
   }

   for (i=0 ; i < N; i++){

       tmp1=0;
       tmp2=0;

       for ( j=0 ; j < N; j++){

           T1 += v2[j+i*N+r2]*v1[i+j*N+r1];
           T2 += v2[j+i*N+r1]*v1[i+j*N+r2];

           tmp1 += v2[i+j*N+r2]*D[j];
           tmp2 += v1[i+j*N+r2]*D[j];
       }

       S1 += tmp1*D[i]-1;
       S2 += tmp2*D[i]-1;
   }

   dist=(S1+S2+T1+T2)/4 ;

   if (dist <0)
       dist=0;

   return dist;
}


static void usage()
{
    fprintf(stderr, "\
Query nquery first vectors from query-file against db-file, write results to result-file or stdout.\n\
Data files are in DISCO format.\n\
\n\
usage: mifquery nref ki ks mpd  db-file-name query-file-name [nquery] [K to query]  [result-file-name]\n\
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


int main (int argc, char *argv[])
{
    const char *outname = NULL, *dbname = NULL, *queryname = NULL;
    int   nref = 0, ks = 0, ki = 0, mpd = 0, nquery = -1, K = 5;
    FILE *outfile;
    disco_file_t dbfile, qfile;
    mif_index_t mif;
    int  *base = NULL;
    int   ndata, ndim;
    clock_t startbuild, startquery;
    float   buildtime, querytime;

    /* args: nref ki ks mpd  db-disco-file-name query-disco-file-name nquery k [result-file-name] */
    switch (argc)
    {	/* fallthrough */
    case 10:	outname	  = argv[9];
    case 9:	K	  = atoi(argv[8]);
    case 8:	nquery	  = atoi(argv[7]);
    case 7:	queryname = argv[6];
    case 6:	dbname	  = argv[5];
    case 5:	mpd	  = atoi(argv[4]);
    case 4:	ks	  = atoi(argv[3]);
    case 3:	ki	  = atoi(argv[2]);
    case 2:	nref	  = atoi(argv[1]);
    break;

    default:	/* wrong number or no args */
    case 1:
	usage();
    break;
    }

    /* open DISCO db file */
    if (dbname)
    {
	disco_file_map(&dbfile, dbname, 0);

	if (dbfile.base == NULL)
	    return -3;
	else
	    fprintf(stderr, "mapped database file '%s' length %ld to pointer %p,\n ndata %d  ndim %d  descr %d\n",
		    dbname, dbfile.len, dbfile.base, dbfile.base->ndata, dbfile.base->ndim, dbfile.base->descrid);
    }
    else
    {
	fprintf(stderr, "no input database file name given\n");
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

    /* open result file */
    if (outname)
	outfile = fopen(outname, "w");
    else
	outfile = stdout;

    base  = dbfile.base;
    ndim  = dbfile.base->ndim;
    ndata = dbfile.base->ndata;

    /* check compatibility*/
    if (ndim != qfile.base->ndim  ||  dbfile.base->descrid != qfile.base->descrid)
    {
	fprintf(stderr, "db and query files are incompatible!\n");
	return -4;
    }


#if 0
    { /* test distance function: voila ce que l'on est sensé obtenir:
	D(O1,O2) = 0.9310
	D(O1,O3) = 1.6384
	D(O2,O3) =  0.2559
	où O1 O2 et O3 sont les trois premiers objets de la base... */
	int i, j;

	for (i = 0; i < 2; i++)
	    for (j = i + 1; j < 3; j++)
	    {
		mif_object_t a = { (void *) base, i };
		mif_object_t b = { (void *) base, j };
		rta_real_t d = disco_KLS(&a, &b);
		fprintf(stderr, "D(o%d, o%d) = %f\n", i+1, j+1, d);
	    }
    }
#endif

    {   /* index every frame (= model) */
	/* check args */
	if (nref <= 0)  /* num. ref. obj. >= 2 * sqrt(nobj) */
	    nref = 2 * sqrt(ndata);
	if (ki <= 0)
	    ki   = nref / 4;
	if (ks <= 0)
	    ks   = nref / 4;
	if (mpd <= 0)
	    mpd   = 5;

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

    {   /* query with first nquery frames from query file */
	int i, j, kfound;
	mif_object_t *obj  = alloca(K * sizeof(*obj));	/* query objects */
	int          *dist = alloca(K * sizeof(*dist));	/* transformed distance to obj */
	mif_object_t  query;
	size_t	      bytesaccessed, bytesaccoopt;

	query.base  = qfile.base;
	startquery  = clock();

	if (nquery < 0)
	    nquery = qfile.base->ndata;

	for (i = 0; i < nquery; i++)
	{
	    query.index = i;
	    
	    fprintf(stderr, "querying %d-NN of query obj %d ", K, i);
	    kfound = mif_search_knn(&mif, &query, K, obj, dist);
	    fprintf(stderr, "--> found %d NN, best is db obj %d with dist %d\n", 
		    kfound, obj[0].index, dist[0]);

	    /* write result file line (outfile as index 1, db file as index 0) */
	    fprintf(outfile, "%d.%d %d ", 1, i, K);
	    for (j = 0; j < kfound; j++)
		fprintf(outfile, "%d.%d %d%c", 0, obj[j].index, dist[j], 
			j < kfound - 1  ?  ' ' : '\n');
	    fflush(outfile);
	}

	querytime = ((float) clock() - startquery) / (float) CLOCKS_PER_SEC;
	mif_profile_print(&mif.profile);
	
	bytesaccessed = mif.profile.placcess * sizeof(mif_postinglist_t) +
			mif.profile.indexaccess * sizeof(mif_pl_entry_t);
	bytesaccoopt  = mif.profile.placcess * sizeof(mif_postinglist_t) +
			mif.profile.indexaccess * 4;
	fprintf(stderr, 
		"#bytes accessed in index:         %f MB = %d blocks of %d\n"
		"#bytes accessed in optimal index: %f MB = %d blocks of %d\n", 
		bytesaccessed / 1e6, (int) ceil(bytesaccessed / BLOCKSIZE), BLOCKSIZE,
		bytesaccoopt  / 1e6, (int) ceil(bytesaccoopt  / BLOCKSIZE), BLOCKSIZE);
	fprintf(stderr, "time for %d queries = %f s, %f s / queryobj\n\n", 
		nquery, querytime, querytime / nquery);

    }

    /* cleanup and close files */
    mif_free(&mif);
    disco_file_unmap(&dbfile);
    disco_file_unmap(&qfile);

    if (outfile != stdout)
	fclose(outfile);

    return 0;
}
