
#include <stdio.h>
#include <math.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/fcntl.h>
#include <unistd.h>

#include "rta.h"
#include "mif.h"

/* calculate pointer to start of frame vector (skipping time) */
#define DISCO_NDIM(o) ((int *) o->base)[1]
#define DISCO_VECTOR(o) (rta_real_t *) (o->base + 4 * sizeof(int) + o->index * DISCO_NDIM(o))

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
/*    return rta_euclidean_distance(a->base + a->index * NDIM, 1,
                 b->base + b->index * NDIM, NDIM);
*/
   rta_real_t* v1 = DISCO_VECTOR(a);
   rta_real_t* v2 = DISCO_VECTOR(b);
   rta_real_t* D;
   rta_real_t dist = 0;
   int i,j, ndim = DISCO_NDIM(a);

   for (i=0;i<10;i++)
       cout<<v1[i]<<" ";

   cout<<endl;

      int r1,r2,N;
   float T1=0;
   float T2=0;
   float S1=0;
   float S2=0;
   float tmp1=0;
   float tmp2=0;

   N=  (int) ((-1 + sqrt(1 + 8 *(ndim - 1))) / 4);
   r1=N;
   r2=N*N+N;

    D =  (rta_real_t*)malloc(sizeof(rta_real_t)*N);

   for (i=0 ; i < N; i++){
       D[i]=v1[i]-v2[i];
   }

   for ( i=0 ; i < N; i++){

       tmp1=0;
       tmp2=0;

       for ( j=0 ; j < N; j++){

           T1=T1 + v2[j+i*N+r2]*v1[i+j*N+r1];
           T2=T2 + v2[j+i*N+r1]*v1[i+j*N+r2];

           tmp1= tmp1 + v2[i+j*N+r2]*D[j];
           tmp2= tmp2 + v1[i+j*N+r2]*D[j];

       }

       S1=S1+tmp1*D[i]-1;
       S2=S2+tmp2*D[i]-1;
   }


   dist=(S1+S2+T1+T2)/4 ;


   if (dist <0)
       dist=0;

   free(D);

   return dist;
}



int main (int argc, char *argv[])
{
    const char *outname = NULL, *inname = NULL;
    int   nref = 0, ks = 0, ki = 0, mpd = 0;
    FILE *outfile;
    int fd;

    int  *base = NULL;
    int   len;
    int   ndata, ndim, descrid;

    mif_index_t mif;

    /* args: nref ki ks mpd  disco-file-name [result-file-name] */
    switch (argc)
    {
    case 7:
	outname = argv[6];
    case 6:
	inname = argv[5];
    case 5:
	mpd = atoi(argv[4]);
    case 4:
	ks = atoi(argv[3]);
    case 3:
	ki = atoi(argv[2]);
    case 2:
	nref = atoi(argv[1]);
    }

    /* open DISCO file */
    if (inname)
    {
	fd = open(inname, O_RDONLY);

	if (fd <= 0)
	{
	    fprintf(stderr, "can't open file '%s'\n", inname);
	    return -3;
	}

	/* get length */
	len = lseek(fd, 0, SEEK_END);
	lseek(fd, 0, SEEK_SET);

	/* map to memory */
	base = mmap(NULL, len, PROT_READ, MAP_FILE | MAP_SHARED, fd, 0);
	
	if (base == MAP_FAILED)
	{
	    fprintf(stderr, "can't map file '%s', errno = %d\n", inname, errno);
	    return -2;
	}

	ndata = base[0];
	ndim  = base[1];
	descrid = base[2];

	fprintf(stderr, "mapped file '%s' length %d to pointer %p,\n ndata %d  ndim %d  descr %d\n",
		inname, len, base, ndata, ndim, descrid);
    }
    else
    {
	fprintf(stderr, "no input file name given\n");
	return -1;
    }

    /* open result file */
    if (outname)
	outfile = fopen(outname, "w");
    else
	outfile = stdout;

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


	mif_init(&mif, disco_euclidean_distance, nref, ki);
	mif_add_data(&mif, 1, (void **) &base, &ndata);
	mif.ks  = ks;
	mif.mpd = mpd;
	mif_print(&mif, 0);
	mif_profile_print(&mif.profile);
	mif_profile_clear(&mif.profile);
    }

    {   /* query with every frame (= model) */
	int K = 5;
	int i, j, kfound;
	mif_object_t *obj  = alloca(K * sizeof(*obj));	/* objects */
	int          *dist = alloca(K * sizeof(*dist));	/* transformed distance to obj */
	mif_object_t  query;
	query.base  = base;

	for (i = 0; i < ndata; i++)
	{
	    query.index = i;
	    
	    kfound = mif_search_knn(&mif, &query, K, obj, dist);

	    //rta_post("--> %d-NN of query obj %d (found %d):  ", K, i, kfound);

	    /* write result file line */
	    fprintf(outfile, "%d.%d %d ", 0, i, K);
	    for (j = 0; j < kfound; j++)
		fprintf(outfile, "%d.%d %d%c", 0, obj[j].index, dist[j], 
			j < kfound - 1  ?  ' ' : '\n');
	}

	mif_profile_print(&mif.profile);
    }

    /* cleanup and close files */
    mif_free(&mif);
    munmap(base, len);
    close(fd);
    if (outfile != stdout)
	fclose(outfile);

    return 0;
}
