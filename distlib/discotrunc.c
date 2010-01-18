
#include <stdio.h>
#include <stdlib.h>
#include "discofile.h"


static void usage()
{
    fprintf(stderr, "\
Write nvec first vectors from input-file to output-file or stdout.\n\
All files are in DISCO format.\n\
\n\
usage: discotrunc nvec in-file-name [out-file-name]\n\
\n");

    exit(1);
}

int main (int argc, char *argv[])
{
    const char *inname = NULL, *outname = NULL;
    int   nvec = 0, nwritten, towrite;
    FILE *outfile;
    disco_file_t infile;

    /* args: nvec disco-file-name out-file-name */
    switch (argc)
    {	/* fallthrough */
    case 4:	outname = argv[3];
    case 3:	inname	= argv[2];
    case 2:	nvec	= atoi(argv[1]);
    break;

    default:	/* wrong number or no args */
    case 1:
	usage();
    break;
    }

    /* open DISCO file */
    if (inname)
    {
	disco_file_map(&infile, inname, 0);

	if (infile.base == NULL)
	    return -3;
	else
	    fprintf(stderr, "mapped database file '%s' length %ld to pointer %p,\n ndata %d  ndim %d  descr %d\n",
		    inname, infile.len, infile.base, infile.base->ndata, infile.base->ndim, infile.base->descrid);
    }
    else
    {
	fprintf(stderr, "no input disco file name given\n");
	usage();
    }

    /* open result file */
    if (outname)
	outfile = fopen(outname, "w");
    else
	outfile = stdout;

    if (nvec > infile.base->ndata)
	nvec = infile.base->ndata;

    /* write new header and nvec vectors to outfile */
    towrite   = sizeof(disco_file_header_t) + infile.base->ndim * sizeof(float) * nvec;
    nwritten  = fwrite(&nvec, 1, sizeof(int), outfile);
    nwritten += fwrite((void *) infile.base + sizeof(int), 1, 
		       towrite - sizeof(int), outfile);
    if (nwritten != towrite)
	fprintf(stderr, "write error: %d bytes written to %s instead of expected %d bytes\n", 
		nwritten, outname ? outname : "stdout", towrite);
    
    /* cleanup and close files */
    disco_file_unmap(&infile);
    if (outfile != stdout)
	fclose(outfile);

    return 0;
}
