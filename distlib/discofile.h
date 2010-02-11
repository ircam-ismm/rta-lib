
/*
 *	DISCO descriptor file format support
 */

#ifndef _DISCOFILE_H_
#define _DISCOFILE_H_

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/* header (12 bytes):
   int32 ndata		number of vectors
   int32 ndim		number of elements of each vector
   int32 descrid	descriptor ID

   data: ndata * ndim float32 values
 */

typedef struct _disco_file_header
{
    int 	ndata;		/* number of vectors */
    int 	ndim;		/* number of elements of each vector */
    int 	descrid;	/* descriptor ID */
    float	data[0];	/* ndata * ndim float32 values */
} disco_file_header_t;


/* open and memory-mapped disco file information */
typedef struct _disco_file
{
    int		fd;		/* file descriptor */
    size_t	len;		/* mapped length in bytes */
    const char *filename;	/* file name opened */
    disco_file_header_t *base;
} disco_file_t;

/*  int disco_file_ndim(disco_file_header_t *db) */
#define disco_file_ndim(file) ((file)->ndim)

/* pointer to start of frame vectors */
#define disco_file_data(file) ((file)->data)


void *disco_file_map (disco_file_t *file, const char *name, int nvec);
void disco_file_unmap (disco_file_t *file);


#ifdef __cplusplus
}
#endif

#endif /* _DISCOFILE_H_ */
