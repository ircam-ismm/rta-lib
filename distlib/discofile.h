/**
@file	discofile.h
@author	Diemo Schwarz
@date	18.1.2010
@brief  DISCO descriptor file format support

Support functions for reading from DISCO format data files.

Copyright (C) 2008 - 2009 by IRCAM-Centre Georges Pompidou, Paris, France.
*/

#ifndef _DISCOFILE_H_
#define _DISCOFILE_H_

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Struct representing a DISCO descriptor data file mapped to memory:
<pre>
   header (12 bytes):
   int32 ndata		number of vectors
   int32 ndim		number of elements of each vector
   int32 descrid	descriptor ID

   data: ndata * ndim float32 values
</pre>
*/
typedef struct _disco_file_header
{
    int 	ndata;		/**< number of vectors */
    int 	ndim;		/**< number of elements of each vector */
    int 	descrid;	/**< descriptor ID */
    float	data[0];	/**< ndata * ndim float32 values */
} disco_file_header_t;


/** open and memory-mapped disco file information */
typedef struct _disco_file
{
    int		fd;		/**< file descriptor */
    size_t	len;		/**< mapped length in bytes */
    const char *filename;	/**< file name opened */
    disco_file_header_t *base;
} disco_file_t;

/**  int disco_file_ndim(disco_file_header_t *db) */
#define disco_file_ndim(file) ((file)->ndim)

/** pointer to start of frame vectors */
#define disco_file_data(file) ((file)->data)

/** open and map file to memory */
void *disco_file_map (disco_file_t *file, const char *name, int nvec);

/** unmap and close file  */
void disco_file_unmap (disco_file_t *file);


#ifdef __cplusplus
}
#endif

#endif /* _DISCOFILE_H_ */
