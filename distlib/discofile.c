
#include <unistd.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/fcntl.h>

#include "discofile.h"

void *disco_file_map (disco_file_t *file, const char *name, int nvec)
{
    file->filename = name;
    file->base = NULL;
    file->fd = open(name, O_RDONLY);

    if (file->fd <= 0)
    {
	fprintf(stderr, "can't open file '%s'\n", name);
	return NULL;
    }

    /* get length (LATER: map only nvec records if >= 0) */
    file->len = lseek(file->fd, 0, SEEK_END);
    lseek(file->fd, 0, SEEK_SET);

    /* map to memory */
    file->base = mmap(NULL, file->len, PROT_READ, MAP_FILE | MAP_SHARED, file->fd, 0);
	
    if (file->base == MAP_FAILED)
    {
	fprintf(stderr, "can't map file '%s', errno = %d\n", name, errno);
	return NULL;
    }

    return (void *) file->base;
}

void disco_file_unmap (disco_file_t *file)
{
    munmap(file->base, file->len);
    close(file->fd);
}
