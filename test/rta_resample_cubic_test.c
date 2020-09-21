/* 

- compile

cc -g ../src/signal/rta_resample.c ../src/signal/rta_cubic.c rta_resample_cubic_test.c -I ../bindings/console/ -I ../src -I ../src/util/ -I ../src/signal/ -o rta_resample_cubic_test

- run

./rta_selection-test

- check

valgrind --error-limit=no ./rta_resample_cubic_test
 
*/


#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rta_configuration.h"
#include "rta_resample.h"

int main (int argc, char *argv[])
{
    for (int nsamples = 1; nsamples < 10000; nsamples *= 2)
    for (int nchannels = 1; nchannels < 8; nchannels ++)
    for (double factor = 0.05; factor <= 5; factor += 0.05)
    {
	int inframes   = nsamples;
	int outframes  = floor((double) (inframes - 1) * (1 / factor)); // shit! (inframes - 1) / factor can be rounded differently
	int insize     = inframes * nchannels;
	int outsize    = outframes * nchannels;
	float *indata  = malloc(insize  * sizeof(float));
	float *outdata = malloc(outsize * sizeof(float));

	// generate data
	for (int i = 0; i < insize; i++)
	    indata[i] = random();
	
	printf("--- channels %d: factor %4g  inframes %d  outframes %d...", nchannels, factor, inframes, outframes);
	int resize = rta_resample_cubic(outdata, indata, inframes, outframes + 1, nchannels, factor);
	printf("real %d\n", resize);

	free(indata);
	free(outdata);
	
	if (inframes > 3)
	    assert(resize == outframes);
    }
}
