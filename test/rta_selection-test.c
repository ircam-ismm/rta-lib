/* 

- compile

cc -g ../src/statistics/rta_selection.c rta_selection-test.c -I ../bindings/console/ -I ../src -I ../src/util/ -I ../src/statistics/ -o rta_selection-test

- run

./rta_selection-test

- check

valgrind --leak-check=yes --track-origins=yes --error-limit=no ./rta_selection-test
 
*/


#include <stdio.h>
#include <stdlib.h>
#include "rta_configuration.h"
#include "rta_selection.h"

int main (int argc, char *argv[])
{
    int dim     = 10;
    int trials = 100;
    
    for (int run = 0; run < trials; run++)
    {
	int numdata = 1000 + run;
	float *tempdata = malloc(numdata * dim * sizeof(float));

	// generate data
	for (int i = 0; i < numdata * dim; i++)
	    tempdata[i] = random();

	
	for (int k = 0; k < dim; k++)
	{
	    float median = rta_selection_stride(tempdata + k, dim, numdata, 0.5 * (float) (numdata - 1));
	    
	    printf(" %g", median);
	}
	free(tempdata);
	
	printf(" --- run %d\n", run);
    }
}
