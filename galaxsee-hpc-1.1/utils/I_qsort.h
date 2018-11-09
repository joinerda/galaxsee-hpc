#ifndef IQSORT
#define IQSORT

#include <stdlib.h>
#include <math.h>
#include <stdio.h>

int * I_qix;
void * I_qval;
size_t I_qsize;
int  (*I_qcompfunc)(const void *,const void *);

void I_qsort(void * values, int * index,
               int n, size_t size, int comp(const void *,const void *));


#endif
