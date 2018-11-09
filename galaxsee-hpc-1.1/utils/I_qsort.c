#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "I_qsort.h"

int I_qcomp(const void * a, const void * b) {
    int ia = *(int *)a;
    int ib = *(int *)b;
    int retval = I_qcompfunc(I_qval+I_qsize*ia,I_qval+I_qsize*ib);
    return retval;
}

void I_qsort(void * values, int * index,
               int n, size_t size, int comp(const void *,const void *)) {
    int i;
    I_qix = index;
    I_qval = values;
    I_qcompfunc = comp;
    I_qsize = size;
    for(i=0;i<n;i++) index[i]=i;
    qsort((void *)I_qix,n,sizeof(int),I_qcomp);
    return;
}


