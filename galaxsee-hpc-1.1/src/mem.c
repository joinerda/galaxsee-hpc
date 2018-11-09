/**********************************************
*  Copyright 2007-2008, David Joiner
*
*  Memory management routines
*
**********************************************/
#include "mem.h"

double length_darray(int n, double *x) {
    int i;
    double sum;

    sum=0.0;
    for(i=0;i<n;i++) {
        sum+=x[i]*x[i];
    }
    return sqrt(sum);
}

double length_farray(int n, float *x) {
    int i;
    float sum;

    sum=0.0;
    for(i=0;i<n;i++) {
        sum+=x[i]*x[i];
    }
    return sqrt(sum);
}

double * alloc_darray(int n) {
    int i;
    double *darray = (double *) malloc(sizeof(double)*n);
    if(darray==NULL) {
        printf("WARNING: failed to allocate memory in alloc_array!\n");
    }
    for (i=0;i<n;i++) {
        darray[i]=0.0;
    }
    return darray;
}

double ** alloc_dmatrix(int n,int m) {
    int i,j;
    double ** dmatrix;
    double * large_array;
    large_array=alloc_darray(n*m);
    dmatrix = (double **) malloc(sizeof(double*)*n); 
    dmatrix[0] = &large_array[0];
    for(i=1;i<n;i++) {
        //dmatrix[i] = &large_array[i*m];
        dmatrix[i] = dmatrix[i-1]+m;
    }
    for(i=0;i<n;i++) {
        for(j=0;j<m;j++) {
            dmatrix[i][j]=0.0;
        }
    }
    return dmatrix;
}

void free_dmatrix(double ** dmatrix) {
    free(*dmatrix);
    free(dmatrix);
}

void free_darray(double * darray) {
    if(darray!=NULL) free(darray);
}

float * alloc_farray(int n) {
    int i;
    float *farray = (float *) malloc(sizeof(float)*n);
    for (i=0;i<n;i++) {
        farray[i]=0.0;
    }
    return farray;
}

float ** alloc_fmatrix(int n,int m) {
    int i,j;
    float ** fmatrix;
    float * large_array;
    large_array=alloc_farray(n*m);
    fmatrix = (float **) malloc(sizeof(float*)*n); 
    fmatrix[0] = &large_array[0];
    for(i=1;i<n;i++) {
        fmatrix[i] = &large_array[i*m];
    }
    for(i=0;i<n;i++) {
        for(j=0;j<m;j++) {
            fmatrix[i][j]=0.0;
        }
    }
    return fmatrix;
}

void free_fmatrix(float ** fmatrix) {
    free(*fmatrix);
    free(fmatrix);
}

void free_farray(float * farray) {
    free(farray);
}

int * alloc_iarray(int n) {
    int i;
    int *iarray = (int *) malloc(sizeof(int)*n);
    for (i=0;i<n;i++) {
        iarray[i]=0.0;
    }
    return iarray;
}

int ** alloc_imatrix(int n,int m) {
    int i,j;
    int ** imatrix;
    int * large_array;
    large_array=alloc_iarray(n*m);
    imatrix = (int **) malloc(sizeof(int*)*n); 
    imatrix[0] = &large_array[0];
    for(i=1;i<n;i++) {
        imatrix[i] = &large_array[i*m];
    }
    for(i=0;i<n;i++) {
        for(j=0;j<m;j++) {
            imatrix[i][j]=0;
        }
    }
    return imatrix;
}

void free_imatrix(int ** imatrix) {
    free(*imatrix);
    free(imatrix);
}

void free_iarray(int * iarray) {
    free(iarray);
}

void print_darray(int n, double * x, const char * string) {
    int i;
    printf("%s",string);
    for(i=0;i<n;i++) printf("[%d] %10.3e\t",i,x[i]);
    printf("\n");
}
