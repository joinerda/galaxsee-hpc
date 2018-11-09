#ifndef MEM_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#undef min
#undef max
#define min(X,Y) ((X) < (Y) ? (X) : (Y))
#define max(X,Y) ((X) > (Y) ? (X) : (Y))

double length_darray(int n,double *x);
double length_farray(int n,float *x);
double * alloc_darray(int n);
double ** alloc_dmatrix(int n,int m);
void free_dmatrix(double ** dmatrix);
void free_darray(double * darray);
float * alloc_farray(int n);
float ** alloc_fmatrix(int n,int m);
void free_fmatrix(float ** fmatrix);
void free_farray(float * farray);
int * alloc_iarray(int n);
int ** alloc_imatrix(int n,int m);
void free_imatrix(int ** imatrix);
void free_iarray(int * iarray);
void print_darray(int n, double * x, const char * string);

#define MEM_H
#endif

