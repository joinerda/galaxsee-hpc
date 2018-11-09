#ifndef RANDTOOLS
#define RANDTOOLS
#include <time.h>
#include <stdlib.h>
#include "mem.h"

void seed_by_time(int offset);
double drand(double min,double max);
double drand_norm(double xbar,double sigma,double alpha);
double inverf(double);
int irand(int min,int max);
void random_direction(int n, double * x);
void random_direction_subset(int n, int m, double * x);

#endif
