/**********************************************
*  Copyright 2007-2008, David Joiner
*
*  Random tools for use with monte carlo methods
*
**********************************************/
#include "rand_tools.h"

void seed_by_time(int offset) {
    time_t the_time;
    time(&the_time);
    srand((int)the_time+offset);
}

double drand_norm(double xbar, double sigma, double alpha) {
    double omalpha = 1.0-alpha;
    return xbar + sigma*inverf(drand(-omalpha,omalpha));
}

double inverf(double y) {
    // inverse error function computed using Newton's method
    double f,fprime,x;

    double con = 1.12837916709551;

    if(y==0) return 0.5;

    x = 0;
    f = erf(x)-y;
    while(fabs(f)>0.0001) {
        fprime = con*exp(-x*x);
        x -= f/fprime;
        f = erf(x)-y;
    }

    return x;
}

int irand(int min, int max) {
    return min+(int)((double)rand()/(double)RAND_MAX*(double)(max-min));
}

double drand(double min,double max) {
    return min+((double)rand()/(double)RAND_MAX)*(max-min);
}

void random_direction_subset(int n, int m, double *x) {
    int * deleted;
    int n_deleted=0;
    int direction,exists,i;
    double sum;

    deleted = alloc_iarray(n);
    random_direction(n,x);
    //eliminate n-m directions
    while(n_deleted<n-m) {
        direction = irand(0,n-1);
        exists=0;
        for(i=0;i<n_deleted&&!exists;i++) {
            if(deleted[i]==direction) exists=1;
        }
        if(!exists) {
            deleted[n_deleted]=direction;
            n_deleted++;
            x[direction]=0.0;
        }
    }
    //renormalize
    sum=0.0;
    for(i=0;i<n;i++) {
        sum+=x[i]*x[i];
    }
    sum=1.0/sqrt(sum);
    for (i=0;i<n;i++) {
        x[i]*=sum;
    }
    
    free(deleted);
}

void random_direction(int n, double * x) {
    double length;
    int i;
    do {
        for (i=0;i<n;i++) {
            x[i]=drand(-1,1);
        }
        length = length_darray(n,x);
        if(length!=0.0) {
            for (i=0;i<n;i++) {
                x[i]/=length;
            }
        }
    } while (length==0.0);
}

