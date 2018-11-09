#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "fcr.h"

#ifndef NAN
#define NAN 0.0/0.0;
#endif

double fcr_guess(double x) {
    if(x>1.0) {
        if(x<1.0e1) return 2.0e0;
        else if(x<1.0e2) return 5.0e0;
        else if(x<1.0e3) return 1.0e1;
        else if(x<1.0e4) return 2.0e1;
        else if(x<1.0e5) return 5.0e1;
        else if(x<1.0e6) return 1.0e2;
        else if(x<1.0e7) return 2.0e2;
        else if(x<1.0e8) return 5.0e2;
        else if(x<1.0e9) return 1.0e3;
        else if(x<1.0e10) return 2.0e3;
        else if(x<1.0e11) return 5.0e3;
        else if(x<1.0e12) return 1.0e4;
        else if(x<1.0e15) return 1.0e5;
        else if(x<1.0e18) return 1.0e6;
        else if(x<1.0e21) return 1.0e7;
        else if(x<1.0e24) return 1.0e8;
        else if(x<1.0e27) return 1.0e9;
        else return 1.0e10;
    } else if(x<1.0) {
        if(x>1.0e-1) return 5.0e-1;
        else if(x>1.0e-2) return 2.0e-1;
        else if(x>1.0e-3) return 1.0e-1;
        else if(x>1.0e-4) return 5.0e-2;
        else if(x>1.0e-5) return 2.0e-2;
        else if(x>1.0e-6) return 1.0e-2;
        else if(x>1.0e-7) return 5.0e-3;
        else if(x>1.0e-8) return 2.0e-3;
        else if(x>1.0e-9) return 1.0e-3;
        else if(x>1.0e-10) return 5.0e-4;
        else if(x>1.0e-11) return 2.0e-4;
        else if(x>1.0e-12) return 1.0e-4;
        else if(x>1.0e-15) return 1.0e-5;
        else if(x>1.0e-18) return 1.0e-6;
        else if(x>1.0e-21) return 1.0e-7;
        else if(x>1.0e-24) return 1.0e-8;
        else if(x>1.0e-27) return 1.0e-9;
        else return 1.0e-10;
    } else if(x==1.0) {
        return 1.0;
    } else {
        return NAN;
    }
}

double fcr(double x, int MAX_ITER) {
    double a,b;
    int done=false;
    int count=0;
    int sign=0;

    if(x==0.0) return 0.0;
    if(x<0.0) { x=-x; sign=1; }
    if(x==1.0) {done=true;}
    a=fcr_guess(x);
    b=a;
    
    while(!done&&count++<MAX_ITER) {
        b = a*(((a*a*a+x)+x)/(a*a*a+(a*a*a+x)));
        if(fabs((b-a)/a)<FCR_EPS) done=true;
        else if(fabs(b-a)<FCR_TINY) done=true;
        a=b;
    }
    if(sign) return -b;
    else return b;
}

