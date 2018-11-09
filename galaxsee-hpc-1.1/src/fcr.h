#ifndef FCR

// see http://en.wikipedia.org/wiki/Cube_root

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define FCR_EPS 1.0e-6
#define FCR_TINY 1.0e-50

#ifndef false
#define false 0
#endif
#ifndef true
#define true 1
#endif

double fcr_guess(double x);
double fcr(double x, int MAX_ITER);

#endif
