#ifndef PPPM_STRUCTS
#define PPPM_STRUCTS


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rand_tools.h"
#ifdef HAS_FFTW3
#include <fftw3.h>
#endif
#include "cubeinterp.h"
#include "nbody.h"

#ifndef OCTTREE_STRUCTS
typedef struct {
    double x;
    double y;
    double z;
} point3d;
#endif

typedef struct {
    int nx,ny,nz,nxy,nxyz,nyz,nzo2,nyzo2,nxyzo2;
    double xmin,xmax,ymin,ymax,zmin,zmax;
    double * x;
    double * y;
    double * z;
    double * kx;
    double * ky;
    double * kz;
    int * n_cell;
    int * max_cell;
    double * density;
    double * dbuffer;
    double * potential;
    double * den_check;
    point3d * force;
#ifdef HAS_FFTW3
    fftw_complex * fft_density;
    fftw_complex * fft_potential;
    fftw_plan pf,pb;
#endif
    int ** cell_contains;
} PPPM;

#endif
