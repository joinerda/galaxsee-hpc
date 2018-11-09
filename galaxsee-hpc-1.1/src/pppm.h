#ifndef PPPM_HEADER
#define PPPM_HEADER

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rand_tools.h"
#ifdef HAS_FFTW3
#include <fftw3.h>
#endif
#include "cubeinterp.h"
#include "nbody.h"
#include "pppm_structs.h"

void freePPPM(PPPM * thePPPM);
void allocPPPM(PPPM * thePPPM,int nx, int ny, int nz);
void pushCellPPPM(PPPM * thePPPM,int i, int j, int k, int l);
double i2rPPPM(int i,double min, double range, int n);
int r2iPPPM(double r,double min, double range, int n);
void populateDensityPPPM(PPPM * thePPPM, NbodyModel * theModel,
        double * x,
        double xmin, double xmax, double ymin, double ymax,
        double zmin, double zmax, double soft);
point3d calculateForcePPPM(int index, double xin, double yin, double zin,
    PPPM* thePPPM, NbodyModel * theModel, double * x,double near);
void setKPPPM(int n,double * k,double range);
void prepPotentialPPPM(PPPM* thePPPM, NbodyModel * theModel);
#endif
