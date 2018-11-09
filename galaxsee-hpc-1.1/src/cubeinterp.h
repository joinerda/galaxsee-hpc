#ifndef CUBEINT
#define CUBEINT
#include <stdio.h>

typedef struct {
    double f000,f001,f010,f011,f100,f101,f110,f111;
    double a,b,c,d,e,f,g,h;
    double xmin,xmax,ymin,ymax,zmin,zmax;
} CubeInterp;

void setCornersCINT(CubeInterp * cint,
        double f000, double f100, double f010, double f110,
        double f001, double f101, double f011, double f111,
        double xmin, double xmax, double ymin, double ymax,
        double zmin, double zmax);
double getValueCINT(CubeInterp cint, double xin, double yin, double zin);
#endif
