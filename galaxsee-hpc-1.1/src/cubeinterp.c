#include <stdio.h>
#include "cubeinterp.h"

void setCornersCINT(CubeInterp * cint,
        double f000, double f100, double f010, double f110,
        double f001, double f101, double f011, double f111,
        double xmin, double xmax, double ymin, double ymax,
        double zmin, double zmax) {
    cint->f000=f000;
    cint->f001=f001;
    cint->f010=f010;
    cint->f011=f011;
    cint->f100=f100;
    cint->f101=f101;
    cint->f110=f110;
    cint->f111=f111;
    cint->xmin=xmin;
    cint->xmax=xmax;
    cint->ymin=ymin;
    cint->ymax=ymax;
    cint->zmin=zmin;
    cint->zmax=zmax;

    cint->a=f000;
    cint->b=f100-cint->a;
    cint->c=f010-cint->a;
    cint->d=f001-cint->a;
    cint->e=f110-cint->a-cint->b-cint->c;
    cint->f=f011-cint->a-cint->b-cint->d;
    cint->g=f101-cint->a-cint->c-cint->d;
    cint->h=f111-cint->a-cint->b-cint->c-cint->d-cint->e-cint->f-cint->g;
}

double getValueCINT(CubeInterp cint, double xin, double yin, double zin) {
    double x,y,z;
    x = (xin-cint.xmin)/(cint.xmax-cint.xmin);
    y = (yin-cint.ymin)/(cint.ymax-cint.ymin);
    z = (zin-cint.zmin)/(cint.zmax-cint.zmin);

    return cint.a+cint.b*x+cint.c*y+cint.d*z+
        cint.e*x*y+cint.f*y*z+cint.g*z*x+cint.h*x*y*z;
}

