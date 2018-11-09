#ifndef TEXT11
#define TEXT11

#include <stdio.h>
#include <stdlib.h>
#include "mem.h"

typedef struct {
    int nx,ny,nb;
    int ** pixels;
    int * boundary;
    char * display;
    double dx_min,dx_max,dy_min,dy_max;
} Text11;

void text11_reset(Text11 * t11p);
void text11_free(Text11 **t11Display);
void text11_initialize(Text11 ** t11Display);
void text11_print(Text11 * t11p);
int  text11_add(Text11 * t11p,double x,double y);
void text11_set_boundary(Text11 * t11p,int max_range);
void text11_set_range(Text11 *t11p,double dx_min,double dy_min,
        double dx_max,double dy_max);

#endif

