#ifndef OCTTREE_STRUCTS	
#define OCTTREE_STRUCTS

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rand_tools.h"

#ifndef PPPM_STRUCTS
typedef struct {
    double x;
    double y;
    double z;
} point3d;
#endif

typedef struct {
    int depth;
    double min_x;
    double max_x;
    double min_y;
    double max_y;
    double min_z;
    double max_z;
    double range2;
    void ** nodes;
    void * parent;
    double com_x;
    double com_y;
    double com_z;
    double mass;
    double srad2;
    int * stars_contained;
    int n_stars_contained;
    int order[8];
} OctTree;

#endif
