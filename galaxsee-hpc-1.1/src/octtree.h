#ifndef OCTTREE
#define OCTTREE

#include "nbody.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rand_tools.h"
#include "octtree_structs.h"

OctTree * allocOctTree();
void initOctTree(OctTree * theTree);
void populateOctTree(OctTree * theTree, OctTree ** theList,
        NbodyModel * theModel,
        int depth, int * stars_contained, int n_stars_contained,
        double min_x, double max_x,
        double min_y, double max_y, double min_z, double max_z);
void resetOctTree(OctTree * theTree, NbodyModel * theModel);
void freeOctTree(OctTree * theTree,int depth);
point3d calculateForceOctTree(int minDepth, int index, double xin, double yin, double zin,
    double G, OctTree * theTree, double soft);
#endif
