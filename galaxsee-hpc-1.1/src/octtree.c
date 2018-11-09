#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rand_tools.h"
#include "nbody.h"
#include "octtree_structs.h"
#include "octtree.h"

OctTree * allocOctTree() {
    OctTree * retval = malloc(sizeof(OctTree));
    initOctTree(retval);
    return retval;
}

void initOctTree(OctTree * theTree) {
    theTree->stars_contained=NULL;
    theTree->nodes=NULL;
    theTree->depth=0;
}

void populateOctTree(OctTree * theTree,  OctTree ** theList,
        NbodyModel * theModel,
        int depth, int * stars_contained, int n_stars_contained,
        double min_x, double max_x,
        double min_y, double max_y, double min_z, double max_z) {
    int i,q,count;
    double q_min_x;
    double q_max_x;
    double ave_x;
    double q_min_y;
    double q_max_y;
    double ave_y;
    double q_min_z;
    double q_max_z;
    double ave_z;
    int * work;
    OctTree * node;
    theTree->min_x=min_x;
    theTree->max_x=max_x;
    theTree->min_y=min_y;
    theTree->max_y=max_y;
    theTree->min_z=min_z;
    theTree->max_z=max_z;
    theTree->range2=theModel->treeRangeCoefficient*
        max(max_x-min_x,max(max_y-min_y,max_z-min_z));
    theTree->range2*=theTree->range2;
    theTree->depth=depth;
    theTree->n_stars_contained=n_stars_contained;
    if(theTree->stars_contained!=NULL) free(theTree->stars_contained);
    theTree->stars_contained = alloc_iarray(n_stars_contained);
    work = alloc_iarray(n_stars_contained);
    theTree->com_x=0.0;
    theTree->com_y=0.0;
    theTree->com_z=0.0;
    theTree->mass=0.0;
    for(i=0;i<n_stars_contained;i++) {
        if(stars_contained!=NULL) {
            theTree->stars_contained[i] = stars_contained[i];
        } else {
            theTree->stars_contained[i] = i;
        }
        theTree->com_x+=theModel->mass[theTree->stars_contained[i]]*
            theModel->x[theTree->stars_contained[i]];
        theTree->com_y+=theModel->mass[theTree->stars_contained[i]]*
            theModel->y[theTree->stars_contained[i]];
        theTree->com_z+=theModel->mass[theTree->stars_contained[i]]*
            theModel->z[theTree->stars_contained[i]];
        theTree->mass+=theModel->mass[theTree->stars_contained[i]];
    }
    theTree->com_x /= theTree->mass;
    theTree->com_y /= theTree->mass;
    theTree->com_z /= theTree->mass;
    theTree->srad2 = computeSoftenedRadius(
        theModel->G*theTree->mass,theModel->tstep*theModel->tstep,
        theModel->srad_factor);
    theTree->srad2*=theTree->srad2;
    if(n_stars_contained==1) {
        theTree->com_x = theModel->x[theTree->stars_contained[0]];
        theTree->com_y = theModel->y[theTree->stars_contained[0]];
        theTree->com_z = theModel->z[theTree->stars_contained[0]];
        theTree->mass = theModel->mass[theTree->stars_contained[0]];
        if(theList != NULL) 
            theList[theTree->stars_contained[0]] = theTree;
        if(theTree->nodes!=NULL) free(theTree->nodes);
        theTree->nodes = NULL;
    } else if(n_stars_contained>1) {
        if(theTree->nodes==NULL)
            theTree->nodes = (void **)malloc(sizeof(OctTree *)*8);
        ave_x = (min_x+max_x)/2.0;
        ave_y = (min_y+max_y)/2.0;
        ave_z = (min_z+max_z)/2.0;
        for(q=0;q<8;q++) {
            switch(q) {
                case 0:
                    q_min_x = ave_x;
                    q_min_y = ave_y;
                    q_min_z = ave_z;
                    q_max_x = max_x;
                    q_max_y = max_y;
                    q_max_z = max_z;
                 break;
                case 1:
                    q_min_x = min_x;
                    q_min_y = ave_y;
                    q_min_z = ave_z;
                    q_max_x = ave_x;
                    q_max_y = max_y;
                    q_max_z = max_z;
                 break;
                case 2:
                    q_min_x = min_x;
                    q_min_y = min_y;
                    q_min_z = ave_z;
                    q_max_x = ave_x;
                    q_max_y = ave_y;
                    q_max_z = max_z;
                 break;
                case 3:
                    q_min_x = ave_x;
                    q_min_y = min_y;
                    q_min_z = ave_z;
                    q_max_x = max_x;
                    q_max_y = ave_y;
                    q_max_z = max_z;
                 break;
                case 4:
                    q_min_x = ave_x;
                    q_min_y = ave_y;
                    q_min_z = min_z;
                    q_max_x = max_x;
                    q_max_y = max_y;
                    q_max_z = ave_z;
                 break;
                case 5:
                    q_min_x = min_x;
                    q_min_y = ave_y;
                    q_min_z = min_z;
                    q_max_x = ave_x;
                    q_max_y = max_y;
                    q_max_z = ave_z;
                 break;
                case 6:
                    q_min_x = min_x;
                    q_min_y = min_y;
                    q_min_z = min_z;
                    q_max_x = ave_x;
                    q_max_y = ave_y;
                    q_max_z = ave_z;
                 break;
                case 7:
                    q_min_x = ave_x;
                    q_min_y = min_y;
                    q_min_z = min_z;
                    q_max_x = max_x;
                    q_max_y = ave_y;
                    q_max_z = ave_z;
                 break;
            }
            count=0;
            for (i=0;i<n_stars_contained;i++) {
                if(theModel->x[theTree->stars_contained[i]]>=q_min_x &&
                   theModel->x[theTree->stars_contained[i]]<q_max_x &&
                   theModel->y[theTree->stars_contained[i]]>=q_min_y &&
                   theModel->y[theTree->stars_contained[i]]<q_max_y &&
                   theModel->z[theTree->stars_contained[i]]>=q_min_z &&
                   theModel->z[theTree->stars_contained[i]]<q_max_z)
                work[count++]=theTree->stars_contained[i];
            }
            if(count==0) {
                theTree->nodes[q]=NULL;
            } else {
                node = allocOctTree();
                node->parent=theTree;
                populateOctTree(node, theList, theModel,
                    depth+1,work,count,
                    q_min_x, q_max_x,
                    q_min_y, q_max_y, q_min_z, q_max_z);
                theTree->nodes[q]=node;
            }
        }
    }
    free(work);
    return;
}

void resetOctTree(OctTree * theTree, NbodyModel * theModel) {
    double sum_x, sum_y, sum_z;
    int q,i;
    if(theTree->nodes!=NULL) {
        sum_x=0.0;
        sum_y=0.0;
        sum_z=0.0;
        for(q=0;q<8;q++) {
            if(theTree->nodes[q]!=NULL) {
                resetOctTree(theTree->nodes[q],theModel);
                sum_x += ((OctTree *)theTree->nodes[q])->com_x*
                         ((OctTree *)theTree->nodes[q])->mass;
                sum_y += ((OctTree *)theTree->nodes[q])->com_y*
                         ((OctTree *)theTree->nodes[q])->mass;
                sum_z += ((OctTree *)theTree->nodes[q])->com_z*
                         ((OctTree *)theTree->nodes[q])->mass;
            }
        }
        sum_x /= theTree->mass;
        sum_y /= theTree->mass;
        sum_z /= theTree->mass;
        theTree->com_x=sum_x;
        theTree->com_y=sum_y;
        theTree->com_z=sum_z;
        sum_x=0.0;
        sum_y=0.0;
        sum_z=0.0;
        for(i=0;i<theTree->n_stars_contained;i++) {
            sum_x += theModel->mass[theTree->stars_contained[i]]*
                     theModel->x[theTree->stars_contained[i]];
            sum_y += theModel->mass[theTree->stars_contained[i]]*
                     theModel->y[theTree->stars_contained[i]];
            sum_z += theModel->mass[theTree->stars_contained[i]]*
                     theModel->z[theTree->stars_contained[i]];
        }
        theTree->com_x = sum_x/theTree->mass;
        theTree->com_y = sum_y/theTree->mass;
        theTree->com_z = sum_z/theTree->mass;
    } 
}


point3d calculateForceOctTree(int minDepth, 
    int index, double xin, double yin, double zin,
    double G, OctTree * theTree, double soft) {
    point3d accel;
    point3d q_accel;
    double dx,dy,dz,r,r2,r3,r3i;
    int q,i_order;

 
    accel.x=0.0;
    accel.y=0.0;
    accel.z=0.0;
    if (theTree->n_stars_contained==1) {
        if (theTree->stars_contained[0]==index) {
            return accel;
        } else {
            // return values based on COM
            dx = xin-theTree->com_x;
            dy = yin-theTree->com_y;
            dz = zin-theTree->com_z;
            r2 = dx*dx+dy*dy+dz*dz+soft*soft;
            if(r2<theTree->srad2) {
                return accel;
            } else {
                r = sqrt(r2);
                r3i = G/(r*r2);
                accel.x = -theTree->mass*dx*r3i;
                accel.y = -theTree->mass*dy*r3i;
                accel.z = -theTree->mass*dz*r3i;
                return accel;
            }
        }
    } else {
        // if depth too low or distance within NRAD
        dx = xin-theTree->com_x;
        dy = yin-theTree->com_y;
        dz = zin-theTree->com_z;
        r2 = dx*dx+dy*dy+dz*dz+soft*soft;
        // order of descent -> not near and not in first
        // near but not in second
        // in third
    
        if(theTree->depth<minDepth||r2<theTree->range2||
                (xin<theTree->max_x&&xin>=theTree->min_x&&
                 yin<theTree->max_y&&yin>=theTree->min_y&&
                 zin<theTree->max_z&&zin>=theTree->min_z)
                 ) {
            for (q=0;q<8;q++) {
                if(theTree->nodes[q]!=NULL) {
                    if(r2>((OctTree*)theTree->nodes[q])->range2) {
                        theTree->order[q]=1;
                    } else {
                        theTree->order[q]=2;
                    }
                }
            }
            for (i_order=1;i_order<=2;i_order++) {
                for (q=0;q<8;q++) {
                    if(theTree->nodes[q]!=NULL&&theTree->order[q]==i_order) {
                        q_accel = calculateForceOctTree(minDepth,
                            index, xin, yin, zin,
                            G, theTree->nodes[q], soft);
                        accel.x += q_accel.x;
                        accel.y += q_accel.y;
                        accel.z += q_accel.z;
                    }
                }
            }
            return accel;
        } else {
            //return value based on COM
            if(r2<theTree->srad2) {
                return accel;
            } else {
                r = sqrt(r2);
                r3 = r*r2;
                accel.x = -theTree->mass*G*dx/r3;
                accel.y = -theTree->mass*G*dy/r3;
                accel.z = -theTree->mass*G*dz/r3;
                return accel;
            }
        }
    }
}

void freeOctTree(OctTree * theTree,int depth) {
    int q;
    if(theTree==NULL) return;
    if(theTree->nodes!=NULL) {
        for(q=0;q<8;q++) {
            if(theTree->nodes[q]!=NULL) {
                freeOctTree(theTree->nodes[q],depth);
            }
        }
    }
    if(theTree->depth>=depth) {
        if(theTree->nodes!=NULL)free(theTree->nodes);
        if(theTree->stars_contained!=NULL)free(theTree->stars_contained);
        free(theTree);
    }
}

