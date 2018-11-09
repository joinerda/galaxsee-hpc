#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "text11.h"
#include "rand_tools.h"
#include "fcr.h"
#include "octtree_structs.h"
#include "octtree.h"
#ifdef USE_PPPM
#include "pppm_structs.h"
#include "pppm.h"
#endif
#include "nbody.h"
#include <time.h>
#include "timer.h"
#ifdef HAS_NETCDF
#include <netcdf.h>
#endif

#ifdef HAS_LIBGD
    #include <gd.h>
    #include <gdfontl.h>
#endif

#ifdef HAS_MPI
    #include <mpi.h>
#endif

#ifdef HAS_X11
    #include <X11/Xlib.h>
    #include <assert.h>
    #include <unistd.h>
    #define NIL (0)
#endif

#ifdef HAS_SDL
#include "sdlwindow.h"
#include <SDL/SDL.h>
#include <SDL/SDL_opengl.h>
#include "quaternion.h"
#endif


NbodyModel * allocateNbodyModel(int n, int ngrid) {
    NbodyModel * theModel;
    theModel = (NbodyModel *)malloc(sizeof(NbodyModel));
    theModel->n=n;
    theModel->mass = (double *)malloc(sizeof(double)*n);
    theModel->x = (double *)malloc(sizeof(double)*n);
    theModel->y = (double *)malloc(sizeof(double)*n);
    theModel->z = (double *)malloc(sizeof(double)*n);
    theModel->draw_x = (double *)malloc(sizeof(double)*n);
    theModel->draw_y = (double *)malloc(sizeof(double)*n);
    theModel->draw_z = (double *)malloc(sizeof(double)*n);
    theModel->srad2 = (double *)malloc(sizeof(double)*n);
    theModel->vx = (double *)malloc(sizeof(double)*n);
    theModel->vy = (double *)malloc(sizeof(double)*n);
    theModel->vz = (double *)malloc(sizeof(double)*n);
    theModel->X = (double *)malloc(sizeof(double)*n*6);
    theModel->dbuffer = (double *)malloc(sizeof(double)*n*6);
    theModel->X2 = (double *)malloc(sizeof(double)*n*6);
    theModel->color = (int *)malloc(sizeof(int)*n);
    theModel->color3 = (double *)malloc(sizeof(double)*n*3);
    theModel->XPRIME = (double *)malloc(sizeof(double)*n*6);
    theModel->XPRIME1 = (double *)malloc(sizeof(double)*n*6);
    theModel->XPRIME2 = (double *)malloc(sizeof(double)*n*6);
    theModel->XPRIME3 = (double *)malloc(sizeof(double)*n*6);
    theModel->XPRIME4 = (double *)malloc(sizeof(double)*n*6);
    theModel->rootTree = allocOctTree();
    (theModel->rootTree)->nodes=NULL;
    (theModel->rootTree)->n_stars_contained=0;
    theModel->treeList = (OctTree **)malloc(sizeof(OctTree *)*n);
#ifdef USE_PPPM
    theModel->thePPPM = (PPPM *)malloc(sizeof(PPPM));
    allocPPPM(theModel->thePPPM,ngrid,ngrid,ngrid);
#endif
    theModel->tstep_last=-1.0;
    theModel->force_total = 0.0;
    theModel->pointsize = 0.02;

    return theModel;
}

void setPointsizeNbodyModel(NbodyModel *theModel,double pointsize) {
    theModel->pointsize = pointsize;
}
void setAnisotropyNbodyModel(NbodyModel *theModel,double anisotropy) {
    theModel->anisotropy = anisotropy;
}
void setDistributionNbodyModel(NbodyModel *theModel,int distribution) {
    theModel->distribution = distribution;
}
void setDistributionSigmaNbodyModel(NbodyModel *theModel,double distribution_sigma) {
    theModel->distribution_sigma = distribution_sigma;
}
void setDistributionXScaleNbodyModel(NbodyModel *theModel,double distribution_x_scale) {
    theModel->distribution_x_scale = distribution_x_scale;
}
void setDistributionYScaleNbodyModel(NbodyModel *theModel,double distribution_y_scale) {
    theModel->distribution_y_scale = distribution_y_scale;
}
void setDistributionZScaleNbodyModel(NbodyModel *theModel,double distribution_z_scale) {
    theModel->distribution_z_scale = distribution_z_scale;
}
void setPPPMCoeffsNbodyModel(NbodyModel *theModel,double ksigma,double knear) {
    theModel->pppm_ksigma = ksigma;
    theModel->pppm_knear = knear;
}
void setExpansionNbodyModel(NbodyModel *theModel,double expansion) {
    theModel->expansion = expansion;
}
void setDragNbodyModel(NbodyModel *theModel,double drag) {
    theModel->drag = drag;
}
void setPrefixNbodyModel(NbodyModel *theModel,const char * prefix) {
    strcpy(theModel->prefix,prefix);
}

void setDefaultsNbodyModel(NbodyModel *theModel) {
    theModel->default_mass=1.0; // solar masses
    theModel->default_scale=10.0; // parsecs
    theModel->default_G=0.0046254; // pc^2/Solar_mass/My^2
    theModel->srad_factor=5.0;
    theModel->softening_factor=0.0;
    theModel->rotation_factor=0.0;
    theModel->initial_v=0.0;
    theModel->tFinal=1000.0;
    theModel->tstep=1.0;
    theModel->tstep_ratio=-1.0;
    theModel->tstep_min=0.01;
    theModel->tstep_max=10.0;
    theModel->int_method=INT_METHOD_RK4;
    theModel->int_method=FORCE_METHOD_DIRECT;
    theModel->treeRangeCoefficient=1.2;
    theModel->drag=0.0;
    theModel->expansion=0.0;
    theModel->distribution=DISTRIBUTION_SPHERICAL_RANDOM_UNIFORM;
    theModel->distribution_sigma=0.2;
    theModel->distribution_x_scale=1.0;
    theModel->distribution_y_scale=1.0;
    theModel->distribution_z_scale=1.0;
    strcpy(theModel->prefix,"out");
}
void setForceMethod(NbodyModel *theModel,int force_method) {
    theModel->force_method=force_method;
}
void setTreeRangeCoefficient(NbodyModel *theModel,double coefficient) {
    theModel->treeRangeCoefficient=coefficient;
}
void setIntMethod(NbodyModel *theModel,int int_method) {
    theModel->int_method=int_method;
}
void setTStep(NbodyModel *theModel,double tstep) {
    theModel->tstep=tstep;
}
void setTStepRatio(NbodyModel *theModel,double tstep_ratio) {
    theModel->tstep_ratio=tstep_ratio;
}
void setTStepMin(NbodyModel *theModel,double tstep_min) {
    theModel->tstep_min=tstep_min;
}
void setTStepMax(NbodyModel *theModel,double tstep_max) {
    theModel->tstep_max=tstep_max;
}
void setRotationFactor(NbodyModel *theModel,double rotation_factor) {
    theModel->rotation_factor=rotation_factor;
}
void setTFinal(NbodyModel *theModel,double tFinal) {
    theModel->tFinal=tFinal;
}
void setInitialV(NbodyModel *theModel,double initial_v) {
    theModel->initial_v=initial_v;
}
void setSofteningNbodyModel(NbodyModel *theModel,double softening_factor) {
    theModel->softening_factor=softening_factor;
}
void setSradNbodyModel(NbodyModel *theModel,double srad_factor) {
    theModel->srad_factor=srad_factor;
}
void setColorNbodyModel(NbodyModel *theModel,int color) {
    theModel->default_color=color;
}
void setMassNbodyModel(NbodyModel *theModel,double mass) {
    theModel->default_mass=mass;
}
void setScaleNbodyModel(NbodyModel *theModel,double scale) {
    theModel->default_scale=scale;
}
void setGNbodyModel(NbodyModel *theModel,double G) {
    theModel->default_G=G;
}

int initializeNbodyModel(NbodyModel *theModel) {
    int i,j,k,l;
    int itest;
    double r2;
    double scale=theModel->default_scale;
    double r2Max=scale*scale;
    double tstep_squared,srad;
    double x,y,z,xmin,ymin,zmin,xstep,ystep,zstep;
    for(i=0;i<theModel->n;i++) {
        theModel->mass[i] = theModel->default_mass;
        if(theModel->distribution==DISTRIBUTION_SPHERICAL_RANDOM_UNIFORM) {
            do {
                theModel->x[i] = scale*(2.0*(double)rand()/
                    (double)RAND_MAX-1.0);
                theModel->y[i] = scale*(2.0*(double)rand()/
                    (double)RAND_MAX-1.0);
                theModel->z[i] = scale*(2.0*(double)rand()/
                    (double)RAND_MAX-1.0);
                r2 = theModel->x[i]*theModel->x[i]+
                    theModel->y[i]*theModel->y[i]+
                    theModel->z[i]*theModel->z[i];
            } while (r2>r2Max) ;
            theModel->x[i] *= theModel->distribution_x_scale;
            theModel->y[i] *= theModel->distribution_y_scale;
            theModel->z[i] *= theModel->distribution_z_scale;
        } else if(theModel->distribution==DISTRIBUTION_RECTANGULAR_RANDOM_UNIFORM) {
            //do {
            theModel->x[i] = theModel->distribution_x_scale*
                scale*(2.0*(double)rand()/
                (double)RAND_MAX-1.0);
            theModel->y[i] = theModel->distribution_y_scale*
                scale*(2.0*(double)rand()/
                (double)RAND_MAX-1.0);
            theModel->z[i] = theModel->distribution_z_scale*
                scale*(2.0*(double)rand()/
                (double)RAND_MAX-1.0);
                r2 = theModel->x[i]*theModel->x[i]+
                    theModel->y[i]*theModel->y[i]+
                    theModel->z[i]*theModel->z[i];
            //} while (r2<r2Max) ;
        } else if(theModel->distribution==DISTRIBUTION_RECTANGULAR_UNIFORM) {
            itest = (int)pow((double)theModel->n+1.0,(1.0/3.0));
            if(itest*itest*itest!=theModel->n) {
                printf("WARNING: RECTANGULAR UNIFORM DIST REQUIRES\n");
                printf("    N BE CUBIC NUMBER\n");
                printf("    %d    %d  \n",itest,theModel->n);
                exit(0);
            }
            xstep = (2.0*theModel->default_scale)/itest;
            ystep = (2.0*theModel->default_scale)/itest;
            zstep = (2.0*theModel->default_scale)/itest;
            xmin = -theModel->default_scale+xstep/2.0;
            ymin = -theModel->default_scale+ystep/2.0;
            zmin = -theModel->default_scale+zstep/2.0;
            for (l=0;l<itest;l++) {
                x = xmin+l*xstep;
                for (j=0;j<itest;j++) {
                    y = ymin+j*ystep;
                    for (k=0;k<itest;k++) {
                        z = zmin+k*zstep;
                        theModel->x[l*itest*itest+j*itest+k]=x+
                            drand(-1.0,1.0)*theModel->anisotropy;
                        theModel->y[l*itest*itest+j*itest+k]=y+
                            drand(-1.0,1.0)*theModel->anisotropy;
                        theModel->z[l*itest*itest+j*itest+k]=z+
                            drand(-1.0,1.0)*theModel->anisotropy;
                    }
                }
            }
            for (l=0;l<itest;l++) {
                for (j=0;j<itest;j++) {
                    for (k=0;k<itest;k++) {
                        theModel->x[l*itest*itest+j*itest+k] *=
                            theModel->distribution_x_scale;
                        theModel->y[l*itest*itest+j*itest+k] *=
                            theModel->distribution_y_scale;
                        theModel->z[l*itest*itest+j*itest+k] *=
                            theModel->distribution_z_scale;
                    }
                }
            }
        } else if(theModel->distribution==DISTRIBUTION_SPHERICAL_RANDOM_NORMAL) {
            theModel->x[i] = theModel->distribution_x_scale*
                drand_norm(0.0,theModel->distribution_sigma*scale,1.0e-3);
            theModel->y[i] = theModel->distribution_y_scale*
                drand_norm(0.0,theModel->distribution_sigma*scale,1.0e-3);
            theModel->z[i] = theModel->distribution_z_scale*
                drand_norm(0.0,theModel->distribution_sigma*scale,1.0e-3);
        } else {
            printf("WARNING: Distribution not understood\n");
            exit(0);
        }
        theModel->vx[i]=0.0;
        theModel->vy[i]=0.0;
        theModel->vz[i]=0.0;
        theModel->color[i]=theModel->default_color;
        for (j=0;j<6;j++) {
            theModel->X[i*6+j]=0.0;
            theModel->XPRIME[i*6+j]=0.0;
        }
        for (j=0;j<3;j++) {
            theModel->color3[i*3+j]=1.0;
        }
    }
    theModel->G=theModel->default_G;
    theModel->t=0.0;
    theModel->iteration=0;
    theModel->abmCounter=-3;
    theModel->bhCounter=0;
    tstep_squared = theModel->tstep*theModel->tstep;
    for (i=0;i<theModel->n;i++) {
        srad = computeSoftenedRadius(theModel->G*theModel->mass[i],
            tstep_squared,theModel->srad_factor);
        theModel->srad2[i]=srad*srad;
    }
/*
    calcStatistics(theModel);
    for(i=0;i<theModel->n;i++) {
        theModel->x[i] -= theModel->comx;
        theModel->y[i] -= theModel->comy;
        theModel->z[i] -= theModel->comz;
        theModel->vx[i] -= theModel->copx;
        theModel->vy[i] -= theModel->copy;
        theModel->vz[i] -= theModel->copz;
    }
*/
/*
    populateOctTree(theModel->rootTree, theModel->treeList, theModel,
        0,
        NULL,
        theModel->n,
        (-theModel->default_scale),
        theModel->default_scale,
        -(theModel->default_scale),
        theModel->default_scale,
        -(theModel->default_scale),
        theModel->default_scale);
*/

    TimerReset(&theModel->step_timer);

    return 1;
}

int freeNbodyModel(NbodyModel *theModel) {
    free(theModel->x);
    free(theModel->y);
    free(theModel->z);
    free(theModel->mass);
    free(theModel->vx);
    free(theModel->vy);
    free(theModel->vz);
    free(theModel->srad2);
    free(theModel->X);
    free(theModel->dbuffer);
    free(theModel->X2);
    free(theModel->draw_x);
    free(theModel->draw_y);
    free(theModel->draw_z);
    free(theModel->color);
    free(theModel->color3);
    free(theModel->XPRIME);
    free(theModel->XPRIME1);
    free(theModel->XPRIME2);
    free(theModel->XPRIME3);
    free(theModel->XPRIME4);
    freeOctTree(theModel->rootTree,0);
    free(theModel->treeList);
#ifdef USE_PPPM
    freePPPM(theModel->thePPPM);
    //free(theModel->thePPPM);
#endif
#ifdef HAS_FFTW3
    fftw_cleanup();
#endif
    free(theModel);
    return 1;
}

void calcDerivs(double * x, double * derivs, double t,double tStep,
        NbodyModel * theModel) {

    TimerReset(&theModel->force_timer);
    TimerStart(&theModel->force_timer);
    switch(theModel->force_method) {
#ifdef USE_PPPM
        case FORCE_METHOD_PPPM:
            calcDerivsPPPM(x,derivs,t,tStep,theModel);
         break;
#endif
        case FORCE_METHOD_TREE:
            calcDerivsBarnesHut(x,derivs,t,tStep,theModel);
         break;
        case FORCE_METHOD_DIRECT:
            calcDerivsDirect(x,derivs,t,tStep,theModel);
         break;
    }
    TimerStop(&theModel->force_timer);
    theModel->force_total += TimerRead(&theModel->force_timer);
}

#ifdef USE_PPPM
void calcDerivsPPPM(double * x, double * derivs, double t,double tStep,
        NbodyModel * theModel) {
    int i,j;
    point3d accel;

    double xmax = theModel->default_scale;
    double xmin = -xmax;
    double ymax = xmax;
    double ymin = xmin;
    double zmax = xmax;
    double zmin = xmin;
    double sigma = theModel->pppm_ksigma*(xmax-xmin)/theModel->thePPPM->nx;
    double near = theModel->pppm_knear*(xmax-xmin)/theModel->thePPPM->nx;

 #ifdef HAS_MPI
    extern int rank;
    extern int size;
 #endif

    //enforce periodic boundaries, ensure that any objects are within range
    for(i=0;i<theModel->n;i++) {
        for(j=0;j<3;j++) {
            while (x[6*i+j]<-theModel->default_scale)
               x[6*i+j] += theModel->default_scale*2.0;
            while (x[6*i+j]>theModel->default_scale)
               x[6*i+j] -= theModel->default_scale*2.0;
        }
    }

    //interpolate masses onto density distribution
    populateDensityPPPM(theModel->thePPPM,theModel,x,xmin,xmax,ymin,ymax,zmin,zmax,sigma);
    prepPotentialPPPM(theModel->thePPPM,theModel);
    //for each object, calculate forces
 #ifdef HAS_MPI
    for(i=0;i<theModel->n*6;i++) {
        derivs[i]=0.0;
    }
 #endif
 #ifdef HAS_MPI
    for(i=rank;i<theModel->n;i+=size) {
 #else
    for(i=0;i<theModel->n;i++) {
 #endif
        derivs[i*6+0] = x[i*6+3];
        derivs[i*6+1] = x[i*6+4];
        derivs[i*6+2] = x[i*6+5];
        near = 0.0;
        accel=calculateForcePPPM(i,x[i*6],x[i*6+1],x[i*6+2],
            theModel->thePPPM,theModel,x,near);
        derivs[i*6+3] = accel.x - theModel->drag*x[i*6+3]/theModel->mass[i];
        derivs[i*6+4] = accel.y - theModel->drag*x[i*6+4]/theModel->mass[i];
        derivs[i*6+5] = accel.z - theModel->drag*x[i*6+5]/theModel->mass[i];
    }
 #ifdef HAS_MPI
    MPI_Allreduce(derivs,theModel->dbuffer,theModel->n*6,
        MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    for(i=0;i<theModel->n*6;i++) derivs[i]=theModel->dbuffer[i];
 #endif
}
#endif

void calcDerivsBarnesHut(double * x, double * derivs, double t,double tStep,
        NbodyModel * theModel) {
    int i;
    point3d accel;
 #ifdef HAS_MPI
    extern int size;
    extern int rank;
 #endif
    for(i=0;i<theModel->n;i++) {
        if(theModel->treeList[i]!=NULL) {
            (theModel->treeList[i])->com_x=x[i*6];
            (theModel->treeList[i])->com_y=x[i*6+1];
            (theModel->treeList[i])->com_z=x[i*6+2];
        }
    }
    resetOctTree(theModel->rootTree,theModel);
 #ifdef HAS_MPI
    for(i=0;i<theModel->n*6;i++) derivs[i]=0.0;
    for(i=rank;i<theModel->n;i+=size) {
 #else
    for(i=0;i<theModel->n;i++) {
 #endif
        derivs[i*6+0] = x[i*6+3];
        derivs[i*6+1] = x[i*6+4];
        derivs[i*6+2] = x[i*6+5];
        accel=calculateForceOctTree(1,i,x[i*6],x[i*6+1],x[i*6+2],
            theModel->G,theModel->rootTree,theModel->softening_factor);
        derivs[i*6+3] = accel.x - theModel->drag*x[i*6+3]/theModel->mass[i];
        derivs[i*6+4] = accel.y - theModel->drag*x[i*6+4]/theModel->mass[i];
        derivs[i*6+5] = accel.z - theModel->drag*x[i*6+5]/theModel->mass[i];
    }
 #ifdef HAS_MPI
    MPI_Allreduce(derivs,theModel->dbuffer,theModel->n*6,MPI_DOUBLE,
        MPI_SUM,MPI_COMM_WORLD);
    for(i=0;i<theModel->n*6;i++) derivs[i]=theModel->dbuffer[i];
 #endif
}

void calcDerivsDirect(double * x, double * derivs, double t,double tStep,
        NbodyModel * theModel) {
    int i,j;
    double r3i;
    double deltaX,deltaY,deltaZ;
    double rad,r2;


#ifdef HAS_MPI
    //extern MPI_Status status;
    extern int size;
    extern int rank;
    double * buffer;
    buffer = (double *) malloc(sizeof(double)*theModel->n*6);
#endif
    // calculate a
    for (i=0;i<theModel->n;i++) {
        derivs[i*6+0] = x[i*6+3];
        derivs[i*6+1] = x[i*6+4];
        derivs[i*6+2] = x[i*6+5];
        derivs[i*6+3] = 0.0;
        derivs[i*6+4] = 0.0;
        derivs[i*6+5] = 0.0;
#ifdef HAS_MPI
       if (i%size==rank) {
#endif
        for (j=0;j<i;j++) {
            deltaX = x[j*6+0]-x[i*6+0];
            deltaY = x[j*6+1]-x[i*6+1];
            deltaZ = x[j*6+2]-x[i*6+2];
            r2 = deltaX*deltaX+deltaY*deltaY+deltaZ*deltaZ+
                theModel->softening_factor*theModel->softening_factor;
            rad=sqrt(r2);
            r3i=theModel->G/(rad*r2);
            if(r2>theModel->srad2[j]) {
                derivs[i*6+3] += theModel->mass[j]*deltaX*r3i;
                derivs[i*6+4] += theModel->mass[j]*deltaY*r3i;
                derivs[i*6+5] += theModel->mass[j]*deltaZ*r3i;
            }
            if(r2>theModel->srad2[i]) {
                derivs[j*6+3] -= theModel->mass[i]*deltaX*r3i;
                derivs[j*6+4] -= theModel->mass[i]*deltaY*r3i;
                derivs[j*6+5] -= theModel->mass[i]*deltaZ*r3i;
            }
        }
#ifdef HAS_MPI
       }
#endif
    }
#ifdef HAS_MPI
    MPI_Allreduce(derivs,buffer,theModel->n*6,
        MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    for(i=0;i<theModel->n;i++) {
        derivs[i*6+3]=buffer[i*6+3];
        derivs[i*6+4]=buffer[i*6+4];
        derivs[i*6+5]=buffer[i*6+5];
    }
    free(buffer);
#endif
    for(i=0;i<theModel->n;i++) {
        derivs[i*6+3]-=theModel->drag*x[i*6+3]/theModel->mass[i];
        derivs[i*6+4]-=theModel->drag*x[i*6+4]/theModel->mass[i];
        derivs[i*6+5]-=theModel->drag*x[i*6+5]/theModel->mass[i];
    }

    return;
}

void copy2X(NbodyModel * theModel) {
    int i;
    for (i=0;i<theModel->n;i++) {
        theModel->X[i*6+0]=theModel->x[i];
        theModel->X[i*6+1]=theModel->y[i];
        theModel->X[i*6+2]=theModel->z[i];
        theModel->X[i*6+3]=theModel->vx[i];
        theModel->X[i*6+4]=theModel->vy[i];
        theModel->X[i*6+5]=theModel->vz[i];
    }
}
void copy2xyz(NbodyModel * theModel) {
    int i;
    for (i=0;i<theModel->n;i++) {
        theModel->x[i]= theModel->X[i*6+0];
        theModel->y[i]= theModel->X[i*6+1];
        theModel->z[i]= theModel->X[i*6+2];
        theModel->vx[i]= theModel->X[i*6+3];
        theModel->vy[i]= theModel->X[i*6+4];
        theModel->vz[i]= theModel->X[i*6+5];
    }
}
int stepNbodyModel(NbodyModel * theModel) {
    int treeSkip=1;
    double treeScaleFactor=4.0;
    int i;
    double srad;
    int success;
    double tstep_squared;
    double rmsqdist,scale;
    double factor;
    TimerStart(&theModel->step_timer);
    if(theModel->tstep_ratio>0.0) {
        theModel->tstep = min(max(theModel->tstep_min,
            theModel->t*theModel->tstep_ratio),theModel->tstep_max);
    }
    if(theModel->tstep!=theModel->tstep_last) {
        tstep_squared = theModel->tstep*theModel->tstep;
        for (i=0;i<theModel->n;i++) {
            srad = computeSoftenedRadius(theModel->G*theModel->mass[i],
                tstep_squared,theModel->srad_factor);
            theModel->srad2[i]=srad*srad;
        }
        // flag to repopulate the tree, reset the repop counter
        theModel->tstep_last=theModel->tstep;
        theModel->abmCounter=-3;
        theModel->bhCounter=0;
    }
    if((theModel->bhCounter==0||theModel->rootTree==NULL)&&theModel->force_method==FORCE_METHOD_TREE) {
        freeOctTree(theModel->rootTree,1);
        scale = theModel->default_scale*treeScaleFactor;
        rmsqdist = 0.0;
        for(i=0;i<theModel->n;i++) {
            theModel->treeList[i]=NULL;
            rmsqdist += theModel->x[i]*theModel->x[i]+
                   theModel->y[i]*theModel->y[i]+
                   theModel->z[i]*theModel->z[i];
        }
        rmsqdist /= (double)theModel->n;
        rmsqdist = sqrt(rmsqdist);
        if(rmsqdist>(theModel->default_scale)) scale = rmsqdist*treeScaleFactor;
        populateOctTree(theModel->rootTree,theModel->treeList,theModel,
        0,NULL,theModel->n,-scale,scale,-scale,scale,-scale,scale);
    } else {
/*
// this is already done in calc derivs no matter what
        for(i=0;i<theModel->n;i++) {
            if(theModel->treeList[i]!=NULL) {
                theModel->treeList[i]->com_x=theModel->x[i];
                theModel->treeList[i]->com_y=theModel->y[i];
                theModel->treeList[i]->com_z=theModel->z[i];
            }
        }
        resetOctTree(theModel->rootTree,theModel);
*/
    }
    theModel->bhCounter=(theModel->bhCounter+1)%treeSkip;
    switch(theModel->int_method) {
        case INT_METHOD_ABM:
            success= stepNbodyModelABM(theModel,theModel->tstep);
         break;
        case INT_METHOD_RK4:
            success= stepNbodyModelRK4(theModel,theModel->tstep);
         break;
        case INT_METHOD_LEAPFROG:
            success= stepNbodyModelLeapfrog(theModel,theModel->tstep);
         break;
        case INT_METHOD_MPEULER:
            success= stepNbodyModelMPEuler(theModel,theModel->tstep);
         break;
        case INT_METHOD_IEULER:
            success= stepNbodyModelIEuler(theModel,theModel->tstep);
         break;
        case INT_METHOD_EULER:
        default:
            success= stepNbodyModelEuler(theModel,theModel->tstep);
         break;
    }
    //update model for periodic boundary conditions
    if(theModel->force_method==FORCE_METHOD_PPPM) {
        for(i=0;i<theModel->n;i++) {
            while (theModel->x[i]<-theModel->default_scale)
                theModel->x[i] += 2.0*theModel->default_scale;
            while (theModel->y[i]<-theModel->default_scale)
                theModel->y[i] += 2.0*theModel->default_scale;
            while (theModel->z[i]<-theModel->default_scale)
                theModel->z[i] += 2.0*theModel->default_scale;
            while (theModel->x[i]>theModel->default_scale)
                theModel->x[i] -= 2.0*theModel->default_scale;
            while (theModel->y[i]>theModel->default_scale)
                theModel->y[i] -= 2.0*theModel->default_scale;
            while (theModel->z[i]>theModel->default_scale)
                theModel->z[i] -= 2.0*theModel->default_scale;
        }
    }
    
    for(i=0;i<theModel->n;i++) {
        theModel->draw_x[i]=theModel->x[i];
        theModel->draw_y[i]=theModel->y[i];
        theModel->draw_z[i]=theModel->z[i];
    }

    // expand universe (not sure this works with ABM method)
    factor = 1.0+theModel->expansion*theModel->tstep/theModel->default_scale;
    theModel->default_scale += theModel->expansion*theModel->tstep;
    for(i=0;i<theModel->n;i++) {
        theModel->x[i]*=factor;
        theModel->y[i]*=factor;
        theModel->z[i]*=factor;
    }
    TimerStop(&theModel->step_timer);
    return success;
}
int stepNbodyModelIEuler(NbodyModel * theModel, double tStep) {
    int i;

    copy2X(theModel);

    calcDerivs(theModel->X,theModel->XPRIME,theModel->t,tStep,theModel);
    // update v,x
    for (i=0;i<theModel->n*6;i++) {
        theModel->X2[i]=theModel->X[i]+theModel->XPRIME[i]*tStep;
    }
    calcDerivs(theModel->X2,theModel->XPRIME2,theModel->t,tStep,theModel);
    for (i=0;i<theModel->n*6;i++) {
        theModel->XPRIME[i]=(theModel->XPRIME[i]+theModel->XPRIME2[i])/2;
    }
    for (i=0;i<theModel->n*6;i++) {
        theModel->X[i]=theModel->X[i]+theModel->XPRIME[i]*tStep;
    }

    copy2xyz(theModel);
    theModel->t += tStep;
    theModel->iteration++;

    return 1;
}
int stepNbodyModelMPEuler(NbodyModel * theModel, double tStep) {
    int i;

    copy2X(theModel);

    calcDerivs(theModel->X,theModel->XPRIME,theModel->t,tStep,theModel);

    // update v,x
    for (i=0;i<theModel->n*6;i++) {
        theModel->X2[i]=theModel->X[i]+theModel->XPRIME[i]*tStep/2;
    }

    calcDerivs(theModel->X2,theModel->XPRIME,theModel->t,tStep,theModel);

    for (i=0;i<theModel->n*6;i++) {
        theModel->X[i]=theModel->X[i]+theModel->XPRIME[i]*tStep;
    }

    copy2xyz(theModel);
    theModel->t += tStep;
    theModel->iteration++;

    return 1;
}
int stepNbodyModelRK4(NbodyModel * theModel, double tStep) {
    int i;

    copy2X(theModel);

    calcDerivs(theModel->X,theModel->XPRIME,theModel->t,tStep,theModel);

    // update v,x
    for (i=0;i<theModel->n*6;i++) {
        theModel->X2[i] = theModel->X[i] + theModel->XPRIME[i]*tStep/2;
    }
    calcDerivs(theModel->X2,theModel->XPRIME2,theModel->t+tStep/2,tStep,theModel);
    for (i=0;i<theModel->n*6;i++) {
        theModel->X2[i] = theModel->X[i] + theModel->XPRIME2[i]*tStep/2;
    }
    calcDerivs(theModel->X2,theModel->XPRIME3,theModel->t+tStep/2,tStep,theModel);
    for (i=0;i<theModel->n*6;i++) {
        theModel->X2[i] = theModel->X[i] + theModel->XPRIME3[i]*tStep;
    }
    calcDerivs(theModel->X2,theModel->XPRIME4,theModel->t+tStep,tStep,theModel);
    for (i=0;i<theModel->n*6;i++) {
        theModel->XPRIME[i] += 2.0*theModel->XPRIME2[i] +
                               2.0*theModel->XPRIME3[i] +
                                   theModel->XPRIME4[i] ;
        theModel->XPRIME[i] /= 6.0;
    }
    for (i=0;i<theModel->n*6;i++) {
        theModel->X[i] += theModel->XPRIME[i]*tStep;
    }

    copy2xyz(theModel);
    theModel->t += tStep;
    theModel->iteration++;

    return 1;
}
int stepNbodyModelLeapfrog(NbodyModel * theModel, double tStep) {
    int i;

    copy2X(theModel);

    calcDerivs(theModel->X,theModel->XPRIME,theModel->t,tStep,theModel);

    if(theModel->t<0.5*tStep) {
        // setup leapfrog on first step, change velocities by a half step
        for (i=0;i<theModel->n;i++) {
            theModel->X[i*6+3] += 0.5*theModel->XPRIME[i*6+3]*tStep;
            theModel->X[i*6+4] += 0.5*theModel->XPRIME[i*6+4]*tStep;
            theModel->X[i*6+5] += 0.5*theModel->XPRIME[i*6+5]*tStep;
        }
    } else {
        // update v,x
        for (i=0;i<theModel->n*6;i++) {
            theModel->X[i] += theModel->XPRIME[i]*tStep;
        }
    }


    copy2xyz(theModel);
    theModel->t += tStep;
    theModel->iteration++;

    return 1;
}
int stepNbodyModelEuler(NbodyModel * theModel, double tStep) {
    int i;

    copy2X(theModel);

    calcDerivs(theModel->X,theModel->XPRIME,theModel->t,tStep,theModel);

    // update v,x
    for (i=0;i<theModel->n*6;i++) {
        theModel->X[i] += theModel->XPRIME[i]*tStep;
    }

    copy2xyz(theModel);
    theModel->t += tStep;
    theModel->iteration++;

    return 1;
}
int stepNbodyModelABM(NbodyModel * theModel, double tStep) {
    //Adams-Bashforth-Moulton Predictor Corrector
    int i;
    double * fk3=NULL;
    double * fk2=NULL;
    double * fk1=NULL;
    double * fk0=NULL;
    double * fkp=NULL;

    // determine if previous steps exist, if not, populate w/ RK4
    if(theModel->abmCounter<0) {
        stepNbodyModelRK4(theModel,tStep);
        if(theModel->abmCounter==-3) {
            for (i=0;i<theModel->n*6;i++)
                theModel->XPRIME4[i]=theModel->XPRIME[i];
        } else if (theModel->abmCounter==-2) {
            for (i=0;i<theModel->n*6;i++)
                theModel->XPRIME3[i]=theModel->XPRIME[i];
        } else {
            for (i=0;i<theModel->n*6;i++)
                theModel->XPRIME2[i]=theModel->XPRIME[i];
        }
    } else {
        copy2X(theModel);
        if(theModel->abmCounter%5==0) {
            fk3=theModel->XPRIME4;
            fk2=theModel->XPRIME3;
            fk1=theModel->XPRIME2;
            fk0=theModel->XPRIME1;
            fkp=theModel->XPRIME;
        } else if (theModel->abmCounter%5==1) {
            fk3=theModel->XPRIME3;
            fk2=theModel->XPRIME2;
            fk1=theModel->XPRIME1;
            fk0=theModel->XPRIME;
            fkp=theModel->XPRIME4;
        } else if (theModel->abmCounter%5==2) {
            fk3=theModel->XPRIME2;
            fk2=theModel->XPRIME1;
            fk1=theModel->XPRIME;
            fk0=theModel->XPRIME4;
            fkp=theModel->XPRIME3;
        } else if (theModel->abmCounter%5==3) {
            fk3=theModel->XPRIME1;
            fk2=theModel->XPRIME;
            fk1=theModel->XPRIME4;
            fk0=theModel->XPRIME3;
            fkp=theModel->XPRIME2;
        } else if (theModel->abmCounter%5==4) {
            fk3=theModel->XPRIME;
            fk2=theModel->XPRIME4;
            fk1=theModel->XPRIME3;
            fk0=theModel->XPRIME2;
            fkp=theModel->XPRIME1;
        }
        calcDerivs(theModel->X,fk0,theModel->t,tStep,theModel);
        for (i=0;i<theModel->n*6;i++) {
            theModel->X2[i] = theModel->X[i] +
                               (tStep/24.0)*(-9.0*fk3[i]+37.0*fk2[i]
                               -59.0*fk1[i]+55.0*fk0[i]);
        }
        calcDerivs(theModel->X2,fkp,theModel->t+tStep,tStep,theModel);
        for (i=0;i<theModel->n*6;i++) {
            theModel->X[i] = theModel->X[i] +
                               (tStep/24.0)*(fk2[i]-5.0*fk1[i]+
                                19.0*fk0[i]+9.0*fkp[i]);
        }
        copy2xyz(theModel);
        theModel->t += tStep;
        theModel->iteration++;
    }

    theModel->abmCounter++;
    return 1;
}

#ifdef HAS_SDL
sdlwindow theWindow;
int window_created_sdl=0;
NbodyModel *sdlModel;

void drawCube(double x, double y, double z, double size) {
    double so2 = size/2;
    x /= sdlModel->default_scale;
    y /= sdlModel->default_scale;
    z /= sdlModel->default_scale;
        glBegin(GL_QUADS);
        glNormal3f(0.0,-1.0,0.0);
        glVertex3f(x-so2,y-so2,z-so2);
        glVertex3f(x-so2,y-so2,z+so2);
        glVertex3f(x+so2,y-so2,z+so2);
        glVertex3f(x+so2,y-so2,z-so2);
        glNormal3f(0.0,+1.0,0.0);
        glVertex3f(x-so2,y+so2,z-so2);
        glVertex3f(x-so2,y+so2,z+so2);
        glVertex3f(x+so2,y+so2,z+so2);
        glVertex3f(x+so2,y+so2,z-so2);
        glNormal3f(-1.0,0.0,0.0);
        glVertex3f(x-so2,y+so2,z-so2);
        glVertex3f(x-so2,y+so2,z+so2);
        glVertex3f(x-so2,y-so2,z+so2);
        glVertex3f(x-so2,y-so2,z-so2);
        glNormal3f(+1.0,0.0,0.0);
        glVertex3f(x+so2,y+so2,z-so2);
        glVertex3f(x+so2,y+so2,z+so2);
        glVertex3f(x+so2,y-so2,z+so2);
        glVertex3f(x+so2,y-so2,z-so2);
        glNormal3f(0.0,0.0,-1.0);
        glVertex3f(x-so2,y+so2,z-so2);
        glVertex3f(x+so2,y+so2,z-so2);
        glVertex3f(x+so2,y-so2,z-so2);
        glVertex3f(x-so2,y-so2,z-so2);
        glNormal3f(0.0,0.0,1.0);
        glVertex3f(x-so2,y+so2,z+so2);
        glVertex3f(x+so2,y+so2,z+so2);
        glVertex3f(x+so2,y-so2,z+so2);
        glVertex3f(x-so2,y-so2,z+so2);
        glEnd();

}

void drawSphere(double r,double xs, double ys, double zs, int lats, int longs) {
    int i, j;
    double lat0 = M_PI * (-0.5 + (double) (0-1) / lats);
    double z0  = r*sin(lat0);
    double zr0 =  r*cos(lat0);
    xs /= sdlModel->default_scale;
    ys /= sdlModel->default_scale;
    zs /= sdlModel->default_scale;
    for(i = 0; i <= lats; i++) {
        double lat1 = M_PI * (-0.5 + (double) i / lats);
        double z1 = r*sin(lat1);
        double zr1 = r*cos(lat1);

        glBegin(GL_QUAD_STRIP);
        for(j = 0; j <= longs; j++) {
            double x0,x1,y0,y1;
            double lng = 2 * M_PI * (double) (j - 1) / longs;
            double x = cos(lng);
            double y = sin(lng);
                
            x0 = x*zr0;
            y0 = y*zr0;
            x1 = x*zr1;
            y1 = y*zr1;
            glNormal3f(x0, y0, z0);
            glVertex3f((x0+xs),
                       (y0+ys),
                       (z0+zs));
            glNormal3f(x1, y1, z1);
            glVertex3f((x1+xs),
                       (y1+ys),
                       (z1+zs));
        }
        glEnd();

        lat0=lat1;
        z0=z1;
        zr0=zr1;
    }
}

void userDrawWorld() {
    int i;
    double scale = 1.0;
    glColor3d(1.0,1.0,1.0);
    glBegin(GL_LINE_STRIP);
    glVertex3d(-scale,-scale,-scale);
    glVertex3d(-scale,-scale,scale);
    glVertex3d(-scale,scale,scale);
    glVertex3d(-scale,scale,-scale);
    glVertex3d(-scale,-scale,-scale);
    glVertex3d(scale,-scale,-scale);
    glVertex3d(scale,-scale,scale);
    glVertex3d(scale,scale,scale);
    glVertex3d(scale,scale,-scale);
    glVertex3d(scale,-scale,-scale);
    glEnd();
    glBegin(GL_LINES);
    glVertex3d(-scale,-scale,scale);
    glVertex3d(scale,-scale,scale);
    glVertex3d(-scale,scale,scale);
    glVertex3d(scale,scale,scale);
    glVertex3d(-scale,scale,-scale);
    glVertex3d(scale,scale,-scale);
    glEnd();
    for(i=0;i<sdlModel->n;i++) {
        glColor3d(sdlModel->color3[i*3],
            sdlModel->color3[i*3+1],sdlModel->color3[i*3+2]);
        //drawCube(sdlModel->x[i],sdlModel->y[i],sdlModel->z[i],0.05);
        drawSphere(sdlModel->pointsize,
            sdlModel->draw_x[i],
            sdlModel->draw_y[i],sdlModel->draw_z[i],7,7);
    }
}
#endif

#ifdef HAS_X11
Display *dpy_x11;
int black_x11;
int white_x11;
Window w_x11;
GC gc_x11;
Pixmap buffer_x11;
Colormap theColormap_x11;
int numXGrayscale_x11=10;
XColor XGrayscale_x11[10];
int numXBluescale_x11=10;
XColor XBluescale_x11[10];
int width_x11=800;
int height_x11=400;
int window_created_x11=0;

void setupWindow(int IMAGE_WIDTH, int IMAGE_HEIGHT) {
    int i,color;

    dpy_x11 = XOpenDisplay(NIL);
    assert(dpy_x11);

    black_x11 = BlackPixel(dpy_x11, DefaultScreen(dpy_x11));
    white_x11 = WhitePixel(dpy_x11, DefaultScreen(dpy_x11));

    w_x11 = XCreateSimpleWindow(dpy_x11, DefaultRootWindow(dpy_x11), 0, 0,
        IMAGE_WIDTH, IMAGE_HEIGHT, 0, black_x11,
        black_x11);
    buffer_x11 = XCreatePixmap(dpy_x11,DefaultRootWindow(dpy_x11),
        IMAGE_WIDTH,IMAGE_HEIGHT,DefaultDepth(dpy_x11,
        DefaultScreen(dpy_x11)));
    theColormap_x11 = XCreateColormap(dpy_x11, DefaultRootWindow(dpy_x11),
        DefaultVisual(dpy_x11,DefaultScreen(dpy_x11)), AllocNone);

    for (i=0;i<numXGrayscale_x11;i++) {
        color = (int)((double)i*35535.0/(double)numXGrayscale_x11)+30000;
        XGrayscale_x11[i].red=0.8*color;
        XGrayscale_x11[i].green=color;
        XGrayscale_x11[i].blue=0.8*color;
        XAllocColor(dpy_x11,theColormap_x11,&(XGrayscale_x11[i]));
        XBluescale_x11[i].red=0.85*color;
        XBluescale_x11[i].green=0.85*color;
        XBluescale_x11[i].blue=color;
        XAllocColor(dpy_x11,theColormap_x11,&(XBluescale_x11[i]));
    }

    XSelectInput(dpy_x11, w_x11, StructureNotifyMask);
    XMapWindow(dpy_x11, w_x11);
    gc_x11 = XCreateGC(dpy_x11, w_x11, 0, NIL);
    XSetForeground(dpy_x11, gc_x11, white_x11);

    for(;;) {
        XEvent e;
        XNextEvent(dpy_x11, &e);
        if (e.type == MapNotify)
        break;
    } 
} 

void make_image(NbodyModel * theModel) {
     
    int i;
    double scale,shift,size;
    int isize;
    int dispX, dispY, dispZ, depthY, depthZ;
     
    XSetForeground(dpy_x11, gc_x11, black_x11);
    XFillRectangle(dpy_x11,buffer_x11,gc_x11,
        0,0,width_x11,height_x11);

    scale=(double)height_x11/(2.0*theModel->default_scale);
    size = scale*theModel->pointsize/2.0;
    isize = (int)size;
    if(isize<2) isize=2;
    shift=(double)width_x11/4.0;
    for (i=0 ; i<theModel->n; i++) {
        dispX = (int)(scale*theModel->x[i]+shift);
        dispY = (int)(scale*theModel->y[i]+shift);
        dispZ = (int)(scale*theModel->z[i]+shift);
        
        depthY = (int)((double)dispY/(double)height_x11*numXBluescale_x11);
        depthZ = (int)((double)dispZ/(double)height_x11*numXGrayscale_x11);
        if (depthY>numXBluescale_x11-1) depthY=numXBluescale_x11-1;
        if (depthZ>numXGrayscale_x11-1) depthZ=numXGrayscale_x11-1;

        if (dispX < width_x11/2) {
            XSetForeground(dpy_x11,gc_x11,XGrayscale_x11[depthZ].pixel);
            XFillRectangle(dpy_x11,buffer_x11,gc_x11,
                dispX-isize/2,dispY+isize/2,isize,isize);
        }
        if (dispX > 0) {
            XSetForeground(dpy_x11,gc_x11,XBluescale_x11[depthY].pixel);
            XFillRectangle(dpy_x11,buffer_x11,gc_x11,
                dispX+width_x11/2-isize/2,dispZ+isize/2,isize,isize);
        }
   }
   XSetForeground(dpy_x11, gc_x11, white_x11);
   XFillRectangle(dpy_x11, buffer_x11, gc_x11,
       width_x11/2-isize/2,0,isize,height_x11);
   XCopyArea(dpy_x11, buffer_x11, w_x11, gc_x11, 0, 0,
       width_x11, height_x11,  0, 0);
   XFlush(dpy_x11);
}
#endif

void speedNbodyModel(NbodyModel *theModel) {
    double pi=atan(1.0)*4.0;
    double theta,phi;
    double v = theModel->initial_v;
    int i;
    for (i=0;i<theModel->n;i++) {
        theta = 2*pi*(double)rand()/(double)RAND_MAX;
        phi = pi*(double)rand()/(double)RAND_MAX;
        theModel->vx[i]+=v*sin(phi)*cos(theta);
        theModel->vy[i]+=v*sin(phi)*sin(theta);
        theModel->vz[i]+=v*cos(phi);
    }
}

double scale_factor(double x, double c,double sigma) {
    return 0.5*(erf((x-c)/sigma)+1.0);
}


void spinNbodyModel(NbodyModel *theModel) {
    int i,j;
    double rij,rij3;
    double ax,ay,az;
    double dx,dy,dz;
    double rplane, rsphere;
    double vr;
    double offset, sigma;

    for(i=0;i<theModel->n;i++) {
        rplane = sqrt(theModel->x[i]*theModel->x[i]+theModel->y[i]*theModel->y[i]);
        rsphere = sqrt(theModel->x[i]*theModel->x[i]+theModel->y[i]*theModel->y[i]+theModel->z[i]*theModel->z[i]);
        ax=0;ay=0;az=0;
        for(j=0;j<theModel->n;j++) {
            if(i!=j) {
                dx = theModel->x[j]-theModel->x[i];
                dy = theModel->y[j]-theModel->y[i];
                dz = theModel->z[j]-theModel->z[i];
                rij = sqrt(dx*dx+dy*dy+dz*dz);
                rij3 = theModel->G*theModel->mass[j]/(rij*rij*rij);
                ax += rij3*dx;
                ay += rij3*dy;
                az += rij3*dz;
            }
        }
        if(rplane>0.0) {
            vr=theModel->rotation_factor*sqrt(sqrt(ax*ax+ay*ay+az*az)*
                rplane*rplane/rsphere)/
                rplane*scale_factor(rplane,
                    0.4*theModel->distribution_sigma*theModel->default_scale,
                    0.4*theModel->distribution_sigma*theModel->default_scale);
        } else {
            vr = 0.0;
        }
        theModel->vx[i] = vr*theModel->y[i];
        theModel->vy[i] = -vr*theModel->x[i];
        theModel->vz[i]=0.0;
    }


    return;
}

int updateNbodyModel(NbodyModel *theModel,int updateMethod) {
    Text11 * t11p;
    int i;
    FILE *dumpout;
    char fname[80];

#ifdef HAS_LIBGD
    gdImagePtr im;
    FILE *pngout;
    int black;
    int white;
    int dx,dy;
    float rx,ry,rz;
    float rmin=-theModel->default_scale;
    float rmax=theModel->default_scale;
    int width=400;
    int height=200;
    int dxmax=width/2;
    int dxmin=0;
    int dymax=0;
    int dymin=height;
    int gd_point_size=2;    // KLUDGE replace with reference to pointsize
#endif

    // output, send to client, send to screen, etc.
    if((updateMethod&UPDATEMETHOD_HASH_TEXT)==UPDATEMETHOD_HASH_TEXT) {
        printf("#");
        fflush(stdout);
    } 

    if((updateMethod&UPDATEMETHOD_BRIEF_TEXT)==UPDATEMETHOD_BRIEF_TEXT) {
        printf("UM%d  Calculated t = %10.3e\n",UPDATEMETHOD_BRIEF_TEXT,theModel->t);
    }

    if((updateMethod&UPDATEMETHOD_VERBOSE_POSITIONS)==UPDATEMETHOD_VERBOSE_POSITIONS) {
        for(i=0;i<theModel->n;i++) {
            printf(
                "UM%d  POS %10.3e\t%d\t%10.3e\t%10.3e\t%10.3e\t%10.3e\t%10.3e\t%10.3e %d\n",UPDATEMETHOD_VERBOSE_POSITIONS,
                theModel->t,i,theModel->x[i],theModel->y[i],
                theModel->z[i],theModel->vx[i],theModel->vy[i],
                theModel->vz[i],theModel->color[i]);
        }
    }
    if((updateMethod&UPDATEMETHOD_VERBOSE_STATISTICS)==UPDATEMETHOD_VERBOSE_STATISTICS) {
        calcStatistics(theModel);
        printf(
            "UM%d  STAT %10.3e\t%10.3e\t%10.3e\t%10.3e\n",UPDATEMETHOD_VERBOSE_STATISTICS,
            theModel->t,theModel->KE,theModel->PE,(theModel->KE+theModel->PE));
    }

    if((updateMethod&UPDATEMETHOD_TEXT11)==UPDATEMETHOD_TEXT11) {
            printf("------------------------------ t=%lf\n",theModel->t);
            text11_initialize(&t11p);
            text11_set_range(t11p,-theModel->default_scale,
                -theModel->default_scale,
                theModel->default_scale,
                theModel->default_scale);
            text11_set_boundary(t11p,theModel->n/2);
            for(i=0;i<theModel->n;i++) {
                text11_add(t11p,theModel->x[i],theModel->z[i]);
            }
            text11_print(t11p);
            printf("==============================\n");
            text11_free(&t11p);
    }

    if((updateMethod&UPDATEMETHOD_DUMP)==UPDATEMETHOD_DUMP) {
        sprintf(fname,"%s%06d.dump",theModel->prefix,theModel->iteration);
        dumpout = fopen(fname, "wb");
        fprintf(dumpout,"%d\t%lf\t%lf\n",
            theModel->n,theModel->t,theModel->default_scale);
        for(i=0;i<theModel->n;i++) {
            double rx,ry,rz,vx,vy,vz,mass;
            rx = theModel->x[i];
            ry = theModel->y[i];
            rz = theModel->z[i];
            vx = theModel->vx[i];
            vy = theModel->vy[i];
            vz = theModel->vz[i];
            mass = theModel->mass[i];
            fprintf(dumpout,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",
                i,rx,ry,rz,vx,vy,vz,mass);
        }
        fclose(dumpout);
    }

#ifdef HAS_NETCDF
    if((updateMethod&UPDATEMETHOD_NETCDF)==UPDATEMETHOD_NETCDF) {
        int ncid,uid,xid,yid,zid,vxid,vyid,vzid,mid,cid;
        int start,count,i,j,k,xgid,ygid,zgid;
        int dimids[3];
        int did,pid,fxid,fyid,fzid;
        double * buffer;


        //if PPPM, enforce periodic boundaries, ensure that any objects are within range
        if(theModel->force_method==FORCE_METHOD_PPPM) {
            for(i=0;i<theModel->n;i++) {
                 while (theModel->x[i]<theModel->default_scale)
                    theModel->x[i] += theModel->default_scale*2.0;
                 while (theModel->x[i]>theModel->default_scale)
                    theModel->x[i] -= theModel->default_scale*2.0;
                 while (theModel->y[i]<theModel->default_scale)
                    theModel->y[i] += theModel->default_scale*2.0;
                 while (theModel->y[i]>theModel->default_scale)
                    theModel->y[i] -= theModel->default_scale*2.0;
                 while (theModel->z[i]<theModel->default_scale)
                    theModel->z[i] += theModel->default_scale*2.0;
                 while (theModel->z[i]>theModel->default_scale)
                    theModel->z[i] -= theModel->default_scale*2.0;
            }
        }




        start=0;
        count=theModel->n;
        sprintf(fname,"%s_part_%06d.nc",theModel->prefix,theModel->iteration);
        nc_create(fname,NC_CLOBBER, &ncid);
        nc_def_dim(ncid,"index",theModel->n, &uid);
        nc_def_var(ncid,"x",NC_DOUBLE,1, &uid, &xid);
        nc_def_var(ncid,"y",NC_DOUBLE,1, &uid, &yid);
        nc_def_var(ncid,"z",NC_DOUBLE,1, &uid, &zid);
        nc_def_var(ncid,"vx",NC_DOUBLE,1, &uid, &vxid);
        nc_def_var(ncid,"vy",NC_DOUBLE,1, &uid, &vyid);
        nc_def_var(ncid,"vz",NC_DOUBLE,1, &uid, &vzid);
        nc_def_var(ncid,"mass",NC_DOUBLE,1, &uid, &mid);
        nc_def_var(ncid,"color",NC_INT,1, &uid, &cid);
        nc_enddef(ncid);
        nc_put_var_double(ncid,xid,&(theModel->x[0]));
        nc_put_var_double(ncid,yid,&(theModel->y[0]));
        nc_put_var_double(ncid,zid,&(theModel->z[0]));
        nc_put_var_double(ncid,vxid,&(theModel->vx[0]));
        nc_put_var_double(ncid,vyid,&(theModel->vy[0]));
        nc_put_var_double(ncid,vzid,&(theModel->vz[0]));
        nc_put_var_double(ncid,mid,&(theModel->mass[0]));
        nc_put_var_int(ncid,cid,&(theModel->color[0]));
        nc_close(ncid);
#ifdef USE_PPPM
        if(theModel->force_method==FORCE_METHOD_PPPM) {
            buffer = (double *)malloc(sizeof(double)*theModel->thePPPM->nxyz);
            sprintf(fname,"%s_mesh_%06d.nc",
                theModel->prefix,theModel->iteration);
            nc_create(fname,NC_CLOBBER, &ncid);
            nc_def_dim(ncid,"x",theModel->thePPPM->nx, &xid);
            nc_def_dim(ncid,"y",theModel->thePPPM->ny, &yid);
            nc_def_dim(ncid,"z",theModel->thePPPM->nz, &zid);
            dimids[0]=zid;
            dimids[1]=yid;
            dimids[2]=xid;
            nc_def_var(ncid,"x",NC_DOUBLE,1, &xid, &xgid);
            nc_def_var(ncid,"y",NC_DOUBLE,1, &yid, &ygid);
            nc_def_var(ncid,"z",NC_DOUBLE,1, &zid, &zgid);
            nc_def_var(ncid,"density",NC_DOUBLE,3, dimids, &did);
            nc_def_var(ncid,"potential",NC_DOUBLE,3, dimids, &pid);
            nc_def_var(ncid,"force_x",NC_DOUBLE,3, dimids, &fxid);
            nc_def_var(ncid,"force_y",NC_DOUBLE,3, dimids, &fyid);
            nc_def_var(ncid,"force_z",NC_DOUBLE,3, dimids, &fzid);
            nc_enddef(ncid);
            nc_put_var_double(ncid,xgid,&(theModel->thePPPM->x[0]));
            nc_put_var_double(ncid,ygid,&(theModel->thePPPM->y[0]));
            nc_put_var_double(ncid,zgid,&(theModel->thePPPM->z[0]));
            for(i=0;i<theModel->thePPPM->nx;i++) 
                for(j=0;j<theModel->thePPPM->ny;j++) 
                    for(k=0;k<theModel->thePPPM->nz;k++) 
                       buffer[k*theModel->thePPPM->nx*theModel->thePPPM->ny+
                              j*theModel->thePPPM->nx+i] =
                       theModel->thePPPM->density[
                              i*theModel->thePPPM->nz*theModel->thePPPM->ny+
                              j*theModel->thePPPM->nz+k];
            nc_put_var_double(ncid,did,buffer);
            for(i=0;i<theModel->thePPPM->nx;i++)
                for(j=0;j<theModel->thePPPM->ny;j++)
                    for(k=0;k<theModel->thePPPM->nz;k++)
                       buffer[k*theModel->thePPPM->nx*theModel->thePPPM->ny+
                              j*theModel->thePPPM->nx+i] =
                       theModel->thePPPM->potential[
                              i*theModel->thePPPM->nz*theModel->thePPPM->ny+
                              j*theModel->thePPPM->nz+k];
            nc_put_var_double(ncid,pid,buffer);
            for(i=0;i<theModel->thePPPM->nx;i++)
                for(j=0;j<theModel->thePPPM->ny;j++)
                    for(k=0;k<theModel->thePPPM->nz;k++)
                       buffer[k*theModel->thePPPM->nx*theModel->thePPPM->ny+
                              j*theModel->thePPPM->nx+i] =
                       theModel->thePPPM->force[
                              i*theModel->thePPPM->nz*theModel->thePPPM->ny+
                              j*theModel->thePPPM->nz+k].x;
            nc_put_var_double(ncid,fxid,buffer);
            for(i=0;i<theModel->thePPPM->nx;i++)
                for(j=0;j<theModel->thePPPM->ny;j++)
                    for(k=0;k<theModel->thePPPM->nz;k++)
                       buffer[k*theModel->thePPPM->nx*theModel->thePPPM->ny+
                              j*theModel->thePPPM->nx+i] =
                       theModel->thePPPM->force[
                              i*theModel->thePPPM->nz*theModel->thePPPM->ny+
                              j*theModel->thePPPM->nz+k].y;
            nc_put_var_double(ncid,fyid,buffer);
            for(i=0;i<theModel->thePPPM->nx;i++)
                for(j=0;j<theModel->thePPPM->ny;j++)
                    for(k=0;k<theModel->thePPPM->nz;k++)
                       buffer[k*theModel->thePPPM->nx*theModel->thePPPM->ny+
                              j*theModel->thePPPM->nx+i] =
                       theModel->thePPPM->force[
                              i*theModel->thePPPM->nz*theModel->thePPPM->ny+
                              j*theModel->thePPPM->nz+k].z;
            nc_put_var_double(ncid,fzid,buffer);
            nc_close(ncid);
            free(buffer);
        }
#endif
    }
#endif

#ifdef HAS_LIBGD
    if((updateMethod&UPDATEMETHOD_GD_IMAGE)==UPDATEMETHOD_GD_IMAGE) {
        // set up GD image
        im = gdImageCreate(width,height);
        black = gdImageColorAllocate(im, 0, 0, 0);
        white = gdImageColorAllocate(im, 255, 255, 255);

        gd_point_size = (int)((double)height * theModel->pointsize);
        if(gd_point_size<1) gd_point_size=1;
    
        gdImageFilledRectangle(im,0,0,width,height,white);
        gdImageFilledRectangle(im,width/2-gd_point_size/2,0,
            width/2+gd_point_size/2+1,height,black);
        for(i=0;i<theModel->n;i++) {
            rx = theModel->x[i];
            ry = theModel->y[i];
            rz = theModel->z[i];
            // plot top view
            dx = dxmin + (int)((rx-rmin)/(rmax-rmin)*(float)(dxmax-dxmin));
            dy = dymin + (int)((ry-rmin)/(rmax-rmin)*(float)(dymax-dymin));
            if(gd_point_size==1)
                gdImageSetPixel(im,dx,dy,black);
            else
                gdImageFilledRectangle(im,dx-gd_point_size/2,
                    dy+gd_point_size/2,
                    dx-gd_point_size/2+gd_point_size,
                    dy+gd_point_size/2-gd_point_size,
                    black);
            // plot side view view
            dx += dxmax;
            dy = dymin + (int)((rz-rmin)/(rmax-rmin)*(float)(dymax-dymin));
            if(gd_point_size==1)
                gdImageSetPixel(im,dx,dy,black);
            else
                gdImageFilledRectangle(im,dx-gd_point_size/2,
                    dy+gd_point_size/2,
                    dx-gd_point_size/2+gd_point_size,
                    dy+gd_point_size/2-gd_point_size,
                    black);
        }
        sprintf(fname,"%s%06d.png",theModel->prefix,theModel->iteration);
        pngout = fopen(fname, "wb");
        gdImagePng(im, pngout);
        fclose(pngout);
    
        //cleanup
        gdImageDestroy(im);
    }
#endif

#ifdef HAS_X11
    if((updateMethod&UPDATEMETHOD_X11)==UPDATEMETHOD_X11) {
        if(!window_created_x11){
            setupWindow(width_x11,height_x11);
            window_created_x11=1;
        }
        make_image(theModel);
    }
#endif

#ifdef HAS_SDL
    if((updateMethod&UPDATEMETHOD_SDL)==UPDATEMETHOD_SDL) {
        if(!window_created_sdl){
            sdlModel = theModel;
            sdlwindow_begin(&theWindow,800,600);
            theWindow.userDrawWorld = userDrawWorld;
            window_created_sdl=1;
        }
        sdlwindow_render(&theWindow);
    }
#endif

    return 1;
}

void nbodyEvents(NbodyModel * theModel, int updateMethod) {
#ifdef HAS_SDL
    if((updateMethod&UPDATEMETHOD_SDL)==UPDATEMETHOD_SDL) {
        if(window_created_sdl) {
            if(sdlwindow_pollAndHandle(&theWindow)) {
                sdlwindow_render(&theWindow);
            }
        }
    }
#endif
}

void printStatistics(NbodyModel * theModel) {
    double total_energy;
    calcStatistics(theModel);

    total_energy = theModel->KE + theModel->PE;
    printf("KE = %le \t",theModel->KE);
    printf("PE = %le \t",theModel->PE);
    printf("TE = %le \n",total_energy);
    printf("COM x = %le (%10.3le) \t COP x = %le (%10.3le) \n",theModel->comx,theModel->rmsd,theModel->copx,theModel->rmsp);
    printf("COM y = %le \t COP y = %le \n",theModel->comy,theModel->copy);
    printf("COM z = %le \t COP z = %le \n",theModel->comz,theModel->copz);
}
void calcStatistics(NbodyModel * theModel) {
    int i,j;
    double v2,r2,r,dx,dy,dz;
    double total_mass;

    theModel->KE=0.0;
    theModel->PE=0.0;
    theModel->comx=0.0;
    theModel->comy=0.0;
    theModel->comz=0.0;
    theModel->copx=0.0;
    theModel->copy=0.0;
    theModel->copz=0.0;
    theModel->rmsd=0.0;
    theModel->rmsp=0.0;
 
    total_mass=0.0;
    for(i=0;i<theModel->n;i++) {
        r2 = theModel->x[i]*theModel->x[i]+
            theModel->y[i]*theModel->y[i]+
            theModel->z[i]*theModel->z[i];
        v2 = theModel->vx[i]*theModel->vx[i]+
            theModel->vy[i]*theModel->vy[i]+
            theModel->vz[i]*theModel->vz[i];
        theModel->KE += 0.5*theModel->mass[i]*v2;
        theModel->comx += theModel->mass[i]*theModel->x[i];
        theModel->comy += theModel->mass[i]*theModel->y[i];
        theModel->comz += theModel->mass[i]*theModel->z[i];
        theModel->copx += theModel->mass[i]*theModel->vx[i];
        theModel->copy += theModel->mass[i]*theModel->vy[i];
        theModel->copz += theModel->mass[i]*theModel->vz[i];
        theModel->rmsd += r2;
        theModel->rmsp += v2*theModel->mass[i]*theModel->mass[i];
        total_mass += theModel->mass[i];
    }
    theModel->comx /= total_mass;
    theModel->comy /= total_mass;
    theModel->comz /= total_mass;
    theModel->rmsd = sqrt(theModel->rmsd/theModel->n);
    theModel->rmsp = sqrt(theModel->rmsp/theModel->n);
    theModel->copx /= theModel->n;
    theModel->copy /= theModel->n;
    theModel->copz /= theModel->n;
    for(i=0;i<theModel->n;i++) {
        for(j=0;j<i;j++) {
            dx = theModel->x[i]-theModel->x[j];
            dy = theModel->y[i]-theModel->y[j];
            dz = theModel->z[i]-theModel->z[j];
            r2 = dx*dx+dy*dy+dz*dz;
            if(r2>theModel->srad2[j]&&r2>theModel->srad2[i]) {
                r = sqrt(r2+theModel->softening_factor*
                    theModel->softening_factor);
                theModel->PE -= theModel->G*
                    theModel->mass[i]*theModel->mass[j]*r2/(r*r*r);
            }
        }
    }
}

void randomizeMassesNbodyModel(NbodyModel * theModel) {
    int i;
    for(i=0;i<theModel->n;i++) {
        //theModel->mass[i]=fabs(drand_norm(theModel->default_mass,theModel->default_mass/5,0.0));
        theModel->mass[i]=drand(0.0,theModel->default_mass*2.0);
    }
}
                
double computeSoftenedRadius(double g_m, double tstep_squared,double srad_factor) {
    // g_m = G*mass;
    if(srad_factor>0.0) {
        return srad_factor*fcr(g_m*tstep_squared,0);
    } else {
        return 0.0;
    }
}
