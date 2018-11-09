#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nbody.h"
#include "text11.h"

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

NbodyModel * allocateNbodyModel(int n) {
    NbodyModel * theModel;
    theModel = (NbodyModel *)malloc(sizeof(NbodyModel));
    theModel->n=n;
    theModel->mass = (double *)malloc(sizeof(double)*n);
    theModel->x = (double *)malloc(sizeof(double)*n);
    theModel->y = (double *)malloc(sizeof(double)*n);
    theModel->z = (double *)malloc(sizeof(double)*n);
    theModel->vx = (double *)malloc(sizeof(double)*n);
    theModel->vy = (double *)malloc(sizeof(double)*n);
    theModel->vz = (double *)malloc(sizeof(double)*n);
    theModel->srad2 = (double *)malloc(sizeof(double)*n);
    theModel->tstep_last=-1.0;
    theModel->X = (double *)malloc(sizeof(double)*n*6);
    theModel->X2 = (double *)malloc(sizeof(double)*n*6);
    theModel->XPRIME = (double *)malloc(sizeof(double)*n*6);
    theModel->XPRIME1 = (double *)malloc(sizeof(double)*n*6);
    theModel->XPRIME2 = (double *)malloc(sizeof(double)*n*6);
    theModel->XPRIME3 = (double *)malloc(sizeof(double)*n*6);
    theModel->XPRIME4 = (double *)malloc(sizeof(double)*n*6);

    return theModel;
}

void setDefaultsNbodyModel(NbodyModel *theModel) {
    theModel->default_mass=1.0; // solar masses
    theModel->default_scale=10.0; // parsecs
    theModel->default_G=0.0046254; // pc^2/Solar_mass/My^2
    theModel->srad_factor=5.0;
    theModel->softening_factor=0.0;
}
void setSofteningNbodyModel(NbodyModel *theModel,double softening_factor) {
    theModel->softening_factor=softening_factor;
}
void setSradNbodyModel(NbodyModel *theModel,double srad_factor) {
    theModel->srad_factor=srad_factor;
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
    int i,j;
    double r2;
    double scale=theModel->default_scale;
    double r2Max=scale*scale;
    for(i=0;i<theModel->n;i++) {
        theModel->mass[i] = theModel->default_mass;
        do {
            theModel->x[i] = scale*(2.0*(double)rand()/(double)RAND_MAX-1.0);
            theModel->y[i] = scale*(2.0*(double)rand()/(double)RAND_MAX-1.0);
            theModel->z[i] = scale*(2.0*(double)rand()/(double)RAND_MAX-1.0);
            r2 = theModel->x[i]*theModel->x[i]+
                theModel->y[i]*theModel->y[i]+
                theModel->z[i]*theModel->z[i];
        } while (r2>r2Max) ;
        theModel->vx[i]=0.0;
        theModel->vy[i]=0.0;
        theModel->vz[i]=0.0;
        for (j=0;j<6;j++) {
            theModel->X[i*6+j]=0.0;
            theModel->XPRIME[i*6+j]=0.0;
        }
    }
    theModel->G=theModel->default_G;
    theModel->t=0.0;
    theModel->iteration=0;
    theModel->abmCounter=-3;

    return 1;
}

int freeNbodyModel(NbodyModel *theModel) {
    free(theModel->x);
    free(theModel->y);
    free(theModel->z);
    free(theModel->vx);
    free(theModel->vy);
    free(theModel->vz);
    free(theModel->srad2);
    free(theModel->X);
    free(theModel->X2);
    free(theModel->XPRIME);
    free(theModel->XPRIME1);
    free(theModel->XPRIME2);
    free(theModel->XPRIME3);
    free(theModel->XPRIME4);
    free(theModel);

    return 1;
}

void calcDerivs(double * x, double * derivs, double t,double tStep,
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
int stepNbodyModel(NbodyModel * theModel, double tStep, int int_method) {
    int i;
    double srad;
    if(tStep!=theModel->tstep_last) {
        for (i=0;i<theModel->n;i++) {
            srad = computeSoftenedRadius(theModel->G*theModel->mass[i],
                tStep,theModel->srad_factor);
            theModel->srad2[i]=srad*srad;
        }
        theModel->tstep_last=tStep;
        theModel->abmCounter=-3;
    }
    switch(int_method) {
        case INT_METHOD_ABM:
            return stepNbodyModelABM(theModel,tStep);
         break;
        case INT_METHOD_RK4:
            return stepNbodyModelRK4(theModel,tStep);
         break;
        case INT_METHOD_LEAPFROG:
            return stepNbodyModelLeapfrog(theModel,tStep);
         break;
        case INT_METHOD_MPEULER:
            return stepNbodyModelMPEuler(theModel,tStep);
         break;
        case INT_METHOD_IEULER:
            return stepNbodyModelIEuler(theModel,tStep);
         break;
        case INT_METHOD_EULER:
        default:
            return stepNbodyModelEuler(theModel,tStep);
         break;
    }
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
            theModel->X[i*6+3] += theModel->XPRIME[i*6+3]*tStep;
            theModel->X[i*6+4] += theModel->XPRIME[i*6+4]*tStep;
            theModel->X[i*6+5] += theModel->XPRIME[i*6+5]*tStep;
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
    double * fk3;
    double * fk2;
    double * fk1;
    double * fk0;
    double * fkp;

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
        XGrayscale_x11[i].red=color;
        XGrayscale_x11[i].green=color;
        XGrayscale_x11[i].blue=color;
        XAllocColor(dpy_x11,theColormap_x11,&(XGrayscale_x11[i]));
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
    double scale,shift;
    int dispX, dispY, dispZ, depthY, depthZ;
     
    XSetForeground(dpy_x11, gc_x11, black_x11);
    XFillRectangle(dpy_x11,buffer_x11,gc_x11,
        0,0,width_x11,height_x11);

    scale=(double)height_x11/(2.0*theModel->default_scale);
    shift=(double)width_x11/4.0;
    for (i=0 ; i<theModel->n; i++) {
        dispX = (int)(scale*theModel->x[i]+shift);
        dispY = (int)(scale*theModel->y[i]+shift);
        dispZ = (int)(scale*theModel->z[i]+shift);
        
        depthY = (int)((double)dispY/(double)height_x11*numXGrayscale_x11);
        depthZ = (int)((double)dispZ/(double)height_x11*numXGrayscale_x11);
        if (depthY>numXGrayscale_x11-1) depthY=numXGrayscale_x11-1;
        if (depthZ>numXGrayscale_x11-1) depthZ=numXGrayscale_x11-1;

        if (dispX < width_x11/2) {
            XSetForeground(dpy_x11,gc_x11,XGrayscale_x11[depthZ].pixel);
            XFillRectangle(dpy_x11,buffer_x11,gc_x11,dispX,dispY,3,3);
        }
        if (dispX > 0) {
            XSetForeground(dpy_x11,gc_x11,XGrayscale_x11[depthY].pixel);
            XFillRectangle(dpy_x11,buffer_x11,gc_x11,
                dispX+width_x11/2,dispZ,3,3);
        }
   }
   XCopyArea(dpy_x11, buffer_x11, w_x11, gc_x11, 0, 0,
       width_x11, height_x11,  0, 0);
   XFlush(dpy_x11);
}
#endif

void speedNbodyModel(NbodyModel *theModel, double v) {
    double pi=atan(1.0)*4.0;
    double theta,phi;
    int i;
    for (i=0;i<theModel->n;i++) {
        theta = 2*pi*(double)rand()/(double)RAND_MAX;
        phi = pi*(double)rand()/(double)RAND_MAX;
        theModel->vx[i]+=v*sin(phi)*cos(theta);
        theModel->vy[i]+=v*sin(phi)*sin(theta);
        theModel->vz[i]+=v*cos(phi);
    }
}

void spinNbodyModel(NbodyModel *theModel, double rotFactor) {
    int i,j;
    double ri,rj,sumM,vTan;

    for (i=0;i<theModel->n;i++) {
        ri = sqrt(pow(theModel->x[i],2.0)+pow(theModel->y[i],2.0));
        sumM=0.0;
        for(j=0;j<theModel->n;j++) {
            if(i!=j) {
                rj = sqrt(pow(theModel->x[j],2.0)+pow(theModel->y[j],2.0));
                if(rj<ri) {
                    sumM+=theModel->mass[j];
                }
            }
        }
        if(sumM>0.0) {
            vTan = sqrt(sumM*theModel->G/ri);
            theModel->vz[i]+=0.0;
            theModel->vx[i]+=-rotFactor*vTan*theModel->y[i]/ri;
            theModel->vy[i]+=rotFactor*vTan*theModel->x[i]/ri;
        } 
    }
    return;
}

int updateNbodyModel(NbodyModel *theModel,int updateMethod) {
    Text11 * t11p;
    int i;

#ifdef HAS_LIBGD
    gdImagePtr im;
    FILE *pngout;
    char fname[80];
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
    int gd_point_size=2;
#endif

    // output, send to client, send to screen, etc.
    if((updateMethod&UPDATEMETHOD_HASH_TEXT)==UPDATEMETHOD_HASH_TEXT) {
        printf("#");
        fflush(stdout);
    } 

    if((updateMethod&UPDATEMETHOD_BRIEF_TEXT)==UPDATEMETHOD_BRIEF_TEXT) {
        printf("Calculated t = %10.3e\n",theModel->t);
    }

    if((updateMethod&UPDATEMETHOD_VERBOSE_POSITIONS)==UPDATEMETHOD_VERBOSE_POSITIONS) {
        for(i=0;i<theModel->n;i++) {
            printf(
                "POS %10.3e\t%d\t%10.3e\t%10.3e\t%10.3e\t%10.3e\t%10.3e\t%10.3e\n",
                theModel->t,i,theModel->x[i],theModel->y[i],
                theModel->z[i],theModel->vx[i],theModel->vy[i],
                theModel->vz[i]);
        }
    }
    if((updateMethod&UPDATEMETHOD_VERBOSE_STATISTICS)==UPDATEMETHOD_VERBOSE_STATISTICS) {
        calcStatistics(theModel);
        printf(
            "STAT %10.3e\t%10.3e\t%10.3e\t%10.3e\n",
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
                text11_add(t11p,theModel->x[i],theModel->y[i]);
            }
            text11_print(t11p);
            printf("==============================\n");
            text11_free(&t11p);
    }

#ifdef HAS_LIBGD
    if((updateMethod&UPDATEMETHOD_GD_IMAGE)==UPDATEMETHOD_GD_IMAGE) {
        // set up GD image
        im = gdImageCreate(width,height);
        black = gdImageColorAllocate(im, 0, 0, 0);
        white = gdImageColorAllocate(im, 255, 255, 255);
    
        gdImageFilledRectangle(im,0,0,width,height,white);
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
        sprintf(fname,"out%03d.png",theModel->iteration);
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

    return 1;
}

void printStatistics(NbodyModel * theModel) {
    double total_energy;
    calcStatistics(theModel);

    total_energy = theModel->KE + theModel->PE;
    printf("KE = %le \t",theModel->KE);
    printf("PE = %le \t",theModel->PE);
    printf("TE = %le \n",total_energy);
}
void calcStatistics(NbodyModel * theModel) {
    int i,j;
    double v2,r2,r,dx,dy,dz;

    theModel->KE=0.0;
    theModel->PE=0.0;
 
    for(i=0;i<theModel->n;i++) {
        v2 = theModel->vx[i]*theModel->vx[i]+
            theModel->vy[i]*theModel->vy[i]+
            theModel->vz[i]*theModel->vz[i];
        theModel->KE += 0.5*theModel->mass[i]*v2;
    }
    for(i=0;i<theModel->n;i++) {
        for(j=0;j<i;j++) {
            dx = theModel->x[i]-theModel->x[j];
            dy = theModel->y[i]-theModel->y[j];
            dz = theModel->z[i]-theModel->z[j];
            r2 = dx*dx+dy*dy+dz*dz;
            r = sqrt(r2+theModel->softening_factor*theModel->softening_factor);
            if(r2>theModel->srad2[j]&&r2>theModel->srad2[i]) {
                theModel->PE -= theModel->G*
                    theModel->mass[i]*theModel->mass[j]*r2/(r*r*r);
            }
        }
    }
}
                
double computeSoftenedRadius(double g_m, double tstep,double srad_factor) {
    // g_m = G*mass;
    if(srad_factor>0.0) {
        return srad_factor*pow(g_m,0.333)*pow(tstep,0.667);
    } else {
        return 0.0;
    }
}
