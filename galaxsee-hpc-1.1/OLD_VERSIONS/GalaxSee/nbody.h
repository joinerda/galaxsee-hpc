#ifndef NBODY
#define NBODY

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define PI 3.14159

#define UPDATEMETHOD_HASH_TEXT 1
#define UPDATEMETHOD_BRIEF_TEXT 2
#define UPDATEMETHOD_VERBOSE_POSITIONS 4
#define UPDATEMETHOD_GD_IMAGE 8
#define UPDATEMETHOD_TEXT11 16
#define UPDATEMETHOD_X11 32
#define UPDATEMETHOD_VERBOSE_STATISTICS 64

#define INT_METHOD_RK4          1
#define INT_METHOD_LEAPFROG     2
#define INT_METHOD_MPEULER      3
#define INT_METHOD_IEULER       4
#define INT_METHOD_EULER        5
#define INT_METHOD_ABM          6


// model structure
typedef struct {
    int n;
    int abmCounter;
    double * mass;
    double default_mass;
    double default_scale;
    double default_G;
    double KE;
    double PE;
    double comx;
    double comy;
    double comz;
    double copx;
    double copy;
    double copz;
    double * x;
    double * y;
    double * z;
    double * vx;
    double * vy;
    double * vz;
    double * srad2;
    double * X;
    double * X2;
    double * XPRIME;
    double * XPRIME1;
    double * XPRIME2;
    double * XPRIME3;
    double * XPRIME4;
    double G;
    double t;
    double tstep_last;
    double srad_factor;
    double softening_factor;
    int iteration;
} NbodyModel;

NbodyModel * allocateNbodyModel(int n);
void setDefaultsNbodyModel(NbodyModel *theModel);
void setMassNbodyModel(NbodyModel *theModel,double mass);
void setScaleNbodyModel(NbodyModel *theModel,double scale);
void setGNbodyModel(NbodyModel *theModel,double G);
int initializeNbodyModel(NbodyModel *theModel);
int freeNbodyModel(NbodyModel *theModel);
void calcDerivs(double * x, double * derivs, double t, double tStep,
    NbodyModel * theModel);
int stepNbodyModelEuler(NbodyModel * theModel, double tStep);
int stepNbodyModelIEuler(NbodyModel * theModel, double tStep);
int stepNbodyModelMPEuler(NbodyModel * theModel, double tStep);
int stepNbodyModelRK4(NbodyModel * theModel, double tStep);
int stepNbodyModelLeapfrog(NbodyModel * theModel, double tStep);
int stepNbodyModel(NbodyModel * theModel, double tStep,int method);
void spinNbodyModel(NbodyModel * theModel, double rotFactor);
void speedNbodyModel(NbodyModel * theModel, double v);
int updateNbodyModel(NbodyModel *theModel,int updateMethod);
double computeSoftenedRadius(double g_m, double tstep,double srad_factor);
void printStatistics(NbodyModel *theModel);
void calcStatistics(NbodyModel *theModel);
void copy2X(NbodyModel *theModel);
void copy2xyz(NbodyModel *theModel);
void setSofteningNbodyModel(NbodyModel *theModel,double softening_factor);
void setSradNbodyModel(NbodyModel *theModel,double srad_factor);



#endif
