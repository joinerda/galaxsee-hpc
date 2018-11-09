#ifndef NBODY
#define NBODY

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "pppm_structs.h"
#include "octtree_structs.h"
#include "readline.h"
#include "timer.h"

#define PI 3.14159

#define UPDATEMETHOD_HASH_TEXT 1
#define UPDATEMETHOD_BRIEF_TEXT 2
#define UPDATEMETHOD_VERBOSE_POSITIONS 4
#define UPDATEMETHOD_GD_IMAGE 8
#define UPDATEMETHOD_TEXT11 16
#define UPDATEMETHOD_X11 32
#define UPDATEMETHOD_VERBOSE_STATISTICS 64
#define UPDATEMETHOD_SDL 128
#define UPDATEMETHOD_DUMP 256
#define UPDATEMETHOD_NETCDF 512

#define INT_METHOD_RK4          1
#define INT_METHOD_LEAPFROG     2
#define INT_METHOD_MPEULER      3
#define INT_METHOD_IEULER       4
#define INT_METHOD_EULER        5
#define INT_METHOD_ABM          6

#define FORCE_METHOD_DIRECT	1
#define FORCE_METHOD_TREE	2
#define FORCE_METHOD_PPPM	3

#define DISTRIBUTION_SPHERICAL_RANDOM_UNIFORM 1
#define DISTRIBUTION_RECTANGULAR_RANDOM_UNIFORM 2
#define DISTRIBUTION_RECTANGULAR_UNIFORM 3
#define DISTRIBUTION_SPHERICAL_RANDOM_NORMAL 4

// model structure
typedef struct {
    int n;			// number of points
    double * mass;		// masses
    double * x;                 // x values
    double * y;			// y values
    double * z;			// z values
    double * vx;		// x velocities
    double * vy;		// y velocities
    double * vz;		// z velocities
    int abmCounter;		// restart counter for ABM integration
    int bhCounter;		// counter for BarnesHut method to determine
				// when to recreate tree
    double default_mass;	// default star mass used for initialization
    double default_scale;	// (1/2) the "box" side length, typical scale
    double default_G;		// Gravitational constant in given units
    double KE;			// Kinetic energy of system
    double PE;			// Potential energy of system
    double comx;		// x center of mass
    double comy;		// y center of mass
    double comz;		// z center of mass
    double copx;		// x system momentum
    double copy;		// y system momentum
    double copz;		// z system momentum
    double rmsd;		// root mean square distance from origin
    double rmsp;		// root mean square momentum
    double * draw_x;		// copy of x for update purposes
    double * draw_y;		// copy of y for update purposes
    double * draw_z;		// copy of z for update purposes
    double * srad2;		// shield radius squared
    double * X;			// Computational scratch space
    double * dbuffer;			// Computational scratch space
    double * X2;		// Computational scratch space
    double * XPRIME;		// Computational scratch space
    double * XPRIME1;		// Computational scratch space
    double * XPRIME2;		// Computational scratch space
    double * XPRIME3;		// Computational scratch space
    double * XPRIME4;		// Computational scratch space
    int * color;		// color code for GalaxSee compatibility
    double * color3;		// Display color for SDL method
    double G;			// Gravitational constant
    int default_color;		// color code for GalaxSee compatibility
    double t;			// model time
    double tFinal;		// model end time
    double tstep;		// model time step
    double tstep_ratio;               // model time step
    double tstep_min;         // model time step
    double tstep_max;         // model time step

    double tstep_last;		// previously used time step
    double srad_factor;		// coefficient for automatic sheild radius
    double softening_factor;	// softening radius
    double rotation_factor;	// initial rotation, scaled so that 1
				// is "equilibrium" rotation
    double treeRangeCoefficient;// multiple of Barnes-Hut node size scale
				// that determine how far away an object
				// must be form the current node to use
				// a tree approximation
    double initial_v;		// initial random velocity
    OctTree * rootTree;		// Barnes-Hut tree structure
    OctTree ** treeList;	// map of each star to its place in the tree
    PPPM * thePPPM;		// PPPM memory structure
    int iteration;		// current iteration
    int int_method;		// integration method
    int force_method;		// force method
    char prefix[READLINE_MAX];  // prefix for output files
    double drag;		// drag coefficient
    double expansion;		// expansion rate (velocity per distance units)
    double pppm_ksigma;		// coefficient for determining PPPM softening
    double pppm_knear; 		// coefficient for determining PPPM nearest
				// neghbor cutoff
    double anisotropy;		// additional anisotropy added to uniform
				//  distributions
    int distribution;		// initial distribution
    double distribution_x_scale;		// distribution x mag
    double distribution_y_scale;		// distribution y mag
    double distribution_z_scale;		// distribution z mag
    double distribution_sigma;          // distribution standard dev for normal
    Timer force_timer;
    Timer step_timer;
    double force_total;
    double pointsize;
} NbodyModel;

NbodyModel * allocateNbodyModel(int n,int ngrid);
int initializeNbodyModel(NbodyModel *theModel);
int freeNbodyModel(NbodyModel *theModel);
void spinNbodyModel(NbodyModel * theModel);
void speedNbodyModel(NbodyModel * theModel);
void randomizeMassesNbodyModel(NbodyModel * theModel);
int updateNbodyModel(NbodyModel *theModel,int updateMethod);
double computeSoftenedRadius(double g_m, double tstep,double srad_factor);
void printStatistics(NbodyModel *theModel);
void calcStatistics(NbodyModel *theModel);
void nbodyEvents();
void setPointsizeNbodyModel(NbodyModel *theModel,double pointsize);
void setAnisotropyNbodyModel(NbodyModel *theModel,double anisotropy);
void setDistributionNbodyModel(NbodyModel *theModel,int distribution);
void setDistributionSigmabodyModel(NbodyModel *theModel,double distribution_sigma);
void setDistributionXScaleNbodyModel(NbodyModel *theModel,double distribution_x_scale);
void setDistributionYScaleNbodyModel(NbodyModel *theModel,double distribution_y_scale);
void setDistributionZScaleNbodyModel(NbodyModel *theModel,double distribution_z_scale);
void setPPPMCoeffsNbodyModel(NbodyModel *theModel,double ksigma,double knear);
void setDragNbodyModel(NbodyModel *theModel,double drag);
void setExpansionNbodyModel(NbodyModel *theModel,double expansion);
void setPrefixNbodyModel(NbodyModel *theModel,const char * prefix);
void setDefaultsNbodyModel(NbodyModel *theModel);
void setMassNbodyModel(NbodyModel *theModel,double mass);
void setColorNbodyModel(NbodyModel *theModel,int color);
void setScaleNbodyModel(NbodyModel *theModel,double scale);
void setGNbodyModel(NbodyModel *theModel,double G);
void setTFinal(NbodyModel *theModel,double tFinal);
void setTStep(NbodyModel *theModel,double tStep);
void setTStepRatio(NbodyModel *theModel,double tStepRatio);
void setTStepMin(NbodyModel *theModel,double tStepMin);
void setTStepMax(NbodyModel *theModel,double tStepMax);
void setIntMethod(NbodyModel *theModel,int int_method);
void setForceMethod(NbodyModel *theModel,int force_method);
void setTreeRangeCoefficient(NbodyModel *theModel,double coefficient);
void setSofteningNbodyModel(NbodyModel *theModel,double softening_factor);
void setSradNbodyModel(NbodyModel *theModel,double srad_factor);
void setRotationFactor(NbodyModel *theModel,double rotation_factor);
void setInitialV(NbodyModel *theModel,double initial_v);
void calcDerivs(double * x, double * derivs, double t, double tStep,
    NbodyModel * theModel);
void calcDerivsDirect(double * x, double * derivs, double t, double tStep,
    NbodyModel * theModel);
void calcDerivsBarnesHut(double * x, double * derivs, double t, double tStep,
    NbodyModel * theModel);
void calcDerivsPPPM(double * x, double * derivs, double t, double tStep,
    NbodyModel * theModel);
int stepNbodyModelEuler(NbodyModel * theModel, double tStep);
int stepNbodyModelIEuler(NbodyModel * theModel, double tStep);
int stepNbodyModelMPEuler(NbodyModel * theModel, double tStep);
int stepNbodyModelRK4(NbodyModel * theModel, double tStep);
int stepNbodyModelLeapfrog(NbodyModel * theModel, double tStep);
int stepNbodyModelABM(NbodyModel * theModel, double tStep);
int stepNbodyModel(NbodyModel * theModel);
void copy2X(NbodyModel *theModel);
void copy2xyz(NbodyModel *theModel);

#endif
