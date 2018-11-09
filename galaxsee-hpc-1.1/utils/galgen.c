#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mem.h"
#include "rand_tools.h"
#include "I_qsort.h"

#define DIST_NORMAL 0
#define DIST_UNIFORM 1
#define USE_INTERIOR_MASS 0
#define USE_CENTRIFUGAL_FORCE 1
#define USE_DISC_APPROX 2
#define USE_SPHERE_APPROX 3

double scale_factor(double x, double c,double sigma) {
    return 0.5*(erf((x-c)/sigma)+1.0);
}

int compare (const void * a, const void * b) {
  double diff;
  diff = ( *(double*)a - *(double*)b );
  if(diff<0.0) return -1;
  else if(diff>0.0) return 1;
  else if(diff==0.0) return 0;
  else {printf("ERROR IN COMPARE\n"); exit(0); }
}

/*
N 5000                          # number of masses
ROTATION_FACTOR 0.8             # unitless rotation factor (equilibrium at 1.0)
SCALE  500.0                    # 1/2 the "box" side length
MASS 200000.0                   # total system mass
G 0.00449                       # Gravitational Constant
*/

int main(int argc, char ** argv)
{
    int n=10;
    int i;

    double *x;
    double *y;
    double *z;
    double *vx;
    double *vy;
    double *vz;
    double px,py,pz;
    double *mass;
    double *rsphere;
    double *rplane;
    double mass_interior;
    double G = 0.00449;
    double scale = 766.0;
    int * ix;
    double mass_total=200000.0;
    double rotation_factor=1.0;
    double frac = 0.0;
    double offset, sigma;
    int j;
    int use_center=0;
    double center_fraction=0.25;
    int dist = DIST_NORMAL;
    int vr_method = USE_CENTRIFUGAL_FORCE;
    double elliptical;
    double xscale=1.0;
    double yscale=1.0;
    double zscale=1.0;

    i=1;
    if(argc>i) sscanf(argv[i++],"%d",&n);
    if(argc>i) sscanf(argv[i++],"%lf",&xscale);
    if(argc>i) sscanf(argv[i++],"%lf",&yscale);
    if(argc>i) sscanf(argv[i++],"%lf",&zscale);
    if(argc>i) sscanf(argv[i++],"%lf",&frac);
    if(argc>i) sscanf(argv[i++],"%d",&dist);
    if(argc>i) sscanf(argv[i++],"%lf",&mass_total);
    if(argc>i) sscanf(argv[i++],"%lf",&scale);
    if(argc>i) sscanf(argv[i++],"%d",&vr_method);
    if(argc>i) sscanf(argv[i++],"%d",&rotation_factor);
    if (argc==1||argc>i) { 
        printf("#### USAGE n xscale yscale zscale frac dist mass_total scale vr_method rotation_factor\n");
        printf("#### DIST_NORMAL 0\n");
        printf("#### DIST_UNIFORM 1\n");
        printf("#### USE_INTERIOR_MASS 0\n");
        printf("#### USE_CENTRIFUGAL_FORCE 1\n");
        printf("#### USE_DISC_APPROX 2\n");
        printf("#### USE_SPHERE_APPROX 3\n");
        exit(0);
    }

    offset = scale*frac;
    sigma = scale*frac;

    printf("N %d\n",n);
    printf("#### XSCALE %lf\n",xscale);
    printf("#### YSCALE %lf\n",yscale);
    printf("#### ZSCALE %lf\n",zscale);
    printf("#### FRAC %lf\n",frac);
    printf("#### DIST %d\n",dist);
    printf("MASS %lf\n",mass_total);
    printf("SCALE %lf\n",scale);
    printf("#### VR_METHOD %d\n",vr_method);
    printf("ROTATION_FACTOR %lf\n",rotation_factor);
    printf("G %lf\n",G);
    printf("#### USE_CENTER %d\n",use_center);
    printf("#### CENTER_FRACTION %lf\n",center_fraction);

    seed_by_time(0);
    ix = (int *)malloc(sizeof(int)*n);
    x = (double *)malloc(sizeof(double)*9*n);
    y = &(x[n]);
    z = &(x[2*n]);
    vx = &(x[3*n]);
    vy = &(x[4*n]);
    vz = &(x[5*n]);
    mass = &(x[6*n]);
    rsphere = &(x[7*n]);
    rplane = &(x[8*n]);

    // initialize values for points
    for(i=0;i<n;i++) {
        if(dist==DIST_NORMAL) {
            x[i] = drand_norm(0.0,0.5*xscale*scale,1.0e-3);
            y[i] = drand_norm(0.0,0.5*yscale*scale,1.0e-3);
            z[i] = drand_norm(0.0,0.5*zscale*scale,1.0e-3);
        } else {
            x[i] = drand(-xscale*scale,xscale*scale);
            y[i] = drand(-yscale*scale,yscale*scale);
            z[i] = drand(-zscale*scale,zscale*scale);
            elliptical = pow(x[i]/(xscale*scale),2.0)+
                         pow(y[i]/(yscale*scale),2.0)+
                         pow(z[i]/(zscale*scale),2.0);
            while(elliptical>1.0) {
                x[i] = drand(-xscale*scale,xscale*scale);
                y[i] = drand(-yscale*scale,yscale*scale);
                z[i] = drand(-zscale*scale,zscale*scale);
                elliptical = pow(x[i]/(xscale*scale),2.0)+
                             pow(y[i]/(yscale*scale),2.0)+
                             pow(z[i]/(zscale*scale),2.0);
            }
        }

        vx[i] = 0.0;
        vy[i] = 0.0;
        vz[i] = 0.0;
        if(use_center) {
            if(i==0) {
                x[i]=0.0;
                y[i]=0.0;
                z[i]=0.0;
                mass[i] = mass_total*center_fraction;
            } else {
                mass[i] = mass_total*(1.0-center_fraction)/(double)(n-1);
            }
        } else {
            mass[i] = mass_total/(double)(n);
        }
        rsphere[i] = sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]);
        rplane[i] = sqrt(x[i]*x[i]+y[i]*y[i]);
    }

    // Set up rotational velocities, so that rotational velocity
    // is equivalent to that needed to counter the central
    // acceleration.

    if(vr_method==USE_INTERIOR_MASS) {
        I_qsort (rplane, ix, n, sizeof(double), compare);
        mass_interior=0.0;
    }
    for(i=0;i<n;i++) {
        double vr,ax,ay,az,dx,dy,dz,rij,rij3;
        if(vr_method==USE_INTERIOR_MASS) {
            int ind=ix[i];
            int indm1=ix[i-1];
            mass_interior += mass[indm1];
            if(rplane[ind]>0.0) {
                vr = rotation_factor*scale_factor(rplane[ind],offset,sigma)*
                    sqrt(G*mass_interior/rsphere[ind]*(rplane[ind]/rsphere[ind]))/rplane[ind];
            } else {
                vr = 0.0;
            }
            vx[ind] = vr*y[ind];
            vy[ind] = -vr*x[ind];
            vz[ind]=0.0;
        } else if(vr_method==USE_SPHERE_APPROX) {
            if(rplane[i]>0.0) {
                vr = rotation_factor*scale_factor(rplane[i],offset,sigma)*
                    sqrt(G*mass_total/scale/scale/scale*
                    rsphere[i]*rplane[i])/rplane[i];
            } else {
                vr = 0.0;
            }
            vx[i] = vr*y[i];
            vy[i] = -vr*x[i];
            vz[i]=0.0;
        } else if(vr_method==USE_DISC_APPROX) {
            if(rplane[i]>0.0) {
                double tfac = scale*scale-rplane[i]*rplane[i];
                tfac = pow(tfac,1.5);
                tfac = scale*scale*scale-tfac;
                tfac = 2.0/3.0*M_PI*tfac;
                tfac = tfac*G*mass_total/scale/scale/scale;
                vr = rotation_factor*scale_factor(rplane[i],offset,sigma)*
                    sqrt(tfac/rplane[i])
                    /rplane[i];
            } else {
                vr = 0.0;
            }
            vx[i] = vr*y[i];
            vy[i] = -vr*x[i];
            vz[i]=0.0;
        } else if(vr_method==USE_CENTRIFUGAL_FORCE) {
            ax=0;ay=0;az=0;
            for(j=0;j<n;j++) {
                if(i!=j) {
                    dx = x[j]-x[i];
                    dy = y[j]-y[i];
                    dz = z[j]-z[i];
                    rij = sqrt(dx*dx+dy*dy+dz*dz);
                    rij3 = G*mass[j]/(rij*rij*rij);
                    ax += rij3*dx;
                    ay += rij3*dy;
                    az += rij3*dz;
                }
            }
            if(rplane[i]>0.0) {
                vr=rotation_factor*sqrt(sqrt(ax*ax+ay*ay+az*az)*
                    rplane[i]*rplane[i]/rsphere[i])/
                    rplane[i]*scale_factor(rplane[i],offset,sigma);
            } else {
                vr = 0.0;
            }
            vx[i] = vr*y[i];
            vy[i] = -vr*x[i];
            vz[i]=0.0;
        }
    }

    // correct for center of momentum, ensure that initial net
    // momentum is zero.
    px = 0.0;
    py = 0.0;
    pz = 0.0;
    for(i=0;i<n;i++) {
        px += mass[i]*vx[i];
        py += mass[i]*vy[i];
        pz += mass[i]*vz[i];
    }
    for(i=0;i<n;i++) {
        vx[i] -= px/mass_total;
        vy[i] -= py/mass_total;
        vz[i] -= pz/mass_total;
    }

    for(i=0;i<n;i++) {
        printf("COORDS %lf %lf %lf %lf %lf %lf %lf ",x[i],y[i],z[i],
            vx[i],vy[i],vz[i],mass[i]);
        printf("%d\n",(int)(rsphere[i]/scale*10.0));
    }

    free(x);
    free(ix);
}

