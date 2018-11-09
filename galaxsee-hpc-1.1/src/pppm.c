#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "rand_tools.h"
#ifdef HAS_FFTW3
#include <fftw3.h>
#endif
#include "cubeinterp.h"
#include "pppm.h"
#include "pppm_structs.h"
 #ifdef HAS_MPI
 #include <mpi.h>
 #endif

void freePPPM(PPPM * thePPPM) {
#ifdef HAS_FFTW3
    int i;
    free(thePPPM->x);
    free(thePPPM->y);
    free(thePPPM->z);
    free(thePPPM->kx);
    free(thePPPM->ky);
    free(thePPPM->kz);
    fftw_free(thePPPM->density);
    fftw_free(thePPPM->dbuffer);
    fftw_free(thePPPM->den_check);
    fftw_free(thePPPM->potential);
    free(thePPPM->force);
    fftw_free(thePPPM->fft_density);
    fftw_free(thePPPM->fft_potential);
    for(i=0;i<thePPPM->nxyz;i++) {
        if(thePPPM->cell_contains[i]!=NULL) {
            free(thePPPM->cell_contains[i]);
        }
    }
    free(thePPPM->cell_contains);
    free(thePPPM->n_cell);
    free(thePPPM->max_cell);
    fftw_destroy_plan(thePPPM->pf);
    fftw_destroy_plan(thePPPM->pb);
    free(thePPPM);
#else
    printf("Warning, pppm methods being called without FFTW3, recompile with FFTW3 and run again\n");
#endif
}

void allocPPPM(PPPM * thePPPM,int nx, int ny, int nz) {
#ifdef HAS_FFTW3
    int i;
    thePPPM->nx=nx;
    thePPPM->ny=ny;
    thePPPM->nxy=nx*ny;
    thePPPM->nyz=ny*nz;
    thePPPM->nz=nz;
    thePPPM->nxyz=thePPPM->nxy*nz;
    thePPPM->nxyzo2=nx*ny*(nz/2+1);
    thePPPM->nyzo2=ny*(nz/2+1);
    thePPPM->nzo2=(nz/2+1);
    thePPPM->x = (double *)malloc(sizeof(double)*thePPPM->nx);
    thePPPM->y = (double *)malloc(sizeof(double)*thePPPM->ny);
    thePPPM->z = (double *)malloc(sizeof(double)*thePPPM->nz);
    thePPPM->kx = (double *)malloc(sizeof(double)*thePPPM->nx);
    thePPPM->ky = (double *)malloc(sizeof(double)*thePPPM->ny);
    thePPPM->kz = (double *)malloc(sizeof(double)*thePPPM->nz);
    thePPPM->density = (double *)fftw_malloc(sizeof(double)*thePPPM->nxyz);
    thePPPM->dbuffer = (double *)fftw_malloc(sizeof(double)*thePPPM->nxyz);
    thePPPM->den_check = (double *)fftw_malloc(sizeof(double)*thePPPM->nxyz);
    thePPPM->potential = (double *)fftw_malloc(sizeof(double)*thePPPM->nxyz);
    thePPPM->force = (point3d *)malloc(sizeof(point3d)*thePPPM->nxyz);
    thePPPM->fft_density = (fftw_complex *)
        fftw_malloc(sizeof(fftw_complex)*thePPPM->nxyzo2);
    thePPPM->fft_potential = (fftw_complex *)
        fftw_malloc(sizeof(fftw_complex)*thePPPM->nxyzo2);
    thePPPM->cell_contains = (int **) malloc(sizeof(int *)*thePPPM->nxyz);
    for(i=0;i<thePPPM->nxyz;i++) thePPPM->cell_contains[i]=NULL;
    thePPPM->n_cell = (int *)malloc(sizeof(int)*thePPPM->nxyz);
    thePPPM->max_cell = (int *)malloc(sizeof(int)*thePPPM->nxyz);
    thePPPM->pf = fftw_plan_dft_r2c_3d(thePPPM->nx,thePPPM->ny,thePPPM->nz,
        thePPPM->density,thePPPM->fft_density,FFTW_ESTIMATE);
    thePPPM->pb = fftw_plan_dft_c2r_3d(thePPPM->nx,thePPPM->ny,thePPPM->nz,
        thePPPM->fft_potential,thePPPM->potential,FFTW_ESTIMATE);
#else
    printf("Warning, pppm methods being called without FFTW3, recompile with FFTW3 and run again\n");
#endif
}

void pushCellPPPM(PPPM * thePPPM,int i, int j, int k, int l) {
#ifdef HAS_FFTW3
    int nyz = thePPPM->nyz;
    int nz = thePPPM->nz;
    int n,m;
    int * copy;
    if(i<0||i>thePPPM->nx-1 ||
       j<0||j>thePPPM->ny-1 ||
       k<0||k>thePPPM->nz-1 ) {
        printf("WARNING!!!!!!! %d %d %d\n",i,j,k);
    }
    if(thePPPM->cell_contains[i*nyz+j*nz+k]==NULL) {
        thePPPM->n_cell[i*nyz+j*nz+k] = 1;
        thePPPM->max_cell[i*nyz+j*nz+k] = 8;
        thePPPM->cell_contains[i*nyz+j*nz+k] =
            (int *)malloc(sizeof(int)*thePPPM->max_cell[i*nyz+j*nz+k]);
        if(thePPPM->cell_contains[i*nyz+j*nz+k]==NULL) {
            printf("OUT OF MEMORY ERROR 1\n");
        }
        thePPPM->cell_contains[i*nyz+j*nz+k][0]=l;
    } else {
        n = thePPPM->n_cell[i*nyz+j*nz+k];
        if(n+1>thePPPM->max_cell[i*nyz+j*nz+k]) {
            thePPPM->max_cell[i*nyz+j*nz+k] *= 2;
            copy =(int*)malloc(sizeof(int)*(thePPPM->max_cell[i*nyz+j*nz+k]));
            if(copy==NULL) {
                printf("OUT OF MEMORY ERROR 1\n");
            }
            for(m=0;m<n;m++) {
                copy[m]=thePPPM->cell_contains[i*nyz+j*nz+k][m];
            }
            free(thePPPM->cell_contains[i*nyz+j*nz+k]);
            thePPPM->cell_contains[i*nyz+j*nz+k]=copy;
        }
        thePPPM->cell_contains[i*nyz+j*nz+k][n]=l;
        thePPPM->n_cell[i*nyz+j*nz+k]=n+1;
    }
#else
    printf("Warning, pppm methods being called without FFTW3, recompile with FFTW3 and run again\n");
#endif
}


double i2rPPPM(int i,double min, double range, int n) {
    return min+range*(double)i/(double)(n);
}
int r2iPPPM(double r,double min, double range, int n) {
    return (int)((r-min)/range*(double)(n));
}

void populateDensityPPPM(PPPM * thePPPM, NbodyModel * theModel,
        double * x,
        double xmin, double xmax, double ymin, double ymax,
        double zmin, double zmax, double sigma) {
#ifdef HAS_FFTW3
    int i,j,k,l,im,jm,km,ip,jp,kp,ir,jr,kr;
    double xl,yl,zl;
    double min,range;
    double s2,r2,svoli,dx,dy,dz;
    double soft=sigma*3.0;
 #ifdef HAS_MPI
    extern int rank;
    extern int size;
 #endif

    for(i=0;i<thePPPM->nxyz;i++) {
        //if(thePPPM->cell_contains[i]!=NULL) {
        //    free(thePPPM->cell_contains[i]);
        //    thePPPM->cell_contains[i]=NULL;
        //}
        thePPPM->n_cell[i]=0;
    }

    thePPPM->xmin=xmin;
    thePPPM->xmax=xmax;
    thePPPM->ymin=ymin;
    thePPPM->ymax=ymax;
    thePPPM->zmin=zmin;
    thePPPM->zmax=zmax;
    min=xmin;
    range=xmax-xmin;
    for(i=0;i<thePPPM->nx;i++) {thePPPM->x[i]=i2rPPPM(i,min,range,thePPPM->nx); }
    min=ymin;
    range=ymax-ymin;
    for(j=0;j<thePPPM->ny;j++) {thePPPM->y[j]=i2rPPPM(j,min,range,thePPPM->ny); }
    min=zmin;
    range=zmax-zmin;
    for(k=0;k<thePPPM->nz;k++) {thePPPM->z[k]=i2rPPPM(k,min,range,thePPPM->nz); }
    for (i=0;i<thePPPM->nx;i++) {
        for(j=0;j<thePPPM->ny;j++) {
            for(k=0;k<thePPPM->nz;k++) {
                    thePPPM->density[i*thePPPM->nyz+j*thePPPM->nz+k] = 0.0;
            }
        }
    }

 #ifdef HAS_MPI
    for(l=rank;l<theModel->n;l+=size) {
 #else
    for(l=0;l<theModel->n;l+=1) {
 #endif
        xl = x[6*l];
        yl = x[6*l+1];
        zl = x[6*l+2];
        i = r2iPPPM(xl,xmin,xmax-xmin,thePPPM->nx);
        j = r2iPPPM(yl,ymin,ymax-ymin,thePPPM->ny);
        k = r2iPPPM(zl,zmin,zmax-zmin,thePPPM->nz);
        i = (i+thePPPM->nx)%thePPPM->nx;
        j = (j+thePPPM->ny)%thePPPM->ny;
        k = (k+thePPPM->nz)%thePPPM->nz;
        pushCellPPPM(thePPPM,i,j,k,l);
        im= r2iPPPM(xl-soft,xmin,xmax-xmin,thePPPM->nx);
        ip= r2iPPPM(xl+soft,xmin,xmax-xmin,thePPPM->nx);
        jm= r2iPPPM(yl-soft,ymin,ymax-ymin,thePPPM->ny);
        jp= r2iPPPM(yl+soft,ymin,ymax-ymin,thePPPM->ny);
        km= r2iPPPM(zl-soft,zmin,zmax-zmin,thePPPM->nz);
        kp= r2iPPPM(zl+soft,zmin,zmax-zmin,thePPPM->nz);
        while(im>i) im -= thePPPM->nx;
        while(ip<i) ip += thePPPM->nx;
        while(jm>j) jm -= thePPPM->ny;
        while(jp<j) jp += thePPPM->ny;
        while(km>k) km -= thePPPM->nz;
        while(kp<k) kp += thePPPM->nz;
        s2=soft*soft;
        svoli = 3.0/(4.0*3.14159*soft*s2);
        for (i=im;i<ip;i++) {
            ir = (i+thePPPM->nx)%thePPPM->nx;
            if(ir<i) { //point wrapped across right boundary
                dx = xl-(thePPPM->x[ir]+(thePPPM->xmax-thePPPM->xmin));
            } else if (ir>i) { //point wrapped across left boundary
                dx = xl-(thePPPM->x[ir]-(thePPPM->xmax-thePPPM->xmin));
            } else {
                dx = xl-thePPPM->x[ir];
            }
            for(j=jm;j<jp;j++) {
                jr = (j+thePPPM->ny)%thePPPM->ny;
                if(jr<j) { //point wrapped across right boundary
                    dy = yl-(thePPPM->y[jr]+(thePPPM->ymax-thePPPM->ymin));
                } else if (jr>j) { //point wrapped across left boundary
                    dy = yl-(thePPPM->y[jr]-(thePPPM->ymax-thePPPM->ymin));
                } else {
                    dy = yl-thePPPM->y[jr];
                }
                for(k=km;k<kp;k++) {
                    kr = (k+thePPPM->nz)%thePPPM->nz;
                    if(kr<k) { //point wrapped across right boundary
                        dz = zl-(thePPPM->z[kr]+(thePPPM->zmax-thePPPM->zmin));
                    } else if (kr>k) { //point wrapped across left boundary
                        dz = zl-(thePPPM->z[kr]-(thePPPM->zmax-thePPPM->zmin));
                    } else {
                        dz = zl-thePPPM->z[kr];
                    }
                    r2 = dx*dx+dy*dy+dz*dz;
                    if(r2<s2) {
                        thePPPM->density[ir*thePPPM->nyz+jr*thePPPM->nz+kr]
                            += theModel->mass[l]*pow(1.0/sigma/sqrt(M_PI),3.0)*
                               exp(-(r2/sigma/sigma));
                    }
                }
            }
        }
    }
 #ifdef HAS_MPI
    MPI_Allreduce(thePPPM->density,thePPPM->dbuffer,thePPPM->nxyz,
        MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    for(i=0;i<thePPPM->nxyz;i++) {
        thePPPM->density[i]=thePPPM->dbuffer[i];
    }
 #endif
#else
    printf("Warning, pppm methods being called without FFTW3, recompile with FFTW3 and run again\n");
#endif
}

point3d calculateForcePPPM(int index, double xin, double yin, double zin,
    PPPM* thePPPM, NbodyModel * theModel,
    double * x,
    double near) {
    point3d retval;
#ifdef HAS_FFTW3
    int i,j,k,ip,jp,kp,im,jm,km,l,ir,jr,kr;
    double cx, cy, cz;
    int nyz,nz;
    CubeInterp cint;
    double dx, dy, dz;
    double xmin,ymin,zmin,xmax,ymax,zmax;
    int * cell_contains;
    int n_cell;
    double r2,r3,r3i;

    while(xin<thePPPM->xmin) xin+=(thePPPM->xmax-thePPPM->xmin);
    while(xin>thePPPM->xmax) xin-=(thePPPM->xmax-thePPPM->xmin);
    while(yin<thePPPM->ymin) yin+=(thePPPM->ymax-thePPPM->ymin);
    while(yin>thePPPM->ymax) yin-=(thePPPM->ymax-thePPPM->ymin);
    while(zin<thePPPM->zmin) zin+=(thePPPM->zmax-thePPPM->zmin);
    while(zin>thePPPM->zmax) zin-=(thePPPM->zmax-thePPPM->zmin);

    nyz=thePPPM->nyz;
    nz=thePPPM->nz;
    xmin = thePPPM->xmin;
    ymin = thePPPM->ymin;
    zmin = thePPPM->zmin;
    xmax = thePPPM->xmax;
    ymax = thePPPM->ymax;
    zmax = thePPPM->zmax;
    retval.x=0.0;
    retval.y=0.0;
    retval.z=0.0;
    i = r2iPPPM(xin,xmin,xmax-xmin,thePPPM->nx);
    j = r2iPPPM(yin,ymin,ymax-ymin,thePPPM->ny);
    k = r2iPPPM(zin,zmin,zmax-zmin,thePPPM->nz);
    i = i%thePPPM->nx;
    j = j%thePPPM->ny;
    k = k%thePPPM->nz;
    ip = (i+1)%thePPPM->nx;
    jp = (j+1)%thePPPM->ny;
    kp = (k+1)%thePPPM->nz;
//need to make sure that this works for periodic boundaries
    xmin = thePPPM->x[i];
    ymin = thePPPM->y[j];
    zmin = thePPPM->z[k];
    xmax = xmin+(thePPPM->xmax-thePPPM->xmin)/(double)thePPPM->nx;
    ymax = ymin+(thePPPM->ymax-thePPPM->ymin)/(double)thePPPM->ny;
    zmax = zmin+(thePPPM->zmax-thePPPM->zmin)/(double)thePPPM->nz;

    // calculate force in x direction at corners of cell
    setCornersCINT(&cint,thePPPM->force[i*nyz+j*nz+k].x,
                        thePPPM->force[ip*nyz+j*nz+k].x,
                        thePPPM->force[i*nyz+jp*nz+k].x,
                        thePPPM->force[ip*nyz+jp*nz+k].x,
                        thePPPM->force[i*nyz+j*nz+kp].x,
                        thePPPM->force[ip*nyz+j*nz+kp].x,
                        thePPPM->force[i*nyz+jp*nz+kp].x,
                        thePPPM->force[ip*nyz+jp*nz+kp].x,
                        xmin,xmax,ymin,ymax,zmin,zmax);
    retval.x = getValueCINT(cint,xin,yin,zin);
    
    // calculate force in y direction at corners of cell
    setCornersCINT(&cint,thePPPM->force[i*nyz+j*nz+k].y,
                        thePPPM->force[ip*nyz+j*nz+k].y,
                        thePPPM->force[i*nyz+jp*nz+k].y,
                        thePPPM->force[ip*nyz+jp*nz+k].y,
                        thePPPM->force[i*nyz+j*nz+kp].y,
                        thePPPM->force[ip*nyz+j*nz+kp].y,
                        thePPPM->force[i*nyz+jp*nz+kp].y,
                        thePPPM->force[ip*nyz+jp*nz+kp].y,
                        xmin,xmax,ymin,ymax,zmin,zmax);
    retval.y = getValueCINT(cint,xin,yin,zin);

    
    // calculate force in z direction at corners of cell
    setCornersCINT(&cint,thePPPM->force[i*nyz+j*nz+k].z,
                        thePPPM->force[ip*nyz+j*nz+k].z,
                        thePPPM->force[i*nyz+jp*nz+k].z,
                        thePPPM->force[ip*nyz+jp*nz+k].z,
                        thePPPM->force[i*nyz+j*nz+kp].z,
                        thePPPM->force[ip*nyz+j*nz+kp].z,
                        thePPPM->force[i*nyz+jp*nz+kp].z,
                        thePPPM->force[ip*nyz+jp*nz+kp].z,
                        xmin,xmax,ymin,ymax,zmin,zmax);
    retval.z = getValueCINT(cint,xin,yin,zin);

    // for nearby cells, check for points in those
    // cells and use them to calculate effect of nearby particles.
    xmin = thePPPM->xmin;
    ymin = thePPPM->ymin;
    zmin = thePPPM->zmin;
    xmax = thePPPM->xmax;
    ymax = thePPPM->ymax;
    zmax = thePPPM->zmax;
    im= r2iPPPM(xin-near,xmin,xmax-xmin,thePPPM->nx);
    ip= r2iPPPM(xin+near,xmin,xmax-xmin,thePPPM->nx)+1;
    jm= r2iPPPM(yin-near,ymin,ymax-ymin,thePPPM->ny);
    jp= r2iPPPM(yin+near,ymin,ymax-ymin,thePPPM->ny)+1;
    km= r2iPPPM(zin-near,zmin,zmax-zmin,thePPPM->nz);
    kp= r2iPPPM(zin+near,zmin,zmax-zmin,thePPPM->nz)+1;
    for (i=im;i<ip;i++) {
        for (j=jm;j<jp;j++) {
            for (k=km;k<kp;k++) {
                ir = (i+thePPPM->nx)%thePPPM->nx;
                jr = (j+thePPPM->ny)%thePPPM->ny;
                kr = (k+thePPPM->nz)%thePPPM->nz;
                if (thePPPM->cell_contains[ir*nyz+jr*nz+kr]!=NULL) {
                    cell_contains=thePPPM->cell_contains[ir*nyz+jr*nz+kr];
                    n_cell = thePPPM->n_cell[ir*nyz+jr*nz+kr];
                    for(l=0;l<n_cell;l++) {
                        cx = x[cell_contains[l]*6];
                        cy = x[cell_contains[l]*6+1];
                        cz = x[cell_contains[l]*6+2];
                        if(i>=thePPPM->nx) cx += thePPPM->xmax-thePPPM->xmin;
                        if(i<0) cx -= thePPPM->xmax-thePPPM->xmin;
                        if(j>=thePPPM->ny) cy += thePPPM->ymax-thePPPM->ymin;
                        if(j<0) cy -= thePPPM->ymax-thePPPM->ymin;
                        if(k>=thePPPM->nz) cz += thePPPM->zmax-thePPPM->zmin;
                        if(k<0) cz -= thePPPM->zmax-thePPPM->zmin;
                        dx = cx-xin;
                        dy = cy-yin;
                        dz = cz-zin;
                        r2 = dx*dx+dy*dy+dz*dz;
                        if(r2>theModel->srad2[cell_contains[l]]&&
                            cell_contains[l]!=index&&r2<near*near) {
                            r3 = sqrt(r2)*r2;
                            r3i = theModel->G/r3;
                            retval.x+=theModel->mass[cell_contains[l]]*
                                      dx*r3i;
                            retval.y+=theModel->mass[cell_contains[l]]*
                                      dy*r3i;
                            retval.z+=theModel->mass[cell_contains[l]]*
                                      dz*r3i;
                        }
                    }
                }
            }
        }
    }

#else
    printf("Warning, pppm methods being called without FFTW3, recompile with FFTW3 and run again\n");
#endif
    return retval;
}
void setKPPPM(int n,double * k,double range) {
    int i;
    double nyquist = (double)n/range/2.0;
    for(i=0;i<n/2+1;i++) {
        k[i] = (double)i/(double)(n/2)*nyquist;
        if(i!=n/2&&i!=0) {
            k[n-i]=k[i];
        }
    }
}

void prepPotentialPPPM(PPPM* thePPPM, NbodyModel * theModel) {
    //calc fft of rho
    // inv fft -G/pi/k/k*fft(rho)
#ifdef HAS_FFTW3
    int i,j,k;
    double k2i;
    int nxyz,nyz,nz;
    double dx,dy,dz;
    int in,inpx,inpy,inpz,inmx,inmy,inmz;
    int nzo2,nyzo2,nxyzo2;

    nxyz = thePPPM->nxyz;
    nyz = thePPPM->nyz;
    nz = thePPPM->nz;
    nzo2 = thePPPM->nzo2;
    nyzo2 = thePPPM->nyzo2;
    nxyzo2 = thePPPM->nxyzo2;

/*
    for(i=0;i<thePPPM->nx;i++) {
        for(j=0;j<thePPPM->ny;j++) {
            for(k=0;k<thePPPM->nz;k++) {
                thePPPM->density[i*nyz+j*nz+k] = 10.0*exp(-1.0*(thePPPM->x[i]*thePPPM->x[i]+thePPPM->y[j]*thePPPM->y[j]+thePPPM->z[k]*thePPPM->z[k]));
            }
        }
    }
*/
    fftw_execute(thePPPM->pf);
    setKPPPM(thePPPM->nx,thePPPM->kx,thePPPM->xmax-thePPPM->xmin);
    setKPPPM(thePPPM->ny,thePPPM->ky,thePPPM->ymax-thePPPM->ymin);
    setKPPPM(thePPPM->nz,thePPPM->kz,thePPPM->zmax-thePPPM->zmin);
    for(i=0;i<thePPPM->nx;i++) {
        for(j=0;j<thePPPM->ny;j++) {
            for(k=0;k<thePPPM->nz/2+1;k++) {
                if(i==0&&j==0&&k==0) {
                    k2i=1.0;
                } else {
                k2i = -theModel->G/M_PI/
                    (thePPPM->kx[i]*thePPPM->kx[i] +
                    thePPPM->ky[j]*thePPPM->ky[j] +
                    thePPPM->kz[k]*thePPPM->kz[k]);
                }

                thePPPM->fft_potential[i*nyzo2+j*nzo2+k][0] =
                    thePPPM->fft_density[i*nyzo2+j*nzo2+k][0]*
                    k2i;
                thePPPM->fft_potential[i*nyzo2+j*nzo2+k][1] =
                    thePPPM->fft_density[i*nyzo2+j*nzo2+k][1]*
                    k2i;
            }
        }
    }
    fftw_execute(thePPPM->pb);
    for(i=0;i<thePPPM->nxyz;i++) {
        thePPPM->potential[i] /= (double)nxyz;
    }
/*
    for(i=1;i<thePPPM->nxyz;i++) {
        thePPPM->potential[i] -= thePPPM->potential[0];
    }
    thePPPM->potential[0]=0.0;
*/

    dx = (thePPPM->x[1]-thePPPM->x[0]);
    dy = (thePPPM->y[1]-thePPPM->y[0]);
    dz = (thePPPM->z[1]-thePPPM->z[0]);
    for(i=0;i<thePPPM->nx;i++) {
        for(j=0;j<thePPPM->ny;j++) {
            for(k=0;k<thePPPM->nz;k++) {
                in = i*nyz+j*nz+k;
                inpx = ((i+1)%thePPPM->nx)*nyz+j*nz+k;
                inmx = ((i-1+thePPPM->nx)%thePPPM->nx)*nyz+j*nz+k;
                inpy = i*nyz+((j+1)%thePPPM->ny)*nz+k;
                inmy = i*nyz+((j-1+thePPPM->ny)%thePPPM->ny)*nz+k;
                inpz = i*nyz+j*nz+(k+1)%thePPPM->nz;
                inmz = i*nyz+j*nz+(k-1+thePPPM->nz)%thePPPM->nz;
                thePPPM->force[in].x=
                    -(thePPPM->potential[inpx]-thePPPM->potential[inmx])/(2.0*dx);
                thePPPM->force[in].y=
                    -(thePPPM->potential[inpy]-thePPPM->potential[inmy])/(2.0*dy);
                thePPPM->force[in].z=
                    -(thePPPM->potential[inpz]-thePPPM->potential[inmz])/(2.0*dz);
/*
                thePPPM->den_check[in]= 1.0/(4.0*M_PI*theModel->G)*
                    ((thePPPM->potential[inpx]-2.0*thePPPM->potential[in]+thePPPM->potential[inmx])/(dx*dx)+
                    (thePPPM->potential[inpy]-2.0*thePPPM->potential[in]+thePPPM->potential[inmy])/(dy*dy)+
                    (thePPPM->potential[inpz]-2.0*thePPPM->potential[in]+thePPPM->potential[inmz])/(dz*dz));
*/
            }
        }
    }
#else
    printf("Warning, pppm methods being called without FFTW3, recompile with FFTW3 and run again\n");
#endif
}
