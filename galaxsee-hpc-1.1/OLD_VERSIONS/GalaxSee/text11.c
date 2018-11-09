#include <stdio.h>
#include <stdlib.h>
#include "text11.h"
#include "mem.h"
/**********************************
 * text11 
 *
 * allows for 2-D binning of x-y
 * point data displayed with a
 * logarithmic iconograph scale
 *********************************/

void text11_set_range(Text11 *t11p,double dx_min,double dy_min,
        double dx_max,double dy_max) {
    t11p->dx_min=dx_min;
    t11p->dy_min=dy_min;
    t11p->dx_max=dx_max;
    t11p->dy_max=dy_max;
}

void text11_reset(Text11 * t11p) {
    int i,j;
    for (i=0;i<t11p->nx;i++) {
        for (j=0;j<t11p->ny;j++) {
            (t11p->pixels)[i][j]=0;
        }
    }
}
void text11_free(Text11 **t11Display) {
    Text11 * t11p = (*t11Display);
    free_imatrix(t11p->pixels);
    free(t11p->boundary);
    free(t11p);
}

void text11_set_boundary(Text11 *t11p,int max_count) {
    double chi;
    t11p->boundary[0]=0;
    t11p->boundary[1]=1;
    t11p->boundary[2]=2;
    chi = pow((double)max_count/2.0,0.25);
    t11p->boundary[3]=(int)(2.0*chi);
    t11p->boundary[4]=(int)(2.0*chi*chi);
    t11p->boundary[5]=(int)(2.0*chi*chi*chi);
    t11p->boundary[6]=(int)(2.0*chi*chi*chi*chi);
}

void text11_initialize(Text11 ** t11Display) {
    Text11 * t11p = (*t11Display);
    (*t11Display)=(Text11 *)malloc(sizeof(Text11));
    t11p=*t11Display;
    
    t11p->nx=48;
    t11p->ny=20;
    t11p->nb=7;
    t11p->display = (char *)malloc(sizeof(char)*(t11p->nb+3));
    t11p->boundary = (int *)malloc(sizeof(int)*(t11p->nb));
    t11p->pixels = alloc_imatrix(t11p->nx,
                                t11p->ny);
    t11p->display = " ,;o*OQ";
    t11p->boundary[0]=0;
    t11p->boundary[1]=1;
    t11p->boundary[2]=2;
    t11p->boundary[3]=10;
    t11p->boundary[4]=50;
    t11p->boundary[5]=100;
    t11p->boundary[5]=500;
    t11p->dx_min=-1.0;
    t11p->dx_max=1.0;
    t11p->dy_min=-1.0;
    t11p->dy_max=1.0;
    text11_reset(t11p);
}
void text11_print(Text11 * t11p) {
    int i,j,k;
    for(j=t11p->ny-1;j>=0;j--) {
        for(i=0;i<t11p->nx;i++) {
            for(k=t11p->nb-1;k>=0;k--) {
                if(t11p->pixels[i][j]>=t11p->boundary[k]) {
                    printf("%c",t11p->display[k]);
                    k=-1;
                }
            }
        }
        printf("\n");
    }
}
int text11_add(Text11 * t11p,double x,double y) {
    double bx_min,bx_max,by_min,by_max;
    int done=0;
    int i,j;
    double xstep,ystep;
    xstep = (t11p->dx_max-t11p->dx_min)/(double)(t11p->nx-1);
    ystep = (t11p->dy_max-t11p->dy_min)/(double)(t11p->ny-1);
    bx_min = t11p->dx_min;
    for(i=0;i<t11p->nx&&!done;i++) {
        bx_max = t11p->dx_min + (double)(i+1)*xstep;
        by_min = t11p->dy_min;
        for(j=0;j<t11p->ny&&!done;j++) {
            by_max = t11p->dy_min + (double)(j+1)*ystep;
            if(x>=bx_min&&x<bx_max&&y>=by_min&&y<by_max) {
                t11p->pixels[i][j]+=1;
                done=1;
            }
            by_min = by_max;
        }
        bx_min=bx_max;
    }
    return done;
}
