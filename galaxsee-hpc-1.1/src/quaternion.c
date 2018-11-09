#include "quaternion.h"
#include <math.h>

void QUAT_multiply(Quaternion quat1,
        Quaternion quat2,Quaternion * product) {
    double dot;
    Quaternion temp;

    dot=quat1.b*quat2.b+quat1.c*quat2.c+quat1.d*quat2.d;
    temp.a=quat1.a*quat2.a-dot;

    temp.b=quat1.a*quat2.b;
    temp.c=quat1.a*quat2.c;
    temp.d=quat1.a*quat2.d;

    temp.b+=quat2.a*quat1.b;
    temp.c+=quat2.a*quat1.c;
    temp.d+=quat2.a*quat1.d;

    temp.b+=quat1.c*quat2.d-quat1.d*quat2.c;
    temp.c+=quat1.d*quat2.b-quat1.b*quat2.d;
    temp.d+=quat1.b*quat2.c-quat1.c*quat2.b;

    product->a=temp.a;
    product->b=temp.b;
    product->c=temp.c;
    product->d=temp.d;
    return;
}

void QUAT_conjugate(Quaternion * quat1) {
    quat1->b=-quat1->b;
    quat1->c=-quat1->c;
    quat1->d=-quat1->d;
}

void QUAT_rotate(Quaternion quat,double x,double y,double z,
    double * rx, double * ry, double * rz) {
    Quaternion vector;
    Quaternion result;
    Quaternion current;

    current.a=quat.a;
    current.b=quat.b;
    current.c=quat.c;
    current.d=quat.d;

    vector.a=0.0;
    vector.b=x;
    vector.c=y;
    vector.d=z;

    QUAT_multiply(current,vector,&result);
    QUAT_conjugate(&current);
    QUAT_multiply(result,current,&result);

    *rx=result.b;
    *ry=result.c;
    *rz=result.d;

    return;
}


void QUAT_mouseRotateSelf(Quaternion * quat,int xcorr, int ycorr,
        int lastx,int lasty,int newx,
        int newy, double rot_con) {
    double tempx,tempy;
    Quaternion rotator;
    rotator.a=1.0; rotator.b=0.0; rotator.c=0.0; rotator.d=0.0;
    tempx=(double)(newx-lastx);
    tempy=(double)(-newy+lasty);
    // check for distance moved
    double rad=sqrt(tempx*tempx+tempy*tempy);
    //Make sure it actually moved (pixel values are assumed for
    // entry
    if(rad<1.0) return;
    // normalize unit vector
    tempx/=rad;
    tempy/=rad;
    // set axis of rotation, increase rotation
    // angle for larger mouse move.
    rotator.a=cos(rad*rot_con);
    if (xcorr==-3) {
        rotator.d=-tempx*sin(rad*rot_con);
    } else if (xcorr==-2) {
        rotator.c=-tempx*sin(rad*rot_con);
    } else if (xcorr==-1) {
        rotator.b=-tempx*sin(rad*rot_con);
    } else if (xcorr==1) {
        rotator.b=tempx*sin(rad*rot_con);
    } else if (xcorr==2) {
        rotator.c=tempx*sin(rad*rot_con);
    } else if (xcorr==3) {
        rotator.d=tempx*sin(rad*rot_con);
    }
    if (ycorr==-3) {
        rotator.d=-tempy*sin(rad*rot_con);
    } else if (ycorr==-2) {
        rotator.c=-tempy*sin(rad*rot_con);
    } else if (ycorr==-1) {
        rotator.b=-tempy*sin(rad*rot_con);
    } else if (ycorr==1) {
        rotator.b=tempy*sin(rad*rot_con);
    } else if (ycorr==2) {
        rotator.c=tempy*sin(rad*rot_con);
    } else if (ycorr==3) {
        rotator.d=tempy*sin(rad*rot_con);
    }
    // update current rotation angle
    QUAT_multiply(rotator,*quat,quat);
}

void QUAT_normalize(Quaternion * quat) {
    double sum=0.0;
    sum = sqrt(quat->a*quat->a+quat->b*quat->b+quat->c*quat->c+quat->d*quat->d);
    quat->a/=sum;
    quat->b/=sum;
    quat->c/=sum;
    quat->d/=sum;
    return;
}

void QUAT_getMatrix(Quaternion quat,fmatrix * retval) {

        QUAT_normalize(&quat);
	float x2 = quat.a * quat.a;
	float y2 = quat.b * quat.b;
	float z2 = quat.c * quat.c;
	float xy = quat.a * quat.b;
	float xz = quat.a * quat.c;
	float yz = quat.b * quat.c;
	float wx = quat.d * quat.a;
	float wy = quat.d * quat.b;
	float wz = quat.d * quat.c;
 
	(*retval)[0] = 1.0f - 2.0f * (y2 + z2);
        (*retval)[1] = 2.0f * (xy - wz);
        (*retval)[2] = 2.0f * (xz + wy);
        (*retval)[3] = 0.0f;
        (*retval)[4] = 2.0f * (xy + wz);
        (*retval)[5] = 1.0f - 2.0f * (x2 + z2);
        (*retval)[6] = 2.0f * (yz - wx);
        (*retval)[7] = 0.0f;
        (*retval)[8] = 2.0f * (xz - wy);
        (*retval)[9] = 2.0f * (yz + wx);
        (*retval)[10] = 1.0f - 2.0f * (x2 + y2);
        (*retval)[11] = 0.0f;
        (*retval)[12] = 0.0f;
        (*retval)[13] = 0.0f;
        (*retval)[14] = 0.0f;
        (*retval)[15] = 1.0f;
        return;
}
