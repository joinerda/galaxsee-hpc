#ifndef QUATERNION
#define QUATERNION

/************************************************************
* Quaternion library
*
*   An important use of quaternions in computer graphics
*   is for rotation of points in 3 dimensional space.
*   You can define a rotation axis and the angle of rotation
*   about that axis by
*     a = cos (theta/2)
*     b = alpha * sin(theta/2)
*     c = beta * sin(theta/2)
*     d = gamma * sin(theta/2)
*   where alpha, beta, and gamma are the direction cosines
*   that define the axis of rotation alpha * i + beta * j + gamma * k
*
*************************************************************
*/

typedef struct {
    double a,b,c,d;
} Quaternion;

typedef float fmatrix[16];

void QUAT_multiply(Quaternion quat1,
            Quaternion quat2,Quaternion * product);
void QUAT_conjugate(Quaternion * quat1);
void QUAT_rotate(Quaternion quat,double x,double y,double z,
        double * rx, double * ry, double * rz);
void QUAT_mouseRotateSelf(Quaternion * quat,int xcorr, int ycorr,
            int lastx,int lasty,int newx,
            int newy, double rot_con);
void QUAT_getMatrix(Quaternion quat,fmatrix *);

#endif
