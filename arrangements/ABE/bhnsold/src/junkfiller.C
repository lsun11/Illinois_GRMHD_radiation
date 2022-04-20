#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <iostream>
#include <fstream>

using namespace std;

#define SQR(x) ((x)*(x))
#define INPUTARRAY_THETASIZE 10000
#define INPUTARRAY_PHISIZE 20000

#define ALLOCATE_2D(array,x,y)			\
  array = new double *[x];			\
  for(int i=0;i<x;i++) {			\
    *(array + i) = new double [y];		\
  }						\


int main() {

  double **inputdata_r1; ALLOCATE_2D(inputdata_r1,INPUTARRAY_THETASIZE,INPUTARRAY_PHISIZE);
  double **inputdata_r2; ALLOCATE_2D(inputdata_r2,INPUTARRAY_THETASIZE,INPUTARRAY_PHISIZE);

  double r0 = 0.2;
  double r1 = 0.4;
  double r2 = 0.5;

  for(int i=0;i<INPUTARRAY_THETASIZE;i++) for(int j=0;j<INPUTARRAY_PHISIZE;j++) {
    double theta_local = (double)i/INPUTARRAY_THETASIZE*M_PI;
    double phi_local = (double)j/INPUTARRAY_PHISIZE*2*M_PI;

    /*
    double x = r1*sin(theta_local)*cos(phi_local);
    double y = r1*sin(theta_local)*sin(phi_local);
    double z = r1*cos(theta_local);
    
    inputdata_r1[i][j] = 1.0;
    inputdata_r2[i][j] = 1.3;
    */

    inputdata_r1[i][j] = sin(r1)*sin(theta_local)*cos(phi_local);
    inputdata_r2[i][j] = sin(r2)*sin(theta_local)*cos(phi_local);
  }

  for(int i=0;i<100;i++) for(int j=0;j<100;j++) {

    double x = (double)i/100.0*0.4;
    double y = (double)j/100.0*0.4;
    double z = 0.0;

    double theta = atan2(sqrt(x*x+y*y),z);
    double phi   = atan2(y,x);
    

    /*
    double theta = (double)i/100.0*M_PI;
    double phi = (double)j/100.0*2*M_PI;

    double x = r1*sin(theta)*cos(phi);
    double y = r1*sin(theta)*sin(phi);
    double z = r1*cos(theta);
    */

    int index_i=(int)(theta*INPUTARRAY_THETASIZE/(M_PI));
    if(index_i==INPUTARRAY_THETASIZE-1) index_i--;
    int index_ip1=index_i+1;

    int index_j=(int)(phi*INPUTARRAY_PHISIZE/(2*M_PI));
    if(index_j==INPUTARRAY_PHISIZE-1) index_j--;
    int index_jp1=index_j+1;

    // Here we set coefficients for interpolation
    //
    // Interpolation coefficients in each dimension are set by solving the following 
    // equation for (a,b,c)[= e.g., (f11,f12,f13)]:
    // a f(x1) + b f(x2) = f(x5), where (x1,x2,x3) are the coordinate points in
    //    the interpolation stencil, and x5 is the interpolated destination point.
    // If f(x) = k = constant, we have
    // a + b = 1
    // If f(x) = k x, we have
    // a x1 + b x2 = x5

    // Using Mathematica, we can solve these equations easily:
    //   CForm[FullSimplify[Solve[{a + b == 1, a*x1 + b*x2 == x5}, {a, b}]]]
    // This has solution: 
    // a = (-x2 + x5)/(x1 - x2);
    // b = (x1 - x5)/(x1 - x2);

    double x1,x2,x5;
    x1 = (double)(index_i)/INPUTARRAY_THETASIZE*M_PI;
    x2 = (double)(index_ip1)/INPUTARRAY_THETASIZE*M_PI;
    x5 = theta;
    //fij indicates direction i, point number j in the interpolation stencil, as in j=1->a, j=2->b
    double f11 = (-x2 + x5)/(x1 - x2);
    double f12 = (x1 - x5)/(x1 - x2);

    if( (fabs(x2-x1) >= fabs(x5-x1)) &&
	(fabs(x2-x1) >= fabs(x5-x2)) ) {
      // do nothing.  this is correct
    } else {
      printf("BAD %e\t%e\t%e\t%d\t%d\n",x1,x5,x2,index_i,index_ip1);
      exit(0);
    }

    x1 = (double)(index_j)/INPUTARRAY_PHISIZE*2.0*M_PI;
    x2 = (double)(index_jp1)/INPUTARRAY_PHISIZE*2.0*M_PI;
    x5 = phi;
    //fij indicates direction i, point number j in the interpolation stencil, as in j=1->a, j=2->b
    double f21 = (-x2 + x5)/(x1 - x2);
    double f22 = (x1 - x5)/(x1 - x2);

    if( (fabs(x2-x1) >= fabs(x5-x1)) &&
	(fabs(x2-x1) >= fabs(x5-x2)) ) {
      // do nothing.  this is correct
    } else {
      printf("BAD %e\t%e\t%e\n",x1,x5,x2);
      exit(0);
    }
  
    double f1 = 
      f11*f21*inputdata_r1[index_i][index_j] +
      f12*f21*inputdata_r1[index_ip1][index_j] +
      f11*f22*inputdata_r1[index_i][index_jp1] +
      f12*f22*inputdata_r1[index_ip1][index_jp1];
    
    double f2 = 
      f11*f21*inputdata_r2[index_i][index_j] +
      f12*f21*inputdata_r2[index_ip1][index_j] +
      f11*f22*inputdata_r2[index_i][index_jp1] +
      f12*f22*inputdata_r2[index_ip1][index_jp1];

    //FullSimplify[Solve[{f1==a*r1^3+b*r1^2+c*r1+d,f2==a*r2^3+b*r2^2+c*r2+d,d==-(a*r0^3+b*r0^2+c*r0),c==-(3*a*r0^2+2*b*r0)},{a,b,c,d}]]
    double a = (-(f2*SQR(r0 - r1)) + f1*SQR(r0 - r2)) / (SQR(r0 - r1)*SQR(r0 - r2)*(r1 - r2));
    double b = (f2*SQR(r0 - r1)*(2*r0 + r1) - f1*SQR(r0 - r2)*(2*r0 + r2)) / (SQR(r0 - r1)*SQR(r0 - r2)*(r1 - r2));
    double c = (r0*(-(f2*SQR(r0 - r1)*(r0 + 2*r1)) + f1*SQR(r0 - r2)*(r0 + 2*r2))) / (SQR(r0 - r1)*SQR(r0 - r2)*(r1 - r2));
    double d = (SQR(r0)*(f2*SQR(r0 - r1)*r1 - f1*SQR(r0 - r2)*r2)) / (SQR(r0 - r1)*SQR(r0 - r2)*(r1 - r2));

    double r = sqrt(x*x+y*y+z*z);
    double outputvalue = a*r*r*r+b*r*r+c*r+d;
  
    //printf("%e\t%e\t%e\t%e\n",a,b,c,d);

    if(r<r0) outputvalue = 0.0;
    //printf("%e\t%e\t%e\t%e\n",b,SQR(r0 - r1)*(2*r0 + r1),SQR(r0 - r2)*(2*r0 + r2),(SQR(r0 - r1)*SQR(r0 - r2)*(r1 - r2)));
    if(r<r2) printf("%e\t%e\t%e\t%e\t%e\n",x,y,z,outputvalue,sin(r)*sin(theta)*cos(phi));
  }
  return 0;
}
