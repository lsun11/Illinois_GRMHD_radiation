#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

using namespace std;

#define SQR(x) ((x)*(x))
#define INDEX3D(SIZEI,SIZEJ, i, j, k)  ((i) + (SIZEI) *		\
					((j) + (SIZEJ) * (k)))

/*
Before this subroutine, we interpolated A^i data onto two spheres, r=r1 and r=r2,
where r=0 is the centroid of the AH.  This subroutine overwrites A^i data inside
r=r1 with data, according to the following algorithm:

To fill in the BH with data at gridpoint (x,y,z), this excision algorithm 
a) first computes the angle theta,phi corresponding to (x,y,z) ,
b) then it performs a 2D linear angular interpolation to theta,phi on spheres
r=r1 and r=r2, where r=0 is the centroid of the AH.
c) Finally, third-order polynomial interpolation is used in the radial direction, 
  fixing the slope and value of the polynomial f(r) at r0 so that f(r0)=f'(r0)=0.  
Thus we only have two remaining unknowns in the interpolating polynomial.
These are set by the values f(r1) and f(r2) on spheres of radius r1 and r2
respectively.  

Thus we have 3 free parameters: r0, r1, and r2.  We can 
*/

void polint_bhns(double xa[], double ya[], int n, double x, double *y, double *dy, double *c, double *d);

extern "C" void CCTK_FCALL CCTK_FNAME(bhns_fill_Ai_BH_interior)
  (const cGH **cctkGH,int *cctk_lsh,int &Symmetry,
   double *Ai, double *xinput,double *yinput,double *zinput,
   double &xbh,double &ybh,double &zbh,
   int &NUM_ZERO_PTS,int &RADIAL_INTERP_ORDER,
   int &INPUTARRAY_THETASIZE,int &INPUTARRAY_PHISIZE,
   double *r_array,
   double *Ai_r);

void bhns_fill_Ai_BH_interior(const cGH *cctkGH,int *cctk_lsh, int &Symmetry,
				  double *Ai, double *xinput,double *yinput,double *zinput,
				  double &xbh,double &ybh,double &zbh,
				  int &NUM_ZERO_PTS,int &RADIAL_INTERP_ORDER,
				  int &INPUTARRAY_THETASIZE,int &INPUTARRAY_PHISIZE,
				  double *r_array,
				  double *Ai_r) {

  double *c = (double *)malloc(sizeof(double)*(RADIAL_INTERP_ORDER+1));
  double *d = (double *)malloc(sizeof(double)*(RADIAL_INTERP_ORDER+1));
  double *f = (double *)malloc(sizeof(double)*(RADIAL_INTERP_ORDER+1));
  double *r_arraynr = (double *)malloc(sizeof(double)*(RADIAL_INTERP_ORDER+1));

  for(int i=1;i<=RADIAL_INTERP_ORDER;i++) {
    r_arraynr[i] = r_array[i-1];
  }

  for(int k=0;k<cctk_lsh[2];k++)
    for(int j=0;j<cctk_lsh[1];j++)
      for(int i=0;i<cctk_lsh[0];i++) {

        int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

	double x = xinput[index]-xbh;
	double y = yinput[index]-ybh;
	double z = zinput[index]-zbh;
	double r = sqrt(x*x+y*y+z*z);

	if(r<r_array[NUM_ZERO_PTS] && z>-1e-8) {

	  double Symmetry_factor;
	  if(Symmetry==0) {
	    Symmetry_factor=1.0;
	  } else if(Symmetry==1) {
	    Symmetry_factor=0.5;
	  } else {
	    printf("Error.  Symmetry=%d not supported in bhns_fill_Ai_BH_interior\n",Symmetry);
	    exit(0);
	  }

	  double theta = acos(z/(r+1e-30));
	  double phi   = atan2(y,x);
	  double xcheck=r*sin(theta)*cos(phi)+xbh-xinput[index];
	  double ycheck=r*sin(theta)*sin(phi)+ybh-yinput[index];
	  if(fabs(xcheck)>1e-8 || fabs(ycheck)>1e-8) { printf("BAD BLAH\n"); exit(0); }

	  // We need the following line, since theta can be >pi/2 when equatorial symmetry is assumed, triggering a "BAD" error below.
	  if(r<1e-8) { theta=0; phi=0; }
	  
	  int index_i=(int)(theta*(INPUTARRAY_THETASIZE-1)/(M_PI*Symmetry_factor));
	  if(index_i==INPUTARRAY_THETASIZE-1) index_i--;
	  int index_ip1=index_i+1;

	  int index_j=(int)((phi+M_PI)*(INPUTARRAY_PHISIZE-1)/(2.0*M_PI));
	  if(index_j==INPUTARRAY_PHISIZE-1) index_j--;
	  int index_jp1=index_j+1;
	  
	  // Here we set coefficients for first-order angular interpolation
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
	  x1 = (double)(index_i)/(INPUTARRAY_THETASIZE-1)*M_PI*Symmetry_factor;
	  x2 = (double)(index_ip1)/(INPUTARRAY_THETASIZE-1)*M_PI*Symmetry_factor;
	  x5 = theta;
	  //fij indicates direction i, point number j in the interpolation stencil, as in j=1->a, j=2->b
	  double f11 = (-x2 + x5)/(x1 - x2);
	  double f12 = (x1 - x5)/(x1 - x2);
	  
	  if( (fabs(x2-x1) >= fabs(x5-x1)) &&
	      (fabs(x2-x1) >= fabs(x5-x2)) ) {
	    // do nothing.  this is correct
	  } else {
	    if(fabs(fabs(x2-x1)-fabs(x5-x1)) < 1e-8 ||
	       fabs(fabs(x2-x1)-fabs(x5-x2)) < 1e-8) {
	      // again do nothing, this is just roundoff error
	    } else {
	      printf("BAD %e\t%e\t%e\t%d\t%d\t%e\t%e\t%e\t%e\n",x1,x5,x2,index_i,index_ip1,x,y,z,r);
	      exit(0);
	    }
	  }

	  x1 = (double)(index_j)/(INPUTARRAY_PHISIZE-1)*2.0*M_PI-M_PI;
	  x2 = (double)(index_jp1)/(INPUTARRAY_PHISIZE-1)*2.0*M_PI-M_PI;
	  x5 = phi;
	  //fij indicates direction i, point number j in the interpolation stencil, as in j=1->a, j=2->b
	  double f21 = (-x2 + x5)/(x1 - x2);
	  double f22 = (x1 - x5)/(x1 - x2);
	  
	  if( (fabs(x2-x1) >= fabs(x5-x1)) &&
	      (fabs(x2-x1) >= fabs(x5-x2)) ) {
	    // do nothing.  this is correct
	  } else {
	    printf("BAD2 coords:%e,%e,%e\t %e\t%e\t%e\t%d\t%d\t%d\t%e\n",x,y,z,x1,x5,x2,index_j,INPUTARRAY_PHISIZE,(int)(phi*INPUTARRAY_PHISIZE/(2.0*M_PI)),M_PI);
	    printf("BAD2 %e\t%e\t%e\t%d\t%d\t%e\t%e\t%e\t%e\n",x1,x5,x2,index_j,index_jp1,x,y,z,r);
	    exit(0);
	  }

          if(NUM_ZERO_PTS>0) {
            for(int whichradius=1;whichradius<=NUM_ZERO_PTS;whichradius++) f[whichradius]=0;
          }
	  for(int whichradius=NUM_ZERO_PTS;whichradius<RADIAL_INTERP_ORDER;whichradius++) {
	    f[whichradius+1] =
	      f11*f21*Ai_r[INDEX3D(RADIAL_INTERP_ORDER-NUM_ZERO_PTS,INPUTARRAY_THETASIZE,whichradius-NUM_ZERO_PTS,index_i,index_j)] +
	      f12*f21*Ai_r[INDEX3D(RADIAL_INTERP_ORDER-NUM_ZERO_PTS,INPUTARRAY_THETASIZE,whichradius-NUM_ZERO_PTS,index_ip1,index_j)] +
	      f11*f22*Ai_r[INDEX3D(RADIAL_INTERP_ORDER-NUM_ZERO_PTS,INPUTARRAY_THETASIZE,whichradius-NUM_ZERO_PTS,index_i,index_jp1)] +
	      f12*f22*Ai_r[INDEX3D(RADIAL_INTERP_ORDER-NUM_ZERO_PTS,INPUTARRAY_THETASIZE,whichradius-NUM_ZERO_PTS,index_ip1,index_jp1)];
	  }

	  double dy,outputvalue;
	  polint_bhns(r_arraynr,f,RADIAL_INTERP_ORDER,r,&outputvalue,&dy,c,d);

	  if(NUM_ZERO_PTS>0) {
	    //Interpolation BCs
	    if(r<r_array[NUM_ZERO_PTS-1]) outputvalue = 0.0;
	  } else {
	    //Extrapolation BCs
	    if(r<r_array[99]) outputvalue = 0.0;
	  }

	  

	  Ai[index] = outputvalue;
	}
      }
  free(c);
  free(d);
  free(f);
  free(r_arraynr);
}

extern "C" void CCTK_FCALL CCTK_FNAME(bhns_fill_Ai_BH_interior)
  (const cGH **cctkGH,int *cctk_lsh,int &Symmetry,
   double *Ai, double *xinput,double *yinput,double *zinput,
   double &xbh,double &ybh,double &zbh,
   int &NUM_ZERO_PTS,int &RADIAL_INTERP_ORDER,
   int &INPUTARRAY_THETASIZE,int &INPUTARRAY_PHISIZE,
   double *r_array,
   double *Ai_r) {
  bhns_fill_Ai_BH_interior(*cctkGH,cctk_lsh,Symmetry,
			       Ai, xinput,yinput,zinput,
			       xbh,ybh,zbh,
			       NUM_ZERO_PTS,RADIAL_INTERP_ORDER,
			       INPUTARRAY_THETASIZE,INPUTARRAY_PHISIZE,
			       r_array,
			       Ai_r);
}

// c and d are temporary arrays here:
void polint_bhns(double xa[], double ya[], int n, double x, double *y, double *dy, double *c, double *d)
{
  int i,m,ns=1;
  double den,dif,dift,ho,hp,w;

  dif=fabs(x-xa[1]);
  for (i=1;i<=n;i++) {
    if ( (dift=fabs(x-xa[i])) < dif) {
      ns=i;
      dif=dift;
    }
    c[i]=ya[i];
    d[i]=ya[i];
  }
  *y=ya[ns--];
  for (m=1;m<n;m++) {
    for (i=1;i<=n-m;i++) {
      ho=xa[i]-x;
      hp=xa[i+m]-x;
      w=c[i+1]-d[i];
      if ( (den=ho-hp) == 0.0) {
	printf("Error in routine polint_bhns");
	exit(1);
      }
      den=w/den;
      d[i]=hp*den;
      c[i]=ho*den;
    }
    *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
  }
}
