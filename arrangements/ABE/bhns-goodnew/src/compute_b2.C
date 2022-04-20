#define F1o4pi 0.07957747154594766788

#include <stdio.h>
#include "cctk.h"
#include <math.h>
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

// Compute b^2 and store it in Pr 

extern "C" void CCTK_FCALL CCTK_FNAME(bhns_compute_b2_cpp)
  (const cGH **cctkGH,int *cctk_lsh, double *phi, double *lapm1, 
   double *shiftx, double *shifty, double *shiftz, double *vx, 
   double *vy, double *vz, double *Bx, double *By, double *Bz, 
   double *gxx, double *gxy, double *gxz, 
   double *gyy, double *gyz, double *gzz, double *Pr);

void bhns_compute_b2_cpp(const cGH *cctkGH,int *cctk_lsh, double *phi, 
   double *lapm1, double *shiftx, double *shifty, double *shiftz, double *vx, 
   double *vy, double *vz, double *Bx, double *By, double *Bz,
   double *gxx, double *gxy, double *gxz, double *gyy, double *gyz, 
   double *gzz, double *Pr) { 
   
#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++)
    for(int j=0;j<cctk_lsh[1];j++)
      for(int i=0;i<cctk_lsh[0];i++) {

        int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	double psi4 = exp(4.0*phi[index]);
	double gxxi = psi4*gxx[index];
        double gxyi = psi4*gxy[index];
        double gxzi = psi4*gxz[index];
        double gyyi = psi4*gyy[index];
        double gyzi = psi4*gyz[index];
        double gzzi = psi4*gzz[index];
        double alpha = 1.0+lapm1[index];
	double vxi = vx[index];
        double vyi = vy[index];
        double vzi = vz[index];
        double Bxi = Bx[index];
        double Byi = By[index];
        double Bzi = Bz[index];
        double shiftxi = shiftx[index];
        double shiftyi = shifty[index];
        double shiftzi = shiftz[index];
	double v2 = (gxxi*(vxi+shiftxi)*(vxi+shiftxi) 
                         + 2.0*gxyi*(vxi+shiftxi)*(vyi+shiftyi) 
                         + 2.0*gxzi*(vxi+shiftxi)*(vzi +shiftzi)
			 + gyyi*(vyi+shiftyi)*(vyi+shiftyi) 
                         + 2.0*gyzi*(vyi+shiftyi)*(vzi +shiftzi) 
		         + gzzi*(vzi +shiftzi)*(vzi +shiftzi) ) / (alpha*alpha);
	// Limit velocity 
	if (v2 > 0.99) {
	   double fac = sqrt(0.99/v2);
	   v2 = 0.99; 
	   vxi = fac*(vxi+shiftxi) - shiftxi;
           vyi = fac*(vyi+shiftyi) - shiftyi;
           vzi = fac*(vzi+shiftzi) - shiftzi;
	}

	double alp_u02 = 1.0/(1.0-v2);
	double u02 = alp_u02 / (alpha*alpha);
        double u0 = sqrt(u02);
        double u_x = u0*(gxxi*(vxi+shiftxi) + gxyi*(vyi+shiftyi) + gxzi*(vzi+shiftzi) );
        double u_y = u0*(gxyi*(vxi+shiftxi) + gyyi*(vyi+shiftyi) + gyzi*(vzi+shiftzi) );
        double u_z = u0*(gxzi*(vxi+shiftxi) + gyzi*(vyi+shiftyi) + gzzi*(vzi+shiftzi) );
	double B2 = gxxi*Bxi*Bxi + 2.0*gxyi*Bxi*Byi + 2.0*gxzi*Bxi*Bzi 
		    + gyyi*Byi*Byi + 2.0*gyzi*Byi*Bzi + gzzi*Bzi*Bzi;
	double udotB = u_x*Bxi + u_y*Byi + u_z*Bzi;
	double b2 = (B2 + udotB*udotB)*F1o4pi/alp_u02;

	if(B2<0) {
	  printf("BAD B2!\n");
	  printf("gij's: %e %e %e %e %e %e\n",gxxi,gxyi,gxzi,gyyi,gyzi,gzzi);
	  printf("Bi's: %e %e %e\n",Bxi,Byi,Bzi);
	  printf("vi's: %e %e %e\n",vxi,vyi,vzi);
	  printf("psi4: %e\n",psi4);
	  printf("alp_u02: %e\n",alp_u02);
	  printf("v2: %e\n",v2);
	  printf("alpha: %e\n",alpha);
	}

        Pr[index] = b2;

  }
}

extern "C" void CCTK_FCALL CCTK_FNAME(bhns_compute_b2_cpp)
  (const cGH **cctkGH,int *cctk_lsh, double *phi, double *lapm1, 
   double *shiftx, double *shifty, double *shiftz, double *vx,
   double *vy, double *vz, double *Bx, double *By, double *Bz,
   double *gxx, double *gxy, double *gxz, double *gyy, double *gyz, double *gzz,   double *Pr) 
{
  bhns_compute_b2_cpp(*cctkGH, cctk_lsh, phi, lapm1, shiftx,shifty,shiftz,
                           vx,vy,vz, Bx,By,Bz, gxx, gxy, gxz, gyy, gyz, gzz, Pr);
}
