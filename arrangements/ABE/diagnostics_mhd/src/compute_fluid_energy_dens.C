#define F1o4pi 0.07957747154594766788
#define v2max  0.99

#include <stdio.h>
#include "cctk.h"
#include <math.h>
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

extern "C" void CCTK_FCALL CCTK_FNAME(compute_fluid_energy_dens_cpp)
  (const cGH **cctkGH,int *cctk_lsh, double *phi, double *lapm1, 
   double *shiftx, double *shifty, double *shiftz, double *vx, 
   double *vy, double *vz, double *rho_star, double *h, double *P, double *gxx, double *gxy, double *gxz, 
   double *gyy, double *gyz, double *gzz, double *Tfluid0_0_sqrtg);

void compute_fluid_energy_dens_cpp(const cGH *cctkGH,int *cctk_lsh, double *phi, 
				   double *lapm1, double *shiftx, double *shifty, double *shiftz, double *vx, 
				   double *vy, double *vz, double *rho_star, double *h, double *P,
				   double *gxx, double *gxy, double *gxz, double *gyy, double *gyz, 
				   double *gzz, double *Tfluid0_0_sqrtg) { 
  
#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++)
    for(int j=0;j<cctk_lsh[1];j++)
      for(int i=0;i<cctk_lsh[0];i++) {

        int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
        double psi2 = exp(2.0*phi[index]);
	double psi4 = psi2*psi2;
        double psi6 = psi2*psi4;
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
	if (v2 > v2max) {
	   double fac = sqrt(v2max/v2);
	   v2 = v2max; 
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
	double shift_xi = gxxi*shiftxi + gxyi*shiftyi + gxzi*shiftzi;
        double shift_yi = gxyi*shiftxi + gyyi*shiftyi + gyzi*shiftzi;
        double shift_zi = gxzi*shiftxi + gyzi*shiftyi + gzzi*shiftzi;
        double beta2 = shiftxi*shift_xi + shiftyi*shift_yi + shiftzi*shift_zi;
        double udotbeta = u_x*shiftxi + u_y*shiftyi + u_z*shiftzi;
	double u_0 = (beta2-alpha*alpha)*u0 + udotbeta;
	Tfluid0_0_sqrtg[index] = -rho_star[index]*(1.0+h[index]*u_0) - P[index]*alpha*psi6; // Tfluid_0^0 sqrt(-g)

  }
}

extern "C" void CCTK_FCALL CCTK_FNAME(compute_fluid_energy_dens_cpp)
  (const cGH **cctkGH,int *cctk_lsh, double *phi, double *lapm1, 
   double *shiftx, double *shifty, double *shiftz, double *vx,
   double *vy, double *vz, double *rho_star, double *h, double *P,
   double *gxx, double *gxy, double *gxz, double *gyy, double *gyz, double *gzz,   
   double *Tfluid0_0_sqrtg) 
{
  compute_fluid_energy_dens_cpp(*cctkGH, cctk_lsh, phi, lapm1, shiftx,shifty,shiftz,
				vx,vy,vz, rho_star,h,P, gxx, gxy, gxz, gyy, gyz, gzz, Tfluid0_0_sqrtg);
}
