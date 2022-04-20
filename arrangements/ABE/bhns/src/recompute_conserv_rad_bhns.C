#define F1o4pi 0.07957747154594766788
#define v2max  0.99

#include <stdio.h>
#include <stdlib.h>
#include "cctk.h"
#include <math.h>
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

#define SQR(x) ((x) * (x)) 

extern "C" void CCTK_FCALL CCTK_FNAME(recompute_conserv_rad_bhns)
  (const cGH **cctkGH, int *cctk_lsh,
   double *tau_rad, double *S_rad_x, double *S_rad_y, double *S_rad_z,
   double *E_rad, double *F_radx, double *F_rady, double *F_radz, double *F_rad0, double *F_rad_scalar, 
   double *P_radxx, double *P_radyy, double *P_radzz, double *P_radxy, double *P_radxz, double *P_radyz,
   double *phi, double *alpha, double *shiftx, double *shifty, double *shiftz, double *gxx, double *gxy, double *gxz, double *gyy, double *gyz, double *gzz,
   double *gupxx, double *gupxy, double *gupxz, double *gupyy, double *gupyz, double *gupzz,
   double *vx, double *vy, double *vz, double *u0, double *chi_rad, double *zeta_rad,
   double &Psi6threshold, int &rad_closure_scheme, double &Erad_atm_cut);

void recompute_conserv_rad_bhns(const cGH *cctkGH, int *cctk_lsh,
				double *tau_rad, double *S_rad_x, double *S_rad_y, double *S_rad_z,
				double *E_rad, double *F_radx, double *F_rady, double *F_radz, double *F_rad0, double *F_rad_scalar,
				double *P_radxx, double *P_radyy, double *P_radzz, double *P_radxy, double *P_radxz, double *P_radyz,
				double *phi, double *alpha, double *shiftx, double *shifty, double *shiftz, double *gxx, double *gxy, double *gxz, double *gyy, double *gyz, double *gzz,
				double *gupxx, double *gupxy, double *gupxz, double *gupyy, double *gupyz, double *gupzz,
				double *vx, double *vy, double *vz, double *u0, double *chi_rad, double *zeta_rad,
				double &Psi6threshold, int &rad_closure_scheme, double &Erad_atm_cut){

#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++)
    for(int j=0;j<cctk_lsh[1];j++)
      for(int i=0;i<cctk_lsh[0];i++) {

        int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	

	double Psi2 = exp(2.0*phi[index]);
	double Psim2= 1.0/Psi2;
	double Psi4 = Psi2*Psi2;
	double Psim4= 1.0/(Psi2*Psi2);
	double Psi6 = Psi2*Psi4;
	double Psim6= 1.0/(Psi2*Psi4);
	double gxxL = gxx[index];
	double gxyL = gxy[index];
	double gxzL = gxz[index];
	double gyyL = gyy[index];
	double gyzL = gyz[index];
	double gzzL = gzz[index];

	double gupxxL = gupxx[index];
        double gupxyL = gupxy[index];
        double gupxzL = gupxz[index];
        double gupyyL = gupyy[index];
        double gupyzL = gupyz[index];
        double gupzzL = gupzz[index];


	double shiftxL = shiftx[index];
	double shiftyL = shifty[index];
	double shiftzL = shiftz[index];

	double shift_xL = Psi4*(shiftxL*gxxL + shiftyL*gxyL + shiftzL*gxzL);
	double shift_yL = Psi4*(shiftxL*gxyL + shiftyL*gyyL + shiftzL*gyzL);
	double shift_zL = Psi4*(shiftxL*gxzL + shiftyL*gyzL + shiftzL*gzzL);

	double u0L = u0[index];
	double uxL = u0L * vx[index];
	double uyL = u0L * vy[index];
	double uzL = u0L * vz[index];

	double u_x = (gxxL*(shiftxL+vx[index]) +
		       gxyL*(shiftyL+vy[index]) +
		       gxzL*(shiftzL+vz[index]))*Psi4*u0L;
	double u_y = (gxyL*(shiftxL+vx[index]) +
		       gyyL*(shiftyL+vy[index]) +
		       gyzL*(shiftzL+vz[index]))*Psi4*u0L;
	double u_z = (gxzL*(shiftxL+vx[index]) +
		       gyzL*(shiftyL+vy[index]) +
		       gzzL*(shiftzL+vz[index]))*Psi4*u0L;


	double alpn1 = alpha[index] + 1.0;
	double alpn1_inv = 1.0/(alpha[index] + 1.0);
	double beta2 = shiftxL*shift_xL + shiftyL*shift_yL + shiftzL*shift_zL;
	double udotbeta = u0L*(vx[index]*shift_xL + vy[index]*shift_yL + vz[index]*shift_zL);
	double g_00L =beta2-alpn1*alpn1;
	double u_0L = g_00L*u0L + udotbeta;

	double E_radl = E_rad[index];
	double F_radxl = F_radx[index];
	double F_radyl = F_rady[index];
	double F_radzl = F_radz[index];
	double F_rad0l = - (F_radxl*u_x + F_radyl*u_y + F_radzl*u_z)/u_0L;
	

	double F_rad_xl = Psi4*(gxxL*F_radxl + gxyL*F_radyl + gxzL*F_radzl) + shift_xL * F_rad0l;
	double F_rad_yl = Psi4*(gxyL*F_radxl + gyyL*F_radyl + gyzL*F_radzl) + shift_yL * F_rad0l;
	double F_rad_zl = Psi4*(gxzL*F_radxl + gyzL*F_radyl + gzzL*F_radzl) + shift_zL * F_rad0l;
	double F_rad_0l = - (F_rad_xl*vx[index] + F_rad_yl*vy[index] + F_rad_zl*vz[index] );


	if (rad_closure_scheme == 0)
	  {
	    double P_radl = E_radl/3.0;
	    double temp_rad = alpn1 * u0L;
	    double temp_rad1 = temp_rad*temp_rad*(E_radl + P_radl) - P_radl + 2.0*alpn1*u0L*F_rad0l;

	    tau_rad[index] =  alpn1*alpn1*Psi6*(E_radl*u0L*u0L + 2.0*F_rad0l*u0L+P_radl*u0L*u0L) - Psi6*P_radl;
	    S_rad_x[index] =  alpn1*Psi6*((E_radl+P_radl)*u0L*u_x + F_rad0l*u_x + F_rad_xl * u0L);
	    S_rad_y[index] =  alpn1*Psi6*((E_radl+P_radl)*u0L*u_y + F_rad0l*u_y + F_rad_yl * u0L);
	    S_rad_z[index] =  alpn1*Psi6*((E_radl+P_radl)*u0L*u_z + F_rad0l*u_z + F_rad_zl * u0L);
	  }
	else{
	  double Fasq = F_rad_0l*F_rad0l + F_rad_xl*F_radxl +  F_rad_yl*F_radyl +  F_rad_zl*F_radzl;
	  double zeta_temp = sqrt(abs(F_rad_0l*F_rad0l + F_rad_xl*F_radxl +  F_rad_yl*F_radyl +  F_rad_zl*F_radzl)/(E_radl*E_radl));

	  double zeta, chi;
	  //	  double zeta_cut = Erad_atm_cut*1.5;
	  double zeta_cut = 1.0e-40;
	  if (E_radl<=zeta_cut){ 
	    zeta = 1.0;
	  }
	  else{
	    zeta = zeta_temp;
	  }
	  
	  if (zeta > 1.0) zeta = 1.0;
				     
	  chi = 1/3.0 + zeta*zeta*(6.0-2.0*zeta+6.0*zeta*zeta)/15.0;

	  double P_radxxl, P_radyyl, P_radzzl, P_radxyl, P_radxzl, P_radyzl;
	  if (E_radl <= Erad_atm_cut)
	    {        
	      P_radxxl = 0.0;
	      P_radyyl = 0.0;
	      P_radzzl = 0.0;
	      P_radxyl = 0.0;
	      P_radxzl = 0.0;
	      P_radyzl = 0.0;
	    }
	  else{
	    if (Fasq <= 0.0){
	      P_radxxl = E_radl*(gupxxL/Psi4 - shiftxL*shiftxL/(alpn1*alpn1) + uxL*uxL)/2.0*(1.0-chi);
	      P_radyyl = E_radl*(gupyyL/Psi4 - shiftyL*shiftyL/(alpn1*alpn1) + uyL*uyL)/2.0*(1.0-chi);
	      P_radzzl = E_radl*(gupzzL/Psi4 - shiftzL*shiftzL/(alpn1*alpn1) + uzL*uzL)/2.0*(1.0-chi);
	      P_radxyl = E_radl*(gupxyL/Psi4 - shiftxL*shiftyL/(alpn1*alpn1) + uxL*uyL)/2.0*(1.0-chi);
	      P_radxzl = E_radl*(gupxzL/Psi4 - shiftxL*shiftzL/(alpn1*alpn1) + uxL*uzL)/2.0*(1.0-chi);
	      P_radyzl = E_radl*(gupyzL/Psi4 - shiftyL*shiftzL/(alpn1*alpn1) + uyL*uzL)/2.0*(1.0-chi);
	    }
	    else{
	      P_radxxl = E_radl*((F_radxl*F_radxl/Fasq)*(3.0*chi -1.0)/2.0 + (gupxxL/Psi4 - shiftxL*shiftxL/(alpn1*alpn1) + uxL*uxL)/2.0*(1.0-chi));
	      P_radyyl = E_radl*((F_radyl*F_radyl/Fasq)*(3.0*chi -1.0)/2.0 + (gupyyL/Psi4 - shiftyL*shiftyL/(alpn1*alpn1) + uyL*uyL)/2.0*(1.0-chi));
	      P_radzzl = E_radl*((F_radzl*F_radzl/Fasq)*(3.0*chi -1.0)/2.0 + (gupzzL/Psi4 - shiftzL*shiftzL/(alpn1*alpn1) + uzL*uzL)/2.0*(1.0-chi));
	      P_radxyl = E_radl*((F_radxl*F_radyl/Fasq)*(3.0*chi -1.0)/2.0 + (gupxyL/Psi4 - shiftxL*shiftyL/(alpn1*alpn1) + uxL*uyL)/2.0*(1.0-chi));
	      P_radxzl = E_radl*((F_radxl*F_radzl/Fasq)*(3.0*chi -1.0)/2.0 + (gupxzL/Psi4 - shiftxL*shiftzL/(alpn1*alpn1) + uxL*uzL)/2.0*(1.0-chi));
	      P_radyzl = E_radl*((F_radyl*F_radzl/Fasq)*(3.0*chi -1.0)/2.0 + (gupyzL/Psi4 - shiftyL*shiftzL/(alpn1*alpn1) + uyL*uzL)/2.0*(1.0-chi));
	    }
	  }
	  
	  double P_rad0xl = - (P_radxxl * u_x + P_radxyl * u_y + P_radxzl * u_z)/u_0L;
	  double P_rad0yl = - (P_radxyl * u_x + P_radyyl * u_y + P_radyzl * u_z)/u_0L;
	  double P_rad0zl = - (P_radxzl * u_x + P_radyzl * u_y + P_radzzl * u_z)/u_0L;
	  double P_rad00l = - (P_rad0xl * u_x + P_rad0yl * u_y + P_rad0zl * u_z)/u_0L;
	  F_rad0[index] = F_rad0l;
	  F_rad_scalar[index] = sqrt(Psi4*(gxxL*SQR(F_radxl)+gyyL*SQR(F_radyl)+gzzL*SQR(F_radzl) + 2.0*( gxyL*F_radxl*F_radyl + gxzL*F_radxl*F_radzl + gyzL*F_radyl*F_radzl))  );
	  P_radxx[index] = P_radxxl;
	  P_radyy[index] = P_radyyl;
	  P_radzz[index] = P_radzzl;
          P_radxy[index] = P_radxyl;
	  P_radxz[index] = P_radxzl;
          P_radyz[index] = P_radyzl;


	  tau_rad[index] = alpn1*alpn1*Psi6*(E_radl*u0L*u0L+2.0*F_rad0l*u0L+P_rad00l);
	  S_rad_x[index] = alpn1*Psi6*(E_radl*u0L*u_x + F_rad0l*u_x + F_rad_xl * u0L + P_rad00l*shift_xL + Psi4*(P_rad0xl*gxxL + P_rad0yl*gxyL + P_rad0zl*gxzL));
	  S_rad_y[index] = alpn1*Psi6*(E_radl*u0L*u_y + F_rad0l*u_y + F_rad_yl * u0L + P_rad00l*shift_yL + Psi4*(P_rad0xl*gxyL + P_rad0yl*gyyL + P_rad0zl*gyzL));
	  S_rad_z[index] = alpn1*Psi6*(E_radl*u0L*u_z + F_rad0l*u_z + F_rad_zl * u0L + P_rad00l*shift_zL + Psi4*(P_rad0xl*gxzL + P_rad0yl*gyzL + P_rad0zl*gzzL));
	}

      }
}

extern "C" void CCTK_FCALL CCTK_FNAME(recompute_conserv_rad_bhns)
  (const cGH **cctkGH, int *cctk_lsh,
   double *tau_rad, double *S_rad_x, double *S_rad_y, double *S_rad_z,
   double *E_rad, double *F_radx, double *F_rady, double *F_radz, double *F_rad0, double *F_rad_scalar,
   double *P_radxx, double *P_radyy, double *P_radzz, double *P_radxy, double *P_radxz, double *P_radyz,
   double *phi, double *alpha, double *shiftx, double *shifty, double *shiftz, double *gxx, double *gxy, double *gxz, double *gyy, double *gyz, double *gzz,
   double *gupxx, double *gupxy, double *gupxz, double *gupyy, double *gupyz, double *gupzz,
   double *vx, double *vy, double *vz, double *u0, double *chi_rad, double *zeta_rad,
   double &Psi6threshold, int &rad_closure_scheme, double &Erad_atm_cut)
   {
  recompute_conserv_rad_bhns(*cctkGH,cctk_lsh,
			     tau_rad, S_rad_x, S_rad_y, S_rad_z, 
			     E_rad, F_radx, F_rady, F_radz, F_rad0, F_rad_scalar,
			     P_radxx, P_radyy, P_radzz, P_radxy, P_radxz, P_radyz,
			     phi, alpha, shiftx, shifty, shiftz, gxx, gxy, gxz, gyy, gyz, gzz,
                             gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,
			     vx, vy, vz, u0, chi_rad,zeta_rad,
			     Psi6threshold, rad_closure_scheme, Erad_atm_cut
			     );

}

