#define F1o4pi 0.07957747154594766788
#define v2max  0.99

#include <stdio.h>
#include "cctk.h"
#include <math.h>
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

extern "C" void CCTK_FCALL bhns_reset_rho_b_atm_
  (const cGH **cctkGH, int *cctk_lsh, 
   double *X, double *Y, double *Z, double *b2,
   double *rho_b,double *P, int &atm_type, double &xNS, 
   double &yNS, double &rhob_fac, double &rad_rhob_fac,
   double &rhob_max,  double &rhobatm_falloff_power, double &rhob_o_b2);

void bhns_reset_rho_b_atm(const cGH *cctkGH, int *cctk_lsh, 
			  double *X, double *Y, double *Z, double *b2,
			  double *rho_b, double *P, int &atm_type, double &xNS, 
			  double &yNS, double &rhob_fac, double &rad_rhob_fac,
			  double &rhob_max,  double &rhobatm_falloff_power, double &rhob_o_b2){
  
#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++)
    for(int j=0;j<cctk_lsh[1];j++)
      for(int i=0;i<cctk_lsh[0];i++) {
	
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	double xc = X[index]-xNS;
	double yc = Y[index]-yNS;
	double zc = Z[index];
	double rad = sqrt(xc*xc+yc*yc+zc*zc);
	
	if (rho_b[index] <= rhob_fac*rhob_max){
	  if (atm_type==0) {
	    rho_b[index] = rhob_fac*rhob_max*pow(rad_rhob_fac/rad,rhobatm_falloff_power);
	    P[index]     = rho_b[index]*rho_b[index];
	  }
	  else if (atm_type==1){
	    rho_b[index] = rhob_o_b2*b2[index];
	  }
	}
      }
}
  
extern "C" void CCTK_FCALL bhns_reset_rho_b_atm_
  (const cGH **cctkGH, int *cctk_lsh, 
   double *X, double *Y, double *Z, double *b2,
   double *rho_b, double *P, int &atm_type, double &xNS, 
   double &yNS, double &rhob_fac, double &rad_rhob_fac,
   double &rhob_max,  double &rhobatm_falloff_power, double &rhob_o_b2)
{
  bhns_reset_rho_b_atm(*cctkGH,cctk_lsh, X,Y,Z, b2,
		       rho_b, P, atm_type, xNS,yNS, rhob_fac, rad_rhob_fac, 
		       rhob_max, rhobatm_falloff_power, rhob_o_b2);
  
}
