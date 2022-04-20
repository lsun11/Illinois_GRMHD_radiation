#include <stdio.h>
#include <stdlib.h>
#include "cctk.h"
#include <math.h>
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

#define KRANC_C
#include "GenericFD.h"

extern "C" void CCTK_FCALL hbpuncture_shift_simple_rhs_
  (const cGH **cctkGH,int *cctk_lsh,int *nghostzones,double *dx,double *dy,double *dz,
   double *eta,
   double *shiftx_old,double *shifty_old,double *shiftz_old,
   double *shiftx_rhs,double *shifty_rhs,double *shiftz_rhs,
   double *Gammax,double *Gammay,double *Gammaz,
   double *X, double *Y,double *Z, double *eta_final_value,int *eta_falloff_enable,double *eta_falloff_radius,double *eta_falloff_dr,
   double *dtshiftx_rhs,double *dtshifty_rhs,double *dtshiftz_rhs);

extern "C" void hbpuncture_shift_SIMPLE_RHS(const cGH *cctkGH,int *cctk_lsh,int *nghostzones,double dx,double dy,double dz,
					    double eta,
					    double *shiftx_old,double *shifty_old,double *shiftz_old,
					    double *shiftx_rhs,double *shifty_rhs,double *shiftz_rhs,
					    double *Gammax,double *Gammay,double *Gammaz,
					    double *X, double *Y,double *Z, double eta_final_value, int eta_falloff_enable, double eta_falloff_radius, double eta_falloff_dr,
					    double *dtshiftx_rhs,double *dtshifty_rhs,double *dtshiftz_rhs) {

  /* Initialise finite differencing variables.  NEED THIS FOR GenericFD.h */
#include "../../GenFD_decl_set_varCPP.h"

  /* Set up variables used in the grid loop for the physical grid points */
  //WARNING: MUST BE CAREFUL WITH 4TH ORDER UPWIND STENCILS: They extend out 3 gridpoints!!
  int istart = nghostzones[0];
  int jstart = nghostzones[1];
  int kstart = nghostzones[2];
  int iend = cctk_lsh[0] - nghostzones[0];
  int jend = cctk_lsh[1] - nghostzones[1];
  int kend = cctk_lsh[2] - nghostzones[2];

  
  if(eta_falloff_enable==0) {
    
#pragma omp parallel for
    for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
	  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  //  dt \beta^i    =3/4 \tilde{Gamma}^i - \eta \beta^i
	  shiftx_rhs[index] = 0.75*Gammax[index] - eta*shiftx_old[index]; // + 0.004 * ( D11gf(shiftx_old,i,j,k) + D22gf(shiftx_old,i,j,k) + D33gf(shiftx_old,i,j,k) );
	  shifty_rhs[index] = 0.75*Gammay[index] - eta*shifty_old[index]; // + 0.004 * ( D11gf(shifty_old,i,j,k) + D22gf(shifty_old,i,j,k) + D33gf(shifty_old,i,j,k) );
	  shiftz_rhs[index] = 0.75*Gammaz[index] - eta*shiftz_old[index]; // + 0.004 * ( D11gf(shiftz_old,i,j,k) + D22gf(shiftz_old,i,j,k) + D33gf(shiftz_old,i,j,k) );

	  dtshiftx_rhs[index] = 0;
	  dtshifty_rhs[index] = 0;
	  dtshiftz_rhs[index] = 0;

	}
  }    
    
  else if(eta_falloff_enable==1) {
#pragma omp parallel for
    for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
	  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  //  dt \beta^i    =3/4 \tilde{Gamma}^i - \eta \beta^i
	  double r=sqrt(X[index]*X[index] + Y[index]*Y[index] + Z[index]*Z[index]);
	  shiftx_rhs[index] = 0.75*Gammax[index] - shiftx_old[index]*(-erf((r-eta_falloff_radius-eta_falloff_dr*2.0)/eta_falloff_dr)*(eta-eta_final_value)*0.5 + (eta-eta_final_value)*0.5 + eta_final_value);
	  shifty_rhs[index] = 0.75*Gammay[index] - shifty_old[index]*(-erf((r-eta_falloff_radius-eta_falloff_dr*2.0)/eta_falloff_dr)*(eta-eta_final_value)*0.5 + (eta-eta_final_value)*0.5 + eta_final_value);
	  shiftz_rhs[index] = 0.75*Gammaz[index] - shiftz_old[index]*(-erf((r-eta_falloff_radius-eta_falloff_dr*2.0)/eta_falloff_dr)*(eta-eta_final_value)*0.5 + (eta-eta_final_value)*0.5 + eta_final_value);
	  
	  dtshiftx_rhs[index] = 0;
	  dtshifty_rhs[index] = 0;
	  dtshiftz_rhs[index] = 0;
	  
	}
  }
  else{
    printf("Stopping: eta_falloff_enable > 1 is not supported with Spatial_Gauge==7\n");
    exit(1);
  }
      
	
//       }


// #pragma omp parallel for
//   for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
// 	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
// 	//  dt \beta^i    =3/4 \tilde{Gamma}^i - \eta \beta^i
// 	shiftx_rhs[index] = 0.75*Gammax[index]; // - eta*shiftx_old[index]; // + 0.004 * ( D11gf(shiftx_old,i,j,k) + D22gf(shiftx_old,i,j,k) + D33gf(shiftx_old,i,j,k) );
// 	shifty_rhs[index] = 0.75*Gammay[index]; // - eta*shifty_old[index]; // + 0.004 * ( D11gf(shifty_old,i,j,k) + D22gf(shifty_old,i,j,k) + D33gf(shifty_old,i,j,k) );
// 	shiftz_rhs[index] = 0.75*Gammaz[index]; // - eta*shiftz_old[index]; // + 0.004 * ( D11gf(shiftz_old,i,j,k) + D22gf(shiftz_old,i,j,k) + D33gf(shiftz_old,i,j,k) );
	
// 	if(eta_falloff_enable==0) {
// 	  shiftx_rhs[index] -= eta*shiftx_old[index]; // + 0.004 * ( D11gf(shiftx_old,i,j,k) + D22gf(shiftx_old,i,j,k) + D33gf(shiftx_old,i,j,k) );
// 	  shifty_rhs[index] -= eta*shifty_old[index]; // + 0.004 * ( D11gf(shifty_old,i,j,k) + D22gf(shifty_old,i,j,k) + D33gf(shifty_old,i,j,k) );
// 	  shiftz_rhs[index] -= eta*shiftz_old[index]; // + 0.004 * ( D11gf(shiftz_old,i,j,k) + D22gf(shiftz_old,i,j,k) + D33gf(shiftz_old,i,j,k) );
// 	}
// 	else if(eta_falloff_enable==1) {
// 	  double r=sqrt(X[index]*X[index] + Y[index]*Y[index] + Z[index]*Z[index]);
// 	  shiftx_rhs[index] -= shiftx_old[index]*(-erf((r-eta_falloff_radius-eta_falloff_dr*2.0)/eta_falloff_dr)*(eta-eta_final_value)*0.5 + (eta-eta_final_value)*0.5 + eta_final_value);
// 	  shifty_rhs[index] -= shifty_old[index]*(-erf((r-eta_falloff_radius-eta_falloff_dr*2.0)/eta_falloff_dr)*(eta-eta_final_value)*0.5 + (eta-eta_final_value)*0.5 + eta_final_value);
// 	  shiftz_rhs[index] -= shiftz_old[index]*(-erf((r-eta_falloff_radius-eta_falloff_dr*2.0)/eta_falloff_dr)*(eta-eta_final_value)*0.5 + (eta-eta_final_value)*0.5 + eta_final_value);
// 	}
// 	else{
// 	  printf("Stopping: eta_falloff_enable > 1 is not supported with Spatial_Gauge==7\n");
// 	  exit(1);
// 	}

	
// 	dtshiftx_rhs[index] = 0;
// 	dtshifty_rhs[index] = 0;
// 	dtshiftz_rhs[index] = 0;
//       }



}

extern "C" void CCTK_FCALL hbpuncture_shift_simple_rhs_
  (const cGH **cctkGH,int *cctk_lsh,int *nghostzones,double *dx,double *dy,double *dz,
   double *eta,
   double *shiftx_old,double *shifty_old,double *shiftz_old,
   double *shiftx_rhs,double *shifty_rhs,double *shiftz_rhs,
   double *Gammax,double *Gammay,double *Gammaz,
   double *X, double *Y,double *Z,double *eta_final_value,int *eta_falloff_enable,double *eta_falloff_radius,double *eta_falloff_dr,
   double *dtshiftx_rhs,double *dtshifty_rhs,double *dtshiftz_rhs)
{
  hbpuncture_shift_SIMPLE_RHS(*cctkGH,cctk_lsh,nghostzones,*dx,*dy,*dz,
			      *eta,
			      shiftx_old,shifty_old,shiftz_old,
			      shiftx_rhs,shifty_rhs,shiftz_rhs,
			      Gammax,Gammay,Gammaz,
			      X,Y,Z,*eta_final_value,*eta_falloff_enable,*eta_falloff_radius,*eta_falloff_dr,
			      dtshiftx_rhs,dtshifty_rhs,dtshiftz_rhs);
}
