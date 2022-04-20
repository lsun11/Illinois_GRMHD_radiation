//-----------------------------------------------------------------------------
// Advect rho_star and e
//-----------------------------------------------------------------------------

#include "math.h"
#include "cctk.h"
#include <stdio.h>

extern "C" void CCTK_FCALL advect_rho_e_cpp_
  (int *flux_direction, const cGH **cctkGH,int *cctk_lsh, int *nghostzones, int *Symmetry,
   double *rho_star_rhs,double *e_rhs,double *frm,double *fem,
   double *dX,double *dY,double *dZ,double *X,double *Y,double *Z, double *drho_b_m_x, double *drho_b_m_xp1);

extern "C" void advect_rho_e_cpp(int flux_direction,const cGH *cctkGH,int *cctk_lsh,int *nghostzones,int Symmetry,
				 double *rho_star_rhs,double *e_rhs,double *frm,double *fem,
				 double dX,double dY,double dZ,double *X,double *Y,double *Z, double *drho_b_m_x, double *drho_b_m_xp1) {
  int AXISYM = 4;

  /* Set up variables used in the grid loop for the physical grid points */
  int istart = 1;
  int jstart = 1;
  int kstart = 1;
  int iend = cctk_lsh[0]-1;
  int jend = cctk_lsh[1]-1;
  int kend = cctk_lsh[2]-1;

  /*
    int istart = nghostzones[0];
    int jstart = nghostzones[1];
    int kstart = nghostzones[2];
    int iend = cctk_lsh[0] - nghostzones[0];
    int jend = cctk_lsh[1] - nghostzones[1];
    int kend = cctk_lsh[2] - nghostzones[2];
  */

  double dxi = 1.0/dX;
  double dyi = 1.0/dY;
  double dzi = 1.0/dZ;


  if (flux_direction==1) {
#pragma omp parallel for
    for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
	  int index   = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  int indexp1 = CCTK_GFINDEX3D(cctkGH,i+1,j,k);
	  rho_star_rhs[index] += (frm[index]-frm[indexp1]) * dxi;
	  e_rhs[index] += (fem[index] - fem[indexp1]) * dxi;

	  if (isnan(e_rhs[index])){
	    printf("inside rho_e.C Z direction tau_rhs is NAN, fem[index]=%e, fem[indexp1]=%e, dxi=%e", fem[index], fem[indexp1], dzi);
	  }

	}
  } else if (flux_direction==2) {
#pragma omp parallel for
      for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
 	int index   = CCTK_GFINDEX3D(cctkGH,i,j,  k);
	int indexp1 = CCTK_GFINDEX3D(cctkGH,i,j+1,k);
	rho_star_rhs[index] += (frm[index]-frm[indexp1]) * dyi;
       	e_rhs[index] += (fem[index] - fem[indexp1]) * dyi;
	  	  
    if (isnan(e_rhs[index])){
      printf("inside rho_e.C Y direction tau_rhs is NAN, fem[index]=%e, fem[indexp1]=%e, dxi=%e", fem[index], fem[indexp1], dyi);
    }
	  }
  }
 else if (flux_direction==3){
#pragma omp parallel for
   for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
	 int index   = CCTK_GFINDEX3D(cctkGH,i,j,k  );
	 int indexp1 = CCTK_GFINDEX3D(cctkGH,i,j,k+1);
	 rho_star_rhs[index] += (frm[index]-frm[indexp1]) * dzi;
	 e_rhs[index] += (fem[index] - fem[indexp1]) * dzi;

	 if (isnan(e_rhs[index])){
	   printf("inside rho_e.C Z direction tau_rhs is NAN, fem[index]=%e, fem[indexp1]=%e, dxi=%e", fem[index], fem[indexp1], dzi);
	 }
     }
 }
}
 

extern "C" void CCTK_FCALL advect_rho_e_cpp_
  (int *flux_direction, const cGH **cctkGH,int *cctk_lsh, int *nghostzones, int *Symmetry,
   double *rho_star_rhs,double *e_rhs,double *frm,double *fem,
   double *dX,double *dY,double *dZ,double *X,double *Y,double *Z, double *drho_b_m_x, double *drho_b_m_xp1) 
{
  advect_rho_e_cpp(*flux_direction,*cctkGH,cctk_lsh, nghostzones, *Symmetry,
		   rho_star_rhs,e_rhs,frm,fem,
		   *dX,*dY,*dZ,X,Y,Z,drho_b_m_x, drho_b_m_xp1);
}
