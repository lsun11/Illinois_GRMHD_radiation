//------------------------------------------------------//
// Advect B^i					       //
//------------------------------------------------------//

#include "cctk.h"
#include "stdlib.h"
#include "stdio.h"
#include "math.h"

#define NO_SYMM 0
#define OCTANT  2
#define AXISYM  4

extern "C" void CCTK_FCALL advect_rad_tau_ct_cpp_
  (int *flux_direction, const cGH **cctkGH,int *ext,double *tau_rad_rhs, 
   double *tau_rad_flux,int &Symmetry, double *dX,double *dY,double *dZ, double *tau_rad_flux_x, double *tau_rad_advect_diag);


//
extern "C" void advect_rad_tau_ct_cpp(int flux_direction, const cGH *cctkGH,int *ext, double *tau_rad_rhs,
				      double *tau_rad_flux,int &Symmetry,double dX,double dY,double dZ, double *tau_rad_flux_x, double *tau_rad_advect_diag) {

  double dxi = 1.0/dX; 
  double dyi = 1.0/dY; 
  double dzi = 1.0/dZ; 

  /* Set up variables used in the grid loop for the physical grid points */
  int imin = 1;
  int jmin = 1;
  int kmin = 1;
  int imax = ext[0] - 1;
  int jmax = ext[1] - 1;
  int kmax = ext[2] - 1;

  //

    if (flux_direction==1){
      #pragma omp parallel for
      for(int k=kmin;k<kmax;k++) for(int j=jmin;j<jmax;j++) for(int i=imin;i<imax;i++) { 
	    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	    int indexp1 = CCTK_GFINDEX3D(cctkGH,i+1,j,k);
	    tau_rad_rhs[index] +=   (tau_rad_flux[index]-tau_rad_flux[indexp1])*dxi;
	    tau_rad_advect_diag[index] +=   (tau_rad_flux[index]-tau_rad_flux[indexp1])*dxi;
	    //	    tau_rad_flux_xp1[index] = tau_rad_flux[indexp1];  
	  }
    }
    else if (flux_direction==2){
       #pragma omp parallel for
      for(int k=kmin;k<kmax;k++) for(int j=jmin;j<jmax;j++) for(int i=imin;i<imax;i++) {
            int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
            int indexp1 = CCTK_GFINDEX3D(cctkGH,i,j+1,k);
            tau_rad_rhs[index] +=   (tau_rad_flux[index]-tau_rad_flux[indexp1])*dyi;
	    tau_rad_advect_diag[index] += (tau_rad_flux[index]-tau_rad_flux[indexp1])*dyi;
	  }
    }
    else if (flux_direction==3){
      #pragma omp parallel for
      for(int k=kmin;k<kmax;k++) for(int j=jmin;j<jmax;j++) for(int i=imin;i<imax;i++) {
            int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
            int indexp1 = CCTK_GFINDEX3D(cctkGH,i,j,k+1);
            tau_rad_rhs[index] +=   (tau_rad_flux[index]-tau_rad_flux[indexp1])*dzi;
	    tau_rad_advect_diag[index] += (tau_rad_flux[index]-tau_rad_flux[indexp1])*dzi;
          }
    }

  //  printf("end advect_rad_tau_ct.C !!!!! \n");
}


extern "C" void CCTK_FCALL advect_rad_tau_ct_cpp_
  (int *flux_direction, const cGH **cctkGH,int *ext,double *tau_rad_rhs,
   double *tau_rad_flux,int &Symmetry,double *dX,double *dY,double *dZ, double *tau_rad_flux_x, double *tau_rad_advect_diag)
{
  advect_rad_tau_ct_cpp(*flux_direction, *cctkGH,ext,tau_rad_rhs,
			tau_rad_flux,Symmetry, *dX, *dY, *dZ, tau_rad_flux_x, tau_rad_advect_diag);
}
