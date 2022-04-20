//-----------------------------------------------------------------------------
// Advect Z_i
//-----------------------------------------------------------------------------
#include "math.h"
#include "cctk.h"
#include "stdio.h"
#include "stdlib.h"

#define SQR(x) ((x) * (x))

extern "C" void CCTK_FCALL CCTK_FNAME(advect_Srad_cpp)
  (int *flux_direction, const cGH **cctkGH,int *cctk_lsh, int *nghostzones, int *Symmetry, 
   double *S_rad_x_rhs, double *S_rad_y_rhs, double *S_rad_z_rhs, 
   double *S_radx_flux, double *S_rady_flux, double *S_radz_flux, 
   double *dX,double *dY,double *dZ);

extern "C" void advect_Srad_cpp(int flux_direction, const cGH *cctkGH,int *cctk_lsh, int *nghostzones, int Symmetry,
			     double *S_rad_x_rhs, double *S_rad_y_rhs, double *S_rad_z_rhs,
			     double *S_radx_flux, double *S_rady_flux, double *S_radz_flux,
			     double dX,double dY,double dZ) {
  double sfpi = sqrt(4.0*M_PI);

  int AXISYM = 4;

  double f1o4pi = 1.0/(4.0*M_PI);

  /* Set up variables used in the grid loop for the physical grid points */

  int istart = 1;
  int jstart = 1;
  int kstart = 1;
  int iend = cctk_lsh[0] - 1;
  int jend = cctk_lsh[1] - 1;
  int kend = cctk_lsh[2] - 1;
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

  if(Symmetry==4) {
    jstart = 0;
    jend = cctk_lsh[1];
    jstart++;
    jend--;
  }
  
  if (flux_direction==1) {
    if (Symmetry==AXISYM) {
#pragma omp parallel for
      for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
	//Lunan: Since axisymmetry is no longer used, we turn this off for radiation. 

	printf("RADIATION FIELD IS NOT SUPPORTED IN AXISYMMETRY RTGHT NOW!\n"); exit(1);
      } 
    } else {
#pragma omp parallel for
      for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
	int index   = CCTK_GFINDEX3D(cctkGH,i,  j,k);  
	int indexp1 = CCTK_GFINDEX3D(cctkGH,i+1,j,k);  
	S_rad_x_rhs[index] += (S_radx_flux[index]-S_radx_flux[indexp1])*dxi;
	S_rad_y_rhs[index] += (S_rady_flux[index]-S_rady_flux[indexp1])*dxi;
	S_rad_z_rhs[index] += (S_radz_flux[index]-S_radz_flux[indexp1])*dxi;
      }
    }
  } else if (flux_direction==2) {
    if (Symmetry != AXISYM) {
#pragma omp parallel for
      for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
	int index   = CCTK_GFINDEX3D(cctkGH,i,j,  k);  
	int indexp1 = CCTK_GFINDEX3D(cctkGH,i,j+1,k);  
	S_rad_x_rhs[index] += (S_radx_flux[index]-S_radx_flux[indexp1])*dyi;
	S_rad_y_rhs[index] += (S_rady_flux[index]-S_rady_flux[indexp1])*dyi;
	S_rad_z_rhs[index] += (S_radz_flux[index]-S_radz_flux[indexp1])*dyi;
      }
    }
  } else if (flux_direction==3) {
#pragma omp parallel for
    for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
      int index   = CCTK_GFINDEX3D(cctkGH,i,j,k  );
      int indexp1 = CCTK_GFINDEX3D(cctkGH,i,j,k+1);
      S_rad_x_rhs[index] += (S_radx_flux[index]-S_radx_flux[indexp1])*dzi;
      S_rad_y_rhs[index] += (S_rady_flux[index]-S_rady_flux[indexp1])*dzi;
      S_rad_z_rhs[index] += (S_radz_flux[index]-S_radz_flux[indexp1])*dzi;
    }
  }
}
extern "C" void CCTK_FCALL CCTK_FNAME(advect_Srad_cpp)
  (int *flux_direction, const cGH **cctkGH,int *cctk_lsh, int *nghostzones, int *Symmetry,
   double *S_rad_x_rhs, double *S_rad_y_rhs, double *S_rad_z_rhs,
   double *S_radx_flux, double *S_rady_flux, double *S_radz_flux,
   double *dX,double *dY,double *dZ)
{
  advect_Srad_cpp(*flux_direction,*cctkGH,cctk_lsh, nghostzones, *Symmetry,
	       S_rad_x_rhs,S_rad_y_rhs, S_rad_z_rhs,
	       S_radx_flux,S_rady_flux,S_radz_flux,	   
	       *dX,*dY,*dZ);
}
