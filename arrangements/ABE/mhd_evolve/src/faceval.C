//-------------------------
// Calculate face averages
//-------------------------
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

#include "stdlib.h"
#include "math.h"

//FORTRAN-callable headers:
extern "C" void CCTK_FCALL CCTK_FNAME(compute_face_avg)
  (int *flux_direction,const cGH **cctkGH,int *cctk_lsh,int *nghostzones,int *Symmetry,
   double *var,double *var_face_avg);

//Compute face value for var
extern "C" void compute_face_avg(int flux_direction,const cGH *cctkGH,int *cctk_lsh,int *nghostzones,int Symmetry,
				 double *var,double *var_face_avg) {
  int istart,jstart,kstart;
  istart=jstart=kstart=0;
  int iend = cctk_lsh[0];
  int jend = cctk_lsh[1];
  int kend = cctk_lsh[2];

  // In axisymmetry, we only update hydro variables on y=0 plane:
  if(Symmetry==4) { jstart=1; jend=2; }

  if(flux_direction==1) {
    istart++;
  } else if(flux_direction==2) {
    jstart++;
  } else if(flux_direction==3){
    kstart++;
  }

#pragma omp parallel for
  for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
    int index,indexm1;
    if(flux_direction==1) {
      index   = CCTK_GFINDEX3D(cctkGH,i,  j,k);
      indexm1 = CCTK_GFINDEX3D(cctkGH,i-1,j,k);
    } else if(flux_direction==2) {
      index   = CCTK_GFINDEX3D(cctkGH,i,j,  k);
      indexm1 = CCTK_GFINDEX3D(cctkGH,i,j-1,k);
    } else if(flux_direction==3){
      index   = CCTK_GFINDEX3D(cctkGH,i,j,k  );
      indexm1 = CCTK_GFINDEX3D(cctkGH,i,j,k-1);
    }
    var_face_avg[index] = 0.5*(var[index] + var[indexm1]);
  }
}

extern "C" void CCTK_FCALL CCTK_FNAME(compute_face_avg)
  (int *flux_direction,const cGH **cctkGH,int *cctk_lsh,int *nghostzones,int *Symmetry,
   double *var,double *var_face_avg)
{
  compute_face_avg(*flux_direction,*cctkGH,cctk_lsh,nghostzones,*Symmetry,var,var_face_avg);
}
