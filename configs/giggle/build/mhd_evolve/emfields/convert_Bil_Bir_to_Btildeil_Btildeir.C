//--------------------------------------------------------------------// 
// Conversion between normal B^i and \tilde{B}^i                      //
//--------------------------------------------------------------------//
//
#include <stdio.h>
#include <math.h>
#include "cctk.h"

extern "C" void CCTK_FCALL convert_bil_bir_to_btildeil_btildeir_
  (const cGH **cctkGH, int *ext,
   double *Bxl,double *Bxr,double *Byl,double *Byr,double *Bzl,double *Bzr,
   double *phi_f);

extern "C" void convert_Bil_Bir_to_Btildeil_Btildeir(const cGH *cctkGH, int *ext,
						     double *Bxl,double *Bxr,double *Byl,double *Byr,double *Bzl,double *Bzr,
						     double *phi_f) {
  //
#pragma omp parallel for
  for(int k=0;k<ext[2];k++) for(int j=0;j<ext[1];j++) for(int i=0;i<ext[0];i++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    double sqrtg = exp(6.0*phi_f[index]);
    Bxl[index] *= sqrtg;
    Bxr[index] *= sqrtg;
    Byl[index] *= sqrtg;
    Byr[index] *= sqrtg;
    Bzl[index] *= sqrtg;
    Bzr[index] *= sqrtg;
  }
}


extern "C" void CCTK_FCALL convert_bil_bir_to_btildeil_btildeir_
  (const cGH **cctkGH, int *ext,
   double *Bxl,double *Bxr,double *Byl,double *Byr,double *Bzl,double *Bzr,
   double *phi_f)
{
  convert_Bil_Bir_to_Btildeil_Btildeir(*cctkGH, ext,
				       Bxl,Bxr,Byl,Byr,Bzl,Bzr,
				       phi_f);
}
