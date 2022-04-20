/*
  Compute gupij's.  You need to do this after every regridding,
  before computing the Ricci tensor,
  and before computing BSSN RHS variables.
  It is possible that one of the above is redundant, but hey, this is a cheap operation.
*/

#include <stdio.h>
#include "cctk.h"
#include <math.h>
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

static char *rcsid="$mew. $";
CCTK_FILEVERSION(BSSN_compute_gupij)


  extern "C" void CCTK_FCALL CCTK_FNAME(BSSN_compute_gupij)
  (const cGH **cctkGH,int *cctk_lsh,
   double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
   double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz);


extern "C" void BSSN_compute_gupij(const cGH *cctkGH,int *cctk_lsh,
				   double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
				   double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz) {
 
  DECLARE_CCTK_PARAMETERS;

  #pragma omp parallel for
    for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
      int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
      double gxxL = gxx[index];
      double gxyL = gxy[index];
      double gxzL = gxz[index];
      double gyyL = gyy[index];
      double gyzL = gyz[index];
      double gzzL = gzz[index];
	
      double gijdet = gxxL * gyyL * gzzL + gxyL * gyzL * gxzL + gxzL * gxyL * gyzL 
	- gxzL * gyyL * gxzL - gxyL * gxyL * gzzL - gxxL * gyzL * gyzL;
      double gijdetinv = 1.0/gijdet;

      gupxx[index] =   ( gyyL * gzzL - gyzL * gyzL )* gijdetinv;
      gupxy[index] = - ( gxyL * gzzL - gyzL * gxzL )* gijdetinv;
      gupxz[index] =   ( gxyL * gyzL - gyyL * gxzL )* gijdetinv;
      gupyy[index] =   ( gxxL * gzzL - gxzL * gxzL )* gijdetinv;
      gupyz[index] = - ( gxxL * gyzL - gxyL * gxzL )* gijdetinv;
      gupzz[index] =   ( gxxL * gyyL - gxyL * gxyL )* gijdetinv;
    }
}


extern "C" void CCTK_FCALL CCTK_FNAME(BSSN_compute_gupij)
  (const cGH **cctkGH,int *cctk_lsh,
   double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
   double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz)
{
  BSSN_compute_gupij(*cctkGH,cctk_lsh,
		     gxx,gxy,gxz,gyy,gyz,gzz,
		     gupxx,gupxy,gupxz,gupyy,gupyz,gupzz);
}
