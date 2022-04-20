//-------------------------------------------------------------------------------------
// Sanity check: Compute determinate of \tilde gamma and trace of A; enforce Tr Aij = 0
//-------------------------------------------------------------------------------------

#include "cctk.h"
#include "math.h"

#define F1o3 0.3333333333333333333333333333

extern "C" void CCTK_FCALL CCTK_FNAME(sanitycheck_restore_Aij)
  (const cGH **cctkGH,int *cctk_lsh,double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz, 
   double *Axx,double *Axy,double *Axz,double *Ayy,double *Ayz,double *Azz);

extern "C" void sanitycheck_restore_Aij(const cGH *cctkGH,int *cctk_lsh,double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz, 
					double *Axx,double *Axy,double *Axz,double *Ayy,double *Ayz,double *Azz) {

  /* Set up variables used in the grid loop for the physical grid points */
  int istart,jstart,kstart;
  istart=jstart=kstart=0;
  int iend = cctk_lsh[0];
  int jend = cctk_lsh[1];
  int kend = cctk_lsh[2];

#pragma omp parallel for
  for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    
    double gxxL = gxx[index];
    double gxyL = gxy[index];
    double gxzL = gxz[index];
    double gyyL = gyy[index];
    double gyzL = gyz[index];
    double gzzL = gzz[index];

    double AxxL = Axx[index];
    double AxyL = Axy[index];
    double AxzL = Axz[index];
    double AyyL = Ayy[index];
    double AyzL = Ayz[index];
    double AzzL = Azz[index];

    // compute and check determinant
    double det =  gxxL * gyyL * gzzL + gxyL * gyzL * gxzL + gxzL * gxyL * gyzL 
      - gxzL * gyyL * gxzL - gxyL * gxyL * gzzL - gxxL * gyzL * gyzL;

    // invert metric
    double gupxx =   ( gyyL * gzzL - gyzL * gyzL ) / det;
    double gupxy = - ( gxyL * gzzL - gyzL * gxzL ) / det;
    double gupxz =   ( gxyL * gyzL - gyyL * gxzL ) / det;
    double gupyy =   ( gxxL * gzzL - gxzL * gxzL ) / det;
    double gupyz = - ( gxxL * gyzL - gxyL * gxzL ) / det;
    double gupzz =   ( gxxL * gyyL - gxyL * gxyL ) / det;

    // Compute trace of A
    double trA =  gupxx * AxxL + gupyy * AyyL + gupzz * AzzL +  
      2.0 * ( gupxy * AxyL + gupxz * AxzL + gupyz * AyzL );

    Axx[index] += - F1o3 * trA * gxxL;
    Axy[index] += - F1o3 * trA * gxyL;
    Axz[index] += - F1o3 * trA * gxzL;
    Ayy[index] += - F1o3 * trA * gyyL;
    Ayz[index] += - F1o3 * trA * gyzL;
    Azz[index] += - F1o3 * trA * gzzL;
  }
}

extern "C" void CCTK_FCALL CCTK_FNAME(sanitycheck_restore_Aij)
  (const cGH **cctkGH,int *cctk_lsh,double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz, 
   double *Axx,double *Axy,double *Axz,double *Ayy,double *Ayz,double *Azz)
{
  sanitycheck_restore_Aij(*cctkGH,cctk_lsh,gxx,gxy,gxz,gyy,gyz,gzz, 
			  Axx,Axy,Axz,Ayy,Ayz,Azz);
}
