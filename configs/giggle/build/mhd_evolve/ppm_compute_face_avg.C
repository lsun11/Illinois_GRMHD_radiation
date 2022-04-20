// Use PPM scheme to reconstruct a smooth function.  This is for interpolating metric quantities onto the cell faces at i-1/2 (x dir'n), j-1/2 (y dir'n), and k-1/2 (z dir'n)
//

// Note: these values are for cell averaged gfs: 
// am2=-1.0/12.0, am1=7.0/12.0, a0=7.0/12.0, a1=-1.0/12.0
// Since the metric gfs store the grid point values instead of the cell average, 
// the following coefficients should be used: 
// am2 = -1/16, am1 = 9/16, a0 = 9/16, a1 = -1/16
#define AM2 -0.0625
#define AM1  0.5625
#define A0   0.5625
#define A1  -0.0625

// Off-center differencing coefficients: 
// (b) compute value at i-1/2 using points at i-3,i-2,i-1,i:
//     bm3=1/16, bm2=-5/16, bm1=15/16, b0 = 5/16
#define bm3  0.0625
#define bm2  -0.3125
#define bm1  0.9375
#define b0   0.3125

// (c) compute value at i-1/2 using points at i-1,i,i+1,i+2:
//     cm1=5/16, c0=15/16, c1=-5/16, c2=1/16
#define cm1  0.3125
#define c0   0.9375
#define c1   -0.3125
#define c2   0.0625

// (d) compute value at i-1/2 using points at i,i+1,i+2,i+3:
//     d0=35/16, d1=-35/16, d2=21/16, d3=-5/16
#define d0   2.1875
#define d1   -2.1875
#define d2   1.3125
#define d3   -0.3125

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "cctk.h"

#define AXISYM 4

extern "C" void CCTK_FCALL compute_face_avg_ppm_cpp_
  (const cGH **cctkGH,int &m,int *ext,int *nghostzones,int &cell_centering_enabled,int &Symmetry,int *sym,double *f,double *ff);

extern "C" void compute_face_avg_ppm_cpp(const cGH *cctkGH,int &m,int *ext,int *nghostzones,int &cell_centering_enabled,int &Symmetry,int *sym,double *f,double *ff) {

  int imin=0,jmin=0,kmin=0;
  int imax=ext[0];
  int jmax=ext[1];
  int kmax=ext[2];
  
  if (Symmetry==AXISYM) { 
    jmin = 1;
    jmax = 1;
  }
  
  //Fill all points except for lower and upper boundaries in m-direction.
  if(m==1) { imin+=2; imax--; }
  if(m==2) { jmin+=2; jmax--; }
  if(m==3) { kmin+=2; kmax--; }
#pragma omp parallel for
  for(int k=kmin;k<kmax;k++) for(int j=jmin;j<jmax;j++) for(int i=imin;i<imax;i++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    int indexm1,indexm2,indexp1;
    if(m==1) {
      indexm1 = CCTK_GFINDEX3D(cctkGH,i-1,j,k);
      indexm2 = CCTK_GFINDEX3D(cctkGH,i-2,j,k);
      indexp1 = CCTK_GFINDEX3D(cctkGH,i+1,j,k);
    } else if(m==2) {
      indexm1 = CCTK_GFINDEX3D(cctkGH,i,j-1,k);
      indexm2 = CCTK_GFINDEX3D(cctkGH,i,j-2,k);
      indexp1 = CCTK_GFINDEX3D(cctkGH,i,j+1,k);
    } else {
      indexm1 = CCTK_GFINDEX3D(cctkGH,i,j,k-1);
      indexm2 = CCTK_GFINDEX3D(cctkGH,i,j,k-2);
      indexp1 = CCTK_GFINDEX3D(cctkGH,i,j,k+1);
    }
    ff[index] = AM2*f[indexm2] + AM1*f[indexm1] + A0*f[index] + A1*f[indexp1];
  }

  //Fill in the upper boundaries 
  if (m==1) {
#pragma omp parallel for     
     for(int k=kmin;k<kmax;k++) for(int j=jmin;j<jmax;j++) {
        int index = CCTK_GFINDEX3D(cctkGH,imax,j,k);
	int indexm1 = CCTK_GFINDEX3D(cctkGH,imax-1,j,k);
	int indexm2 = CCTK_GFINDEX3D(cctkGH,imax-2,j,k);
	int indexm3 = CCTK_GFINDEX3D(cctkGH,imax-3,j,k);
	ff[index] = bm3*f[indexm3] + bm2*f[indexm2] + bm1*f[indexm1] + b0*f[index];
     }
  } else if (m==2) {
#pragma omp parallel for
     for(int k=kmin;k<kmax;k++) for(int i=imin;i<imax;i++) {
        int index = CCTK_GFINDEX3D(cctkGH,i,jmax,k);
        int indexm1 = CCTK_GFINDEX3D(cctkGH,i,jmax-1,k);
        int indexm2 = CCTK_GFINDEX3D(cctkGH,i,jmax-2,k);
        int indexm3 = CCTK_GFINDEX3D(cctkGH,i,jmax-3,k);
        ff[index] = bm3*f[indexm3] + bm2*f[indexm2] + bm1*f[indexm1] + b0*f[index];
     }
  } else {
#pragma omp parallel for
     for(int j=jmin;j<jmax;j++) for(int i=imin;i<imax;i++) {
        int index = CCTK_GFINDEX3D(cctkGH,i,j,kmax);
        int indexm1 = CCTK_GFINDEX3D(cctkGH,i,j,kmax-1);
        int indexm2 = CCTK_GFINDEX3D(cctkGH,i,j,kmax-2);
        int indexm3 = CCTK_GFINDEX3D(cctkGH,i,j,kmax-3);
        ff[index] = bm3*f[indexm3] + bm2*f[indexm2] + bm1*f[indexm1] + b0*f[index];
     }
  }
 
  //Fill in the lower boundaries
  if (m==1) {
     imin -=2;
#pragma omp parallel for
     for(int k=kmin;k<kmax;k++) for(int j=jmin;j<jmax;j++) {
        int index = CCTK_GFINDEX3D(cctkGH,imin,j,k);
        int index1=CCTK_GFINDEX3D(cctkGH,imin+1,j,k);
        int index2=CCTK_GFINDEX3D(cctkGH,imin+2,j,k);
        int index3=CCTK_GFINDEX3D(cctkGH,imin+3,j,k);
        ff[index] = d0*f[index] + d1*f[index1] + d2*f[index2] + d3*f[index3];
        ff[index1] = cm1*f[index] + c0*f[index1] + c1*f[index2] + c2*f[index3];
     }
  } else if (m==2) {
     jmin -=2;
#pragma omp parallel for
     for(int k=kmin;k<kmax;k++) for(int i=imin;i<imax;i++) {
        int index = CCTK_GFINDEX3D(cctkGH,i,jmin,k);
        int index1=CCTK_GFINDEX3D(cctkGH,i,jmin+1,k);
        int index2=CCTK_GFINDEX3D(cctkGH,i,jmin+2,k);
        int index3=CCTK_GFINDEX3D(cctkGH,i,jmin+3,k);
        ff[index] = d0*f[index] + d1*f[index1] + d2*f[index2] + d3*f[index3];
        ff[index1] = cm1*f[index] + c0*f[index1] + c1*f[index2] + c2*f[index3];
     }
  } else {
     kmin -=2;
#pragma omp parallel for
     for(int j=jmin;j<jmax;j++) for(int i=imin;i<imax;i++) {
        int index = CCTK_GFINDEX3D(cctkGH,i,j,kmin);
        int index1=CCTK_GFINDEX3D(cctkGH,i,j,kmin+1);
        int index2=CCTK_GFINDEX3D(cctkGH,i,j,kmin+2);
        int index3=CCTK_GFINDEX3D(cctkGH,i,j,kmin+3);
        ff[index] = d0*f[index] + d1*f[index1] + d2*f[index2] + d3*f[index3];
        ff[index1] = cm1*f[index] + c0*f[index1] + c1*f[index2] + c2*f[index3];
     }
  }
}

extern "C" void CCTK_FCALL compute_face_avg_ppm_cpp_
  (const cGH **cctkGH,int &m,int *ext,int *nghostzones,int &cell_centering_enabled,int &Symmetry,int *sym,double *f,double *ff)
{
  compute_face_avg_ppm_cpp(*cctkGH,m,ext,nghostzones,cell_centering_enabled,Symmetry,sym,f,ff);
}
