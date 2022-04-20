#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "cctk.h"

#define AXISYM 4

extern "C" void CCTK_FCALL ppm_shift_v2_cpp_
  (const cGH **cctkGH,int *ext,double *ar,double *al,double *temp,int &m,int &Symmetry);

extern "C" void ppm_shift_v2_cpp(const cGH *cctkGH,int *ext,double *ar,double *al,double *temp,int &m,int &Symmetry) {
  int imin=0,jmin=0,kmin=0;
  int imax=ext[0];
  int jmax=ext[1];
  int kmax=ext[2];
  
  if (Symmetry==AXISYM) { 
    jmin = 1;
    jmax = 1;
  }
  
  //Fill all points except for lower and upper boundaries in m-direction.
  //It should not be necessary to fill lower and upper boundaries, since the values will be overwritten at the next sync call.
#pragma omp parallel for
  for(int k=kmin;k<kmax;k++) for(int j=jmin;j<jmax;j++) for(int i=imin;i<imax;i++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    temp[index] = ar[index];
    ar[index] = al[index];
  }

  if(m==1) { imin++; }
  if(m==2) { jmin++; }
  if(m==3) { kmin++; }
#pragma omp parallel for
  for(int k=kmin;k<kmax;k++) for(int j=jmin;j<jmax;j++) for(int i=imin;i<imax;i++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    int indexm1,indexm2,indexp1;
    if(m==1) {
      indexm1 = CCTK_GFINDEX3D(cctkGH,i-1,j,k);
    } else if(m==2) {
      indexm1 = CCTK_GFINDEX3D(cctkGH,i,j-1,k);
    } else {
      indexm1 = CCTK_GFINDEX3D(cctkGH,i,j,k-1);
    }
    al[index] = temp[indexm1];
  }
}

extern "C" void CCTK_FCALL ppm_shift_v2_cpp_
  (const cGH **cctkGH,int *ext,double *ar,double *al,double *temp,int &m,int &Symmetry)
{
  ppm_shift_v2_cpp(*cctkGH,ext,ar,al,temp,m,Symmetry);
}
