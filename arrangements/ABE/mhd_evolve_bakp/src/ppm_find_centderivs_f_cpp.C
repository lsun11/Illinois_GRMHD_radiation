//-----------------------------------------------------------------//
// Calculate center differences df = (f_{j+1} - f_{j-1} )/2
//                           d^2 f = f_{j+1} - 2 f_j + f_{j-1}
//-----------------------------------------------------------------//

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "cctk.h"

#define AXISYM 4

extern "C" void CCTK_FCALL CCTK_FNAME(find_centderivs_f_cpp)
  (const cGH **cctkGH,int *ext,int *nghostzones,double *f,double *df,double *d2f,int &m,double *sym,  
   int &Symmetry,int &cell_centering_enabled);


extern "C" void find_centderivs_f_cpp(const cGH *cctkGH,int *ext,int *nghostzones,double *f,double *df,double *d2f,int &m,double *sym,  
				      int &Symmetry,int &cell_centering_enabled) {
  // 
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
  if(m==1) { imin++; imax--; }
  if(m==2) { jmin++; jmax--; }
  if(m==3) { kmin++; kmax--; }
#pragma omp parallel for
  for(int k=kmin;k<kmax;k++) for(int j=jmin;j<jmax;j++) for(int i=imin;i<imax;i++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    int indexm1,indexp1;
    if(m==1) {
      indexm1 = CCTK_GFINDEX3D(cctkGH,i-1,j,k);
      indexp1 = CCTK_GFINDEX3D(cctkGH,i+1,j,k);
    } else if(m==2) {
      indexm1 = CCTK_GFINDEX3D(cctkGH,i,j-1,k);
      indexp1 = CCTK_GFINDEX3D(cctkGH,i,j+1,k);
    } else {
      indexm1 = CCTK_GFINDEX3D(cctkGH,i,j,k-1);
      indexp1 = CCTK_GFINDEX3D(cctkGH,i,j,k+1);
    }
    df[index] = 0.5*(f[indexp1] - f[indexm1]);
    d2f[index] = f[indexp1] - 2.0*f[index] + f[indexm1];
  }
}

extern "C" void CCTK_FCALL CCTK_FNAME(find_centderivs_f_cpp)
  (const cGH **cctkGH,int *ext,int *nghostzones,double *f,double *df,double *d2f,int &m,double *sym,  
   int &Symmetry,int &cell_centering_enabled)
{
  find_centderivs_f_cpp(*cctkGH,ext,nghostzones,f,df,d2f,m,sym,Symmetry,cell_centering_enabled);
}
