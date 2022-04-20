//-----------------------------------------------------------------//
// ensure monotonicity of the interpolating polynomial (remove the PPM plus stuff)
//-----------------------------------------------------------------//
//

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "cctk.h"

#define AXISYM 4
#define ONESIXTH 0.166666666666666666666666666666666666666666666
#define SQR(x) ((x) * (x))

extern "C" void CCTK_FCALL ppm_monotonize_v2_cpp_
  (const cGH **cctkGH,int *ext,double *a,double *ar,double *al,int &Symmetry);

extern "C" void ppm_monotonize_v2_cpp (const cGH *cctkGH,int *ext,double *a,double *ar,double *al,int &Symmetry) {
  int imin=0,jmin=0,kmin=0;
  int imax=ext[0];
  int jmax=ext[1];
  int kmax=ext[2];
  
  if (Symmetry==AXISYM) { 
    jmin = 1;
    jmax = 1;
  }
  
  //Fill all points.
#pragma omp parallel for
  for(int k=kmin;k<kmax;k++) for(int j=jmin;j<jmax;j++) for(int i=imin;i<imax;i++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

    double arijk=ar[index];
    double alijk=al[index];
    double aijk = a[index];
    double daijk = arijk - alijk;
    double maijk = 0.5*(arijk+alijk);

    if ( (arijk-aijk)*(aijk-alijk) <= 0.0) { 
      ar[index] = aijk;
      al[index] = aijk;
    } else if ( daijk*(aijk-maijk) > ONESIXTH*SQR(daijk)) { 
      al[index] = 3.0*aijk - 2.0*arijk;
    } else if ( daijk*(aijk-maijk) < -ONESIXTH*SQR(daijk)) {
      ar[index] = 3.0*aijk - 2.0*alijk;
    }
  }
}

extern "C" void CCTK_FCALL ppm_monotonize_v2_cpp_
  (const cGH **cctkGH,int *ext,double *a,double *ar,double *al,int &Symmetry)
{
  ppm_monotonize_v2_cpp(*cctkGH,ext,a,ar,al,Symmetry);
}
