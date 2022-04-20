//----------------------------------------------------------------------!
// Calculate a_{j+1/2}, set a_R = a_{j+1/2}, a_L  = a_{j-1/2} = a_R(j-1)
//----------------------------------------------------------------------!

#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "cctk.h"

#define AXISYM 4

#define ONESIXTH 0.1666666666666666666666666666666666666666666666

extern "C" void CCTK_FCALL CCTK_FNAME(ppm_find_face_vals_v2_cpp)
  (const cGH **cctkGH,int *ext,int *nghostzones,double *a,double *delta_a,double *ar,double *al,int &m,int *sym,  
   int &Symmetry,int &cell_centering_enabled);

extern "C" void ppm_find_face_vals_v2_cpp(const cGH *cctkGH,int *ext,int *nghostzones,double *a,double *delta_a,double *ar,double *al,int &m,int *sym,
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
  
  //Fill all points except for upper boundaries in m-direction.
  //Below is an old comment.  I'm not sure if it is correct:
  ///     It should not be necessary to fill lower and upper boundaries, since the values will be overwritten at the next sync call.  Besides, the values at the boundaries aren't correct.

  if(m==1) { imax--; }
  if(m==2) { jmax--; }
  if(m==3) { kmax--; }

#pragma omp parallel for
  for(int k=kmin;k<kmax;k++) for(int j=jmin;j<jmax;j++) for(int i=imin;i<imax;i++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    int indexm1,indexp1;
    if(m==1) {
      indexp1 = CCTK_GFINDEX3D(cctkGH,i+1,j,k);
    } else if(m==2) {
      indexp1 = CCTK_GFINDEX3D(cctkGH,i,j+1,k);
    } else {
      indexp1 = CCTK_GFINDEX3D(cctkGH,i,j,k+1);
    }
    ar[index] = a[index] + 0.5*(a[indexp1]-a[index]) + 
      ONESIXTH*(delta_a[index]-delta_a[indexp1]);
  }

  // Fill in the upper boundary points 
  if (m==1) { 
#pragma omp parallel for
     for(int k=kmin;k<kmax;k++) for(int j=jmin;j<jmax;j++) {
	int index = CCTK_GFINDEX3D(cctkGH,imax,j,k);
	ar[index] = a[index] + 0.5*delta_a[index];
     }
  } else if (m==2) {
#pragma omp parallel for
    for(int k=kmin;k<kmax;k++) for(int i=imin;i<imax;i++) {
        int index = CCTK_GFINDEX3D(cctkGH,i,jmax,k);
        ar[index] = a[index] + 0.5*delta_a[index];
    }
  } else {
#pragma omp parallel for
    for(int j=jmin;j<jmax;j++) for(int i=imin;i<imax;i++) {
        int index = CCTK_GFINDEX3D(cctkGH,i,j,kmax);
        ar[index] = a[index] + 0.5*delta_a[index];
    }
  }


  imax=ext[0];
  jmax=ext[1];
  kmax=ext[2];

  if(m==1) { imin++; }
  if(m==2) { jmin++; }
  if(m==3) { kmin++; }
#pragma omp parallel for
  for(int k=kmin;k<kmax;k++) for(int j=jmin;j<jmax;j++) for(int i=imin;i<imax;i++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    int indexm1;
    if(m==1) {
      indexm1 = CCTK_GFINDEX3D(cctkGH,i-1,j,k);
    } else if(m==2) {
      indexm1 = CCTK_GFINDEX3D(cctkGH,i,j-1,k);
    } else {
      indexm1 = CCTK_GFINDEX3D(cctkGH,i,j,k-1);
    }
    al[index] = ar[indexm1];
  }

  // Fill in the lower boundary points 
  if (m==1) {
     imin--;
#pragma omp parallel for
     for(int k=kmin;k<kmax;k++) for(int j=jmin;j<jmax;j++) {
	int index = CCTK_GFINDEX3D(cctkGH,imin,j,k);
        al[index] = a[index] - 0.5*delta_a[index];
     }
  } else if (m==2) {
     jmin--;
#pragma omp parallel for
     for(int k=kmin;k<kmax;k++) for(int i=imin;i<imax;i++) {
        int index = CCTK_GFINDEX3D(cctkGH,i,jmin,k);
        al[index] = a[index] - 0.5*delta_a[index];
     }
  } else {
     kmin--;
#pragma omp parallel for
     for(int j=jmin;j<jmax;j++) for(int i=imin;i<imax;i++) {
        int index = CCTK_GFINDEX3D(cctkGH,i,j,kmin);
        al[index] = a[index] - 0.5*delta_a[index];
     }
  }
}


extern "C" void CCTK_FCALL CCTK_FNAME(ppm_find_face_vals_v2_cpp)
  (const cGH **cctkGH,int *ext,int *nghostzones,double *a,double *delta_a,double *ar,double *al,int &m,int *sym,  
   int &Symmetry,int &cell_centering_enabled) 
{
  ppm_find_face_vals_v2_cpp(*cctkGH,ext,nghostzones,a,delta_a,ar,al,m,sym,  
			Symmetry,cell_centering_enabled);
}
