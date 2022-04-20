//-----------------------------------------------------------------//
// apply flattening                                                //
//-----------------------------------------------------------------//
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "cctk.h"

#define AXISYM 4

#define ONESIXTH 0.1666666666666666666666666666666666666666666666

double max_ppmflatten(double a,double b) {
  if(a>b) return a;
  return b;
}

extern "C" void CCTK_FCALL CCTK_FNAME(ppm_f_cpp)   
  (const cGH **cctkGH,int *ext,int *m,double *ftilde, double *f, double *P);

extern "C" void ppm_f_cpp(const cGH *cctkGH,int *ext,int m,double *ftilde,double *f, double *P) {
  // 
  int imin=0,jmin=0,kmin=0;
  int imax=ext[0];
  int jmax=ext[1];
  int kmax=ext[2];

  //Fill all points except for upper boundaries in m-direction.
  //It should not be necessary to fill lower and upper boundaries, since the values will be overwritten at the next sync call.
  // Besides, the values at the boundaries aren't correct.
  if(m==1) { imin++; imax--; }
  if(m==2) { jmin++; jmax--; }
  if(m==3) { kmin++; kmax--; }

#pragma omp parallel for
    for(int k=kmin;k<kmax;k++) for(int j=jmin;j<jmax;j++) for(int i=imin;i<imax;i++) {
          int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
          int indexm1,indexp1,indexps;
          if(m==1) {
            indexm1 = CCTK_GFINDEX3D(cctkGH,i-1,j,k);
            indexp1 = CCTK_GFINDEX3D(cctkGH,i+1,j,k);
            // Calculate s
            int s=-1;
            if ((P[indexp1] - P[indexm1]) <= 0.0) s=1;
            indexps = CCTK_GFINDEX3D(cctkGH,i+s,j,k);
          } else if(m==2) {
            indexm1 = CCTK_GFINDEX3D(cctkGH,i,j-1,k);
            indexp1 = CCTK_GFINDEX3D(cctkGH,i,j+1,k);
            // Calculate s
            int s=-1;
            if ((P[indexp1] - P[indexm1]) <= 0.0) s=1;
            indexps = CCTK_GFINDEX3D(cctkGH,i,j+s,k);
          } else {
            indexm1 = CCTK_GFINDEX3D(cctkGH,i,j,k-1);
            indexp1 = CCTK_GFINDEX3D(cctkGH,i,j,k+1);
            // Calculate s
            int s=-1;
            if ((P[indexp1] - P[indexm1]) <= 0.0) s=1;
            indexps = CCTK_GFINDEX3D(cctkGH,i,j,k+s);
          }
          // Calculate f
          f[index] = max_ppmflatten(ftilde[index],ftilde[indexps]);
        }
}


extern "C" void CCTK_FCALL CCTK_FNAME(ppm_flatten_1var_cpp)   
  (const cGH **cctkGH,int *ext, int &m, double *f, double *u, double *ur, double *ul);

extern "C" void ppm_flatten_1var_cpp(const cGH *cctkGH,int *ext, int &m, double *f,
                          double *u, double *ur, double *ul) {
  //
  int imin=0,jmin=0,kmin=0;
  int imax=ext[0];
  int jmax=ext[1];
  int kmax=ext[2];

  //Fill all points except for upper boundaries in m-direction.
  //It should not be necessary to fill lower and upper boundaries, since the values will be overwritten at the next sync call.
  // Besides, the values at the boundaries aren't correct.
  if(m==1) { imin++; imax--; }
  if(m==2) { jmin++; jmax--; }
  if(m==3) { kmin++; kmax--; }
#pragma omp parallel for
  for(int k=kmin;k<kmax;k++) for(int j=jmin;j<jmax;j++) for(int i=imin;i<imax;i++) 
    {    
     int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
     double fi = f[index];
     double omf = 1.0-fi;
     ur[index]   = u[index]*fi + ur[index]*omf;
     ul[index]   = u[index]*fi + ul[index]*omf;
    }
}

extern "C" void CCTK_FCALL CCTK_FNAME(ppm_f_cpp)
  (const cGH **cctkGH,int *ext,int *m,double *ftilde, double *f, double *P)
{
  ppm_f_cpp(*cctkGH,ext,*m,ftilde,f,P);
}

extern "C" void CCTK_FCALL CCTK_FNAME(ppm_flatten_1var_cpp)
  (const cGH **cctkGH,int *ext, int &m, double *f, double *u, double *ur, double *ul)
{
  ppm_flatten_1var_cpp(*cctkGH,ext,m,f,u,ur,ul);
}
