//-----------------------------------------------------------------//
// calculate \tilde{f}_j                                           //
//-----------------------------------------------------------------//
//
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "cctk.h"

#define AXISYM 4

#define OMEGA1   0.75
#define OMEGA2  10.0
#define EPSILON2 0.33
// The following parameter choices are Font-like.
//#define OMEGA1   0.52
//#define OMEGA2  10.0
//#define EPSILON2 0.5

double min_ppmftilde(double x,double y);
double max_ppmftilde(double x,double y);

extern "C" void CCTK_FCALL ppm_ftilde_v2_cpp_ 
  (const cGH **cctkGH,int *ext,int *nghostzones,double *ftilde,
   double *P,double *vx,double *vy,double *vz,int *symP,int *symvx,int *symvy,int *symvz,  
   int &m,int &Symmetry,int &cell_centering_enabled);

extern "C" void ppm_ftilde_v2_cpp(const cGH *cctkGH,int *ext,int *nghostzones,double *ftilde,
			      double *P,double *vx,double *vy,double *vz,int *symP,int *symvx,int *symvy,int *symvz,  
			      int &m,int &Symmetry,int &cell_centering_enabled) {
  int imin=0,jmin=0,kmin=0;
  int imax=ext[0];
  int jmax=ext[1];
  int kmax=ext[2];
  
  if (Symmetry==AXISYM) { 
    jmin = 1;
    jmax = 1;
  }

  for(int k=kmin;k<kmax;k++) for(int j=jmin;j<jmax;j++) for(int i=imin;i<imax;i++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    ftilde[index] = 1.0;
  }

  //Fill all points except for lower and upper boundaries in m-direction.
  //It should not be necessary to fill lower and upper boundaries, since the values will be overwritten at the next sync call.
  // Now compute ftilde[index] for the rest of i
  if(m==1) { imin+=2; imax-=2; }
  if(m==2) { jmin+=2; jmax-=2; }
  if(m==3) { kmin+=2; kmax-=2; }
  //make more efficient by adding openmp support
#pragma omp parallel for
  for(int k=kmin;k<kmax;k++) for(int j=jmin;j<jmax;j++) for(int i=imin;i<imax;i++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    int indexm1,indexm2,indexp1,indexp2;
    if(m==1) {
      indexm1 = CCTK_GFINDEX3D(cctkGH,i-1,j,k);
      indexp1 = CCTK_GFINDEX3D(cctkGH,i+1,j,k);
      indexm2 = CCTK_GFINDEX3D(cctkGH,i-2,j,k);
      indexp2 = CCTK_GFINDEX3D(cctkGH,i+2,j,k);
    } else if(m==2) {
      indexm1 = CCTK_GFINDEX3D(cctkGH,i,j-1,k);
      indexp1 = CCTK_GFINDEX3D(cctkGH,i,j+1,k);
      indexm2 = CCTK_GFINDEX3D(cctkGH,i,j-2,k);
      indexp2 = CCTK_GFINDEX3D(cctkGH,i,j+2,k);
    } else {
      indexm1 = CCTK_GFINDEX3D(cctkGH,i,j,k-1);
      indexp1 = CCTK_GFINDEX3D(cctkGH,i,j,k+1);
      indexm2 = CCTK_GFINDEX3D(cctkGH,i,j,k-2);
      indexp2 = CCTK_GFINDEX3D(cctkGH,i,j,k+2);
    }

    double dP1 = P[indexp1]-P[indexm1];
    double dP2 = P[indexp2]-P[indexm2];
    double dP1odP2=1.0;
    if (dP2 != 0.0) dP1odP2 = dP1/dP2;

    double q1 = (dP1odP2-OMEGA1)*OMEGA2;
    double q2 = dP1/min_ppmftilde(P[indexp1],P[indexm1]);
    double w=0.0;
    if(m==1) {
      if (q2 > EPSILON2 &&  q2*(vx[indexm1]-vx[indexp1]) > 0.0) w = 1.0; 
    } if(m==2) {
      if (q2 > EPSILON2 &&  q2*(vy[indexm1]-vy[indexp1]) > 0.0) w = 1.0; 
    } else {
      if (q2 > EPSILON2 &&  q2*(vz[indexm1]-vz[indexp1]) > 0.0) w = 1.0; 
    }
    ftilde[index] = min_ppmftilde(1.0, w*max_ppmftilde(0.0,q1) );
  }

  /*
    double dP1 = P[indexp1]-P[indexm1];
    double dP2 = P[indexp2]-P[indexm2];
    double dP1odP2 = 1.0;
    if (dP2 != 0.0) dP1odP2 = dP1/dP2;

    double q1 = (dP1odP2-OMEGA1)*OMEGA2;
    double q2 = dP1/min_ppmftilde(P[indexp1],P[indexm1]);

    double w=0.0;
    if (q2 > EPSILON2 && q2*(vx[indexm1]-vx[indexp1]) > 0.0) {
      w = 1.0;
    }

    ftilde[index] = min_ppmftilde(1.0, w*max_ppmftilde(0.0,q1) );
  }
  */
}

double min_ppmftilde(double x,double y) {
  if(x<y) return x;
  return y;
}

double max_ppmftilde(double x,double y) {
  if(x>y) return x;
  return y;
}


extern "C" void CCTK_FCALL ppm_ftilde_v2_cpp_ 
  (const cGH **cctkGH,int *ext,int *nghostzones,double *ftilde,
   double *P,double *vx,double *vy,double *vz,int *symP,int *symvx,int *symvy,int *symvz,  
   int &m,int &Symmetry,int &cell_centering_enabled)
{
  ppm_ftilde_v2_cpp(*cctkGH,ext,nghostzones,ftilde,
		P,vx,vy,vz,symP,symvx,symvy,symvz,  
		m,Symmetry,cell_centering_enabled);
}
