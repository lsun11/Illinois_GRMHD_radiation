//-----------------------------------------------------------------------------
// Advect Z_i
//-----------------------------------------------------------------------------

#include "math.h"
#include "cctk.h"
#include "stdio.h"

#define SQR(x) ((x) * (x))

extern "C" void CCTK_FCALL CCTK_FNAME(advect_z_cpp)
  (int *flux_direction, const cGH **cctkGH,int *cctk_lsh, int *nghostzones, int *Symmetry,
   double *X, double *Z_y,double *vy, double *P,double *alpha,double *phi, 
   double *betax,double *betay,double *betaz, 
   double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz, 
   double *sbt,double *sbx,double *sby,double *sbz, 
   double *st_x_rhs, double *st_y_rhs, double *st_z_rhs, 
   double *fmx, double *fmy, double *fmz, 
   double *dX,double *dY,double *dZ);

extern "C" void advect_z_cpp(int flux_direction, const cGH *cctkGH,int *cctk_lsh, int *nghostzones, int Symmetry,
			     double *X, double *Z_y,double *vy, double *P,double *alpha,double *phi, 
			     double *betax,double *betay,double *betaz, 
			     double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz, 
			     double *sbt,double *sbx,double *sby,double *sbz, 
			     double *st_x_rhs, double *st_y_rhs, double *st_z_rhs, 
			     double *fmx, double *fmy, double *fmz, 
			     double dX,double dY,double dZ) {
  double sfpi = sqrt(4.0*M_PI);

  int AXISYM = 4;

  double f1o4pi = 1.0/(4.0*M_PI);

  /* Set up variables used in the grid loop for the physical grid points */

  int istart = 1;
  int jstart = 1;
  int kstart = 1;
  int iend = cctk_lsh[0] - 1;
  int jend = cctk_lsh[1] - 1;
  int kend = cctk_lsh[2] - 1;

  double dxi = 1.0/dX;
  double dyi = 1.0/dY;
  double dzi = 1.0/dZ;

  if(Symmetry==4) {
    jstart = 0;
    jend = cctk_lsh[1];
    jstart++;
    jend--;
  }
  
  if (flux_direction==1) {
    if (Symmetry==AXISYM) {
#pragma omp parallel for
      for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
	int index   = CCTK_GFINDEX3D(cctkGH,i,  j,k);  
	int indexp1 = CCTK_GFINDEX3D(cctkGH,i+1,j,k);  
	double al = 1.0 + alpha[index];
	double betaxL = betax[index];
	double betayL = betay[index];
	double betazL = betaz[index];
	double Psi = exp(phi[index]);
	double Psi4 = Psi*Psi*Psi*Psi;
	double Psi6 = Psi4*Psi*Psi;
	double gxxL = gxx[index];
	double gxyL = gxy[index];
	double gxzL = gxz[index];
	double gyyL = gyy[index];
	double gyzL = gyz[index];
	double gzzL = gzz[index];
  
	double sbtL = sbt[index];
	double sbxL = sbx[index];
	double sbyL = sby[index];
	double sbzL = sbz[index];

	double xLi = 1.0/X[index];
	double sb_y = Psi4*( gxyL*(sbtL*betaxL + sbxL) + gyyL*(sbtL*betayL + sbyL) + 
			     gyzL*(sbtL*betazL + sbzL) );
	double b2 = - SQR(al*sbtL) + Psi4*( gxxL*SQR(sbxL + betaxL*sbtL) + 
					    2.0*gxyL*(sbxL + betaxL*sbtL)*(sbyL + betayL*sbtL) + 
					    2.0*gxzL*(sbxL + betaxL*sbtL)*(sbzL + betazL*sbtL) + 
					    gyyL*SQR(sbyL + betayL*sbtL) + 
					    2.0*gyzL*(sbyL + betayL*sbtL)*(sbzL + betazL*sbtL) + 
					    gzzL*SQR(sbzL + betazL*sbtL) );

	st_x_rhs[index] += ( (fmx[index]-fmx[indexp1])*dxi + Z_y[index]*vy[index] + 
			     al*Psi6*(P[index] + f1o4pi/SQR(al) *(0.5*b2 + sb_y * (sbtL*vy[index] - sbyL) ) ))*xLi;
	st_y_rhs[index] += (fmy[index]-fmy[indexp1])*dxi*xLi*xLi;
	st_z_rhs[index] += (fmz[index]-fmz[indexp1])*dxi*xLi;
      }
    } else {
#pragma omp parallel for
      for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
	int index   = CCTK_GFINDEX3D(cctkGH,i,  j,k);  
	int indexp1 = CCTK_GFINDEX3D(cctkGH,i+1,j,k);  
	st_x_rhs[index] += (fmx[index]-fmx[indexp1])*dxi;
	st_y_rhs[index] += (fmy[index]-fmy[indexp1])*dxi;
	st_z_rhs[index] += (fmz[index]-fmz[indexp1])*dxi;
      }
    }
  } else if (flux_direction==2) {
    if (Symmetry != AXISYM) {
#pragma omp parallel for
      for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
	int index   = CCTK_GFINDEX3D(cctkGH,i,j,  k);  
	int indexp1 = CCTK_GFINDEX3D(cctkGH,i,j+1,k);  
	st_x_rhs[index] += (fmx[index]-fmx[indexp1])*dyi;
	st_y_rhs[index] += (fmy[index]-fmy[indexp1])*dyi;
	st_z_rhs[index] += (fmz[index]-fmz[indexp1])*dyi;
      }
    }
  } else if (flux_direction==3) {
#pragma omp parallel for
    for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
      int index   = CCTK_GFINDEX3D(cctkGH,i,j,k  );
      int indexp1 = CCTK_GFINDEX3D(cctkGH,i,j,k+1);
      st_x_rhs[index] += (fmx[index]-fmx[indexp1])*dzi;
      st_y_rhs[index] += (fmy[index]-fmy[indexp1])*dzi;
      st_z_rhs[index] += (fmz[index]-fmz[indexp1])*dzi;
    }
  }
}
extern "C" void CCTK_FCALL CCTK_FNAME(advect_z_cpp)
  (int *flux_direction, const cGH **cctkGH,int *cctk_lsh, int *nghostzones, int *Symmetry,
   double *X, double *Z_y,double *vy, double *P,double *alpha,double *phi, 
   double *betax,double *betay,double *betaz, 
   double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz, 
   double *sbt,double *sbx,double *sby,double *sbz, 
   double *st_x_rhs, double *st_y_rhs, double *st_z_rhs, 
   double *fmx, double *fmy, double *fmz, 
   double *dX,double *dY,double *dZ)
{
  advect_z_cpp(*flux_direction,*cctkGH,cctk_lsh, nghostzones, *Symmetry,
	       X, Z_y,vy, P,alpha,phi, 
	       betax,betay,betaz, 
	       gxx,gxy,gxz,gyy,gyz,gzz, 
	       sbt,sbx,sby,sbz, 
	       st_x_rhs, st_y_rhs, st_z_rhs, 
	       fmx, fmy, fmz, 
	       *dX,*dY,*dZ);
}
