//--------------------------------------------------------------------
// Add the extrinsic curvature terms to tau_rad_rhs 
//--------------------------------------------------------------------

#define THIRD 0.333333333333333333333333333
#define SQR(x) ((x) * (x))

#include "math.h"
#include "cctk.h"

extern "C" void CCTK_FCALL CCTK_FNAME(tau_rad_curvature_cpp)
  (const cGH **cctkGH,int *cctk_lsh, int *nghostzones, int *Symmetry,
   double *tau_rad_rhs, double *E_rad,double *F_rad0,double *F_radx,
   double *F_rady,double *F_radz,double *P_rad,
   double *u0,double *vx,double *vy,double *vz,
   double *alpha,double *betax,double *betay,double *betaz,
   double *phi,double *trK,
   double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
   double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,
   double *Axx,double *Axy,double *Axz,double *Ayy,double *Ayz,double *Azz,
   double *dX,double *dY,double *dZ);

extern "C" void tau_rad_curvature_cpp(const cGH *cctkGH,int *cctk_lsh, int *nghostzones, int Symmetry,
				      double *tau_rad_rhs, double *E_rad,double *F_rad0,double *F_radx,
				      double *F_rady,double *F_radz,double *P_rad,
				      double *u0,double *vx,double *vy,double *vz,
				      double *alpha,double *betax,double *betay,double *betaz,
				      double *phi,double *trK,
				      double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
				      double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,
				      double *Axx,double *Axy,double *Axz,double *Ayy,double *Ayz,double *Azz,
				      double dX,double dY,double dZ) {

  double sqrtmfourpi = 1.0/sqrt(4.0*M_PI);

  int AXISYM = 4;

  /* Set up variables used in the grid loop for the physical grid points */
  int istart = 0;
  int jstart = 0;
  int kstart = 0;
  int iend = cctk_lsh[0];
  int jend = cctk_lsh[1];
  int kend = cctk_lsh[2];
  /*
    int istart = nghostzones[0];
    int jstart = nghostzones[1];
    int kstart = nghostzones[2];
    int iend = cctk_lsh[0] - nghostzones[0];
    int jend = cctk_lsh[1] - nghostzones[1];
    int kend = cctk_lsh[2] - nghostzones[2];
  */

  /*
  //Following lines needed since nghostzones[0] = ORDER, and 
  //   not ORDER-1 in axisymmetry 
  //   (so that rotation can be done on multiprocessor runs)
  if(Symmetry==4) {
    istart--;
    iend++;
  }
  */

  if(Symmetry==4) {
    jstart = 1; jend = 2;
  }

#pragma omp parallel for
  for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    double alphaL = 1.0+alpha[index];
    double alphaLm1 = 1.0/alphaL;
    double Psi4 = exp(4.0*phi[index]);
    double sqrtg = alphaL *exp(6.0*phi[index]); // Psi4 * sqrt(Psi4); // == alphaL*exp(6.0*phi)
    
    double betaxL = betax[index];
    double betayL = betay[index];
    double betazL = betaz[index];

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

    double trKL = trK[index];
    double u0L  = u0[index];

    double E_radL  = E_rad[index];
    double F_rad0L = F_rad0[index];
    double F_radxL = F_radx[index];
    double F_radyL = F_rady[index];
    double F_radzL = F_radz[index];
    double P_radL  = P_rad[index];

    double vxL = vx[index];
    double vyL = vy[index];
    double vzL = vz[index];


    // Now add extrinsic curvature terms to tau_rad_rhs
    double tau_rad_rhsL = tau_rad_rhs[index] + (sqrtg*(E_radL + P_radL)*SQR(u0L))*Psi4* ( 
								  (AxxL+THIRD*gxxL*trKL)*SQR(vxL+betaxL) + 
								  2.0*(AxyL+THIRD*gxyL*trKL)*(vxL+betaxL)*(vyL+betayL) + 
								  2.0*(AxzL+THIRD*gxzL*trKL)*(vxL+betaxL)*(vzL+betazL) + 
								  (AyyL+THIRD*gyyL*trKL)*SQR(vyL+betayL) + 
								  2.0*(AyzL+THIRD*gyzL*trKL)*(vyL+betayL)*(vzL+betazL) + 
								  (AzzL+THIRD*gzzL*trKL)*SQR(vzL+betazL) );

    tau_rad_rhsL += 2.0* (sqrtg*u0L*F_rad0L)*Psi4* ( 
					   (AxxL+THIRD*gxxL*trKL)*(SQR(betaxL)+betaxL*vxL) + 
					   (AyyL+THIRD*gyyL*trKL)*(SQR(betayL)+betayL*vyL) +
					   (AzzL+THIRD*gzzL*trKL)*(SQR(betazL)+betazL*vzL) +
					   (AxyL+THIRD*gxyL*trKL)*(2.0*betaxL*betayL + vxL*betayL + vyL*betaxL) + 
					   (AxyL+THIRD*gxzL*trKL)*(2.0*betaxL*betazL + vxL*betazL + vzL*betayL) +
					   (AxyL+THIRD*gyzL*trKL)*(2.0*betayL*betazL + vyL*betazL + vzL*betayL) );

    tau_rad_rhs[index] = tau_rad_rhsL + sqrtg*(P_radL)*trKL;

  }
}

extern "C" void CCTK_FCALL CCTK_FNAME(tau_rad_curvature_cpp)
  (const cGH **cctkGH,int *cctk_lsh, int *nghostzones, int *Symmetry,
   double *tau_rad_rhs, double *E_rad,double *F_rad0,double *F_radx,
   double *F_rady,double *F_radz,double *P_rad,
   double *u0,double *vx,double *vy,double *vz,
   double *alpha,double *betax,double *betay,double *betaz,
   double *phi,double *trK,
   double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
   double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,
   double *Axx,double *Axy,double *Axz,double *Ayy,double *Ayz,double *Azz,
   double *dX,double *dY,double *dZ)
{
  tau_rad_curvature_cpp(*cctkGH,cctk_lsh, nghostzones, *Symmetry,
			tau_rad_rhs, E_rad, 
			F_rad0, F_radx,F_rady,F_radz,P_rad,
			u0,vx,vy,vz,
			alpha,betax,betay,betaz,
			phi,trK,
			gxx,gxy,gxz,gyy,gyz,gzz,
			gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,
			Axx,Axy,Axz,Ayy,Ayz,Azz,
			*dX,*dY,*dZ);
}
