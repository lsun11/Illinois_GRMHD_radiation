//--------------------------------------------------------------------
// Add the extrinsic curvature terms to tau_rhs 
//--------------------------------------------------------------------

#define THIRD 0.333333333333333333333333333
#define SQR(x) ((x) * (x))

#include "math.h"
#include "cctk.h"

extern "C" void CCTK_FCALL mhd_tau_curvature_cpp_
  (const cGH **cctkGH,int *cctk_lsh, int *nghostzones, int *Symmetry,
   double *tau_rhs,double *rho_star,double *h,double *P,
   double *sbt,double *sbx,double *sby,double *sbz,
   double *u0,double *vx,double *vy,double *vz,
   double *alpha,double *betax,double *betay,double *betaz,
   double *phi,double *trK,
   double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
   double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,
   double *Axx,double *Axy,double *Axz,double *Ayy,double *Ayz,double *Azz,
   double *dX,double *dY,double *dZ);

extern "C" void mhd_tau_curvature_cpp(const cGH *cctkGH,int *cctk_lsh, int *nghostzones, int Symmetry,
				      double *tau_rhs,double *rho_star,double *h,double *P,
				      double *sbt,double *sbx,double *sby,double *sbz,
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
    
    // Divide b^a by alpha sqrt(4 pi)
    double sbtL = sbt[index]*alphaLm1*sqrtmfourpi;
    double sbxL = sbx[index]*alphaLm1*sqrtmfourpi;
    double sbyL = sby[index]*alphaLm1*sqrtmfourpi;
    double sbzL = sbz[index]*alphaLm1*sqrtmfourpi;

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
 
    double rho_starL = rho_star[index];
    double u0L = u0[index];
    double hL = h[index];
    double PL = P[index];

    double vxL = vx[index];
    double vyL = vy[index];
    double vzL = vz[index];

    // Compute b^2 
    double b2 = -SQR(alphaL*sbtL) + Psi4*( gxxL*SQR(sbxL+betaxL*sbtL) + 
					   2.0*gxyL*(sbxL+betaxL*sbtL)*(sbyL+betayL*sbtL) + 
					   2.0*gxzL*(sbxL+betaxL*sbtL)*(sbzL+betazL*sbtL) + 
					   gyyL*SQR(sbyL+betayL*sbtL) + 
					   2.0*gyzL*(sbyL+betayL*sbtL)*(sbzL+betazL*sbtL) + 
					   gzzL*SQR(sbzL+betazL*sbtL) );

    // Now add extrinsic curvature terms to tau_rhs
    double tau_rhsL = tau_rhs[index] + (rho_starL*hL*u0L + 
					sqrtg*b2*SQR(u0L))*Psi4* ( 
								  (AxxL+THIRD*gxxL*trKL)*SQR(vxL+betaxL) + 
								  2.0*(AxyL+THIRD*gxyL*trKL)*(vxL+betaxL)*(vyL+betayL) + 
								  2.0*(AxzL+THIRD*gxzL*trKL)*(vxL+betaxL)*(vzL+betazL) + 
								  (AyyL+THIRD*gyyL*trKL)*SQR(vyL+betayL) + 
								  2.0*(AyzL+THIRD*gyzL*trKL)*(vyL+betayL)*(vzL+betazL) + 
								  (AzzL+THIRD*gzzL*trKL)*SQR(vzL+betazL) );

    tau_rhsL += - sqrtg*Psi4* ( 
			       (AxxL+THIRD*gxxL*trKL)*SQR(sbxL+sbtL*betaxL) + 
			       2.0*(AxyL+THIRD*gxyL*trKL)*(sbxL+sbtL*betaxL)*(sbyL+sbtL*betayL) + 
			       2.0*(AxzL+THIRD*gxzL*trKL)*(sbxL+sbtL*betaxL)*(sbzL+sbtL*betazL) + 
			       (AyyL+THIRD*gyyL*trKL)*SQR(sbyL+sbtL*betayL) + 
			       2.0*(AyzL+THIRD*gyzL*trKL)*(sbyL+sbtL*betayL)*(sbzL+sbtL*betazL) + 
			       (AzzL+THIRD*gzzL*trKL)*SQR(sbzL+sbtL*betazL) );

    tau_rhs[index] = tau_rhsL + sqrtg*(PL+0.5*b2)*trKL;

  }
}

extern "C" void CCTK_FCALL mhd_tau_curvature_cpp_
  (const cGH **cctkGH,int *cctk_lsh, int *nghostzones, int *Symmetry,
   double *tau_rhs,double *rho_star,double *h,double *P,
   double *sbt,double *sbx,double *sby,double *sbz,
   double *u0,double *vx,double *vy,double *vz,
   double *alpha,double *betax,double *betay,double *betaz,
   double *phi,double *trK,
   double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
   double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,
   double *Axx,double *Axy,double *Axz,double *Ayy,double *Ayz,double *Azz,
   double *dX,double *dY,double *dZ)
{
  mhd_tau_curvature_cpp(*cctkGH,cctk_lsh, nghostzones, *Symmetry,
			tau_rhs,rho_star,h,P,
			sbt,sbx,sby,sbz,
			u0,vx,vy,vz,
			alpha,betax,betay,betaz,
			phi,trK,
			gxx,gxy,gxz,gyy,gyz,gzz,
			gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,
			Axx,Axy,Axz,Ayy,Ayz,Azz,
			*dX,*dY,*dZ);
}
