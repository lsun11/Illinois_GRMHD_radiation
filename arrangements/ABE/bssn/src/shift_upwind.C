
//----------------------------------------------------------------------
//
// Compute advection terms in right hand sides of field equations
// See Shibata in Prog. Theor. Phys., 101, 1199 (1999) (also gr-qc/9905058)
// Note:  centered advective term has already been included in ()_rhs by 
// compute_rhs.f90.
//
//
// I optimized the crap out of this code.  -Zach
//
//-----------------------------------------------------------------------------

#define EPSILON 1.0e-8
#define F2o3 0.666666666666666666666666666666666666
#define F1o3 0.333333333333333333333333333333333333

/* Define macros used in calculations */
#define INITVALUE  (42)
#define INV(x) ((1) / (x))
#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * (x) * (x))
#define QAD(x) ((x) * (x) * (x) * (x))

#include <stdio.h>
#include "cctk.h"
#include <math.h>
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

extern "C" void CCTK_FCALL CCTK_FNAME(shift_upwind2)
  (const cGH **cctkGH,int *cctk_lsh,int *nghostzones,double *dT, double *dx, double *dy, double *dz,
   double *gxx, double *gxy, double *gxz, double *gyy, double *gyz, double *gzz, 
   double *gxx_rhs, double *gxy_rhs, double *gxz_rhs, double *gyy_rhs, double *gyz_rhs, double *gzz_rhs, 
   double *Axx, double *Axy, double *Axz, double *Ayy, double *Ayz, double *Azz, 
   double *Axx_rhs, double *Axy_rhs, double *Axz_rhs, double *Ayy_rhs, double *Ayz_rhs, double *Azz_rhs, 
   double *phi, double *phi_rhs, 
   double *trK, double *trK_rhs, 
   double *Gamx, double *Gamy, double *Gamz,
   double *Gamx_rhs, double *Gamy_rhs, double *Gamz_rhs, 
   double *betax, double *betay, double *betaz, int *Symmetry);


void shift_upwind2(const cGH *cctkGH, int *cctk_lsh,int *nghostzones,double dT, double dx, double dy, double dz, 
		   double *gxx, double *gxy, double *gxz, double *gyy, double *gyz, double *gzz, 
		   double *gxx_rhs, double *gxy_rhs, double *gxz_rhs, double *gyy_rhs, double *gyz_rhs, double *gzz_rhs, 
		   double *Axx, double *Axy, double *Axz, double *Ayy, double *Ayz, double *Azz, 
		   double *Axx_rhs, double *Axy_rhs, double *Axz_rhs, double *Ayy_rhs, double *Ayz_rhs, double *Azz_rhs, 
		   double *phi, double *phi_rhs, 
		   double *trK, double *trK_rhs, 
		   double *Gamx, double *Gamy, double *Gamz,
		   double *Gamx_rhs, double *Gamy_rhs, double *Gamz_rhs, 
		   double *betax, double *betay, double *betaz, int Symmetry) {
  /* Initialise finite differencing variables.  NEED THIS FOR GenericFD.h */
#include "../../GenFD_decl_set_varCPP.h"

  /* Set up variables used in the grid loop for the physical grid points */
  //WARNING: MUST BE CAREFUL WITH 4TH ORDER UPWIND STENCILS: They extend out 3 gridpoints!!
  int istart = nghostzones[0];
  int jstart = nghostzones[1];
  int kstart = nghostzones[2];
  int iend = cctk_lsh[0] - nghostzones[0];
  int jend = cctk_lsh[1] - nghostzones[1];
  int kend = cctk_lsh[2] - nghostzones[2];

  //Following lines needed since nghostzones[0] = ORDER, and 
  //   not ORDER-1 in axisymmetry 
  //   (so that rotation can be done on multiprocessor runs)
  if(Symmetry==4) {
    istart--;
    iend++;
  }

  //  printf("NGHOSTZONES = %d %d %d\n",nghostzones[0],nghostzones[1],nghostzones[2]);
  
#pragma omp parallel for
  for(int k=kstart;k<kend;k++)
    for(int j=jstart;j<jend;j++)
      for(int i=istart;i<iend;i++) {
	double fac1,fac2;
	double gxxs,gxys,gxzs,gyys,gyzs,gzzs;
	double Axxs,Axys,Axzs,Ayys,Ayzs,Azzs;
	double phis,trKs;
	double Gamxs,Gamys,Gamzs;

	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

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

	double betaxL = betax[index];
	double betayL = betay[index];
	double betazL = betaz[index];

	double betaxL2fac = betaxL*betaxL*dT;
	double betayL2fac = betayL*betayL*dT;
	double betazL2fac = betazL*betazL*dT;

	double absbetaxLfac = fabs(betaxL)*dx;
	double absbetayLfac = fabs(betayL)*dy;
	double absbetazLfac = fabs(betazL)*dz;

	double GamxL = Gamx[index];
	double GamyL = Gamy[index];
	double GamzL = Gamz[index];

	double phiL = phi[index];
	double trKL = trK[index];

	double gxx_rhsL = gxx_rhs[index];
	double gxy_rhsL = gxy_rhs[index];
	double gxz_rhsL = gxz_rhs[index];
	double gyy_rhsL = gyy_rhs[index];
	double gyz_rhsL = gyz_rhs[index];
	double gzz_rhsL = gzz_rhs[index];

	double Axx_rhsL = Axx_rhs[index];
	double Axy_rhsL = Axy_rhs[index];
	double Axz_rhsL = Axz_rhs[index];
	double Ayy_rhsL = Ayy_rhs[index];
	double Ayz_rhsL = Ayz_rhs[index];
	double Azz_rhsL = Azz_rhs[index];

	double phi_rhsL = phi_rhs[index];
	double trK_rhsL = trK_rhs[index];

	double Gamx_rhsL = Gamx_rhs[index];
	double Gamy_rhsL = Gamy_rhs[index];
	double Gamz_rhsL = Gamz_rhs[index];


	//Advect gxx :	
	// x-direction
	fac1 = (gxx[CCTK_GFINDEX3D(cctkGH,i+1,j,k)] - gxxL);
	fac2 = (gxxL - gxx[CCTK_GFINDEX3D(cctkGH,i-1,j,k)]);
	gxxs=0;
	if(fac1*fac2 >= 0.0) gxxs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	gxx_rhsL = gxx_rhsL  
          + ( (1.0 - gxxs)*absbetaxLfac + gxxs*betaxL2fac )*hdxi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order)!
	// y-direction
	fac1 = (gxx[CCTK_GFINDEX3D(cctkGH,i,j+1,k)] - gxxL);
	fac2 = (gxxL - gxx[CCTK_GFINDEX3D(cctkGH,i,j-1,k)]);
	gxxs=0;
	if(fac1*fac2 >= 0.0) gxxs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	gxx_rhsL = gxx_rhsL 
          + ( (1.0 - gxxs)*absbetayLfac + gxxs*betayL2fac )*hdyi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order)!
	// z-direction
	fac1 = (gxx[CCTK_GFINDEX3D(cctkGH,i,j,k+1)] - gxxL);
	fac2 = (gxxL - gxx[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]);
	gxxs=0;
	if(fac1*fac2 >= 0.0) gxxs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	gxx_rhs[index] = gxx_rhsL  
          + ( (1.0 - gxxs)*absbetazLfac + gxxs*betazL2fac)*hdzi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order)!

	//Advect gxy :	
	// x-direction
	fac1 = (gxy[CCTK_GFINDEX3D(cctkGH,i+1,j,k)] - gxyL);
	fac2 = (gxyL - gxy[CCTK_GFINDEX3D(cctkGH,i-1,j,k)]);
	gxys=0;
	if(fac1*fac2 >= 0.0) gxys = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	gxy_rhsL = gxy_rhsL  
          + ( (1.0 - gxys)*absbetaxLfac + gxys*betaxL2fac )*hdxi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order)! 
	// y-direction
	fac1 = (gxy[CCTK_GFINDEX3D(cctkGH,i,j+1,k)] - gxyL);
	fac2 = (gxyL - gxy[CCTK_GFINDEX3D(cctkGH,i,j-1,k)]);
	gxys=0;
	if(fac1*fac2 >= 0.0) gxys = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	gxy_rhsL = gxy_rhsL  
          + ( (1.0 - gxys)*absbetayLfac + gxys*betayL2fac )*hdyi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order)!
	// z-direction
	fac1 = (gxy[CCTK_GFINDEX3D(cctkGH,i,j,k+1)] - gxyL);
	fac2 = (gxyL - gxy[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]);
	gxys=0;
	if(fac1*fac2 >= 0.0) gxys = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	gxy_rhs[index] = gxy_rhsL  
          + ( (1.0 - gxys)*absbetazLfac + gxys*betazL2fac)*hdzi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order)!
	
	//Advect gxz :	
	// x-direction
	fac1 = (gxz[CCTK_GFINDEX3D(cctkGH,i+1,j,k)] - gxzL);
	fac2 = (gxzL - gxz[CCTK_GFINDEX3D(cctkGH,i-1,j,k)]);
	gxzs=0;
	if(fac1*fac2 >= 0.0) gxzs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	gxz_rhsL = gxz_rhsL  
          + ( (1.0 - gxzs)*absbetaxLfac + gxzs*betaxL2fac )*hdxi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order)!
	// y-direction
	fac1 = (gxz[CCTK_GFINDEX3D(cctkGH,i,j+1,k)] - gxzL);
	fac2 = (gxzL - gxz[CCTK_GFINDEX3D(cctkGH,i,j-1,k)]);
	gxzs=0;
	if(fac1*fac2 >= 0.0) gxzs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	gxz_rhsL = gxz_rhsL  
          + ( (1.0 - gxzs)*absbetayLfac + gxzs*betayL2fac )*hdyi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order)!
	// z-direction
	fac1 = (gxz[CCTK_GFINDEX3D(cctkGH,i,j,k+1)] - gxzL);
	fac2 = (gxzL - gxz[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]);
	gxzs=0;
	if(fac1*fac2 >= 0.0) gxzs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	gxz_rhs[index] = gxz_rhsL  
          + ( (1.0 - gxzs)*absbetazLfac + gxzs*betazL2fac)*hdzi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order)!
	
	//Advect gyy :	
	// x-direction
	fac1 = (gyy[CCTK_GFINDEX3D(cctkGH,i+1,j,k)] - gyyL);
	fac2 = (gyyL - gyy[CCTK_GFINDEX3D(cctkGH,i-1,j,k)]);
	gyys=0;
	if(fac1*fac2 >= 0.0) gyys = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	gyy_rhsL = gyy_rhsL  
          + ( (1.0 - gyys)*absbetaxLfac + gyys*betaxL2fac )*hdxi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order)!
	// y-direction
	fac1 = (gyy[CCTK_GFINDEX3D(cctkGH,i,j+1,k)] - gyyL);
	fac2 = (gyyL - gyy[CCTK_GFINDEX3D(cctkGH,i,j-1,k)]);
	gyys=0;
	if(fac1*fac2 >= 0.0) gyys = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	gyy_rhsL = gyy_rhsL  
          + ( (1.0 - gyys)*absbetayLfac + gyys*betayL2fac )*hdyi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order)!
	// z-direction
	fac1 = (gyy[CCTK_GFINDEX3D(cctkGH,i,j,k+1)] - gyyL);
	fac2 = (gyyL - gyy[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]);
	gyys=0;
	if(fac1*fac2 >= 0.0) gyys = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	gyy_rhs[index] = gyy_rhsL  
          + ( (1.0 - gyys)*absbetazLfac + gyys*betazL2fac)*hdzi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order)!
	
	//Advect gyz :	
	// x-direction
	fac1 = (gyz[CCTK_GFINDEX3D(cctkGH,i+1,j,k)] - gyzL);
	fac2 = (gyzL - gyz[CCTK_GFINDEX3D(cctkGH,i-1,j,k)]);
	gyzs=0;
	if(fac1*fac2 >= 0.0) gyzs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	gyz_rhsL = gyz_rhsL  
          + ( (1.0 - gyzs)*absbetaxLfac + gyzs*betaxL2fac )*hdxi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order)!
	// y-direction
	fac1 = (gyz[CCTK_GFINDEX3D(cctkGH,i,j+1,k)] - gyzL);
	fac2 = (gyzL - gyz[CCTK_GFINDEX3D(cctkGH,i,j-1,k)]);
	gyzs=0;
	if(fac1*fac2 >= 0.0) gyzs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	gyz_rhsL = gyz_rhsL  
          + ( (1.0 - gyzs)*absbetayLfac + gyzs*betayL2fac )*hdyi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order)!
	// z-direction
	fac1 = (gyz[CCTK_GFINDEX3D(cctkGH,i,j,k+1)] - gyzL);
	fac2 = (gyzL - gyz[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]);
	gyzs=0;
	if(fac1*fac2 >= 0.0) gyzs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	gyz_rhs[index] = gyz_rhsL  
          + ( (1.0 - gyzs)*absbetazLfac + gyzs*betazL2fac)*hdzi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order): D33gf(gyz,i,j,k);

	//Advect gzz :	
	// x-direction
	fac1 = (gzz[CCTK_GFINDEX3D(cctkGH,i+1,j,k)] - gzzL);
	fac2 = (gzzL - gzz[CCTK_GFINDEX3D(cctkGH,i-1,j,k)]);
	gzzs=0;
	if(fac1*fac2 >= 0.0) gzzs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	gzz_rhsL = gzz_rhsL 
          + ( (1.0 - gzzs)*absbetaxLfac + gzzs*betaxL2fac )*hdxi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order): D11gf(gzz,i,j,k);
	// y-direction
	fac1 = (gzz[CCTK_GFINDEX3D(cctkGH,i,j+1,k)] - gzzL);
	fac2 = (gzzL - gzz[CCTK_GFINDEX3D(cctkGH,i,j-1,k)]);
	gzzs=0;
	if(fac1*fac2 >= 0.0) gzzs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	gzz_rhsL = gzz_rhsL  
          + ( (1.0 - gzzs)*absbetayLfac + gzzs*betayL2fac )*hdyi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order): D22gf(gzz,i,j,k);
	// z-direction
	fac1 = (gzz[CCTK_GFINDEX3D(cctkGH,i,j,k+1)] - gzzL);
	fac2 = (gzzL - gzz[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]);
	gzzs=0;
	if(fac1*fac2 >= 0.0) gzzs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	gzz_rhs[index] = gzz_rhsL
          + ( (1.0 - gzzs)*absbetazLfac + gzzs*betazL2fac)*hdzi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order): D33gf(gzz,i,j,k);
  
	//Advect Axx :	
	// x-direction
	fac1 = (Axx[CCTK_GFINDEX3D(cctkGH,i+1,j,k)] - AxxL);
	fac2 = (AxxL - Axx[CCTK_GFINDEX3D(cctkGH,i-1,j,k)]);
	Axxs=0;
	if(fac1*fac2 >= 0.0) Axxs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	Axx_rhsL = Axx_rhsL  
          + ( (1.0 - Axxs)*absbetaxLfac + Axxs*betaxL2fac )*hdxi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order): D11gf(Axx,i,j,k);
	// y-direction
	fac1 = (Axx[CCTK_GFINDEX3D(cctkGH,i,j+1,k)] - AxxL);
	fac2 = (AxxL - Axx[CCTK_GFINDEX3D(cctkGH,i,j-1,k)]);
	Axxs=0;
	if(fac1*fac2 >= 0.0) Axxs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	Axx_rhsL = Axx_rhsL 
          + ( (1.0 - Axxs)*absbetayLfac + Axxs*betayL2fac )*hdyi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order): D22gf(Axx,i,j,k);
	// z-direction
	fac1 = (Axx[CCTK_GFINDEX3D(cctkGH,i,j,k+1)] - AxxL);
	fac2 = (AxxL - Axx[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]);
	Axxs=0;
	if(fac1*fac2 >= 0.0) Axxs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	Axx_rhs[index] = Axx_rhsL  
          + ( (1.0 - Axxs)*absbetazLfac + Axxs*betazL2fac)*hdzi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order): D33gf(Axx,i,j,k);

	//Advect Axy :	
	// x-direction
	fac1 = (Axy[CCTK_GFINDEX3D(cctkGH,i+1,j,k)] - AxyL);
	fac2 = (AxyL - Axy[CCTK_GFINDEX3D(cctkGH,i-1,j,k)]);
	Axys=0;
	if(fac1*fac2 >= 0.0) Axys = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	Axy_rhsL = Axy_rhsL  
          + ( (1.0 - Axys)*absbetaxLfac + Axys*betaxL2fac )*hdxi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order): D11gf(Axy,i,j,k);
	// y-direction
	fac1 = (Axy[CCTK_GFINDEX3D(cctkGH,i,j+1,k)] - AxyL);
	fac2 = (AxyL - Axy[CCTK_GFINDEX3D(cctkGH,i,j-1,k)]);
	Axys=0;
	if(fac1*fac2 >= 0.0) Axys = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	Axy_rhsL = Axy_rhsL  
          + ( (1.0 - Axys)*absbetayLfac + Axys*betayL2fac )*hdyi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order): D22gf(Axy,i,j,k);
	// z-direction
	fac1 = (Axy[CCTK_GFINDEX3D(cctkGH,i,j,k+1)] - AxyL);
	fac2 = (AxyL - Axy[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]);
	Axys=0;
	if(fac1*fac2 >= 0.0) Axys = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	Axy_rhs[index] = Axy_rhsL  
          + ( (1.0 - Axys)*absbetazLfac + Axys*betazL2fac)*hdzi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order): D33gf(Axy,i,j,k);
	
	//Advect Axz :	
	// x-direction
	fac1 = (Axz[CCTK_GFINDEX3D(cctkGH,i+1,j,k)] - AxzL);
	fac2 = (AxzL - Axz[CCTK_GFINDEX3D(cctkGH,i-1,j,k)]);
	Axzs=0;
	if(fac1*fac2 >= 0.0) Axzs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	Axz_rhsL = Axz_rhsL  
          + ( (1.0 - Axzs)*absbetaxLfac + Axzs*betaxL2fac )*hdxi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order): D11gf(Axz,i,j,k);
	// y-direction
	fac1 = (Axz[CCTK_GFINDEX3D(cctkGH,i,j+1,k)] - AxzL);
	fac2 = (AxzL - Axz[CCTK_GFINDEX3D(cctkGH,i,j-1,k)]);
	Axzs=0;
	if(fac1*fac2 >= 0.0) Axzs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	Axz_rhsL = Axz_rhsL  
          + ( (1.0 - Axzs)*absbetayLfac + Axzs*betayL2fac )*hdyi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order): D22gf(Axz,i,j,k);
	// z-direction
	fac1 = (Axz[CCTK_GFINDEX3D(cctkGH,i,j,k+1)] - AxzL);
	fac2 = (AxzL - Axz[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]);
	Axzs=0;
	if(fac1*fac2 >= 0.0) Axzs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	Axz_rhs[index] = Axz_rhsL  
          + ( (1.0 - Axzs)*absbetazLfac + Axzs*betazL2fac)*hdzi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order): D33gf(Axz,i,j,k);
	
	//Advect Ayy :	
	// x-direction
	fac1 = (Ayy[CCTK_GFINDEX3D(cctkGH,i+1,j,k)] - AyyL);
	fac2 = (AyyL - Ayy[CCTK_GFINDEX3D(cctkGH,i-1,j,k)]);
	Ayys=0;
	if(fac1*fac2 >= 0.0) Ayys = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	Ayy_rhsL = Ayy_rhsL  
          + ( (1.0 - Ayys)*absbetaxLfac + Ayys*betaxL2fac )*hdxi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order): D11gf(Ayy,i,j,k);
	// y-direction
	fac1 = (Ayy[CCTK_GFINDEX3D(cctkGH,i,j+1,k)] - AyyL);
	fac2 = (AyyL - Ayy[CCTK_GFINDEX3D(cctkGH,i,j-1,k)]);
	Ayys=0;
	if(fac1*fac2 >= 0.0) Ayys = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	Ayy_rhsL = Ayy_rhsL  
          + ( (1.0 - Ayys)*absbetayLfac + Ayys*betayL2fac )*hdyi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order): D22gf(Ayy,i,j,k);
	// z-direction
	fac1 = (Ayy[CCTK_GFINDEX3D(cctkGH,i,j,k+1)] - AyyL);
	fac2 = (AyyL - Ayy[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]);
	Ayys=0;
	if(fac1*fac2 >= 0.0) Ayys = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	Ayy_rhs[index] = Ayy_rhsL  
          + ( (1.0 - Ayys)*absbetazLfac + Ayys*betazL2fac)*hdzi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order): D33gf(Ayy,i,j,k);
	
	//Advect Ayz :	
	// x-direction
	fac1 = (Ayz[CCTK_GFINDEX3D(cctkGH,i+1,j,k)] - AyzL);
	fac2 = (AyzL - Ayz[CCTK_GFINDEX3D(cctkGH,i-1,j,k)]);
	Ayzs=0;
	if(fac1*fac2 >= 0.0) Ayzs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	Ayz_rhsL = Ayz_rhsL  
          + ( (1.0 - Ayzs)*absbetaxLfac + Ayzs*betaxL2fac )*hdxi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order): D11gf(Ayz,i,j,k);
	// y-direction
	fac1 = (Ayz[CCTK_GFINDEX3D(cctkGH,i,j+1,k)] - AyzL);
	fac2 = (AyzL - Ayz[CCTK_GFINDEX3D(cctkGH,i,j-1,k)]);
	Ayzs=0;
	if(fac1*fac2 >= 0.0) Ayzs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	Ayz_rhsL = Ayz_rhsL  
          + ( (1.0 - Ayzs)*absbetayLfac + Ayzs*betayL2fac )*hdyi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order): D22gf(Ayz,i,j,k);
	// z-direction
	fac1 = (Ayz[CCTK_GFINDEX3D(cctkGH,i,j,k+1)] - AyzL);
	fac2 = (AyzL - Ayz[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]);
	Ayzs=0;
	if(fac1*fac2 >= 0.0) Ayzs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	Ayz_rhs[index] = Ayz_rhsL  
          + ( (1.0 - Ayzs)*absbetazLfac + Ayzs*betazL2fac)*hdzi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order): D33gf(Ayz,i,j,k);

	//Advect Azz :	
	// x-direction
	fac1 = (Azz[CCTK_GFINDEX3D(cctkGH,i+1,j,k)] - AzzL);
	fac2 = (AzzL - Azz[CCTK_GFINDEX3D(cctkGH,i-1,j,k)]);
	Azzs=0;
	if(fac1*fac2 >= 0.0) Azzs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	Azz_rhsL = Azz_rhsL 
          + ( (1.0 - Azzs)*absbetaxLfac + Azzs*betaxL2fac )*hdxi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order): D11gf(Azz,i,j,k);
	// y-direction
	fac1 = (Azz[CCTK_GFINDEX3D(cctkGH,i,j+1,k)] - AzzL);
	fac2 = (AzzL - Azz[CCTK_GFINDEX3D(cctkGH,i,j-1,k)]);
	Azzs=0;
	if(fac1*fac2 >= 0.0) Azzs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	Azz_rhsL = Azz_rhsL  
          + ( (1.0 - Azzs)*absbetayLfac + Azzs*betayL2fac )*hdyi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order): D22gf(Azz,i,j,k);
	// z-direction
	fac1 = (Azz[CCTK_GFINDEX3D(cctkGH,i,j,k+1)] - AzzL);
	fac2 = (AzzL - Azz[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]);
	Azzs=0;
	if(fac1*fac2 >= 0.0) Azzs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	Azz_rhs[index] = Azz_rhsL
          + ( (1.0 - Azzs)*absbetazLfac + Azzs*betazL2fac)*hdzi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order): D33gf(Azz,i,j,k);

	//Advect phi :	
	// x-direction
	fac1 = (phi[CCTK_GFINDEX3D(cctkGH,i+1,j,k)] - phiL);
	fac2 = (phiL - phi[CCTK_GFINDEX3D(cctkGH,i-1,j,k)]);
	phis=0;
	if(fac1*fac2 >= 0.0) phis = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	phi_rhsL = phi_rhsL 
          + ( (1.0 - phis)*absbetaxLfac + phis*betaxL2fac )*hdxi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order): D11gf(phi,i,j,k);
	// y-direction
	fac1 = (phi[CCTK_GFINDEX3D(cctkGH,i,j+1,k)] - phiL);
	fac2 = (phiL - phi[CCTK_GFINDEX3D(cctkGH,i,j-1,k)]);
	phis=0;
	if(fac1*fac2 >= 0.0) phis = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	phi_rhsL = phi_rhsL  
          + ( (1.0 - phis)*absbetayLfac + phis*betayL2fac )*hdyi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order): D22gf(phi,i,j,k);
	// z-direction
	fac1 = (phi[CCTK_GFINDEX3D(cctkGH,i,j,k+1)] - phiL);
	fac2 = (phiL - phi[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]);
	phis=0;
	if(fac1*fac2 >= 0.0) phis = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	phi_rhs[index] = phi_rhsL
          + ( (1.0 - phis)*absbetazLfac + phis*betazL2fac)*hdzi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order): D33gf(phi,i,j,k);


	//Advect trK :	
	// x-direction
	fac1 = (trK[CCTK_GFINDEX3D(cctkGH,i+1,j,k)] - trKL);
	fac2 = (trKL - trK[CCTK_GFINDEX3D(cctkGH,i-1,j,k)]);
	trKs=0;
	if(fac1*fac2 >= 0.0) trKs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	trK_rhsL = trK_rhsL 
          + ( (1.0 - trKs)*absbetaxLfac + trKs*betaxL2fac )*hdxi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order): D11gf(trK,i,j,k);
	// y-direction
	fac1 = (trK[CCTK_GFINDEX3D(cctkGH,i,j+1,k)] - trKL);
	fac2 = (trKL - trK[CCTK_GFINDEX3D(cctkGH,i,j-1,k)]);
	trKs=0;
	if(fac1*fac2 >= 0.0) trKs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	trK_rhsL = trK_rhsL  
          + ( (1.0 - trKs)*absbetayLfac + trKs*betayL2fac )*hdyi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order): D22gf(trK,i,j,k);
	// z-direction
	fac1 = (trK[CCTK_GFINDEX3D(cctkGH,i,j,k+1)] - trKL);
	fac2 = (trKL - trK[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]);
	trKs=0;
	if(fac1*fac2 >= 0.0) trKs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	trK_rhs[index] = trK_rhsL
          + ( (1.0 - trKs)*absbetazLfac + trKs*betazL2fac)*hdzi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order): D33gf(trK,i,j,k);



	//Advect Gamx :	
	// x-direction
	fac1 = (Gamx[CCTK_GFINDEX3D(cctkGH,i+1,j,k)] - GamxL);
	fac2 = (GamxL - Gamx[CCTK_GFINDEX3D(cctkGH,i-1,j,k)]);
	Gamxs=0;
	if(fac1*fac2 >= 0.0) Gamxs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	Gamx_rhsL = Gamx_rhsL  
          + ( (1.0 - Gamxs)*absbetaxLfac + Gamxs*betaxL2fac )*hdxi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order): D11gf(Gamx,i,j,k);
	// y-direction
	fac1 = (Gamx[CCTK_GFINDEX3D(cctkGH,i,j+1,k)] - GamxL);
	fac2 = (GamxL - Gamx[CCTK_GFINDEX3D(cctkGH,i,j-1,k)]);
	Gamxs=0;
	if(fac1*fac2 >= 0.0) Gamxs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	Gamx_rhsL = Gamx_rhsL 
          + ( (1.0 - Gamxs)*absbetayLfac + Gamxs*betayL2fac )*hdyi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order): D22gf(Gamx,i,j,k);
	// z-direction
	fac1 = (Gamx[CCTK_GFINDEX3D(cctkGH,i,j,k+1)] - GamxL);
	fac2 = (GamxL - Gamx[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]);
	Gamxs=0;
	if(fac1*fac2 >= 0.0) Gamxs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	Gamx_rhs[index] = Gamx_rhsL  
          + ( (1.0 - Gamxs)*absbetazLfac + Gamxs*betazL2fac)*hdzi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order): D33gf(Gamx,i,j,k);

	//Advect Gamy :	
	// x-direction
	fac1 = (Gamy[CCTK_GFINDEX3D(cctkGH,i+1,j,k)] - GamyL);
	fac2 = (GamyL - Gamy[CCTK_GFINDEX3D(cctkGH,i-1,j,k)]);
	Gamys=0;
	if(fac1*fac2 >= 0.0) Gamys = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	Gamy_rhsL = Gamy_rhsL  
          + ( (1.0 - Gamys)*absbetaxLfac + Gamys*betaxL2fac )*hdxi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order): D11gf(Gamy,i,j,k);
	// y-direction
	fac1 = (Gamy[CCTK_GFINDEX3D(cctkGH,i,j+1,k)] - GamyL);
	fac2 = (GamyL - Gamy[CCTK_GFINDEX3D(cctkGH,i,j-1,k)]);
	Gamys=0;
	if(fac1*fac2 >= 0.0) Gamys = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	Gamy_rhsL = Gamy_rhsL  
          + ( (1.0 - Gamys)*absbetayLfac + Gamys*betayL2fac )*hdyi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order): D22gf(Gamy,i,j,k);
	// z-direction
	fac1 = (Gamy[CCTK_GFINDEX3D(cctkGH,i,j,k+1)] - GamyL);
	fac2 = (GamyL - Gamy[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]);
	Gamys=0;
	if(fac1*fac2 >= 0.0) Gamys = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	Gamy_rhs[index] = Gamy_rhsL  
          + ( (1.0 - Gamys)*absbetazLfac + Gamys*betazL2fac)*hdzi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order): D33gf(Gamy,i,j,k);
	
	//Advect Gamz :	
	// x-direction
	fac1 = (Gamz[CCTK_GFINDEX3D(cctkGH,i+1,j,k)] - GamzL);
	fac2 = (GamzL - Gamz[CCTK_GFINDEX3D(cctkGH,i-1,j,k)]);
	Gamzs=0;
	if(fac1*fac2 >= 0.0) Gamzs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	Gamz_rhsL = Gamz_rhsL
          + ( (1.0 - Gamzs)*absbetaxLfac + Gamzs*betaxL2fac )*hdxi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order): D11gf(Gamz,i,j,k);
	// y-direction
	fac1 = (Gamz[CCTK_GFINDEX3D(cctkGH,i,j+1,k)] - GamzL);
	fac2 = (GamzL - Gamz[CCTK_GFINDEX3D(cctkGH,i,j-1,k)]);
	Gamzs=0;
	if(fac1*fac2 >= 0.0) Gamzs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	Gamz_rhsL = Gamz_rhsL  
          + ( (1.0 - Gamzs)*absbetayLfac + Gamzs*betayL2fac )*hdyi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order): D22gf(Gamz,i,j,k);
	// z-direction
	fac1 = (Gamz[CCTK_GFINDEX3D(cctkGH,i,j,k+1)] - GamzL);
	fac2 = (GamzL - Gamz[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]);
	Gamzs=0;
	if(fac1*fac2 >= 0.0) Gamzs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	Gamz_rhs[index] = Gamz_rhsL
          + ( (1.0 - Gamzs)*absbetazLfac + Gamzs*betazL2fac)*hdzi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order): D33gf(Gamz,i,j,k);
	

      }
}

extern "C" void CCTK_FCALL CCTK_FNAME(shift_upwind2)
  (const cGH **cctkGH,int *cctk_lsh,int *nghostzones,double *dT, double *dx, double *dy, double *dz,
   double *gxx, double *gxy, double *gxz, double *gyy, double *gyz, double *gzz, 
   double *gxx_rhs, double *gxy_rhs, double *gxz_rhs, double *gyy_rhs, double *gyz_rhs, double *gzz_rhs, 
   double *Axx, double *Axy, double *Axz, double *Ayy, double *Ayz, double *Azz, 
   double *Axx_rhs, double *Axy_rhs, double *Axz_rhs, double *Ayy_rhs, double *Ayz_rhs, double *Azz_rhs, 
   double *phi, double *phi_rhs, 
   double *trK, double *trK_rhs, 
   double *Gamx, double *Gamy, double *Gamz,
   double *Gamx_rhs, double *Gamy_rhs, double *Gamz_rhs, 
   double *betax, double *betay, double *betaz, int *Symmetry)
{
  shift_upwind2(*cctkGH, cctk_lsh, nghostzones, *dT,  *dx,  *dy,  *dz,
		gxx, gxy, gxz, gyy, gyz, gzz, 
		gxx_rhs, gxy_rhs, gxz_rhs, gyy_rhs, gyz_rhs, gzz_rhs, 
		Axx, Axy, Axz, Ayy, Ayz, Azz, 
		Axx_rhs, Axy_rhs, Axz_rhs, Ayy_rhs, Ayz_rhs, Azz_rhs, 
		phi, phi_rhs, 
		trK, trK_rhs, 
		Gamx, Gamy, Gamz,
		Gamx_rhs, Gamy_rhs, Gamz_rhs, 
		betax, betay, betaz, *Symmetry);
}
