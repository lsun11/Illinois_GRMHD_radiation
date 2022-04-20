//-----------------------------------------------------------------------------
//
// Carry out predictor step
//
//-----------------------------------------------------------------------------

/* Define macros used in calculations */
#define INITVALUE  (42)
#define EPSILON 1.0e-8
#define INV(x) ((1) / (x))
#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * (x) * (x))
#define QAD(x) ((x) * (x) * (x) * (x))

#define F1o3 0.333333333333333333333333333333333
#define F1o6 0.166666666666666666666666666666666
#define F2o3 0.666666666666666666666666666666666

#include <stdio.h>
#include <stdlib.h>
#include "cctk.h"
#include <math.h>
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

#define KRANC_C
#include "GenericFD.h"

extern "C" void CCTK_FCALL CCTK_FNAME(hbpuncture_shiftRHS)
  (const cGH **cctkGH,int *nghostzones,int *cctk_lsh,
   double *dT,double *dx,double *dy,double *dz,int *enable_shift_advection,double *eta,
   double *RadiusDerivative,
   double *shiftx_old,double *shifty_old,double *shiftz_old,
   double *shiftx_rhs,double *shifty_rhs,double *shiftz_rhs,
   double *dtshiftx_old,double *dtshifty_old,double *dtshiftz_old,
   double *dtshiftx_rhs,double *dtshifty_rhs,double *dtshiftz_rhs,
   double *Gammax,double *Gammay,double *Gammaz,
   double *Gammax_rhs,double *Gammay_rhs,double *Gammaz_rhs,
   double *tempx,double *tempy,double *tempz,int *Symmetry, 
   int *hbpuncture_shift_convert_Gammai_fisheye_to_physical,
   double *r,double *eta_final_value,int *eta_falloff_enable,double *eta_falloff_radius,double *eta_falloff_dr,
   double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,double *phi);

extern "C" void hbpuncture_shiftRHS(const cGH *cctkGH,int *nghostzones,int *cctk_lsh,
				    double dT,double dx,double dy,double dz,int enable_shift_advection,double eta,
				    double *RadiusDerivative,
				    double *shiftx_old,double *shifty_old,double *shiftz_old,
				    double *shiftx_rhs,double *shifty_rhs,double *shiftz_rhs,
				    double *dtshiftx_old,double *dtshifty_old,double *dtshiftz_old,
				    double *dtshiftx_rhs,double *dtshifty_rhs,double *dtshiftz_rhs,
				    double *Gammax,double *Gammay,double *Gammaz,
				    double *Gammax_rhs,double *Gammay_rhs,double *Gammaz_rhs,
				    double *tempx,double *tempy,double *tempz,int Symmetry, 
				    int hbpuncture_shift_convert_Gammai_fisheye_to_physical,
				    double *r, double eta_final_value, int eta_falloff_enable, double eta_falloff_radius, double eta_falloff_dr,
				    double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,double *phi) {
  /* Initialise finite differencing variables.  NEED THIS FOR GenericFD.h */
#include "../../GenFD_decl_set_varCPP.h"

  /* Set up variables used in the grid loop for the physical grid points */
  //We're dealing with centered derivatives in the main loop below!
#ifdef FD_C2
  int ghostsize=1;
#endif
#ifdef FD_C4
  int ghostsize=2;
#endif
#ifdef FD_C6
  int ghostsize=3;
#endif

  int istart = ghostsize;
  int jstart = ghostsize;
  int kstart = ghostsize;
  int iend = cctk_lsh[0] - ghostsize;
  int jend = cctk_lsh[1] - ghostsize;
  int kend = cctk_lsh[2] - ghostsize;

  //Following lines needed since nghostzones[0] = ORDER, and
  //   not ORDER-1 in axisymmetry
  //   (so that rotation can be done on multiprocessor runs)
  if(Symmetry==4) {
    istart--;
    iend++;
  }

#pragma omp parallel for
  for(int k=kstart;k<kend;k++)
    for(int j=jstart;j<jend;j++)
      for(int i=istart;i<iend;i++) {
	double fac1,fac2;
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

	double dtshiftx_oldL = dtshiftx_old[index];
	double dtshifty_oldL = dtshifty_old[index];
	double dtshiftz_oldL = dtshiftz_old[index];

	double shiftx_oldL = shiftx_old[index];
	double shifty_oldL = shifty_old[index];
	double shiftz_oldL = shiftz_old[index];

	double Gammax_rhsL = Gammax_rhs[index];
	double Gammay_rhsL = Gammay_rhs[index];
	double Gammaz_rhsL = Gammaz_rhs[index];

	double Gammaxrhs_conv,Gammayrhs_conv,Gammazrhs_conv;

	double RadiusDerivativeL2inv = 1.0/SQR(RadiusDerivative[index]);
	
	if(hbpuncture_shift_convert_Gammai_fisheye_to_physical==0) {
	  Gammaxrhs_conv = Gammax_rhsL;
	  Gammayrhs_conv = Gammay_rhsL;
	  Gammazrhs_conv = Gammaz_rhsL;
	} else if(hbpuncture_shift_convert_Gammai_fisheye_to_physical==1) {
	  Gammaxrhs_conv = Gammax_rhsL*RadiusDerivativeL2inv;
	  Gammayrhs_conv = Gammay_rhsL*RadiusDerivativeL2inv;
	  Gammazrhs_conv = Gammaz_rhsL*RadiusDerivativeL2inv;
	}
  
	//  dt B^i      = dt gam^i   -eta * B_i
	/*
	  dtshiftx_rhs[index]=Gammaxrhs_conv - eta*dtshiftx_oldL;
	  dtshifty_rhs[index]=Gammayrhs_conv - eta*dtshifty_oldL; 
	  dtshiftz_rhs[index]=Gammazrhs_conv - eta*dtshiftz_oldL;
	*/

	//etaevolve_rhs[index] = 0.0;

	if(eta_falloff_enable==0) {
	  dtshiftx_rhs[index]=Gammaxrhs_conv - eta*dtshiftx_oldL;
	  dtshifty_rhs[index]=Gammayrhs_conv - eta*dtshifty_oldL; 
	  dtshiftz_rhs[index]=Gammazrhs_conv - eta*dtshiftz_oldL;
	} else if(eta_falloff_enable==1) {
	  //-erf((x-20-5.0*2)/5.0)*(5-0.25)*0.5+(5-0.25)*0.5+0.25
	  dtshiftx_rhs[index]=Gammaxrhs_conv - dtshiftx_oldL*(-erf((r[index]-eta_falloff_radius-eta_falloff_dr*2.0)/eta_falloff_dr)*(eta-eta_final_value)*0.5 + (eta-eta_final_value)*0.5 + eta_final_value);
	  dtshifty_rhs[index]=Gammayrhs_conv - dtshifty_oldL*(-erf((r[index]-eta_falloff_radius-eta_falloff_dr*2.0)/eta_falloff_dr)*(eta-eta_final_value)*0.5 + (eta-eta_final_value)*0.5 + eta_final_value);
	  dtshiftz_rhs[index]=Gammazrhs_conv - dtshiftz_oldL*(-erf((r[index]-eta_falloff_radius-eta_falloff_dr*2.0)/eta_falloff_dr)*(eta-eta_final_value)*0.5 + (eta-eta_final_value)*0.5 + eta_final_value);
	} else if(eta_falloff_enable==2) {
	  // Jena & RIT groups' temporal & spatially-varying eta (arXiv:1009.0292)

	  printf("ERROR.  We have currently implemented a shifting-shift condition, which is incompatible with Jena & RIT groups' temporal & spatially-varying eta.\n"); exit(1);

	  double R0 = 1.31;

	  double WL=exp(-2.0*phi[index]);
	  double phix=D1gf(phi,i,j,k);
	  double phiy=D2gf(phi,i,j,k);
	  double phiz=D3gf(phi,i,j,k);
	  double Wx = -2.0*WL*phix;
	  double Wy = -2.0*WL*phiy;
	  double Wz = -2.0*WL*phiz;

	  //Note that gupij was set inside the BSSN RHS function.
	  double gupxxL = gupxx[index];
	  double gupxyL = gupxy[index];
	  double gupxzL = gupxz[index];
	  double gupyyL = gupyy[index];
	  double gupyzL = gupyz[index];
	  double gupzzL = gupzz[index];
	  
	  //assume a=2:
	  double onemWa = 1.0-WL*WL;
	  
	  //assume b=2:
	  double denom = onemWa*onemWa;

	  double etanew = R0 * sqrt(Wx*Wx*gupxxL +
				    2.0*Wx*Wy*gupxyL +
				    2.0*Wx*Wz*gupxzL +
				    Wy*Wy*gupyyL +
				    2.0*Wy*Wz*gupyzL +
				    Wz*Wz*gupzzL)/denom;

	  //if(etanew>5.0) etanew=5.0;
	  //if(etanew<0.125) etanew=0.125;

	  dtshiftx_rhs[index]=Gammaxrhs_conv - etanew*dtshiftx_oldL;
	  dtshifty_rhs[index]=Gammayrhs_conv - etanew*dtshifty_oldL; 
	  dtshiftz_rhs[index]=Gammazrhs_conv - etanew*dtshiftz_oldL;
	} else if(eta_falloff_enable==3) {	
	  /*
	  // AEI's time-evolved eta prescription:

	  double R0 = 1.31;

	  double WL=exp(-2.0*phi[index]);
	  double phix=D1gf(phi,i,j,k);
	  double phiy=D2gf(phi,i,j,k);
	  double phiz=D3gf(phi,i,j,k);
	  double Wx = -2.0*WL*phix;
	  double Wy = -2.0*WL*phiy;
	  double Wz = -2.0*WL*phiz;

	  //Note that gupij was set inside the BSSN RHS function.
	  double gupxxL = gupxx[index];
	  double gupxyL = gupxy[index];
	  double gupxzL = gupxz[index];
	  double gupyyL = gupyy[index];
	  double gupyzL = gupyz[index];
	  double gupzzL = gupzz[index];
	  
	  //assume a=1 FOR AEI DAMPING EVOLUTION:
	  double onemWa = 1.0-WL;
	  
	  //assume b=2:
	  double denom = onemWa*onemWa;

	  double etanew = R0 * sqrt(Wx*Wx*gupxxL +
				    2.0*Wx*Wy*gupxyL +
				    2.0*Wx*Wz*gupxzL +
				    Wy*Wy*gupyyL +
				    2.0*Wy*Wz*gupyzL +
				    Wz*Wz*gupzzL)/denom;


	  dtshiftx_rhs[index]=Gammaxrhs_conv - etaevolve[index]*dtshiftx_oldL;
	  dtshifty_rhs[index]=Gammayrhs_conv - etaevolve[index]*dtshifty_oldL; 
	  dtshiftz_rhs[index]=Gammazrhs_conv - etaevolve[index]*dtshiftz_oldL;
	  */
	}


	//  dt beta^i  =          3/4 B^i
	shiftx_rhs[index] = 0.75*dtshiftx_oldL;
	shifty_rhs[index] = 0.75*dtshifty_oldL;
	shiftz_rhs[index] = 0.75*dtshiftz_oldL;

	if(enable_shift_advection == 1) {
	  double dxshiftx = D1gf(shiftx_old,i,j,k);
	  double dyshiftx = D2gf(shiftx_old,i,j,k);
	  double dzshiftx = D3gf(shiftx_old,i,j,k);
	
	  double dxshifty = D1gf(shifty_old,i,j,k);
	  double dyshifty = D2gf(shifty_old,i,j,k);
	  double dzshifty = D3gf(shifty_old,i,j,k);
	
	  double dxshiftz = D1gf(shiftz_old,i,j,k);
	  double dyshiftz = D2gf(shiftz_old,i,j,k);
	  double dzshiftz = D3gf(shiftz_old,i,j,k);

	  double dxgamx = D1gf(Gammax,i,j,k);
	  double dygamx = D2gf(Gammax,i,j,k);
	  double dzgamx = D3gf(Gammax,i,j,k);
	
	  double dxgamy = D1gf(Gammay,i,j,k);
	  double dygamy = D2gf(Gammay,i,j,k);
	  double dzgamy = D3gf(Gammay,i,j,k);
	
	  double dxgamz = D1gf(Gammaz,i,j,k);
	  double dygamz = D2gf(Gammaz,i,j,k);
	  double dzgamz = D3gf(Gammaz,i,j,k);

	  double dxbx = D1gf(dtshiftx_old,i,j,k);
	  double dybx = D2gf(dtshiftx_old,i,j,k);
	  double dzbx = D3gf(dtshiftx_old,i,j,k);
	
	  double dxby = D1gf(dtshifty_old,i,j,k);
	  double dyby = D2gf(dtshifty_old,i,j,k);
	  double dzby = D3gf(dtshifty_old,i,j,k);
	
	  double dxbz = D1gf(dtshiftz_old,i,j,k);
	  double dybz = D2gf(dtshiftz_old,i,j,k);
	  double dzbz = D3gf(dtshiftz_old,i,j,k);

	  //                        +  beta^j dj beta^i
	  shiftx_rhs[index] += (shiftx_oldL*dxshiftx + shifty_oldL*dyshiftx + shiftz_oldL*dzshiftx);
	  shifty_rhs[index] += (shiftx_oldL*dxshifty + shifty_oldL*dyshifty + shiftz_oldL*dzshifty);
	  shiftz_rhs[index] += (shiftx_oldL*dxshiftz + shifty_oldL*dyshiftz + shiftz_oldL*dzshiftz);
	  //                        - beta^k dk gam^i + beta^j djB^i
	  dtshiftx_rhs[index] += -RadiusDerivativeL2inv*(shiftx_oldL*dxgamx+shifty_oldL*dygamx+shiftz_oldL*dzgamx)+ 
	    (shiftx_oldL*dxbx+shifty_oldL*dybx+shiftz_oldL*dzbx);
	  dtshifty_rhs[index] += -RadiusDerivativeL2inv*(shiftx_oldL*dxgamy+shifty_oldL*dygamy+shiftz_oldL*dzgamy)+ 
	    (shiftx_oldL*dxby+shifty_oldL*dyby+shiftz_oldL*dzby);
	  dtshiftz_rhs[index] += -RadiusDerivativeL2inv*(shiftx_oldL*dxgamz+shifty_oldL*dygamz+shiftz_oldL*dzgamz)+ 
	    (shiftx_oldL*dxbz+shifty_oldL*dybz+shiftz_oldL*dzbz);
	}
      }
}

extern "C" void CCTK_FCALL CCTK_FNAME(hbpuncture_shiftRHS)
  (const cGH **cctkGH,int *nghostzones,int *cctk_lsh,
   double *dT,double *dx,double *dy,double *dz,int *enable_shift_advection,double *eta,
   double *RadiusDerivative,
   double *shiftx_old,double *shifty_old,double *shiftz_old,
   double *shiftx_rhs,double *shifty_rhs,double *shiftz_rhs,
   double *dtshiftx_old,double *dtshifty_old,double *dtshiftz_old,
   double *dtshiftx_rhs,double *dtshifty_rhs,double *dtshiftz_rhs,
   double *Gammax,double *Gammay,double *Gammaz,
   double *Gammax_rhs,double *Gammay_rhs,double *Gammaz_rhs,
   double *tempx,double *tempy,double *tempz,int *Symmetry, 
   int *hbpuncture_shift_convert_Gammai_fisheye_to_physical,
   double *r,double *eta_final_value,int *eta_falloff_enable,double *eta_falloff_radius,double *eta_falloff_dr,
   double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,double *phi)
{
  hbpuncture_shiftRHS(*cctkGH,nghostzones,cctk_lsh,
		      *dT,*dx,*dy,*dz,*enable_shift_advection,*eta,
		      RadiusDerivative,
		      shiftx_old,shifty_old,shiftz_old,
		      shiftx_rhs,shifty_rhs,shiftz_rhs,
		      dtshiftx_old,dtshifty_old,dtshiftz_old,
		      dtshiftx_rhs,dtshifty_rhs,dtshiftz_rhs,
		      Gammax,Gammay,Gammaz,
		      Gammax_rhs,Gammay_rhs,Gammaz_rhs,
		      tempx,tempy,tempz,*Symmetry, 
		      *hbpuncture_shift_convert_Gammai_fisheye_to_physical,
		      r,*eta_final_value,*eta_falloff_enable,*eta_falloff_radius,*eta_falloff_dr,
		      gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,phi);
}
