//-----------------------------------------------------------------------------
//
// RHS for first order shift (for puncture calculations, typically):
// \partial_t beta^i = 3/4 \tilde{Gamma}^i - eta beta^i
//
//-----------------------------------------------------------------------------

/* Define macros used in calculations */
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

extern "C" void CCTK_FCALL CCTK_FNAME(puncture_firstorder_shiftRHS)
  (const cGH **cctkGH,int *nghostzones,int *cctk_lsh, int *Symmetry,
   double *dx,double *dy,double *dz,
   double *x,double *y,double *z,double *r,
   double *PhysicalRadius,double *RadiusDerivative,double *RadiusDerivative2,
   double *Gammax,double *Gammay,double *Gammaz,
   double *gupxx,double *gupxy,double *gupxz,
   double *gupyy,double *gupyz,double *gupzz, 
   double *betax,double *betay,double *betaz,
   double *betax_rhs,double *betay_rhs,double *betaz_rhs,
   double *eta,int *firstorder_shift_enable_convert_fisheye_to_physical);

extern "C" void puncture_firstorder_shiftRHS(const cGH *cctkGH,int *nghostzones,int *cctk_lsh, int Symmetry,
					     double dx,double dy,double dz,
					     double *x,double *y,double *z,double *r,
					     double *PhysicalRadius,double *RadiusDerivative,double *RadiusDerivative2,
					     double *Gammax,double *Gammay,double *Gammaz,
					     double *gupxx,double *gupxy,double *gupxz,
					     double *gupyy,double *gupyz,double *gupzz, 
					     double *betax,double *betay,double *betaz,
					     double *betax_rhs,double *betay_rhs,double *betaz_rhs,
					     double eta,int firstorder_shift_enable_convert_fisheye_to_physical) {
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

  for(int k=kstart;k<kend;k++)
    for(int j=jstart;j<jend;j++)
      for(int i=istart;i<iend;i++) {
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	double GammaxL,GammayL,GammazL;
	if(firstorder_shift_enable_convert_fisheye_to_physical==1) {
	  double xL=x[index];
	  double yL=y[index];
	  double zL=z[index];
	  double rL=r[index];
	  double rLinv = 1.0/rL;
	  
	  double nxL=xL*rLinv;
	  double nyL=yL*rLinv;
	  double nzL=zL*rLinv;
	  
	  double prL=PhysicalRadius[index];
	  double rdL=RadiusDerivative[index];
	  double rd2L=RadiusDerivative2[index];
	  double rbrL=prL*rLinv;
	  
	  // pow() is expensive!
	  double fac1 = pow(SQR(rbrL) * rdL,F1o6);
	  double fac2 = F1o3*(rdL/prL-rLinv) + F1o6*rd2L/rdL;
	  
	  //   we need to calculate \tilde{Gamma}^i from 
	  //   \bar{\tilde{Gamma}}^i, via jacobian and derivatives
	  //   f4 is 'F' in Josh's notes, and df4*nx gives you dF/dx
	  double f4=1.0/QAD(fac1);
	  double df4=F1o3*f4*(-4.0/prL*(rdL-rbrL)-2.0*rd2L/rdL);
	  
	  double df4x=df4*nxL;
	  double df4y=df4*nyL;
	  double df4z=df4*nzL;
	  double jxx=rbrL+nxL*nxL*(rdL-rbrL);
	  double jxy=nxL*nyL*(rdL-rbrL);
	  double jxz=nxL*nzL*(rdL-rbrL);
	  double jyy=rbrL+nyL*nyL*(rdL-rbrL);
	  double jyz=nyL*nzL*(rdL-rbrL);
	  double jzz=rbrL+nzL*nzL*(rdL-rbrL);
	  
	  double jxxx=CUB(nxL)*rd2L+(rdL-rbrL)*rLinv*(3.0*nxL-3.0*CUB(nxL));
	  double jyyy=CUB(nyL)*rd2L+(rdL-rbrL)*rLinv*(3.0*nyL-3.0*CUB(nyL));
	  double jzzz=CUB(nzL)*rd2L+(rdL-rbrL)*rLinv*(3.0*nzL-3.0*CUB(nzL));
	  
	  double jxxy=SQR(nxL)*nyL*rd2L+(rdL-rbrL)*rLinv*(nyL-3.0*SQR(nxL)*nyL);
	  double jxyy=nxL*SQR(nyL)*rd2L+(rdL-rbrL)*rLinv*(nxL-3.0*nxL*SQR(nyL));
	  double jxxz=SQR(nxL)*nzL*rd2L+(rdL-rbrL)*rLinv*(nzL-3.0*SQR(nxL)*nzL);
	  double jxzz=nxL*SQR(nzL)*rd2L+(rdL-rbrL)*rLinv*(nxL-3.0*nxL*SQR(nzL));
	  double jyyz=SQR(nyL)*nzL*rd2L+(rdL-rbrL)*rLinv*(nzL-3.0*SQR(nyL)*nzL);
	  double jyzz=nyL*SQR(nzL)*rd2L+(rdL-rbrL)*rLinv*(nyL-3.0*nyL*SQR(nzL));
	  
	  double jxyz=nxL*nyL*nzL*(rd2L-3.0*(rdL-rbrL)*rLinv);

	  double gamxL=Gammax[index];
	  double gamyL=Gammay[index];
	  double gamzL=Gammaz[index];
	  double gupxxL=gupxx[index];
	  double gupxyL=gupxy[index];
	  double gupxzL=gupxz[index];
	  double gupyyL=gupyy[index];
	  double gupyzL=gupyz[index];
	  double gupzzL=gupzz[index];
	  
	  double gamxpL=f4*(jxx*gamxL+jxy*gamyL+jxz*gamzL)+ 
	    gupxxL*(0.5*df4x*jxx-f4*jxxx)+ 
	    gupyyL*(0.5*df4y*jxy-f4*jxyy)+ 
	    gupzzL*(0.5*df4z*jxz-f4*jxzz)+ 
	    gupxyL*(0.5*(df4x*jxy+df4y*jxx)-2.0*f4*jxxy)+ 
	    gupxzL*(0.5*(df4x*jxz+df4z*jxx)-2.0*f4*jxxz)+ 
	    gupyzL*(0.5*(df4y*jxz+df4z*jxy)-2.0*f4*jxyz);
	  double gamypL=f4*(jxy*gamxL+jyy*gamyL+jyz*gamzL)+ 
	    gupxxL*(0.5*df4x*jxy-f4*jxxy)+ 
	    gupyyL*(0.5*df4y*jyy-f4*jyyy)+ 
	    gupzzL*(0.5*df4z*jyz-f4*jyzz)+ 
	    gupxyL*(0.5*(df4x*jyy+df4y*jxy)-2.0*f4*jxyy)+ 
	    gupxzL*(0.5*(df4x*jyz+df4z*jxy)-2.0*f4*jxyz)+ 
	    gupyzL*(0.5*(df4y*jyz+df4z*jyy)-2.0*f4*jyyz);
	  double gamzpL=f4*(jxz*gamxL+jyz*gamyL+jzz*gamzL)+ 
	    gupxxL*(0.5*df4x*jxz-f4*jxxz)+ 
	    gupyyL*(0.5*df4y*jyz-f4*jyyz)+ 
	    gupzzL*(0.5*df4z*jzz-f4*jzzz)+ 
	    gupxyL*(0.5*(df4x*jyz+df4y*jxz)-2.0*f4*jxyz)+ 
	    gupxzL*(0.5*(df4x*jzz+df4z*jxz)-2.0*f4*jxzz)+ 
	    gupyzL*(0.5*(df4y*jzz+df4z*jyz)-2.0*f4*jyzz);
	  GammaxL = gamxpL;
	  GammayL = gamypL;
	  GammazL = gamzpL;
	} else {
	  GammaxL = Gammax[index];
	  GammayL = Gammay[index];
	  GammazL = Gammaz[index];
	}
	
	// Okay.  Now we have Gamma^i's in physical coordinates.  Now we update der shift:
	betax_rhs[index] = 0.75 * GammaxL - eta * betax[index];
	betay_rhs[index] = 0.75 * GammayL - eta * betay[index];
	betaz_rhs[index] = 0.75 * GammazL - eta * betaz[index];
      }
}



extern "C" void CCTK_FCALL CCTK_FNAME(puncture_firstorder_shiftRHS)
  (const cGH **cctkGH,int *nghostzones,int *cctk_lsh, int *Symmetry,
   double *dx,double *dy,double *dz,
   double *x,double *y,double *z,double *r,
   double *PhysicalRadius,double *RadiusDerivative,double *RadiusDerivative2,
   double *Gammax,double *Gammay,double *Gammaz,
   double *gupxx,double *gupxy,double *gupxz,
   double *gupyy,double *gupyz,double *gupzz, 
   double *betax,double *betay,double *betaz,
   double *betax_rhs,double *betay_rhs,double *betaz_rhs,
   double *eta,int *firstorder_shift_enable_convert_fisheye_to_physical)
{
  puncture_firstorder_shiftRHS(*cctkGH,nghostzones,cctk_lsh, *Symmetry,
			       *dx, *dy, *dz,
			       x, y, z, r,
			       PhysicalRadius, RadiusDerivative, RadiusDerivative2,
			       Gammax, Gammay, Gammaz,
			       gupxx, gupxy, gupxz,
			       gupyy, gupyz, gupzz, 
			       betax, betay, betaz,
			       betax_rhs, betay_rhs, betaz_rhs,
			       *eta,*firstorder_shift_enable_convert_fisheye_to_physical);
}
