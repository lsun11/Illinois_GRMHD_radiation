//-----------------------------------------------------------------------------
// Compute 1+log lapse RHS: \partial_t \alpha = -2 \alpha trK
//-----------------------------------------------------------------------------
/* Define macros used in calculations */
#define INV(x) ((1) / (x))
#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * (x) * (x))
#define QAD(x) ((x) * (x) * (x) * (x))

#define EPSILON 1.e-8

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
#define FD_SET_BY_USER
#define FD_C4
#include "GenericFD.h"

extern "C" void CCTK_FCALL CCTK_FNAME(opl_lapse_upwind_4)
  (const cGH **cctkGH,int *nghostzones,int *cctk_lsh, int *Symmetry,int *opl_advect_enable,
   double *dT,double *dx, double *dy, double *dz,
   double *lapsex,double *lapsey,double *lapsez,
   double *betax,double *betay,double *betaz,
   double *lapse_old, double *lapse_rhs);

extern "C" void opl_lapse_upwind_4(const cGH *cctkGH,int *nghostzones,int *cctk_lsh, int Symmetry, int opl_advect_enable,
				   double dT,double dx, double dy, double dz,
				   double *lapsex,double *lapsey,double *lapsez,
				   double *betax,double *betay,double *betaz,
				   double *lapse_old, double *lapse_rhs) {
  /* Initialise finite differencing variables.  NEED THIS FOR GenericFD.h */
#include "../../GenFD_decl_set_varCPP.h"

  /* Set up variables used in the grid loop for the physical grid points */
  //WARNING: MUST BE CAREFUL WITH 4TH ORDER UPWIND STENCILS: They extend out 3 gridpoints!!
  int istart = 2;
  int jstart = 2;
  int kstart = 2;
  int iend = cctk_lsh[0]-2;
  int jend = cctk_lsh[1]-2;
  int kend = cctk_lsh[2]-2;

  //Following lines needed since nghostzones[0] = ORDER, and
  //   not ORDER-1 in axisymmetry
  //   (so that rotation can be done on multiprocessor runs)
  if(Symmetry==4) {
    istart--;
    iend++;
  }
  if(opl_advect_enable == 2) {
#pragma omp parallel for
    for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
      if(i==2 || j==2 || k==2 || i==cctk_lsh[0]-3 || j==cctk_lsh[1]-3 || k==cctk_lsh[2]-3) {
	// The following if() statement is quite complex and takes a while to compute, 
	//     so we don't want to evaluate it for all i,j,k!  Thus we have the above if()
	//     statement to reduce the number of evaluations.
	if((i==2 && j>2 && k>2 && i<cctk_lsh[0]-3 && j<cctk_lsh[1]-3 && k<cctk_lsh[2]-3) || 
	   (j==2 && i>2 && k>2 && i<cctk_lsh[0]-3 && j<cctk_lsh[1]-3 && k<cctk_lsh[2]-3) || 
	   (k==2 && i>2 && j>2 && i<cctk_lsh[0]-3 && j<cctk_lsh[1]-3 && k<cctk_lsh[2]-3) || 
	   (i==cctk_lsh[0]-3 && i>2 && j>2 && k>2 && j<cctk_lsh[1]-3 && k<cctk_lsh[2]-3) || 
	   (j==cctk_lsh[1]-3 && i>2 && j>2 && k>2 && i<cctk_lsh[0]-3 && k<cctk_lsh[2]-3) || 
	   (k==cctk_lsh[2]-3 && i>2 && j>2 && k>2 && i<cctk_lsh[0]-3 && j<cctk_lsh[1]-3)) {
	  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

	  double betaxL = betax[index];
	  double betayL = betay[index];
	  double betazL = betaz[index];
      
	  double lapse_rhsL = lapse_rhs[index];
      
	  if(betaxL > 0.0) {
	    if(i==cctk_lsh[0]-3) {
	      lapse_rhsL += betaxL*D1gf(lapse_old,i,j,k);
	    } else {
	      lapse_rhsL += betaxL*D1_up_gt_gf(lapse_old,i,j,k);
	    }
	  } else {
	    if(i==2) {
	      lapse_rhsL += betaxL*D1gf(lapse_old,i,j,k);
	    } else {
	      lapse_rhsL += betaxL*D1_up_lt_gf(lapse_old,i,j,k);
	    }
	  }
      
	  if(betayL > 0.0) {
	    if(j==cctk_lsh[1]-3) {
	      lapse_rhsL += betayL*D2gf(lapse_old,i,j,k);
	    } else {
	      lapse_rhsL += betayL*D2_up_gt_gf(lapse_old,i,j,k);
	    }
	  } else {
	    if(j==2) {
	      lapse_rhsL += betayL*D2gf(lapse_old,i,j,k);
	    } else {
	      lapse_rhsL += betayL*D2_up_lt_gf(lapse_old,i,j,k);
	    }
	  }
      
	  if(betazL > 0.0) {
	    if(k==cctk_lsh[2]-3) {
	      lapse_rhsL += betazL*D3gf(lapse_old,i,j,k);
	    } else {
	      lapse_rhsL += betazL*D3_up_gt_gf(lapse_old,i,j,k);
	    }
	  } else {
	    if(k==2) {
	      lapse_rhsL += betazL*D3gf(lapse_old,i,j,k);
	    } else {
	      lapse_rhsL += betazL*D3_up_lt_gf(lapse_old,i,j,k);
	    }
	  }
	  lapse_rhs[index] = lapse_rhsL;
	}
      }
    }
  } else if(opl_advect_enable == 1) {
    kstart=jstart=istart=2;
    kend=cctk_lsh[2]-2;
    jend=cctk_lsh[1]-2;
    iend=cctk_lsh[0]-2;
#pragma omp parallel for
    for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
      // Use mixed 1st/2nd order upwind advection!
      int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
      double fac1,fac2;

      double lapse_rhsL = lapse_rhs[index];
      double lapse_oldL = lapse_old[index];
      
      double betaxL = betax[index];
      double betayL = betay[index];
      double betazL = betaz[index];

      double betaxL2fac = betaxL*betaxL*dT;
      double betayL2fac = betayL*betayL*dT;
      double betazL2fac = betazL*betazL*dT;

      double absbetaxLfac = fabs(betaxL)*dx;
      double absbetayLfac = fabs(betayL)*dy;
      double absbetazLfac = fabs(betazL)*dz;
      
      //Add shift^i partial_i lapse
      lapse_rhsL += opl_advect_enable*(betaxL*lapsex[index] + betayL*lapsey[index] + betazL*lapsez[index]);
      
      //Advect lapse_rhs
      // x-direction
      fac1 = (lapse_old[CCTK_GFINDEX3D(cctkGH,i+1,j,k)] - lapse_oldL);
      fac2 = (lapse_oldL - lapse_old[CCTK_GFINDEX3D(cctkGH,i-1,j,k)]);
      double lapse_olds=0;
      if(fac1*fac2 >= 0.0) lapse_olds = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
      lapse_rhsL = lapse_rhsL  
	+ ( (1.0 - lapse_olds)*absbetaxLfac + lapse_olds*betaxL2fac )*hdxi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order)!
      // y-direction
      fac1 = (lapse_old[CCTK_GFINDEX3D(cctkGH,i,j+1,k)] - lapse_oldL);
      fac2 = (lapse_oldL - lapse_old[CCTK_GFINDEX3D(cctkGH,i,j-1,k)]);
      lapse_olds=0;
      if(fac1*fac2 >= 0.0) lapse_olds = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
      lapse_rhsL = lapse_rhsL
	+ ( (1.0 - lapse_olds)*absbetayLfac + lapse_olds*betayL2fac )*hdyi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order)!
      // z-direction
      fac1 = (lapse_old[CCTK_GFINDEX3D(cctkGH,i,j,k+1)] - lapse_oldL);
      fac2 = (lapse_oldL - lapse_old[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]);
      lapse_olds=0;
      if(fac1*fac2 >= 0.0) lapse_olds = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
      lapse_rhs[index] = lapse_rhsL
	+ ( (1.0 - lapse_olds)*absbetazLfac + lapse_olds*betazL2fac)*hdzi2*(fac1-fac2); //fac1-fac2 is usual second derivative numerator (to second order)!
    }
  } else {
    printf("UPWINDING TYPE %d IS NOT SUPPORTED!\n",opl_advect_enable);
  }
}
extern "C" void CCTK_FCALL CCTK_FNAME(opl_lapse_upwind_4)
  (const cGH **cctkGH,int *nghostzones,int *cctk_lsh,int *Symmetry,int *opl_advect_enable,
   double *dT,double *dx, double *dy, double *dz,
   double *lapsex,double *lapsey,double *lapsez,
   double *betax,double *betay,double *betaz,
   double *lapse_old, double *lapse_rhs)
{
  opl_lapse_upwind_4(*cctkGH,nghostzones,cctk_lsh,*Symmetry,*opl_advect_enable,
		     *dT,*dx,*dy,*dz,
		     lapsex,lapsey,lapsez,
		     betax,betay,betaz,
		     lapse_old, lapse_rhs);
}
