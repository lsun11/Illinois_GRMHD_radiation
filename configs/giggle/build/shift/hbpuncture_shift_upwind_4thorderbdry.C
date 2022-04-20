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

extern "C" void CCTK_FCALL hbpuncture_upwind_4_
  (const cGH **cctkGH,int *nghostzones,int *cctk_lsh,double *dT, double *dx, double *dy, double *dz,
   double *RadiusDerivative,int *shift_advect_enable,int *bssn_enable_shift_upwind,
   double *shiftx_old,double *shifty_old,double *shiftz_old,
   double *shiftx_rhs,double *shifty_rhs,double *shiftz_rhs,
   double *dtshiftx_old,double *dtshifty_old,double *dtshiftz_old,
   double *dtshiftx_rhs,double *dtshifty_rhs,double *dtshiftz_rhs,
   double *Gammax,double *Gammay,double *Gammaz,
   double *tempx,double *tempy,double *tempz,int *Symmetry);

extern "C" void hbpuncture_upwind_4(const cGH *cctkGH,int *nghostzones,int *cctk_lsh,
				    double dT,double dx,double dy,double dz,
				    double *RadiusDerivative, int shift_advect_enable,int bssn_enable_shift_upwind,
				    double *shiftx_old,double *shifty_old,double *shiftz_old,
				    double *shiftx_rhs,double *shifty_rhs,double *shiftz_rhs,
				    double *dtshiftx_old,double *dtshifty_old,double *dtshiftz_old,
				    double *dtshiftx_rhs,double *dtshifty_rhs,double *dtshiftz_rhs,
				    double *Gammax,double *Gammay,double *Gammaz,
				    double *tempx,double *tempy,double *tempz,int Symmetry) {
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

  if(shift_advect_enable==2) {
    //Note that with this shift advection choice, we do BOTH advection + upwinding here!
    for(int k=kstart;k<kend;k++)
      for(int j=jstart;j<jend;j++)
	for(int i=istart;i<iend;i++) {
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

	      double RadiusDerivativeL2inv = 1.0/SQR(RadiusDerivative[index]);
	      double betaxL = shiftx_old[index];
	      double betayL = shifty_old[index];
	      double betazL = shiftz_old[index];

	      double dtshiftxL = dtshiftx_old[index];
	      double dtshiftyL = dtshifty_old[index];
	      double dtshiftzL = dtshiftz_old[index];

	      double betax_rhsL = shiftx_rhs[index];
	      double betay_rhsL = shifty_rhs[index];
	      double betaz_rhsL = shiftz_rhs[index];

	      double dtshiftx_rhsL = dtshiftx_rhs[index];
	      double dtshifty_rhsL = dtshifty_rhs[index];
	      double dtshiftz_rhsL = dtshiftz_rhs[index];
	
	      if(betaxL > 0.0) {
		if(i==cctk_lsh[0]-3) {
		  betax_rhsL += betaxL*D1gf(shiftx_old,i,j,k);
		  betay_rhsL += betaxL*D1gf(shifty_old,i,j,k);
		  betaz_rhsL += betaxL*D1gf(shiftz_old,i,j,k);

		  dtshiftx_rhsL += betaxL*D1gf(dtshiftx_old,i,j,k);
		  dtshifty_rhsL += betaxL*D1gf(dtshifty_old,i,j,k);
		  dtshiftz_rhsL += betaxL*D1gf(dtshiftz_old,i,j,k);
		} else {
		  betax_rhsL += betaxL*D1_up_gt_gf(shiftx_old,i,j,k);
		  betay_rhsL += betaxL*D1_up_gt_gf(shifty_old,i,j,k);
		  betaz_rhsL += betaxL*D1_up_gt_gf(shiftz_old,i,j,k);

		  dtshiftx_rhsL += betaxL*D1_up_gt_gf(dtshiftx_old,i,j,k);
		  dtshifty_rhsL += betaxL*D1_up_gt_gf(dtshifty_old,i,j,k);
		  dtshiftz_rhsL += betaxL*D1_up_gt_gf(dtshiftz_old,i,j,k);
		}
	      } else {
		if(i==2) {
		  betax_rhsL += betaxL*D1gf(shiftx_old,i,j,k);
		  betay_rhsL += betaxL*D1gf(shifty_old,i,j,k);
		  betaz_rhsL += betaxL*D1gf(shiftz_old,i,j,k);

		  dtshiftx_rhsL += betaxL*D1gf(dtshiftx_old,i,j,k);
		  dtshifty_rhsL += betaxL*D1gf(dtshifty_old,i,j,k);
		  dtshiftz_rhsL += betaxL*D1gf(dtshiftz_old,i,j,k);
		} else {
		  betax_rhsL += betaxL*D1_up_lt_gf(shiftx_old,i,j,k);
		  betay_rhsL += betaxL*D1_up_lt_gf(shifty_old,i,j,k);
		  betaz_rhsL += betaxL*D1_up_lt_gf(shiftz_old,i,j,k);

		  dtshiftx_rhsL += betaxL*D1_up_lt_gf(dtshiftx_old,i,j,k);
		  dtshifty_rhsL += betaxL*D1_up_lt_gf(dtshifty_old,i,j,k);
		  dtshiftz_rhsL += betaxL*D1_up_lt_gf(dtshiftz_old,i,j,k);
		}
	      }

	      if(betayL > 0.0) {
		if(j==cctk_lsh[1]-3) {
		  betax_rhsL += betayL*D2gf(shiftx_old,i,j,k);
		  betay_rhsL += betayL*D2gf(shifty_old,i,j,k);
		  betaz_rhsL += betayL*D2gf(shiftz_old,i,j,k);

		  dtshiftx_rhsL += betayL*D2gf(dtshiftx_old,i,j,k);
		  dtshifty_rhsL += betayL*D2gf(dtshifty_old,i,j,k);
		  dtshiftz_rhsL += betayL*D2gf(dtshiftz_old,i,j,k);
		} else {
		  betax_rhsL += betayL*D2_up_gt_gf(shiftx_old,i,j,k);
		  betay_rhsL += betayL*D2_up_gt_gf(shifty_old,i,j,k);
		  betaz_rhsL += betayL*D2_up_gt_gf(shiftz_old,i,j,k);

		  dtshiftx_rhsL += betayL*D2_up_gt_gf(dtshiftx_old,i,j,k);
		  dtshifty_rhsL += betayL*D2_up_gt_gf(dtshifty_old,i,j,k);
		  dtshiftz_rhsL += betayL*D2_up_gt_gf(dtshiftz_old,i,j,k);
		}
	      } else {
		if(j==2) {
		  betax_rhsL += betayL*D2gf(shiftx_old,i,j,k);
		  betay_rhsL += betayL*D2gf(shifty_old,i,j,k);
		  betaz_rhsL += betayL*D2gf(shiftz_old,i,j,k);

		  dtshiftx_rhsL += betayL*D2gf(dtshiftx_old,i,j,k);
		  dtshifty_rhsL += betayL*D2gf(dtshifty_old,i,j,k);
		  dtshiftz_rhsL += betayL*D2gf(dtshiftz_old,i,j,k);
		} else {
		  betax_rhsL += betayL*D2_up_lt_gf(shiftx_old,i,j,k);
		  betay_rhsL += betayL*D2_up_lt_gf(shifty_old,i,j,k);
		  betaz_rhsL += betayL*D2_up_lt_gf(shiftz_old,i,j,k);

		  dtshiftx_rhsL += betayL*D2_up_lt_gf(dtshiftx_old,i,j,k);
		  dtshifty_rhsL += betayL*D2_up_lt_gf(dtshifty_old,i,j,k);
		  dtshiftz_rhsL += betayL*D2_up_lt_gf(dtshiftz_old,i,j,k);

		}
	      }

	      if(betazL > 0.0) {
		if(k==cctk_lsh[2]-3) {
		  betax_rhsL += betazL*D3gf(shiftx_old,i,j,k);
		  betay_rhsL += betazL*D3gf(shifty_old,i,j,k);
		  betaz_rhsL += betazL*D3gf(shiftz_old,i,j,k);

		  dtshiftx_rhsL += betazL*D3gf(dtshiftx_old,i,j,k);
		  dtshifty_rhsL += betazL*D3gf(dtshifty_old,i,j,k);
		  dtshiftz_rhsL += betazL*D3gf(dtshiftz_old,i,j,k);
		} else {
		  betax_rhsL += betazL*D3_up_gt_gf(shiftx_old,i,j,k);
		  betay_rhsL += betazL*D3_up_gt_gf(shifty_old,i,j,k);
		  betaz_rhsL += betazL*D3_up_gt_gf(shiftz_old,i,j,k);

		  dtshiftx_rhsL += betazL*D3_up_gt_gf(dtshiftx_old,i,j,k);
		  dtshifty_rhsL += betazL*D3_up_gt_gf(dtshifty_old,i,j,k);
		  dtshiftz_rhsL += betazL*D3_up_gt_gf(dtshiftz_old,i,j,k);
		}
	      } else {
		if(k==2) {
		  betax_rhsL += betazL*D3gf(shiftx_old,i,j,k);
		  betay_rhsL += betazL*D3gf(shifty_old,i,j,k);
		  betaz_rhsL += betazL*D3gf(shiftz_old,i,j,k);

		  dtshiftx_rhsL += betazL*D3gf(dtshiftx_old,i,j,k);
		  dtshifty_rhsL += betazL*D3gf(dtshifty_old,i,j,k);
		  dtshiftz_rhsL += betazL*D3gf(dtshiftz_old,i,j,k);
		} else {
		  betax_rhsL += betazL*D3_up_lt_gf(shiftx_old,i,j,k);
		  betay_rhsL += betazL*D3_up_lt_gf(shifty_old,i,j,k);
		  betaz_rhsL += betazL*D3_up_lt_gf(shiftz_old,i,j,k);

		  dtshiftx_rhsL += betazL*D3_up_lt_gf(dtshiftx_old,i,j,k);
		  dtshifty_rhsL += betazL*D3_up_lt_gf(dtshifty_old,i,j,k);
		  dtshiftz_rhsL += betazL*D3_up_lt_gf(dtshiftz_old,i,j,k);
		}
	      }

	      //Note: Gamma^i's on RHS of dtshift equations must NOT contain shift advection terms!
	      //  ... so here we subtract them off:
	      if(bssn_enable_shift_upwind == 2) {
		if(betaxL > 0.0) {
		  if(i==cctk_lsh[0]-3) {
		    dtshiftx_rhsL -= betaxL*D1gf(Gammax,i,j,k)*RadiusDerivativeL2inv;
		    dtshifty_rhsL -= betaxL*D1gf(Gammay,i,j,k)*RadiusDerivativeL2inv;
		    dtshiftz_rhsL -= betaxL*D1gf(Gammaz,i,j,k)*RadiusDerivativeL2inv;
		  } else {
		    dtshiftx_rhsL -= betaxL*D1_up_gt_gf(Gammax,i,j,k)*RadiusDerivativeL2inv;
		    dtshifty_rhsL -= betaxL*D1_up_gt_gf(Gammay,i,j,k)*RadiusDerivativeL2inv;
		    dtshiftz_rhsL -= betaxL*D1_up_gt_gf(Gammaz,i,j,k)*RadiusDerivativeL2inv;
		  }
		} else {
		  if(i==2) {
		    dtshiftx_rhsL -= betaxL*D1gf(Gammax,i,j,k)*RadiusDerivativeL2inv;
		    dtshifty_rhsL -= betaxL*D1gf(Gammay,i,j,k)*RadiusDerivativeL2inv;
		    dtshiftz_rhsL -= betaxL*D1gf(Gammaz,i,j,k)*RadiusDerivativeL2inv;
		  } else {
		    dtshiftx_rhsL -= betaxL*D1_up_lt_gf(Gammax,i,j,k)*RadiusDerivativeL2inv;
		    dtshifty_rhsL -= betaxL*D1_up_lt_gf(Gammay,i,j,k)*RadiusDerivativeL2inv;
		    dtshiftz_rhsL -= betaxL*D1_up_lt_gf(Gammaz,i,j,k)*RadiusDerivativeL2inv;
		  }
		}
		if(betayL > 0.0) {
		  if(j==cctk_lsh[1]-3) {
		    dtshiftx_rhsL -= betayL*D2gf(Gammax,i,j,k)*RadiusDerivativeL2inv;
		    dtshifty_rhsL -= betayL*D2gf(Gammay,i,j,k)*RadiusDerivativeL2inv;
		    dtshiftz_rhsL -= betayL*D2gf(Gammaz,i,j,k)*RadiusDerivativeL2inv;
		  } else {
		    dtshiftx_rhsL -= betayL*D2_up_gt_gf(Gammax,i,j,k)*RadiusDerivativeL2inv;
		    dtshifty_rhsL -= betayL*D2_up_gt_gf(Gammay,i,j,k)*RadiusDerivativeL2inv;
		    dtshiftz_rhsL -= betayL*D2_up_gt_gf(Gammaz,i,j,k)*RadiusDerivativeL2inv;
		  }
		} else {
		  if(j==2) {
		    dtshiftx_rhsL -= betayL*D2gf(Gammax,i,j,k)*RadiusDerivativeL2inv;
		    dtshifty_rhsL -= betayL*D2gf(Gammay,i,j,k)*RadiusDerivativeL2inv;
		    dtshiftz_rhsL -= betayL*D2gf(Gammaz,i,j,k)*RadiusDerivativeL2inv;
		  } else {
		    dtshiftx_rhsL -= betayL*D2_up_lt_gf(Gammax,i,j,k)*RadiusDerivativeL2inv;
		    dtshifty_rhsL -= betayL*D2_up_lt_gf(Gammay,i,j,k)*RadiusDerivativeL2inv;
		    dtshiftz_rhsL -= betayL*D2_up_lt_gf(Gammaz,i,j,k)*RadiusDerivativeL2inv;
		  }
		}
		if(betazL > 0.0) {
		  if(k==cctk_lsh[2]-3) {
		    dtshiftx_rhsL -= betazL*D3gf(Gammax,i,j,k)*RadiusDerivativeL2inv;
		    dtshifty_rhsL -= betazL*D3gf(Gammay,i,j,k)*RadiusDerivativeL2inv;
		    dtshiftz_rhsL -= betazL*D3gf(Gammaz,i,j,k)*RadiusDerivativeL2inv;
		  } else {
		    dtshiftx_rhsL -= betazL*D3_up_gt_gf(Gammax,i,j,k)*RadiusDerivativeL2inv;
		    dtshifty_rhsL -= betazL*D3_up_gt_gf(Gammay,i,j,k)*RadiusDerivativeL2inv;
		    dtshiftz_rhsL -= betazL*D3_up_gt_gf(Gammaz,i,j,k)*RadiusDerivativeL2inv;
		  }
		} else {
		  if(k==2) {
		    dtshiftx_rhsL -= betazL*D3gf(Gammax,i,j,k)*RadiusDerivativeL2inv;
		    dtshifty_rhsL -= betazL*D3gf(Gammay,i,j,k)*RadiusDerivativeL2inv;
		    dtshiftz_rhsL -= betazL*D3gf(Gammaz,i,j,k)*RadiusDerivativeL2inv;
		  } else {
		    dtshiftx_rhsL -= betazL*D3_up_lt_gf(Gammax,i,j,k)*RadiusDerivativeL2inv;
		    dtshifty_rhsL -= betazL*D3_up_lt_gf(Gammay,i,j,k)*RadiusDerivativeL2inv;
		    dtshiftz_rhsL -= betazL*D3_up_lt_gf(Gammaz,i,j,k)*RadiusDerivativeL2inv;
		  }
		}
	      } else if(bssn_enable_shift_upwind == 0) {
		dtshiftx_rhsL -= RadiusDerivativeL2inv*(betaxL*D1gf(Gammax,i,j,k) + betayL*D2gf(Gammax,i,j,k) + betazL*D3gf(Gammax,i,j,k));
		dtshifty_rhsL -= RadiusDerivativeL2inv*(betaxL*D1gf(Gammay,i,j,k) + betayL*D2gf(Gammay,i,j,k) + betazL*D3gf(Gammay,i,j,k));	
		dtshiftz_rhsL -= RadiusDerivativeL2inv*(betaxL*D1gf(Gammaz,i,j,k) + betayL*D2gf(Gammaz,i,j,k) + betazL*D3gf(Gammaz,i,j,k));
	      } else if(bssn_enable_shift_upwind == 1) {
		printf("Shibata upwinding in BSSN sector (bssn_enable_shift_upwind == 1) and shift_advect_enable==2 NOT compatible!\n");
		printf("Please choose another set of parameters.\n");
	      }
	  
	      shiftx_rhs[index] = betax_rhsL;
	      shifty_rhs[index] = betay_rhsL;
	      shiftz_rhs[index] = betaz_rhsL;

	      dtshiftx_rhs[index] = dtshiftx_rhsL;
	      dtshifty_rhs[index] = dtshifty_rhsL;
	      dtshiftz_rhs[index] = dtshiftz_rhsL;
	    }
	  }
	}
  } else if(shift_advect_enable==1) {
    //Note that with this shift advection choice, we have already computed \beta^i partial_i [blah] in the RHS function.
    for(int k=0;k<cctk_lsh[2];k++)
      for(int j=0;j<cctk_lsh[1];j++)
	for(int i=0;i<cctk_lsh[0];i++) {
	  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  double RadiusDerivativeL2inv = 1.0/SQR(RadiusDerivative[index]);
	
	  tempx[index] = (dtshiftx_old[index]-Gammax[index]*RadiusDerivativeL2inv);
	  tempy[index] = (dtshifty_old[index]-Gammay[index]*RadiusDerivativeL2inv);
	  tempz[index] = (dtshiftz_old[index]-Gammaz[index]*RadiusDerivativeL2inv);
	}

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

	  double tempxL = tempx[index];
	  double tempyL = tempy[index];
	  double tempzL = tempz[index];
	
	  double dtshiftx_rhsL = dtshiftx_rhs[index];
	  double dtshifty_rhsL = dtshifty_rhs[index];
	  double dtshiftz_rhsL = dtshiftz_rhs[index];
	
	  //Advect dtshiftx (recall tempx[] array is set above):  
	  // x-direction
	  fac1 = (tempx[CCTK_GFINDEX3D(cctkGH,i+1,j,k)] - tempxL);
	  fac2 = (tempxL - tempx[CCTK_GFINDEX3D(cctkGH,i-1,j,k)]);
	  double tempxs=0;
	  if(fac1*fac2 >= 0.0) tempxs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	  dtshiftx_rhsL = dtshiftx_rhsL  + 
	    ((1.0-tempxs)*fabs(shiftx_oldL)*dx+tempxs*dT*SQR(shiftx_oldL))*hdxi2*(fac1-fac2); //(fac1-fac2)/(dx^i)^2 is usual second derivative (to second order)!
	  // y-direction
	  fac1 = (tempx[CCTK_GFINDEX3D(cctkGH,i,j+1,k)] - tempxL);
	  fac2 = (tempxL - tempx[CCTK_GFINDEX3D(cctkGH,i,j-1,k)]);
	  tempxs=0;
	  if(fac1*fac2 >= 0.0) tempxs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	  dtshiftx_rhsL = dtshiftx_rhsL  + 
	    ((1.0-tempxs)*fabs(shifty_oldL)*dy+tempxs*dT*SQR(shifty_oldL))*hdyi2*(fac1-fac2); //(fac1-fac2)/(dx^i)^2 is usual second derivative (to second order)!
	  // z-direction
	  fac1 = (tempx[CCTK_GFINDEX3D(cctkGH,i,j,k+1)] - tempxL);
	  fac2 = (tempxL - tempx[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]);
	  tempxs=0;
	  if(fac1*fac2 >= 0.0) tempxs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	  dtshiftx_rhs[index] = dtshiftx_rhsL  + 
	    ((1.0-tempxs)*fabs(shiftz_oldL)*dz+tempxs*dT*SQR(shiftz_oldL))*hdzi2*(fac1-fac2); //(fac1-fac2)/(dx^i)^2 is usual second derivative (to second order)!
	  
	  //Advect dtshifty (recall tempy[] array is set above):  
	  // x-direction
	  fac1 = (tempy[CCTK_GFINDEX3D(cctkGH,i+1,j,k)] - tempyL);
	  fac2 = (tempyL - tempy[CCTK_GFINDEX3D(cctkGH,i-1,j,k)]);
	  double tempys=0;
	  if(fac1*fac2 >= 0.0) tempys = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	  dtshifty_rhsL = dtshifty_rhsL  + 
	    ((1.0-tempys)*fabs(shiftx_oldL)*dx+tempys*dT*SQR(shiftx_oldL))*hdxi2*(fac1-fac2); //(fac1-fac2)/(dx^i)^2 is usual second derivative (to second order)!
	  // y-direction
	  fac1 = (tempy[CCTK_GFINDEX3D(cctkGH,i,j+1,k)] - tempyL);
	  fac2 = (tempyL - tempy[CCTK_GFINDEX3D(cctkGH,i,j-1,k)]);
	  tempys=0;
	  if(fac1*fac2 >= 0.0) tempys = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	  dtshifty_rhsL = dtshifty_rhsL  + 
	    ((1.0-tempys)*fabs(shifty_oldL)*dy+tempys*dT*SQR(shifty_oldL))*hdyi2*(fac1-fac2); //(fac1-fac2)/(dx^i)^2 is usual second derivative (to second order)!
	  // z-direction
	  fac1 = (tempy[CCTK_GFINDEX3D(cctkGH,i,j,k+1)] - tempyL);
	  fac2 = (tempyL - tempy[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]);
	  tempys=0;
	  if(fac1*fac2 >= 0.0) tempys = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	  dtshifty_rhs[index] = dtshifty_rhsL  + 
	    ((1.0-tempys)*fabs(shiftz_oldL)*dz+tempys*dT*SQR(shiftz_oldL))*hdzi2*(fac1-fac2); //(fac1-fac2)/(dx^i)^2 is usual second derivative (to second order)!
	  
	  //Advect dtshiftz (recall tempz[] array is set above):  
	  // x-direction
	  fac1 = (tempz[CCTK_GFINDEX3D(cctkGH,i+1,j,k)] - tempzL);
	  fac2 = (tempzL - tempz[CCTK_GFINDEX3D(cctkGH,i-1,j,k)]);
	  double tempzs=0;
	  if(fac1*fac2 >= 0.0) tempzs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	  dtshiftz_rhsL = dtshiftz_rhsL  + 
	    ((1.0-tempzs)*fabs(shiftx_oldL)*dx+tempzs*dT*SQR(shiftx_oldL))*hdxi2*(fac1-fac2); //(fac1-fac2)/(dx^i)^2 is usual second derivative (to second order)!
	  // y-direction
	  fac1 = (tempz[CCTK_GFINDEX3D(cctkGH,i,j+1,k)] - tempzL);
	  fac2 = (tempzL - tempz[CCTK_GFINDEX3D(cctkGH,i,j-1,k)]);
	  tempzs=0;
	  if(fac1*fac2 >= 0.0) tempzs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	  dtshiftz_rhsL = dtshiftz_rhsL  + 
	    ((1.0-tempzs)*fabs(shifty_oldL)*dy+tempzs*dT*SQR(shifty_oldL))*hdyi2*(fac1-fac2); //(fac1-fac2)/(dx^i)^2 is usual second derivative (to second order)!
	  // z-direction
	  fac1 = (tempz[CCTK_GFINDEX3D(cctkGH,i,j,k+1)] - tempzL);
	  fac2 = (tempzL - tempz[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]);
	  tempzs=0;
	  if(fac1*fac2 >= 0.0) tempzs = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	  dtshiftz_rhs[index] = dtshiftz_rhsL  + 
	    ((1.0-tempzs)*fabs(shiftz_oldL)*dz+tempzs*dT*SQR(shiftz_oldL))*hdzi2*(fac1-fac2); //(fac1-fac2)/(dx^i)^2 is usual second derivative (to second order)!
	  //
	  // done advecting dtshift
	  //

	  double shiftx_rhsL = shiftx_rhs[index];
	  double shifty_rhsL = shifty_rhs[index];
	  double shiftz_rhsL = shiftz_rhs[index];
	
	  //Advect shiftx:
	  // x-direction
	  fac1 = (shiftx_old[CCTK_GFINDEX3D(cctkGH,i+1,j,k)] - shiftx_oldL);
	  fac2 = (shiftx_oldL - shiftx_old[CCTK_GFINDEX3D(cctkGH,i-1,j,k)]);
	  double shiftx_olds=0;
	  if(fac1*fac2 >= 0.0) shiftx_olds = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	  shiftx_rhsL = shiftx_rhsL  +  
	    ((1.0-shiftx_olds)*fabs(shiftx_oldL)*dx+shiftx_olds*dT*SQR(shiftx_oldL))*hdxi2*(fac1-fac2); //(fac1-fac2)/(dx^i)^2 is usual second derivative (to second order)!
	  // y-direction
	  fac1 = (shiftx_old[CCTK_GFINDEX3D(cctkGH,i,j+1,k)] - shiftx_oldL);
	  fac2 = (shiftx_oldL - shiftx_old[CCTK_GFINDEX3D(cctkGH,i,j-1,k)]);
	  shiftx_olds=0;
	  if(fac1*fac2 >= 0.0) shiftx_olds = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	  shiftx_rhsL = shiftx_rhsL  +  
	    ((1.0-shiftx_olds)*fabs(shifty_oldL)*dy+shiftx_olds*dT*SQR(shifty_oldL))*hdyi2*(fac1-fac2); //(fac1-fac2)/(dx^i)^2 is usual second derivative (to second order)!
	  // z-direction
	  fac1 = (shiftx_old[CCTK_GFINDEX3D(cctkGH,i,j,k+1)] - shiftx_oldL);
	  fac2 = (shiftx_oldL - shiftx_old[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]);
	  shiftx_olds=0;
	  if(fac1*fac2 >= 0.0) shiftx_olds = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	  shiftx_rhs[index] = shiftx_rhsL  +  
	    ((1.0-shiftx_olds)*fabs(shiftz_oldL)*dz+shiftx_olds*dT*SQR(shiftz_oldL))*hdzi2*(fac1-fac2); //(fac1-fac2)/(dx^i)^2 is usual second derivative (to second order)!

	  //Advect shifty:
	  // x-direction
	  fac1 = (shifty_old[CCTK_GFINDEX3D(cctkGH,i+1,j,k)] - shifty_oldL);
	  fac2 = (shifty_oldL - shifty_old[CCTK_GFINDEX3D(cctkGH,i-1,j,k)]);
	  double shifty_olds=0;
	  if(fac1*fac2 >= 0.0) shifty_olds = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	  shifty_rhsL = shifty_rhsL  +  
	    ((1.0-shifty_olds)*fabs(shiftx_oldL)*dx+shifty_olds*dT*SQR(shiftx_oldL))*hdxi2*(fac1-fac2); //(fac1-fac2)/(dx^i)^2 is usual second derivative (to second order)!
	  //if(i==30 && j==31 && k==28) printf("shifty1 %.15e %.15e %.15e\n",shifty_rhsL,shifty_oldL,shifty_old[CCTK_GFINDEX3D(cctkGH,i-1,j,k)],shifty_olds);
	  // y-direction
	  fac1 = (shifty_old[CCTK_GFINDEX3D(cctkGH,i,j+1,k)] - shifty_oldL);
	  fac2 = (shifty_oldL - shifty_old[CCTK_GFINDEX3D(cctkGH,i,j-1,k)]);
	  shifty_olds=0;
	  if(fac1*fac2 >= 0.0) shifty_olds = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	  shifty_rhsL = shifty_rhsL  +  
	    ((1.0-shifty_olds)*fabs(shifty_oldL)*dy+shifty_olds*dT*SQR(shifty_oldL))*hdyi2*(fac1-fac2); //(fac1-fac2)/(dx^i)^2 is usual second derivative (to second order)!
	  //if(i==30 && j==31 && k==28) printf("shifty2 %.15e\n",shifty_rhsL);
	  // z-direction
	  fac1 = (shifty_old[CCTK_GFINDEX3D(cctkGH,i,j,k+1)] - shifty_oldL);
	  fac2 = (shifty_oldL - shifty_old[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]);
	  shifty_olds=0;
	  if(fac1*fac2 >= 0.0) shifty_olds = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	  shifty_rhs[index] = shifty_rhsL  +  
	    ((1.0-shifty_olds)*fabs(shiftz_oldL)*dz+shifty_olds*dT*SQR(shiftz_oldL))*hdzi2*(fac1-fac2); //(fac1-fac2)/(dx^i)^2 is usual second derivative (to second order)!
	  //if(i==30 && j==31 && k==28) printf("shifty3 %.15e\n",shifty_rhs[index]);

	  //Advect shiftz:
	  // x-direction
	  fac1 = (shiftz_old[CCTK_GFINDEX3D(cctkGH,i+1,j,k)] - shiftz_oldL);
	  fac2 = (shiftz_oldL - shiftz_old[CCTK_GFINDEX3D(cctkGH,i-1,j,k)]);
	  double shiftz_olds=0;
	  if(fac1*fac2 >= 0.0) shiftz_olds = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	  shiftz_rhsL = shiftz_rhsL  +  
	    ((1.0-shiftz_olds)*fabs(shiftx_oldL)*dx+shiftz_olds*dT*SQR(shiftx_oldL))*hdxi2*(fac1-fac2); //(fac1-fac2)/(dx^i)^2 is usual second derivative (to second order)!
	  // y-direction
	  fac1 = (shiftz_old[CCTK_GFINDEX3D(cctkGH,i,j+1,k)] - shiftz_oldL);
	  fac2 = (shiftz_oldL - shiftz_old[CCTK_GFINDEX3D(cctkGH,i,j-1,k)]);
	  shiftz_olds=0;
	  if(fac1*fac2 >= 0.0) shiftz_olds = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	  shiftz_rhsL = shiftz_rhsL  +  
	    ((1.0-shiftz_olds)*fabs(shifty_oldL)*dy+shiftz_olds*dT*SQR(shifty_oldL))*hdyi2*(fac1-fac2); //(fac1-fac2)/(dx^i)^2 is usual second derivative (to second order)!
	  // z-direction
	  fac1 = (shiftz_old[CCTK_GFINDEX3D(cctkGH,i,j,k+1)] - shiftz_oldL);
	  fac2 = (shiftz_oldL - shiftz_old[CCTK_GFINDEX3D(cctkGH,i,j,k-1)]);
	  shiftz_olds=0;
	  if(fac1*fac2 >= 0.0) shiftz_olds = (2.0*fac1*fac2 + EPSILON) / (fac1*fac1 + fac2*fac2 + EPSILON);
	  shiftz_rhs[index] = shiftz_rhsL  +  
	    ((1.0-shiftz_olds)*fabs(shiftz_oldL)*dz+shiftz_olds*dT*SQR(shiftz_oldL))*hdzi2*(fac1-fac2); //(fac1-fac2)/(dx^i)^2 is usual second derivative (to second order)!

	}
  }
}

extern "C" void CCTK_FCALL hbpuncture_upwind_4_
  (const cGH **cctkGH,int *nghostzones,int *cctk_lsh,double *dT, double *dx, double *dy, double *dz,
   double *RadiusDerivative, int *shift_advect_enable,int *bssn_enable_shift_upwind,
   double *shiftx_old,double *shifty_old,double *shiftz_old,
   double *shiftx_rhs,double *shifty_rhs,double *shiftz_rhs,
   double *dtshiftx_old,double *dtshifty_old,double *dtshiftz_old,
   double *dtshiftx_rhs,double *dtshifty_rhs,double *dtshiftz_rhs,
   double *Gammax,double *Gammay,double *Gammaz,
   double *tempx,double *tempy,double *tempz,int *Symmetry) 
{
  hbpuncture_upwind_4(*cctkGH,nghostzones,cctk_lsh, *dT,  *dx,  *dy,  *dz,
		      RadiusDerivative, *shift_advect_enable, *bssn_enable_shift_upwind,
		      shiftx_old,shifty_old,shiftz_old,
		      shiftx_rhs,shifty_rhs,shiftz_rhs,
		      dtshiftx_old,dtshifty_old,dtshiftz_old,
		      dtshiftx_rhs,dtshifty_rhs,dtshiftz_rhs,
		      Gammax,Gammay,Gammaz,
		      tempx,tempy,tempz,*Symmetry);
}
