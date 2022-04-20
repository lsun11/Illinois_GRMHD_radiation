//-----------------------------------------------------------------------------
//
// Upwinding advection term for first order shift:
// \partial_t beta^i += [ \beta^j \partial_j \beta^i ]
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
//#define FD_SET_BY_USER
//#define FD_C2

#include "GenericFD.h"

extern "C" void CCTK_FCALL puncture_firstorder_shift_upwind_
  (const cGH **cctkGH,int *nghostzones,int *cctk_lsh, int *Symmetry,
   double *dx,double *dy,double *dz,
   double *betax,double *betay,double *betaz,
   double *betax_rhs,double *betay_rhs,double *betaz_rhs);

extern "C" void puncture_firstorder_shift_upwind(const cGH *cctkGH,int *nghostzones,int *cctk_lsh, int Symmetry,
					     double dx,double dy,double dz,
					     double *betax,double *betay,double *betaz,
					     double *betax_rhs,double *betay_rhs,double *betaz_rhs) {
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
	double betaxL = betax[index];
	double betayL = betay[index];
	double betazL = betaz[index];

	double betax_rhsL = betax_rhs[index];
	double betay_rhsL = betay_rhs[index];
	double betaz_rhsL = betaz_rhs[index];

	if(betaxL > 0.0) {
	  betax_rhsL += betaxL*D1_up_gt_gf(betax,i,j,k);
	  betay_rhsL += betaxL*D1_up_gt_gf(betay,i,j,k);
	  betaz_rhsL += betaxL*D1_up_gt_gf(betaz,i,j,k);
	} else {
	  betax_rhsL += betaxL*D1_up_lt_gf(betax,i,j,k);
	  betay_rhsL += betaxL*D1_up_lt_gf(betay,i,j,k);
	  betaz_rhsL += betaxL*D1_up_lt_gf(betaz,i,j,k);
	}

	if(betayL > 0.0) {
	  betax_rhsL += betayL*D2_up_gt_gf(betax,i,j,k);
	  betay_rhsL += betayL*D2_up_gt_gf(betay,i,j,k);
	  betaz_rhsL += betayL*D2_up_gt_gf(betaz,i,j,k);
	} else {
	  betax_rhsL += betayL*D2_up_lt_gf(betax,i,j,k);
	  betay_rhsL += betayL*D2_up_lt_gf(betay,i,j,k);
	  betaz_rhsL += betayL*D2_up_lt_gf(betaz,i,j,k);
	}

	if(betazL > 0.0) {
	  betax_rhsL += betazL*D3_up_gt_gf(betax,i,j,k);
	  betay_rhsL += betazL*D3_up_gt_gf(betay,i,j,k);
	  betaz_rhsL += betazL*D3_up_gt_gf(betaz,i,j,k);
	} else {
	  betax_rhsL += betazL*D3_up_lt_gf(betax,i,j,k);
	  betay_rhsL += betazL*D3_up_lt_gf(betay,i,j,k);
	  betaz_rhsL += betazL*D3_up_lt_gf(betaz,i,j,k);
	}
	betax_rhs[index] = betax_rhsL;
	betay_rhs[index] = betay_rhsL;
	betaz_rhs[index] = betaz_rhsL;

      }
}



extern "C" void CCTK_FCALL puncture_firstorder_shift_upwind_
  (const cGH **cctkGH,int *nghostzones,int *cctk_lsh, int *Symmetry,
   double *dx,double *dy,double *dz,
   double *betax,double *betay,double *betaz,
   double *betax_rhs,double *betay_rhs,double *betaz_rhs)
{
  puncture_firstorder_shift_upwind(*cctkGH,nghostzones,cctk_lsh, *Symmetry,
				   *dx, *dy, *dz,
				   betax, betay, betaz,
				   betax_rhs, betay_rhs, betaz_rhs);
}
