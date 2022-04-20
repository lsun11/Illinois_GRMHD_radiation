//-----------------------------------------------------------------------------
// Compute 1+log lapse RHS: \partial_t \alpha = -2 \alpha trK
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

extern "C" void CCTK_FCALL CCTK_FNAME(opl_lapse_rhs)
  (const cGH **cctkGH,int *nghostzones,int *cctk_lsh,int *Symmetry,int *alap,
   double *r,double *PhysicalRadius,double *RadiusDerivative,
   double *lapse_old, double *lapse_rhs,double *phi, double *trK,double *dx,double *dy,double *dz);

extern "C" void opl_lapse_rhs(const cGH *cctkGH,int *nghostzones,int *cctk_lsh,int Symmetry,int alap,
			      double *r,double *PhysicalRadius,double *RadiusDerivative,
			      double *lapse_old, double *lapse_rhs, 
			      double *phi,double *trK,double dx,double dy,double dz) {
  
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

  if (alap==1) {
    printf("YOU CRAZY!\n"); exit(1);
  }

#pragma omp parallel for
  for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    double lapse_rhsL = lapse_rhs[index];
    double lapse_oldL = lapse_old[index];
    double trKL = trK[index];

    //if (alap!=1) {
    lapse_rhsL = - 2.0 * ( (lapse_oldL + 1.0) * trKL);
    /*
      } else {
      double phiL = phi[index];
      double RadiusDerivativeL = RadiusDerivative[index];
      double rL=r[index];
      double PhysicalRadiusL = PhysicalRadius[index];
      lapse_rhsL=pow(RadiusDerivativeL,(alap/6.0)) * 
      pow(PhysicalRadiusL/rL,(alap/3.0));
      lapse_rhsL=- 2.0 *lapse_rhsL*((lapse_oldL+1.0)*trKL)*exp(-1.0*alap*phiL);
      }
    */
    
    //ZACH SAYS: TESTING ONLY:
    //lapse_rhsL -= 0.2 * lapse_oldL;
    lapse_rhsL += 0.002 * ( D11gf(lapse_old,i,j,k) + D22gf(lapse_old,i,j,k) + D33gf(lapse_old,i,j,k)); // - 0.2*(lapse_oldL+1.0)*exp(-2.0*phi[index]);
    
    lapse_rhs[index] = lapse_rhsL;
  }
}

extern "C" void CCTK_FCALL CCTK_FNAME(opl_lapse_rhs)
  (const cGH **cctkGH,int *nghostzones,int *cctk_lsh, int *Symmetry, int *alap,
   double *r,double *PhysicalRadius,double *RadiusDerivative,
   double *lapse_old, double *lapse_rhs,double *phi, double *trK,double *dx,double *dy,double *dz)
{
  opl_lapse_rhs(*cctkGH,nghostzones,cctk_lsh,*Symmetry,*alap,
		r,PhysicalRadius,RadiusDerivative,
		lapse_old, lapse_rhs, phi, trK,*dx,*dy,*dz);
}
