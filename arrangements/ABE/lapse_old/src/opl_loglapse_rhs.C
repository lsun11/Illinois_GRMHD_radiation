//-----------------------------------------------------------------------------
// Compute 1+log lapse RHS: \partial_t log \alpha = -2 \alpha trK
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

extern "C" void CCTK_FCALL CCTK_FNAME(opl_loglapse_rhs)
  (const cGH **cctkGH,int *nghostzones,int *cctk_lsh,int *Symmetry,
   double *lapse_rhs, double *trK);

extern "C" void opl_loglapse_rhs(const cGH *cctkGH,int *nghostzones,int *cctk_lsh,int Symmetry,
				 double *lapse_rhs, double *trK) {
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
#pragma omp parallel for
  for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    double lapse_rhsL = lapse_rhs[index];
    double trKL = trK[index];

    lapse_rhsL = - 2.0 * trKL;

    lapse_rhs[index] = lapse_rhsL;
  }
}

extern "C" void CCTK_FCALL CCTK_FNAME(opl_loglapse_rhs)
  (const cGH **cctkGH,int *nghostzones,int *cctk_lsh, int *Symmetry,
   double *lapse_rhs, double *trK)
{
  opl_loglapse_rhs(*cctkGH,nghostzones,cctk_lsh,*Symmetry,
		   lapse_rhs,trK);
}
