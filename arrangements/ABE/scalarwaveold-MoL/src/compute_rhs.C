//-----------------------------------------------------------------------------
// Laplace Operator of Function f (for interior, and inner boundaries)
//-----------------------------------------------------------------------------

#include <stdio.h>
#include "cctk.h"
#include <math.h>
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

extern "C" void compute_Laplacian(const cGH *cctkGH,int *cctk_lsh,int *nghostzones,int Symmetry,double dx,double dy,double dz,double *phi,double *phidot_rhs);


extern "C" void scalarwaveMoL_compute_rhs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  double dx = CCTK_DELTA_SPACE(0);
  double dy = CCTK_DELTA_SPACE(1);
  double dz = CCTK_DELTA_SPACE(2);
  
  /* Set up variables used in the grid loop for the physical grid points */
  int istart = cctk_nghostzones[0];
  int jstart = cctk_nghostzones[1];
  int kstart = cctk_nghostzones[2];
  int iend = cctk_lsh[0] - cctk_nghostzones[0];
  int jend = cctk_lsh[1] - cctk_nghostzones[1];
  int kend = cctk_lsh[2] - cctk_nghostzones[2];

  printf("HI! Updating rhs's!  dt = %e\n",CCTK_DELTA_TIME);

  //Following lines needed since cctk_nghostzones[0] = ORDER, and 
  //   not ORDER-1 in axisymmetry 
  //   (so that rotation can be done on multiprocessor runs)
  if(scalarwave_Symmetry==4) {
    istart--;
    iend++;
  }

  //First compute phidot_rhs = Laplacian[phi]
  compute_Laplacian(cctkGH,cctk_lsh,cctk_nghostzones,Symmetry,dx,dy,dz,phi,phidot_rhs);

  // Next compute phi_rhs = phidot everywhere except:
  // 1) Interprocessor ghost zones.
  // 2) Symmetry ghost zones. (i.e., across symmetry axes)
  for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

    phi_rhs[index] = phidot[index];
  }
}
