//-----------------------------------------------------------------------------
// Laplace Operator of Function f (for interior, and inner boundaries)
//-----------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>
#include "cctk.h"
#include <math.h>
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

//Following line declares finite differentation (FD) definitions found in arrangements/ABE/GenericFD.h
#define KRANC_C
#include "GenericFD.h"

extern "C" void compute_Laplacian(const cGH *cctkGH,int *cctk_lsh,int *nghostzones,int Symmetry,double dx,double dy,double dz,double *phi,double *phidot_rhs);

extern "C" void scalarwaveMoL_compute_rhs(CCTK_ARGUMENTS) {
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  /* Set up variables used in the grid loop for the physical grid points */
  int imin[3],imax[3];
  for(int i=0;i<3;i++) {
    imin[i] = cctk_nghostzones[i]; //FD_nghostzones;
    imax[i] = cctk_lsh[i] - cctk_nghostzones[i]; //-FD_nghostzones;
  }

  //Following lines needed since cctk_nghostzones[0] = ORDER, and 
  //   not ORDER-1 in axisymmetry 
  //   (so that rotation can be done on multiprocessor runs)
  if(scalarwave_Symmetry==4) {
    imin[0]--;
    imax[0]++;
  }

  //First compute phidot_rhs = Laplacian[phi]
  double dx = CCTK_DELTA_SPACE(0);
  double dy = CCTK_DELTA_SPACE(1);
  double dz = CCTK_DELTA_SPACE(2);
  compute_Laplacian(cctkGH,cctk_lsh,cctk_nghostzones,scalarwave_Symmetry,dx,dy,dz,phi,phidot_rhs);
  //compute_Laplacian(cctkGH,cctk_lsh,cctk_nghostzones,scalarwave_Symmetry,dx,dy,dz,phi_stagger,phidot_stagger_rhs);

  // Next compute phi_rhs = phidot everywhere except:
  // 1) Interprocessor ghost zones.
  // 2) Symmetry ghost zones. (i.e., across symmetry axes)
  //for(int k=imin[2];k<imax[2];k++) for(int j=imin[1];j<imax[1];j++) for(int i=imin[0];i<imax[0];i++) {
  for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    phi_rhs[index] = phidot[index]; 

   // phi_stagger_rhs[index] = phidot_stagger[index];


    // *** TEST ***
    phi_stagger_rhs[index] = phidot_stagger[index];
    //double xl = x[index];
    //double yl = y[index] + 0.5*dy;
    //double zl = z[index] + 0.5*dz;
    phidot_stagger_rhs[index] = 0.0;
    // ************
  }
}
