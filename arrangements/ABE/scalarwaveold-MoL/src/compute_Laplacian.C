//-----------------------------------------------------------------------------
// Laplace Operator of Function f (for interior, and inner boundaries)
//-----------------------------------------------------------------------------


#include <stdio.h>
#include "cctk.h"
#include <math.h>
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

//Following line needed to enable GenericFD.h
#define KRANC_C
#include "GenericFD.h"

extern "C" void compute_Laplacian(const cGH *cctkGH,int *cctk_lsh,int *nghostzones,int scalarwave_Symmetry,double dx,double dy,double dz,
		       double *f,double *DDf) {
  
  /* Initialise finite differencing variables */
  double dxi = 1 / dx;
  double dyi = 1 / dy;
  double dzi = 1 / dz;
  double dxi2 = dxi*dxi;
  double dyi2 = dyi*dyi;
  double dzi2 = dzi*dzi;
  double hdxi = 0.5 * dxi;
  double hdyi = 0.5 * dyi;
  double hdzi = 0.5 * dzi;

  /* Set up variables used in the grid loop for the physical grid points */
  int istart = nghostzones[0];
  int jstart = nghostzones[1];
  int kstart = nghostzones[2];
  int iend = cctk_lsh[0] - nghostzones[0];
  int jend = cctk_lsh[1] - nghostzones[1];
  int kend = cctk_lsh[2] - nghostzones[2];

  //Following lines needed since nghostzones[0] = ORDER, and 
  //   not ORDER-1 in axisymmetry 
  //   (so that rotation can be done on multiprocessor runs)
  if(scalarwave_Symmetry==4) {
    istart--;
    iend++;
  }

  // Following fills in values for DDf = [Laplacian of f] everywhere except:
  // 1) Interprocessor ghost zones.
  // 2) scalarwave_Symmetry ghost zones. (i.e., across symmetry axes)
  for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    // You'll find D11gf() stuff inside GenericFD.h
    double fxx = D11gf(f,i,j,k);;
    double fyy = D22gf(f,i,j,k);;
    double fzz = D33gf(f,i,j,k);;

    DDf[index] = fxx + fyy + fzz;
  }
}
