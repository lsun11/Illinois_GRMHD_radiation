/*
  Differences between this function and Derivs():
  1) Derivs() is strictly second-order; Derivs_interior is the same as your spatial order
  2) Derivs() computes derivatives everywhere, Derivs_interior compute derivatives everywhere but the ghostzones
*/

#define KRANC_C
#include <stdio.h>
#include "cctk.h"
#include <math.h>
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

#include "GenericFD.h"


extern "C" void CCTK_FCALL CCTK_FNAME(Derivs_interior)
     (const cGH **cctkGH, double *dx, double *dy, double *dz,
      int *nghostzones,int *cctk_lsh, int *Symmetry,
      double *f, double *fx, double *fy, double *fz);

extern "C" void Derivs_interior(const cGH *cctkGH, double dx, double dy, double dz,
				int *nghostzones,int *cctk_lsh, int Symmetry,
				double *f, double *fx, double *fy, double *fz) {
  
  //  DECLARE_CCTK_PARAMETERS;

  /* Initialise finite differencing variables.  NEED THIS FOR GenericFD.h */
#include "../../GenFD_decl_set_varCPP.h"

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

  /* Set up variables used in the grid loop for the physical grid points */
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
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

	fx[index] = D1gf(f,i,j,k);
	fy[index] = D2gf(f,i,j,k);
	fz[index] = D3gf(f,i,j,k);
      }
}

extern "C" void CCTK_FCALL CCTK_FNAME(Derivs_interior)
     (const cGH **cctkGH, double *dx, double *dy, double *dz,
      int *nghostzones,int *cctk_lsh, int *Symmetry,
      double *f, double *fx, double *fy, double *fz)
{
  Derivs_interior(*cctkGH, *dx,  *dy,  *dz,
		  nghostzones, cctk_lsh, *Symmetry,
		  f,fx,fy,fz);
}
