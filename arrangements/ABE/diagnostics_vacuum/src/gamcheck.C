#define KRANC_C

#include <stdio.h>
#include "cctk.h"
#include <math.h>
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

#include "GenericFD.h"

extern "C" void CCTK_FCALL CCTK_FNAME(gamcheck)
  (const cGH **cctkGH, int *cctk_lsh, int *nghostzones, int *Symmetry,
   double *dx, double *dy, double *dz, 
   double *gconx, double *gcony, double *gconz,
   double *Gammax, double *Gammay, double *Gammaz, 
   double *gupxx, double *gupxy, double *gupxz, double *gupyy, double *gupyz, double *gupzz);


//-----------------------------------------------------
// compares Gammai to the divergence of gupij
//-----------------------------------------------------
extern "C" void gamcheck(const cGH *cctkGH, int *cctk_lsh, int *nghostzones, int Symmetry,
			 double dx, double dy, double dz, 
			 double *gconx, double *gcony, double *gconz,
			 double *Gammax, double *Gammay, double *Gammaz, 
			 double *gupxx, double *gupxy, double *gupxz, double *gupyy, double *gupyz, double *gupzz) {

  /* Initialise finite differencing variables.  NEED THIS FOR GenericFD.h */
#include "../../GenFD_decl_set_varCPP.h"

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
  if(Symmetry==4) {
    istart--;
    iend++;
  }

#pragma omp parallel for
  for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    gconx[index] = -(D1gf(gupxx,i,j,k) + D2gf(gupxy,i,j,k) + D3gf(gupxz,i,j,k)) - Gammax[index];
    gcony[index] = -(D1gf(gupxy,i,j,k) + D2gf(gupyy,i,j,k) + D3gf(gupyz,i,j,k)) - Gammay[index];
    gconz[index] = -(D1gf(gupxz,i,j,k) + D2gf(gupyz,i,j,k) + D3gf(gupzz,i,j,k)) - Gammaz[index];
  }
}

extern "C" void CCTK_FCALL CCTK_FNAME(gamcheck)
  (const cGH **cctkGH, int *cctk_lsh, int *nghostzones, int *Symmetry,
   double *dx, double *dy, double *dz, 
   double *gconx, double *gcony, double *gconz,
   double *Gammax, double *Gammay, double *Gammaz, 
   double *gupxx, double *gupxy, double *gupxz, double *gupyy, double *gupyz, double *gupzz)
{
  gamcheck(*cctkGH, cctk_lsh, nghostzones, *Symmetry,
	    *dx,  *dy,  *dz, 
	   gconx, gcony, gconz,
	   Gammax, Gammay, Gammaz, 
	   gupxx, gupxy, gupxz, gupyy, gupyz, gupzz);
}
