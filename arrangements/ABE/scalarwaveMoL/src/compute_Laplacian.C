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

extern "C" void compute_Laplacian(const cGH *cctkGH,int *cctk_lsh,int *nghostzones,int scalarwave_Symmetry,double dx,double dy,double dz,
				  double *f,double *DDf) {

  int imin[3],imax[3];
  for(int i=0;i<3;i++) {
    imin[i] = nghostzones[i]; //FD_nghostzones;
    imax[i] = cctk_lsh[i] - nghostzones[i]; //-FD_nghostzones;
  }

  /* Initialise finite differencing variables.  NEED THIS FOR GenericFD.h */
#include "../../GenFD_decl_set_varCPP.h"

  //Following lines needed since nghostzones[0] = ORDER, and 
  //   not ORDER-1 in axisymmetry 
  //   (so that rotation can be done on multiprocessor runs)
  if(scalarwave_Symmetry==4) {
    imin[0]--;
    imax[0]++;
  }

  // Following fills in values for DDf = [Laplacian of f] everywhere, as described in arXiv:0706.0740v1, pgs 4-5 (Jena group)
  // ... also described in arXiv:0711.1165, pg 4 (RIT group)
#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

    int local_spatial_order=4;
    if       (i==0 || j==0 || k==0 || i==cctk_lsh[0]-1 || j==cctk_lsh[1]-1 || k==cctk_lsh[2]-1) {
      local_spatial_order=0;
    } else if(i==1 || j==1 || k==1 || i==cctk_lsh[0]-2 || j==cctk_lsh[1]-2 || k==cctk_lsh[2]-2) {
      local_spatial_order=2;
    }

    double fxx,fyy,fzz;

    if(local_spatial_order==4) {
      // You'll find D11gf() stuff inside GenericFD.h
      fxx = D11_c4(f,i,j,k);
      fyy = D22_c4(f,i,j,k);
      fzz = D33_c4(f,i,j,k);
    } else if(local_spatial_order==2) {
      fxx = D11_c2(f,i,j,k);
      fyy = D22_c2(f,i,j,k);
      fzz = D33_c2(f,i,j,k);
    } else {
      fxx = 0.0;
      fyy = 0.0;
      fzz = 0.0;
    }

    DDf[index] = fxx + fyy + fzz;
  }
  
  /* OLD WAY OF DOING IT:
    for(int k=imin[2];k<imax[2];k++) for(int j=imin[1];j<imax[1];j++) for(int i=imin[0];i<imax[0];i++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

    int local_spatial_order=4;
    if       (i==0 || j==0 || k==0 || i==cctk_lsh[0]-1 || j==cctk_lsh[1]-1 || k==cctk_lsh[2]-1) {
    local_spatial_order=0;
    } else if(i==1 || j==1 || k==1 || i==cctk_lsh[0]-2 || j==cctk_lsh[1]-2 || k==cctk_lsh[2]-2) {
    local_spatial_order=2;
    }

    double fxx,fyy,fzz;

    if(local_spatial_order==4) {
    // You'll find D11gf() stuff inside GenericFD.h
    fxx = D11_c4(f,i,j,k);
    fyy = D22_c4(f,i,j,k);
    fzz = D33_c4(f,i,j,k);
    } else if(local_spatial_order==2) {
    fxx = D11_c2(f,i,j,k);
    fyy = D22_c2(f,i,j,k);
    fzz = D33_c2(f,i,j,k);
    } else {
    fxx = 0.0;
    fyy = 0.0;
    fzz = 0.0;
    }

    DDf[index] = fxx + fyy + fzz;
    }
  */
}
