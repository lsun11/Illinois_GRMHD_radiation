/*
  Compute gupij's.  You need to do this after every regridding,
  before computing the Ricci tensor,
  and before computing BSSN RHS variables.
  It is possible that one of the above is redundant, but hey, this is a cheap operation.
*/

#include <stdio.h>
#include "cctk.h"
#include <math.h>
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

static char *rcsid="$mew. $";
CCTK_FILEVERSION(compute_h17_rhs)


  extern "C" void CCTK_FCALL CCTK_FNAME(compute_h17_rhs)
  (const cGH **cctkGH,int *cctk_lsh,
   double *X,double *Y,double *Z,
   double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,
   double *lapm1,double *phi,double *h17_rhs);


extern "C" void compute_h17_rhs(const cGH *cctkGH,int *cctk_lsh,
				double *X,double *Y,double *Z, 
				double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,
				double *lapm1,double *phi,double *h17_rhs) {
 
  DECLARE_CCTK_PARAMETERS;

#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	double xL = X[index];
	double yL = Y[index];
	double zL = Z[index];
	double gupxxL = gupxx[index];
	double gupxyL = gupxy[index];
	double gupxzL = gupxz[index];
	double gupyyL = gupyy[index];
	double gupyzL = gupyz[index];
	double gupzzL = gupzz[index];
	
	double lapm1L = lapm1[index];
	double phiL = phi[index];
	double r_coord = sqrt(xL*xL+yL*yL+zL*zL); 
	double guprrL = exp(-4.0*phiL)*(gupxxL*xL*xL+gupyyL*yL*yL+gupzzL*zL*zL+2.0*(gupxyL*xL*yL+gupxzL*xL*zL+gupyzL*yL*zL))/(r_coord*r_coord);
	double r_areal = sqrt(exp(6.0*phiL)*sqrt(guprrL))*r_coord;
	double C = 1.24672;
	h17_rhs[index] = 1-2.0/r_areal + C*C*exp(lapm1L+1.0)/pow(r_areal,4);
      }
}


extern "C" void CCTK_FCALL CCTK_FNAME(compute_h17_rhs)
  (const cGH **cctkGH,int *cctk_lsh,
   double *X,double *Y,double *Z, 
   double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,
   double *lapm1,double *phi,double *h17_rhs)
{
  compute_h17_rhs(*cctkGH,cctk_lsh,
		  X,Y,Z,
		  gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,
		  lapm1,phi,h17_rhs);
}
