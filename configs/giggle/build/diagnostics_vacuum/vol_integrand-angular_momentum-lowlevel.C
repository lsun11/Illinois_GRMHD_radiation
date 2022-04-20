#include <stdio.h>
#include "cctk.h"
#include <math.h>

#define F1o3 0.333333333333333333333333333333

#define KRANC_C
#include "GenericFD.h"

extern "C" void CCTK_FCALL jdens_cpp_
  (const cGH **cctkGH,int *cctk_lsh,int *nghostzones, int *Symmetry,
   double *Jigrnd, 
   double *dx,double *dy,double *dz, 
   double *X,double *Y,double *Z, 
   double *phi, double *Sx, double *Sy, double *trK, 
   double *Axx, double *Axy, double *Axz, double *Ayy, double *Ayz, double *Azz, 
   double *gupxx, double *gupxy, double *gupxz, double *gupyy, double *gupyz, double *gupzz);



//-----------------------------------------------------------------------------
// Compute angular momentum as in Shibata PRD60 104052
//-----------------------------------------------------------------------------
extern "C" void jdens_cpp(const cGH *cctkGH,int *cctk_lsh,int *nghostzones,int Symmetry,
			  double *Jigrnd, 
			  double dx,double dy,double dz, 
			  double *X,double *Y,double *Z, 
			  double *phi, double *Sx, double *Sy, double *trK, 
			  double *Axx, double *Axy, double *Axz, double *Ayy, double *Ayz, double *Azz, 
			  double *gupxx, double *gupxy, double *gupxz, double *gupyy, double *gupyz, double *gupzz) {

  double f1o8p = 1.0/(8.0*M_PI);

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

    double AxxL = Axx[index];
    double AxyL = Axy[index];
    double AxzL = Axz[index];
    double AyyL = Ayy[index];
    double AyzL = Ayz[index];
    double AzzL = Azz[index];
    
    double gupxxL = gupxx[index];
    double gupxyL = gupxy[index];
    double gupxzL = gupxz[index];
    double gupyyL = gupyy[index];
    double gupyzL = gupyz[index];
    double gupzzL = gupzz[index];
    
    //-----------------------------------------------------------------------------
    // Compute first derivatives of the shift and physical metric
    //-----------------------------------------------------------------------------
    double trKx = D1gf(trK,i,j,k);
    double trKy = D2gf(trK,i,j,k);
    double trKz = D3gf(trK,i,j,k);

    double gupxxx = D1gf(gupxx,i,j,k);
    double gupxxy = D2gf(gupxx,i,j,k);
    double gupxxz = D3gf(gupxx,i,j,k);

    double gupxyx = D1gf(gupxy,i,j,k);
    double gupxyy = D2gf(gupxy,i,j,k);
    double gupxyz = D3gf(gupxy,i,j,k);

    double gupxzx = D1gf(gupxz,i,j,k);
    double gupxzy = D2gf(gupxz,i,j,k);
    double gupxzz = D3gf(gupxz,i,j,k);

    double gupyyx = D1gf(gupyy,i,j,k);
    double gupyyy = D2gf(gupyy,i,j,k);
    double gupyyz = D3gf(gupyy,i,j,k);

    double gupyzx = D1gf(gupyz,i,j,k);
    double gupyzy = D2gf(gupyz,i,j,k);
    double gupyzz = D3gf(gupyz,i,j,k);

    double gupzzx = D1gf(gupzz,i,j,k);
    double gupzzy = D2gf(gupzz,i,j,k);
    double gupzzz = D3gf(gupzz,i,j,k);

    double kxy = AxyL * gupxxL + AyyL * gupxyL + AyzL * gupxzL;
    double kyx = AxxL * gupxyL + AxyL * gupyyL + AxzL * gupyzL;

    double xL = X[index];
    double yL = Y[index];
    double zL = Z[index];


    //-----------------------------------------------------------------------------
    // compute integrand
    //-----------------------------------------------------------------------------
    Jigrnd[index] = exp(6.0*phi[index]) * (
					    xL*Sy[index] - yL*Sx[index] + 
					    f1o8p*(kxy - kyx) 
					    + f1o8p*(2.0/3.0) * (xL*trKy - yL*trKx) 
					    - 0.50*f1o8p*(AxxL*(xL*gupxxy - yL*gupxxx) 
							  + 2.0*AxyL*(xL*gupxyy - yL*gupxyx) 
							  + 2.0*AxzL*(xL*gupxzy - yL*gupxzx) 
							  + AyyL*(xL*gupyyy - yL*gupyyx) 
							  + 2.0*AyzL*(xL*gupyzy - yL*gupyzx) 
							  + AzzL*(xL*gupzzy - yL*gupzzx) ) );
  }
}

extern "C" void CCTK_FCALL jdens_cpp_
  (const cGH **cctkGH,int *cctk_lsh,int *nghostzones, int *Symmetry,
   double *Jigrnd, 
   double *dx,double *dy,double *dz, 
   double *X,double *Y,double *Z, 
   double *phi, double *Sx, double *Sy, double *trK, 
   double *Axx, double *Axy, double *Axz, double *Ayy, double *Ayz, double *Azz, 
   double *gupxx, double *gupxy, double *gupxz, double *gupyy, double *gupyz, double *gupzz)
{
  jdens_cpp(*cctkGH,cctk_lsh,nghostzones, *Symmetry,
	    Jigrnd, 
	    *dx,*dy,*dz,
	    X,Y,Z, 
	    phi, Sx, Sy, trK, 
	    Axx, Axy, Axz, Ayy, Ayz, Azz, 
	    gupxx, gupxy, gupxz, gupyy, gupyz, gupzz);
}
