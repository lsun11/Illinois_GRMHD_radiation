/*
  compute Ricci tensor and Hamiltonian+momentum constraints
  Optimized from FORTRAN version by Zach Etienne: Apr, 2007
  (3x speed increase with better scalability)
  Speed increase due primarily to
  1) Fewer calls to main memory 
  -> can do calculations entirely in (uber-fast) processor cache!
  2) Replace derivative function calls (eww) 
  with "#define" finite differencing (yay).
  3) No more allocating entire temporary grid functions!
*/

#define KRANC_C
#define FD_SET_BY_USER
#define FD_C4

/* Define macros used in calculations */
#define INITVALUE  (42)
#define INV(x) ((1) / (x))
#define SQR(x) ((x) * (x))
#define CUB(x) ((x) * (x) * (x))
#define QAD(x) ((x) * (x) * (x) * (x))

#define F2o3 0.666666666666666666666666666666666

#include <stdio.h>
#include "cctk.h"
#include <math.h>
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Symmetry.h"

#include "GenericFD.h"

static char *rcsid="$mew. $";
CCTK_FILEVERSION(spatial_metric_derivs_4)


  extern "C" void CCTK_FCALL spatial_metric_derivs_4_
  (const cGH **cctkGH,double *dT, double *dx, double *dy, double *dz,
   int *nghostzones,int *cctk_lsh,
   double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
   double *lapm1,double *shiftx,double *shifty,double *shiftz,
   double *lapsex,double *lapsey, double *lapsez,
   double *shiftxx, double *shiftxy, double *shiftxz,
   double *shiftyx, double *shiftyy, double *shiftyz,
   double *shiftzx, double *shiftzy, double *shiftzz,
   double *gxxx, double *gxxy, double *gxxz,
   double *gxyx, double *gxyy, double *gxyz,
   double *gxzx, double *gxzy, double *gxzz,
   double *gyyx, double *gyyy, double *gyyz,
   double *gyzx, double *gyzy, double *gyzz,
   double *gzzx, double *gzzy, double *gzzz);
  
extern "C" void spatial_metric_derivs_4(const cGH *cctkGH, double dT, double dx, double dy, double dz,
					int *nghostzones,int *cctk_lsh,
					double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
					double *lapm1,double *shiftx,double *shifty,double *shiftz,
					double *lapsex,double *lapsey, double *lapsez,
					double *shiftxx, double *shiftxy, double *shiftxz,
					double *shiftyx, double *shiftyy, double *shiftyz,
					double *shiftzx, double *shiftzy, double *shiftzz,
					double *gxxx, double *gxxy, double *gxxz,
					double *gxyx, double *gxyy, double *gxyz,
					double *gxzx, double *gxzy, double *gxzz,
					double *gyyx, double *gyyy, double *gyyz,
					double *gyzx, double *gyzy, double *gyzz,
					double *gzzx, double *gzzy, double *gzzz) {
  
  DECLARE_CCTK_PARAMETERS

  /* Initialise finite differencing variables.  NEED THIS FOR GenericFD.h */
#include "../../GenFD_decl_set_varCPP.h"

  /* Set up variables used in the grid loop for the physical grid points */
  int istart = 2;
  int jstart = 2;
  int kstart = 2;
  int iend = cctk_lsh[0]-2;
  int jend = cctk_lsh[1]-2;
  int kend = cctk_lsh[2]-2;

  //Following lines needed since nghostzones[0] = ORDER, and
  //   not ORDER-1 in axisymmetry
  //   (so that rotation can be done on multiprocessor runs)
  if(Symmetry==4) {
    istart--;
    iend++;
  }

  long w;
  //#pragma omp parallel for  private(w)
#pragma omp parallel 
  {
#pragma omp for
    for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
      if(i==2 || j==2 || k==2 || i==cctk_lsh[0]-3 || j==cctk_lsh[1]-3 || k==cctk_lsh[2]-3) {
	// The following if() statement is quite complex and takes a while to compute, 
	//     so we don't want to evaluate it for all i,j,k!  Thus we have the above if()
	//     statement to reduce the number of evaluations.
	if((i==2 && j>2 && k>2 && i<cctk_lsh[0]-3 && j<cctk_lsh[1]-3 && k<cctk_lsh[2]-3) || 
	   (j==2 && i>2 && k>2 && i<cctk_lsh[0]-3 && j<cctk_lsh[1]-3 && k<cctk_lsh[2]-3) || 
	   (k==2 && i>2 && j>2 && i<cctk_lsh[0]-3 && j<cctk_lsh[1]-3 && k<cctk_lsh[2]-3) || 
	   (i==cctk_lsh[0]-3 && i>2 && j>2 && k>2 && j<cctk_lsh[1]-3 && k<cctk_lsh[2]-3) || 
	   (j==cctk_lsh[1]-3 && i>2 && j>2 && k>2 && i<cctk_lsh[0]-3 && k<cctk_lsh[2]-3) || 
	   (k==cctk_lsh[2]-3 && i>2 && j>2 && k>2 && i<cctk_lsh[0]-3 && j<cctk_lsh[1]-3)) {
	  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

	  double gxxL = gxx[index];
	  double gxyL = gxy[index];
	  double gxzL = gxz[index];
	  double gyyL = gyy[index];
	  double gyzL = gyz[index];
	  double gzzL = gzz[index];

	  double gxxxL = D1gf(gxx,i,j,k);
	  double gxyxL = D1gf(gxy,i,j,k);
	  double gxzxL = D1gf(gxz,i,j,k);
	  double gyyxL = D1gf(gyy,i,j,k);
	  double gyzxL = D1gf(gyz,i,j,k);
	  double gzzxL = D1gf(gzz,i,j,k);
	  double lapsexL = D1gf(lapm1,i,j,k);
	  double shiftxxL = D1gf(shiftx,i,j,k);
	  double shiftyxL = D1gf(shifty,i,j,k);
	  double shiftzxL = D1gf(shiftz,i,j,k);

	  double gxxyL = D2gf(gxx,i,j,k);
	  double gxyyL = D2gf(gxy,i,j,k);
	  double gxzyL = D2gf(gxz,i,j,k);
	  double gyyyL = D2gf(gyy,i,j,k);
	  double gyzyL = D2gf(gyz,i,j,k);
	  double gzzyL = D2gf(gzz,i,j,k);
	  double lapseyL = D2gf(lapm1,i,j,k);
	  double shiftxyL = D2gf(shiftx,i,j,k);
	  double shiftyyL = D2gf(shifty,i,j,k);
	  double shiftzyL = D2gf(shiftz,i,j,k);
	  
	  double gxxzL = D3gf(gxx,i,j,k);
	  double gxyzL = D3gf(gxy,i,j,k);
	  double gxzzL = D3gf(gxz,i,j,k);
	  double gyyzL = D3gf(gyy,i,j,k);
	  double gyzzL = D3gf(gyz,i,j,k);
	  double gzzzL = D3gf(gzz,i,j,k);
	  double lapsezL = D3gf(lapm1,i,j,k);
	  double shiftxzL = D3gf(shiftx,i,j,k);
	  double shiftyzL = D3gf(shifty,i,j,k);
	  double shiftzzL = D3gf(shiftz,i,j,k);
	  

	  gxxx[index] = gxxxL;
	  gxyx[index] = gxyxL;
	  gxzx[index] = gxzxL;
	  gyyx[index] = gyyxL;
	  gyzx[index] = gyzxL;
	  gzzx[index] = gzzxL;

	  gxxy[index] = gxxyL;
	  gxyy[index] = gxyyL;
	  gxzy[index] = gxzyL;
	  gyyy[index] = gyyyL;
	  gyzy[index] = gyzyL;
	  gzzy[index] = gzzyL;

	  gxxz[index] = gxxzL;
	  gxyz[index] = gxyzL;
	  gxzz[index] = gxzzL;
	  gyyz[index] = gyyzL;
	  gyzz[index] = gyzzL;
	  gzzz[index] = gzzzL;
	  
	  lapsex[index] = lapsexL;
	  lapsey[index] = lapseyL;
	  lapsez[index] = lapsezL;
	  
	  shiftxx[index] = shiftxxL;
	  shiftxy[index] = shiftxyL;
	  shiftxz[index] = shiftxzL;
	  shiftyx[index] = shiftyxL;
	  shiftyy[index] = shiftyyL;
	  shiftyz[index] = shiftyzL;
	  shiftzx[index] = shiftzxL;
	  shiftzy[index] = shiftzyL;
	  shiftzz[index] = shiftzzL;
	  
	}
      }
    }
  }
}



extern "C" void CCTK_FCALL spatial_metric_derivs_4_
  (const cGH **cctkGH,double *dT, double *dx, double *dy, double *dz,
   int *nghostzones,int *cctk_lsh,
   double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
   double *lapm1,double *shiftx,double *shifty,double *shiftz,
   double *lapsex,double *lapsey, double *lapsez,
   double *shiftxx, double *shiftxy, double *shiftxz,
   double *shiftyx, double *shiftyy, double *shiftyz,
   double *shiftzx, double *shiftzy, double *shiftzz,
   double *gxxx, double *gxxy, double *gxxz,
   double *gxyx, double *gxyy, double *gxyz,
   double *gxzx, double *gxzy, double *gxzz,
   double *gyyx, double *gyyy, double *gyyz,
   double *gyzx, double *gyzy, double *gyzz,
   double *gzzx, double *gzzy, double *gzzz)
{
  spatial_metric_derivs_4(*cctkGH,  *dT,  *dx,  *dy,  *dz,
			  nghostzones, cctk_lsh,
			  gxx, gxy, gxz, gyy, gyz, gzz,
			  lapm1,shiftx,shifty,shiftz,
			  lapsex,lapsey, lapsez,
			  shiftxx, shiftxy, shiftxz,
			  shiftyx, shiftyy, shiftyz,
			  shiftzx, shiftzy, shiftzz,
			  gxxx, gxxy, gxxz,
			  gxyx, gxyy, gxyz,
			  gxzx, gxzy, gxzz,
			  gyyx, gyyy, gyyz,
			  gyzx, gyzy, gyzz,
			  gzzx, gzzy, gzzz);
}
