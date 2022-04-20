/* 
Compute \partial_x g_{\mu \nu} and calculate the boosted metric
*/

#define KRANC_C

/* Define macros used in calculations */ #define INITVALUE  (42)
#define INV(x) ((1) / (x)) 
#define SQR(x) ((x) * (x))

#include <stdio.h> 
#include "cctk.h"
#include <math.h>
#include "cctk_Arguments.h" 
#include "cctk_Parameters.h"
#include "Symmetry.h"
#include "GenericFD.h"

static char *rcsid="$mew. $";

CCTK_FILEVERSION(metric_derivatives_and_boost)

  extern "C" void CCTK_FCALL CCTK_FNAME(metric_derivatives_and_boost)
  (const cGH **cctkGH,double *dx, double *dy, double *dz,
   int *nghostzones, int *cctk_lsh, double *phi,  double *lapm1, double *shiftx,
   double *shifty, double *shiftz, double *gtt, double *gtx, double *gty,
   double *gtz, double *gxx, double *gxy, double *gxz, double *gyy, 
   double *gyz, double *gzz, double *gttx, double *gtxx, double *gtyx, 
   double *gtzx, double *gxxx, double *gxyx, double *gxzx,
   double *gyyx, double *gyzx, double *gzzx, double *vx_boost); 

  extern "C" void metric_derivatives_and_boost(const cGH *cctkGH,
   double dx, double dy, double dz,
   int *nghostzones, int *cctk_lsh, double *phi, double *lapm1, double *shiftx,
   double *shifty, double *shiftz, double *gtt, double *gtx, double *gty,
   double *gtz, double *gxx, double *gxy, double *gxz, double *gyy,
   double *gyz, double *gzz, double *gttx, double *gtxx, double *gtyx,
   double *gtzx, double *gxxx, double *gxyx, double *gxzx,
   double *gyyx, double *gyzx, double *gzzx, double vx_boost){ 

  // DECLARE_CCTK_PARAMETERS;

  /* Initialise finite differencing variables.  NEED THIS FOR GenericFD.h */
#include "../../GenFD_decl_set_varCPP.h"

  /* Set up variables used in the grid loop for the physical grid points */
  int istart = nghostzones[0];
  int jstart = nghostzones[1];
  int kstart = nghostzones[2];
  int iend = cctk_lsh[0] - nghostzones[0];
  int jend = cctk_lsh[1] - nghostzones[1];
  int kend = cctk_lsh[2] - nghostzones[2];

  double lorentz_gamma = 1.0/sqrt(1.0-SQR(vx_boost));

// Setup gtt and gti
#pragma omp for
    for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
      int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
      double psi4L = exp(4.0*phi[index]);
      double gxxL = gxx[index]*psi4L;
      double gxyL = gxy[index]*psi4L;
      double gxzL = gxz[index]*psi4L;
      double gyyL = gyy[index]*psi4L;
      double gyzL = gyz[index]*psi4L;
      double gzzL = gzz[index]*psi4L;
      double betax = shiftx[index];
      double betay = shifty[index];
      double betaz = shiftz[index];
      double beta_x = gxxL*betax + gxyL*betay + gxzL*betaz;
      double beta_y = gxyL*betax + gyyL*betay + gyzL*betaz;
      double beta_z = gxzL*betax + gyzL*betay + gzzL*betaz;
      double beta2 = betax*beta_x + betay*beta_y + betaz*beta_z;
      gtt[index] = beta2 - SQR(1.0+lapm1[index]);
      gtx[index] = beta_x;
      gty[index] = beta_y;
      gtz[index] = beta_z;
      gxx[index] = gxxL;
      gxy[index] = gxyL;
      gxz[index] = gxzL;
      gyy[index] = gyyL;
      gyz[index] = gyzL;
      gzz[index] = gzzL;
  }

// Compute derivatives. These are used to compute the time derivative of 
// the boosted metric.
#pragma omp for
    for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
      int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
      // The 1/lorentz_gamma comes from the fact that 
      // Delta x = lorentz_gamma * Delta x_boost 
      gttx[index] = D1gf(gtt,i,j,k)/lorentz_gamma;
      gtxx[index] = D1gf(gtx,i,j,k)/lorentz_gamma;
      gtyx[index] = D1gf(gty,i,j,k)/lorentz_gamma;
      gtzx[index] = D1gf(gtz,i,j,k)/lorentz_gamma;
      gxxx[index] = D1gf(gxx,i,j,k)/lorentz_gamma;
      gxyx[index] = D1gf(gxy,i,j,k)/lorentz_gamma;
      gxzx[index] = D1gf(gxz,i,j,k)/lorentz_gamma;
      gyyx[index] = D1gf(gyy,i,j,k)/lorentz_gamma;
      gyzx[index] = D1gf(gyz,i,j,k)/lorentz_gamma;
      gzzx[index] = D1gf(gzz,i,j,k)/lorentz_gamma;

// *** TEST ***
if (i==44 && j==29 && k==3 && dx < 0.026) {
  printf("gxxx, lorentz_gamma: %e %e\n",gxxx[index],lorentz_gamma);
  printf("gxx-, gxx0, gxx+: %e %e %e\n",gxx[CCTK_GFINDEX3D(cctkGH,i-1,j,k)],gxx[index],gxx[CCTK_GFINDEX3D(cctkGH,i+1,j,k)]);
}
// ***************
    }

// Boost the metric
#pragma omp for
    for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
      int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
      double gttL = gtt[index];
      double gtxL = gtx[index];
      double gtyL = gty[index];
      double gtzL = gtz[index];
      double gxxL = gxx[index];
      double gxyL = gxy[index];
      double gxzL = gxz[index];
      double gyyL = gyy[index];
      double gyzL = gyz[index];
      double gzzL = gzz[index];
      double alphaL = lapm1[index] + 1.0;
      double f1oalp2 = 1.0/(alphaL*alphaL);
      double betaxL = shiftx[index];
      double betayL = shifty[index];
      double betazL = shiftz[index];
      double gijdet = gxxL * gyyL * gzzL + gxyL * gyzL * gxzL 
                      + gxzL * gxyL * gyzL - gxzL * gyyL * gxzL 
                      - gxyL * gxyL * gzzL - gxxL * gyzL * gyzL;
      double gijdetinv = 1.0/gijdet;
      double gupttL = -f1oalp2;
      double guptxL = betaxL*f1oalp2;
      double guptyL = betayL*f1oalp2;
      double guptzL = betazL*f1oalp2;
      double gupxxL =   ( gyyL * gzzL - gyzL * gyzL )* gijdetinv - betaxL*betaxL*f1oalp2;
      double gupxyL = - ( gxyL * gzzL - gyzL * gxzL )* gijdetinv - betaxL*betayL*f1oalp2;
      double gupxzL =   ( gxyL * gyzL - gyyL * gxzL )* gijdetinv - betaxL*betazL*f1oalp2;
      double gupyyL =   ( gxxL * gzzL - gxzL * gxzL )* gijdetinv - betayL*betayL*f1oalp2;
      double gupyzL = - ( gxxL * gyzL - gxyL * gxzL )* gijdetinv - betayL*betazL*f1oalp2;
      double gupzzL =   ( gxxL * gyyL - gxyL * gxyL )* gijdetinv - betazL*betazL*f1oalp2;
      // compute L^a_b = \partial x_boost^{a}/\partial x^b
      double Ltt = lorentz_gamma;
      double Ltx = lorentz_gamma*vx_boost;
      double Lxt = lorentz_gamma*vx_boost;
      double Lxx = lorentz_gamma;

      // Boost g^{ta} and then calculate the new lapse and shift
      double guptt_boost = Ltt*Ltt*gupttL + 2.0*Ltt*Ltx*guptxL 
			   + Lxx*Lxx*gupxxL;
      double guptx_boost = Ltt*Lxt*gupttL + (Ltt*Lxx + Ltx*Lxt)*guptxL 
			   + Ltx*Lxx*gupxxL;
      double gupty_boost = Ltt*guptyL + Ltx*gupxyL;
      double guptz_boost = Ltt*guptzL + Ltx*gupxzL;
      double alpha_boost = 1.0/sqrt(-guptt_boost);
      double betax_boost = -guptx_boost/guptt_boost;
      double betay_boost = -gupty_boost/guptt_boost;
      double betaz_boost = -guptz_boost/guptt_boost;

      // Boost g_{ij}
      Ltx = -lorentz_gamma*vx_boost;
      Lxx = lorentz_gamma;
      double gxx_boost = Ltx*Ltx*gttL + 2.0*Ltx*Lxx*gtxL + Lxx*Lxx*gxxL;
      double gxy_boost = Ltx*gtyL + Lxx*gxyL;
      double gxz_boost = Ltx*gtzL + Lxx*gxzL;
      double gyy_boost = gyyL;
      double gyz_boost = gyzL;
      double gzz_boost = gzzL;

      lapm1[index] = alpha_boost - 1.0;
      shiftx[index] = betax_boost;
      shifty[index] = betay_boost;
      shiftz[index] = betaz_boost;
      // Note that gij store the boosted physical 3-metric for the moment.
      // They will be used to compute the Christoffel symbol.
      gxx[index] = gxx_boost;
      gxy[index] = gxy_boost;
      gxz[index] = gxz_boost;
      gyy[index] = gyy_boost;
      gyz[index] = gyz_boost;
      gzz[index] = gzz_boost;
    }
}

extern "C" void CCTK_FCALL CCTK_FNAME(metric_derivatives_and_boost)
(const cGH **cctkGH,double *dx, double *dy, double *dz,
   int *nghostzones, int *cctk_lsh, double *phi,  double *lapm1, double *shiftx,
   double *shifty, double *shiftz, double *gtt, double *gtx, double *gty,
   double *gtz, double *gxx, double *gxy, double *gxz, double *gyy,
   double *gyz, double *gzz, double *gttx, double *gtxx, double *gtyx,
   double *gtzx, double *gxxx, double *gxyx, double *gxzx,
   double *gyyx, double *gyzx, double *gzzx, double *vx_boost) 
{
  metric_derivatives_and_boost(*cctkGH, *dx, *dy, *dz,
   nghostzones, cctk_lsh, phi,  lapm1, shiftx,
   shifty, shiftz, gtt, gtx, gty, gtz, gxx, gxy, gxz, gyy,
   gyz, gzz, gttx, gtxx, gtyx, gtzx, gxxx, gxyx, gxzx,
   gyyx, gyzx, gzzx, *vx_boost);
}
