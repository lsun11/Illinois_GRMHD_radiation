/* 
Compute the boosted Aij
*/

#define KRANC_C

/* Define macros used in calculations */ #define INITVALUE  (42)
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

CCTK_FILEVERSION(compute_Aij_boost)

  extern "C" void CCTK_FCALL CCTK_FNAME(compute_Aij_boost)
  (const cGH **cctkGH,double *dx, double *dy, double *dz,
   double *x, double *y, double *z,
   int *nghostzones, int *cctk_lsh, double *lapm1, double *phi,
   double *shiftx, double *shifty, double *shiftz,
   double *gxx, double *gxy, double *gxz, double *gyy, 
   double *gyz, double *gzz, 
   double *gupxx, double *gupxy, double *gupxz, double *gupyy,
   double *gupyz, double *gupzz,
   double *gttx0, double *gtxx0, 
   double *gtyx0, double *gtzx0, double *gxxx0, double *gxyx0, 
   double *gxzx0, double *gyyx0, double *gyzx0, double *gzzx0,
   //double *Gammaxxx, double *Gammaxxy, double *Gammaxxz, double *Gammaxyy, double *Gammaxyz, double *Gammaxzz,    
   //double *Gammayxx, double *Gammayxy, double *Gammayxz, double *Gammayyy, double *Gammayyz, double *Gammayzz,    
   //double *Gammazxx, double *Gammazxy, double *Gammazxz, double *Gammazyy, double *Gammazyz, double *Gammazzz, 
   double *vx_boost, int *coordinate_type, double *sam,
   double *Axx, double *Axy, double *Axz, double *Ayy, double *Ayz, 
   double *Azz, double *trK);

  extern "C" void compute_Aij_boost(const cGH *cctkGH,
   double dx, double dy, double dz, 
   double *x, double *y, double *z,
   int *nghostzones, int *cctk_lsh, double *lapm1, double *phi,
   double *shiftx, double *shifty, double *shiftz,
   double *gxx, double *gxy, double *gxz, double *gyy,
   double *gyz, double *gzz, 
   double *gupxx, double *gupxy, double *gupxz, double *gupyy,
   double *gupyz, double *gupzz,
   double *gttx0, double *gtxx0,
   double *gtyx0, double *gtzx0, double *gxxx0, double *gxyx0,
   double *gxzx0, double *gyyx0, double *gyzx0, double *gzzx0,
   //double *Gammaxxx, double *Gammaxxy, double *Gammaxxz, double *Gammaxyy, double *Gammaxyz, double *Gammaxzz,    
   //double *Gammayxx, double *Gammayxy, double *Gammayxz, double *Gammayyy, double *Gammayyz, double *Gammayzz,    
   //double *Gammazxx, double *Gammazxy, double *Gammazxz, double *Gammazyy, double *Gammazyz, double *Gammazzz, 
   double vx_boost, int coordinate_type, double sam,
   double *Axx, double *Axy, double *Axz, double *Ayy, double *Ayz,    
   double *Azz, double *trK){ 

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
  double gamma_v = lorentz_gamma*vx_boost;

  // Compute beta_i and temporarily store in shifti. 
  // Also, initialize Aij and trK
#pragma omp for
   for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
      int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
      double gxxL = gxx[index];
      double gxyL = gxy[index];
      double gxzL = gxz[index];
      double gyyL = gyy[index];
      double gyzL = gyz[index];
      double gzzL = gzz[index];
      double betaxL = shiftx[index];
      double betayL = shifty[index];
      double betazL = shiftz[index];
      shiftx[index] = gxxL*betaxL + gxyL*betayL + gxzL*betazL;
      shifty[index] = gxyL*betaxL + gyyL*betayL + gyzL*betazL;
      shiftz[index] = gxzL*betaxL + gyzL*betayL + gzzL*betazL;
      Axx[index] = 0.0;
      Axy[index] = 0.0;
      Axz[index] = 0.0;
      Ayy[index] = 0.0;
      Ayz[index] = 0.0;
      Azz[index] = 0.0;
      trK[index] = 0.0;
   }

#pragma omp for
    for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
      int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
      double gxxL = gxx[index];
      double gxyL = gxy[index];
      double gxzL = gxz[index];
      double gyyL = gyy[index];
      double gyzL = gyz[index];
      double gzzL = gzz[index];
      
      double gijdet = gxxL * gyyL * gzzL + gxyL * gyzL * gxzL 
                      + gxzL * gxyL * gyzL - gxzL * gyyL * gxzL 
                      - gxyL * gxyL * gzzL - gxxL * gyzL * gyzL;
      double gijdetinv = 1.0/gijdet;
      double gupxxL =   ( gyyL * gzzL - gyzL * gyzL )* gijdetinv;
      double gupxyL = - ( gxyL * gzzL - gyzL * gxzL )* gijdetinv;
      double gupxzL =   ( gxyL * gyzL - gyyL * gxzL )* gijdetinv;
      double gupyyL =   ( gxxL * gzzL - gxzL * gxzL )* gijdetinv;
      double gupyzL = - ( gxxL * gyzL - gxyL * gxzL )* gijdetinv;
      double gupzzL =   ( gxxL * gyyL - gxyL * gxyL )* gijdetinv;

      double gxxxL = D1gf(gxx,i,j,k);
      double gxyxL = D1gf(gxy,i,j,k);
      double gxzxL = D1gf(gxz,i,j,k);
      double gyyxL = D1gf(gyy,i,j,k);
      double gyzxL = D1gf(gyz,i,j,k);
      double gzzxL = D1gf(gzz,i,j,k);

      double gxxyL = D2gf(gxx,i,j,k);
      double gxyyL = D2gf(gxy,i,j,k);
      double gxzyL = D2gf(gxz,i,j,k);
      double gyyyL = D2gf(gyy,i,j,k);
      double gyzyL = D2gf(gyz,i,j,k);
      double gzzyL = D2gf(gzz,i,j,k);

      double gxxzL = D3gf(gxx,i,j,k);
      double gxyzL = D3gf(gxy,i,j,k);
      double gxzzL = D3gf(gxz,i,j,k);
      double gyyzL = D3gf(gyy,i,j,k);
      double gyzzL = D3gf(gyz,i,j,k);
      double gzzzL = D3gf(gzz,i,j,k);
     
      double beta_xxL = D1gf(shiftx,i,j,k);
      double beta_yxL = D1gf(shifty,i,j,k);
      double beta_zxL = D1gf(shiftz,i,j,k);
      double beta_xyL = D2gf(shiftx,i,j,k);
      double beta_yyL = D2gf(shifty,i,j,k);
      double beta_zyL = D2gf(shiftz,i,j,k);
      double beta_xzL = D3gf(shiftx,i,j,k);
      double beta_yzL = D3gf(shifty,i,j,k);
      double beta_zzL = D3gf(shiftz,i,j,k);
    
      double gttx0L = gttx0[index]; 
      double gtxx0L = gtxx0[index];
      double gtyx0L = gtyx0[index];
      double gtzx0L = gtzx0[index];
      double gxxx0L = gxxx0[index];
      double gxyx0L = gxyx0[index];
      double gxzx0L = gxzx0[index];
      double gyyx0L = gyyx0[index];
      double gyzx0L = gyzx0[index];
      double gzzx0L = gzzx0[index];
// *** TEST &***
if ( fabs(x[index]-0.38789746925)<1.e-5 && fabs(y[index]-0.01289746925)<1.e-5 && fabs(z[index])<1.e-5 && dx < 0.026)
{
  printf("gxxx, delta gxxx: %e %e\n",gxxxL,gxxx0L);
  printf("gxx-, gxx0, gxx+: %e %e %e\n",gxx[CCTK_GFINDEX3D(cctkGH,i-1,j,k)],gxx[index],gxx[CCTK_GFINDEX3D(cctkGH,i+1,j,k)]);
  printf("i,j,k: %d %d %d\n",i,j,k);
}
// **********

      double Ltx = -gamma_v;
      double Lxx = lorentz_gamma;
    
      // Compute the time derivative of the boosted metric:
      // partial_t g_{ij} = -gamma v (partial x^mu/partial x_boost^i)
      //                     (partial x^nu/partial x_boost^j) 
      //                      partial_x g0_{mu nu} 
      // where g0_{mu nu} is the unboosted metric
      double gxxtL = -gamma_v*(Ltx*Ltx*gttx0L + 2.0*Ltx*Lxx*gtxx0L 
				+ Lxx*Lxx*gxxx0L);
      double gxytL = -gamma_v*(Ltx*gtyx0L + Lxx*gxyx0L);
      double gxztL = -gamma_v*(Ltx*gtzx0L + Lxx*gxzx0L);
      double gyytL = -gamma_v*gyyx0L;
      double gyztL = -gamma_v*gyzx0L;
      double gzztL = -gamma_v*gzzx0L;

      // Now compute Gamma^i_{jk} associated with the boosted 3-metric
      double GamxxxL = 0.5*(gupxxL *gxxxL +gupxyL *(2.0*gxyxL -gxxyL )+gupxzL *(2.0*gxzxL -gxxzL ));
      double GamyxxL = 0.5*(gupxyL *gxxxL +gupyyL *(2.0*gxyxL -gxxyL )+gupyzL *(2.0*gxzxL -gxxzL ));
      double GamzxxL = 0.5*(gupxzL *gxxxL +gupyzL *(2.0*gxyxL -gxxyL )+gupzzL *(2.0*gxzxL -gxxzL ));

      double GamxyyL = 0.5*(gupxxL *(2.0*gxyyL -gyyxL )+gupxyL *gyyyL +gupxzL *(2.0*gyzyL -gyyzL ));
      double GamyyyL = 0.5*(gupxyL *(2.0*gxyyL -gyyxL )+gupyyL *gyyyL +gupyzL *(2.0*gyzyL -gyyzL ));
      double GamzyyL = 0.5*(gupxzL *(2.0*gxyyL -gyyxL )+gupyzL *gyyyL +gupzzL *(2.0*gyzyL -gyyzL ));

      double GamxzzL = 0.5*(gupxxL *(2.0*gxzzL -gzzxL )+gupxyL *(2.0*gyzzL -gzzyL )+gupxzL *gzzzL );
      double GamyzzL = 0.5*(gupxyL *(2.0*gxzzL -gzzxL )+gupyyL *(2.0*gyzzL -gzzyL )+gupyzL *gzzzL );
      double GamzzzL = 0.5*(gupxzL *(2.0*gxzzL -gzzxL )+gupyzL *(2.0*gyzzL -gzzyL )+gupzzL *gzzzL );

      double GamxxyL = 0.5*( gupxxL * gxxyL + gupxyL * gyyxL + gupxzL *( gxzyL + gyzxL - gxyzL ) );
      double GamyxyL = 0.5*( gupxyL * gxxyL + gupyyL * gyyxL + gupyzL *( gxzyL + gyzxL - gxyzL ) );
      double GamzxyL = 0.5*( gupxzL * gxxyL + gupyzL * gyyxL + gupzzL *( gxzyL + gyzxL - gxyzL ) );

      double GamxxzL = 0.5*( gupxxL * gxxzL + gupxyL *( gxyzL + gyzxL - gxzyL ) + gupxzL * gzzxL );
      double GamyxzL = 0.5*( gupxyL * gxxzL + gupyyL *( gxyzL + gyzxL - gxzyL ) + gupyzL * gzzxL );
      double GamzxzL = 0.5*( gupxzL * gxxzL + gupyzL *( gxyzL + gyzxL - gxzyL ) + gupzzL * gzzxL );

      double GamxyzL = 0.5*( gupxxL *( gxyzL + gxzyL - gyzxL ) + gupxyL * gyyzL + gupxzL * gzzyL );
      double GamyyzL = 0.5*( gupxyL *( gxyzL + gxzyL - gyzxL ) + gupyyL * gyyzL + gupyzL * gzzyL );
      double GamzyzL = 0.5*( gupxzL *( gxyzL + gxzyL - gyzxL ) + gupyzL * gyyzL + gupzzL * gzzyL );

      // Now we can compute the boosted K_ij
      double beta_xL = shiftx[index];
      double beta_yL = shifty[index];
      double beta_zL = shiftz[index];
      double f1o2alp = 0.5/(1.0 + lapm1[index]);

      // Note: we need to flip the sign of K_ij inside the horizon for 
      // "puncture" initial data
      if (coordinate_type==3) {
         double riso = sqrt(SQR(lorentz_gamma*x[index]) + SQR(y[index]) + SQR(z[index]));
         double m0 = 1.0;
         double a = sam*m0;
         double rpluso4 = 0.25*( m0 + sqrt(m0*m0-a*a) );
         if (riso < rpluso4) f1o2alp *=-1.0;
      }

      double KxxL = (2.0*beta_xxL - 2.0*GamxxxL*beta_xL - 2.0*GamyxxL*beta_yL 
			- 2.0*GamzxxL*beta_zL - gxxtL)*f1o2alp;
      double KxyL = (beta_xyL + beta_yxL - 2.0*GamxxyL*beta_xL - 2.0*GamyxyL*beta_yL
                        - 2.0*GamzxyL*beta_zL - gxytL)*f1o2alp;
      double KxzL = (beta_xzL + beta_zxL - 2.0*GamxxzL*beta_xL - 2.0*GamyxzL*beta_yL
                        - 2.0*GamzxzL*beta_zL - gxztL)*f1o2alp;
      double KyyL = (2.0*beta_yyL - 2.0*GamxyyL*beta_xL - 2.0*GamyyyL*beta_yL
                        - 2.0*GamzyyL*beta_zL - gyytL)*f1o2alp;
      double KyzL = (beta_yzL + beta_zyL - 2.0*GamxyzL*beta_xL - 2.0*GamyyzL*beta_yL
                        - 2.0*GamzyzL*beta_zL - gyztL)*f1o2alp;
      double KzzL = (2.0*beta_zzL - 2.0*GamxzzL*beta_xL - 2.0*GamyzzL*beta_yL
                        - 2.0*GamzzzL*beta_zL - gzztL)*f1o2alp;

      // Trace of K_ij:
      double trKL = gupxxL*KxxL + 2.0*gupxyL*KxyL + 2.0*gupxzL*KxzL 
			+ gupyyL*KyyL + 2.0*gupyzL*KyzL + gupzzL*KzzL;
      trK[index] = trKL;

      // Finally, A_ij:
      double f1o3 = 1.0/3.0;
      Axx[index] = KxxL - f1o3*gxxL*trKL;
      Axy[index] = KxyL - f1o3*gxyL*trKL;
      Axz[index] = KxzL - f1o3*gxzL*trKL;
      Ayy[index] = KyyL - f1o3*gyyL*trKL;
      Ayz[index] = KyzL - f1o3*gyzL*trKL;
      Azz[index] = KzzL - f1o3*gzzL*trKL;
    }

 // Convert shifti back and compute the tilde metric and Aij
#pragma omp for    
    for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
      int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

      double gxxL = gxx[index];
      double gxyL = gxy[index];
      double gxzL = gxz[index];
      double gyyL = gyy[index];
      double gyzL = gyz[index];
      double gzzL = gzz[index];

      double gijdet = gxxL * gyyL * gzzL + gxyL * gyzL * gxzL
                      + gxzL * gxyL * gyzL - gxzL * gyyL * gxzL
                      - gxyL * gxyL * gzzL - gxxL * gyzL * gyzL;       
      double gijdetinv = 1.0/gijdet;
      double gupxxL =   ( gyyL * gzzL - gyzL * gyzL )* gijdetinv;
      double gupxyL = - ( gxyL * gzzL - gyzL * gxzL )* gijdetinv;
      double gupxzL =   ( gxyL * gyzL - gyyL * gxzL )* gijdetinv;
      double gupyyL =   ( gxxL * gzzL - gxzL * gxzL )* gijdetinv;
      double gupyzL = - ( gxxL * gyzL - gxyL * gxzL )* gijdetinv;
      double gupzzL =   ( gxxL * gyyL - gxyL * gxyL )* gijdetinv;

      double AxxL = Axx[index];
      double AxyL = Axy[index];
      double AxzL = Axz[index];
      double AyyL = Ayy[index];
      double AyzL = Ayz[index];
      double AzzL = Azz[index];
 
      double beta_xL = shiftx[index];
      double beta_yL = shifty[index];
      double beta_zL = shiftz[index];    

      shiftx[index] = gupxxL*beta_xL + gupxyL*beta_yL + gupxzL*beta_zL;
      shifty[index] = gupxyL*beta_xL + gupyyL*beta_yL + gupyzL*beta_zL;
      shiftz[index] = gupxzL*beta_xL + gupyzL*beta_yL + gupzzL*beta_zL;

      double psi4 = pow(gijdet, 1.0/3.0);
      double psim4 = 1.0/psi4;
      phi[index] = log(gijdet)/12.0;
      gxx[index] = psim4*gxxL;
      gxy[index] = psim4*gxyL;
      gxz[index] = psim4*gxzL;
      gyy[index] = psim4*gyyL;
      gyz[index] = psim4*gyzL;
      gzz[index] = psim4*gzzL;
      gupxx[index] = psi4*gupxxL;
      gupxy[index] = psi4*gupxyL;
      gupxz[index] = psi4*gupxzL;
      gupyy[index] = psi4*gupyyL;
      gupyz[index] = psi4*gupyzL;
      gupzz[index] = psi4*gupzzL;
      Axx[index] = psim4*AxxL;
      Axy[index] = psim4*AxyL;
      Axz[index] = psim4*AxzL;
      Ayy[index] = psim4*AyyL;
      Ayz[index] = psim4*AyzL;
      Azz[index] = psim4*AzzL;

  }

}

extern "C" void CCTK_FCALL CCTK_FNAME(compute_Aij_boost)
  (const cGH **cctkGH,double *dx, double *dy, double *dz,
   double *x, double *y, double *z,
   int *nghostzones, int *cctk_lsh, double *lapm1, double *phi,
   double *shiftx, double *shifty, double *shiftz,
   double *gxx, double *gxy, double *gxz, double *gyy,
   double *gyz, double *gzz,
   double *gupxx, double *gupxy, double *gupxz, double *gupyy,
   double *gupyz, double *gupzz,
   double *gttx0, double *gtxx0,
   double *gtyx0, double *gtzx0, double *gxxx0, double *gxyx0,
   double *gxzx0, double *gyyx0, double *gyzx0, double *gzzx0,
   //double *Gammaxxx, double *Gammaxxy, double *Gammaxxz, double *Gammaxyy, double *Gammaxyz, double *Gammaxzz,
   //double *Gammayxx, double *Gammayxy, double *Gammayxz, double *Gammayyy, double *Gammayyz, double *Gammayzz,
   //double *Gammazxx, double *Gammazxy, double *Gammazxz, double *Gammazyy, double *Gammazyz, double *Gammazzzi, 
   double *vx_boost, int *coordinate_type, double *sam,
   double *Axx, double *Axy, double *Axz, double *Ayy, double *Ayz,
   double *Azz, double *trK)
{ 
   compute_Aij_boost(*cctkGH, *dx, *dy, *dz, x,y,z,
   nghostzones, cctk_lsh, lapm1, phi, shiftx, shifty, shiftz,
   gxx, gxy, gxz, gyy, gyz, gzz, gupxx, gupxy, gupxz, gupyy,
   gupyz, gupzz, gttx0, gtxx0, gtyx0, gtzx0, gxxx0, gxyx0,
   gxzx0, gyyx0, gyzx0, gzzx0,
   //Gammaxxx, Gammaxxy, Gammaxxz, Gammaxyy, Gammaxyz, Gammaxzz,
   //Gammayxx, Gammayxy, Gammayxz, Gammayyy, Gammayyz, Gammayzz,
   //Gammazxx, Gammazxy, Gammazxz, Gammazyy, Gammazyz, Gammazzzi, 
   *vx_boost, *coordinate_type, *sam,
   Axx, Axy, Axz, Ayy, Ayz, Azz, trK);
}
