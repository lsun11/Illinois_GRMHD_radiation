#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "PsiKadelia_constants.h"
#include "nPF_decl.h"
#include "stdio.h"
#include "stdlib.h"

#define KRANC_C 
#include "GenericFD.h"

#define SQR(x) ((x) * (x))

extern "C" void CCTK_FCALL psikadelia_lowlevel_
(const cGH **cctkGH,int *cctk_lsh,int *cctk_nghostzones,
 double &dx,double &dy,double &dz,double *r,double *x,double *y,double *z,
 double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
 double *gxx_phys,double *gxy_phys,double *gxz_phys,double *gyy_phys,double *gyz_phys,double *gzz_phys,
 double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,
 double *phi,double *psi,
 double *Kxx,double *Kxy,double *Kxz,double *Kyy,double *Kyz,double *Kzz,
 double *Rxx,double *Rxy,double *Rxz,double *Ryy,double *Ryz,double *Rzz,
 double *Gammaxxx,double *Gammayxx,double *Gammazxx,double *Gammaxyy,double *Gammayyy,double *Gammazyy,
 double *Gammaxzz,double *Gammayzz,double *Gammazzz,double *Gammaxxy,double *Gammayxy,double *Gammazxy,
 double *Gammaxxz,double *Gammayxz,double *Gammazxz,double *Gammaxyz,double *Gammayyz,double *Gammazyz,
 double *psi0re,double *psi0im);


extern "C" void PsiKadelia_lowlevel(const cGH *cctkGH,int *cctk_lsh,int *cctk_nghostzones,
				    double &dx,double &dy,double &dz,double *r,double *x,double *y,double *z,
				    double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
				    double *gxx_phys,double *gxy_phys,double *gxz_phys,double *gyy_phys,double *gyz_phys,double *gzz_phys,
				    double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,
				    double *phi,double *psi,
				    double *kxx,double *kxy,double *kxz,double *kyy,double *kyz,double *kzz,
				    double *Rxx,double *Rxy,double *Rxz,double *Ryy,double *Ryz,double *Rzz,
				    double *Gammaxxx,double *Gammayxx,double *Gammazxx,double *Gammaxyy,double *Gammayyy,double *Gammazyy,
				    double *Gammaxzz,double *Gammayzz,double *Gammazzz,double *Gammaxxy,double *Gammayxy,double *Gammazxy,
				    double *Gammaxxz,double *Gammayxz,double *Gammazxz,double *Gammaxyz,double *Gammayyz,double *Gammazyz,
				    double *psi0re,double *psi0im) {

  /* Initialise finite differencing variables.  NEED THIS FOR GenericFD.h */
#include "../../../GenFD_decl_set_varCPP.h"

  DECLARE_CCTK_PARAMETERS

  int whichvec = 0;
  if (CCTK_Equals(psif_vec,"radial")==1) {
    whichvec = NEWAGEPSIFRIENDS_RADIAL_VEC;
  } else if (CCTK_Equals(psif_vec,"cartesian")==1) {
    whichvec = NEWAGEPSIFRIENDS_CARTESIAN_VEC;
  } else if (CCTK_Equals(psif_vec,"metric_diag")==1) { 
    whichvec = NEWAGEPSIFRIENDS_METDIAG_VEC;
  } else if (CCTK_Equals(psif_vec,"shock")==1) { 
    whichvec = NEWAGEPSIFRIENDS_SHOCK_VEC;
  }
  
#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    if(r[index] > compute_Psi4_min_radius && r[index] < compute_Psi4_max_radius) {

      double fac;

      double v1xL,v2xL,v3xL;
      double v1yL,v2yL,v3yL;
      double v1zL,v2zL,v3zL; 
      double ricxx_loc,ricxy_loc,ricxz_loc,ricyy_loc,ricyz_loc,riczz_loc,ricscal_loc;
      double elecxx,elecxy,elecxz,elecyy,elecyz,eleczz;
      double magxx,magxy,magxz,magyy,magyz,magzz;
      double psi0re_loc,psi0im_loc,psi1re_loc,psi1im_loc,psi2re_loc;
      double psi2im_loc,psi3re_loc,psi3im_loc,psi4re_loc,psi4im_loc;

#include "ricci_temps.h"

      double PsiL = psi[index];
      double psi_4thpower = PsiL*PsiL*PsiL*PsiL;
      double two = 2.0;
      double one = 1.0;

      
      int local_spatial_order=2;

      if(i==1 || j==1 || k==1 || i==cctk_lsh[0]-2 || j==cctk_lsh[1]-2 || k==cctk_lsh[2]-2) {
	local_spatial_order=0;
      }
	
#ifdef FD_C4
      if(i==1 || j==1 || k==1 || i==cctk_lsh[0]-2 || j==cctk_lsh[1]-2 || k==cctk_lsh[2]-2) {
	local_spatial_order=2;
      } else {
	local_spatial_order=4;
      }
#endif
#ifdef FD_C6
      if(i==1 || j==1 || k==1 || i==cctk_lsh[0]-2 || j==cctk_lsh[1]-2 || k==cctk_lsh[2]-2) {
	local_spatial_order=2;
      } else if(i==2 || j==2 || k==2 || i==cctk_lsh[0]-3 || j==cctk_lsh[1]-3 || k==cctk_lsh[2]-3) {
	local_spatial_order=4;
      } else {
	local_spatial_order=6;
      }
#endif
      if(i==0 || j==0 || k==0 || i==cctk_lsh[0]-1 || j==cctk_lsh[1]-1 || k==cctk_lsh[2]-1) {
	local_spatial_order=0;
      } 



      //     Define local k,g variables to pass to subroutines
      double kkxx=kxx[index];
      double kkxy=kxy[index];
      double kkxz=kxz[index];
      double kkyy=kyy[index];
      double kkyz=kyz[index];
      double kkzz=kzz[index];
      double ggxx=gxx_phys[index];
      double ggxy=gxy_phys[index];
      double ggxz=gxz_phys[index];
      double ggyy=gyy_phys[index];
      double ggyz=gyz_phys[index];
      double ggzz=gzz_phys[index];
      double xx=x[index];
      double yy=y[index];
      double zz=z[index];

      // FIRST: Set up the tetrad vectors
#include "nPF_vecs.C"

      //     Find determinant.
      double det= -(SQR(gxz_phys[index])*gyy_phys[index]) + 
	2*gxy_phys[index]*gxz_phys[index]*gyz_phys[index] - 
	gxx_phys[index]*SQR(gyz_phys[index])  -
	SQR(gxy_phys[index])*gzz_phys[index] + 
	gxx_phys[index]*gyy_phys[index]*gzz_phys[index];
      
      //     Invert metric. This is the conformal upper metric.
      double uxx=(-SQR(gyz_phys[index]) + gyy_phys[index]*gzz_phys[index])/det;
      double uxy=(gxz_phys[index]*gyz_phys[index] - gxy_phys[index]*gzz_phys[index])/det;
      double uyy=(-SQR(gxz_phys[index]) + gxx_phys[index]*gzz_phys[index])/det;
      double uxz=(-gxz_phys[index]*gyy_phys[index] + gxy_phys[index]*gyz_phys[index])/det;
      double uyz=(gxy_phys[index]*gxz_phys[index] - gxx_phys[index]*gyz_phys[index])/det;
      double uzz=(-SQR(gxy_phys[index]) + gxx_phys[index]*gyy_phys[index])/det;

      //     Form traceK.
      double traceK =  (uxx*kxx[index]+
			uyy*kyy[index]+
			uzz*kzz[index]+
			two*uxy*kxy[index]+
			two*uxz*kxz[index]+
			two*uyz*kyz[index]);
      
      double phixL,phiyL,phizL,func;
      double phixx,phixy,phixz,phiyy,phiyz,phizz;
      double gupxxL,gupxyL,gupxzL,gupyyL,gupyzL,gupzzL;
      double gxxL,gxyL,gxzL,gyyL,gyzL,gzzL;
      double GamxxxL, GamxxyL, GamxxzL, GamxyyL, GamxyzL, GamxzzL;
      double GamyxxL, GamyxyL, GamyxzL, GamyyyL, GamyyzL, GamyzzL;
      double GamzxxL, GamzxyL, GamzxzL, GamzyyL, GamzyzL, GamzzzL;
      double gxxxL,gxxyL,gxxzL;
      double gxyxL,gxyyL,gxyzL;
      double gxzxL,gxzyL,gxzzL;
      double gyyxL,gyyyL,gyyzL;
      double gyzxL,gyzyL,gyzzL;
      double gzzxL,gzzyL,gzzzL;
      double psixL,psiyL,psizL;
      double dxkxx, dxgxx, cdxkxx;
      double dxkxy, dxgxy, cdxkxy;
      double dxkxz, dxgxz, cdxkxz;
      double dxkyy, dxgyy, cdxkyy;
      double dxkyz, dxgyz, cdxkyz;
      double dxkzz, dxgzz, cdxkzz;
      double dykxx, dygxx, cdykxx;
      double dykxy, dygxy, cdykxy;
      double dykxz, dygxz, cdykxz;
      double dykyy, dygyy, cdykyy;
      double dykyz, dygyz, cdykyz;
      double dykzz, dygzz, cdykzz;
      double dzkxx, dzgxx, cdzkxx;
      double dzkxy, dzgxy, cdzkxy;
      double dzkxz, dzgxz, cdzkxz;
      double dzkyy, dzgyy, cdzkyy;
      double dzkyz, dzgyz, cdzkyz;
      double dzkzz, dzgzz, cdzkzz;
      double d2xx_gxx,d2xx_gxy,d2xx_gxz,d2xx_gyy,d2xx_gyz,d2xx_gzz;
      double d2xy_gxx,d2xy_gxy,d2xy_gxz,d2xy_gyy,d2xy_gyz,d2xy_gzz;
      double d2xz_gxx,d2xz_gxy,d2xz_gxz,d2xz_gyy,d2xz_gyz,d2xz_gzz;
      double d2yy_gxx,d2yy_gxy,d2yy_gxz,d2yy_gyy,d2yy_gyz,d2yy_gzz;
      double d2yz_gxx,d2yz_gxy,d2yz_gxz,d2yz_gyy,d2yz_gyz,d2yz_gzz;
      double d2zz_gxx,d2zz_gxy,d2zz_gxz,d2zz_gyy,d2zz_gyz,d2zz_gzz;

      if( whichvec != NEWAGEPSIFRIENDS_RADIAL_VEC || xx*yy!=0) {
	//     Calculate only if not on the z-axis for the "radial" tetrad case
	
	//     Form derivatives of K and G  and second derivatives of G.
	if(local_spatial_order==6) {
	  dxkxx = D1_c6(kxx,i,j,k);
	  dxkxy = D1_c6(kxy,i,j,k);
	  dxkxz = D1_c6(kxz,i,j,k);
	  dxkyy = D1_c6(kyy,i,j,k);
	  dxkyz = D1_c6(kyz,i,j,k);
	  dxkzz = D1_c6(kzz,i,j,k);
	  dykxx = D2_c6(kxx,i,j,k);
	  dykxy = D2_c6(kxy,i,j,k);
	  dykxz = D2_c6(kxz,i,j,k);
	  dykyy = D2_c6(kyy,i,j,k); 
	  dykyz = D2_c6(kyz,i,j,k); 
	  dykzz = D2_c6(kzz,i,j,k);
	  dzkxx = D3_c6(kxx,i,j,k);
	  dzkxy = D3_c6(kxy,i,j,k);
	  dzkxz = D3_c6(kxz,i,j,k);
	  dzkyy = D3_c6(kyy,i,j,k);
	  dzkyz = D3_c6(kyz,i,j,k);
	  dzkzz = D3_c6(kzz,i,j,k);

	  gxxxL = D1_c6(gxx,i,j,k);
	  gxxyL = D2_c6(gxx,i,j,k);
	  gxxzL = D3_c6(gxx,i,j,k);
	  gxyxL = D1_c6(gxy,i,j,k);
	  gxyyL = D2_c6(gxy,i,j,k);
	  gxyzL = D3_c6(gxy,i,j,k);
	  gxzxL = D1_c6(gxz,i,j,k);
	  gxzyL = D2_c6(gxz,i,j,k);
	  gxzzL = D3_c6(gxz,i,j,k);
	  gyyxL = D1_c6(gyy,i,j,k);
	  gyyyL = D2_c6(gyy,i,j,k);
	  gyyzL = D3_c6(gyy,i,j,k);
	  gyzxL = D1_c6(gyz,i,j,k);
	  gyzyL = D2_c6(gyz,i,j,k);
	  gyzzL = D3_c6(gyz,i,j,k);
	  gzzxL = D1_c6(gzz,i,j,k);
	  gzzyL = D2_c6(gzz,i,j,k);
	  gzzzL = D3_c6(gzz,i,j,k);

	  d2xx_gxx = D11_c6(gxx_phys,i,j,k);
	  d2xx_gxy = D11_c6(gxy_phys,i,j,k);
	  d2xx_gxz = D11_c6(gxz_phys,i,j,k);
	  d2xx_gyy = D11_c6(gyy_phys,i,j,k);
	  d2xx_gyz = D11_c6(gyz_phys,i,j,k);
	  d2xx_gzz = D11_c6(gzz_phys,i,j,k);

	  d2xy_gxx = D21_c6(gxx_phys,i,j,k);
	  d2xy_gxy = D21_c6(gxy_phys,i,j,k);
	  d2xy_gxz = D21_c6(gxz_phys,i,j,k);
	  d2xy_gyy = D21_c6(gyy_phys,i,j,k);
	  d2xy_gyz = D21_c6(gyz_phys,i,j,k);
	  d2xy_gzz = D21_c6(gzz_phys,i,j,k);

	  d2xz_gxx = D31_c6(gxx_phys,i,j,k);
	  d2xz_gxy = D31_c6(gxy_phys,i,j,k);
	  d2xz_gxz = D31_c6(gxz_phys,i,j,k);
	  d2xz_gyy = D31_c6(gyy_phys,i,j,k);
	  d2xz_gyz = D31_c6(gyz_phys,i,j,k);
	  d2xz_gzz = D31_c6(gzz_phys,i,j,k);

	  d2yy_gxx = D22_c6(gxx_phys,i,j,k);
	  d2yy_gxy = D22_c6(gxy_phys,i,j,k);
	  d2yy_gxz = D22_c6(gxz_phys,i,j,k);
	  d2yy_gyy = D22_c6(gyy_phys,i,j,k);
	  d2yy_gyz = D22_c6(gyz_phys,i,j,k);
	  d2yy_gzz = D22_c6(gzz_phys,i,j,k);

	  d2yz_gxx = D32_c6(gxx_phys,i,j,k);
	  d2yz_gxy = D32_c6(gxy_phys,i,j,k);
	  d2yz_gxz = D32_c6(gxz_phys,i,j,k);
	  d2yz_gyy = D32_c6(gyy_phys,i,j,k);
	  d2yz_gyz = D32_c6(gyz_phys,i,j,k);
	  d2yz_gzz = D32_c6(gzz_phys,i,j,k);

	  d2zz_gxx = D33_c6(gxx_phys,i,j,k);
	  d2zz_gxy = D33_c6(gxy_phys,i,j,k);
	  d2zz_gxz = D33_c6(gxz_phys,i,j,k);
	  d2zz_gyy = D33_c6(gyy_phys,i,j,k);
	  d2zz_gyz = D33_c6(gyz_phys,i,j,k);
	  d2zz_gzz = D33_c6(gzz_phys,i,j,k);

	  psixL = D1_c6(psi,i,j,k);
	  psiyL = D2_c6(psi,i,j,k);
	  psizL = D3_c6(psi,i,j,k);

	//     Calculate only if not on the z-axis for the "radial" tetrad case
	
	//     Form derivatives of K and G  and second derivatives of G.
	} else if(local_spatial_order==4) {
	  dxkxx = D1_c4(kxx,i,j,k);
	  dxkxy = D1_c4(kxy,i,j,k);
	  dxkxz = D1_c4(kxz,i,j,k);
	  dxkyy = D1_c4(kyy,i,j,k);
	  dxkyz = D1_c4(kyz,i,j,k);
	  dxkzz = D1_c4(kzz,i,j,k);
	  dykxx = D2_c4(kxx,i,j,k);
	  dykxy = D2_c4(kxy,i,j,k);
	  dykxz = D2_c4(kxz,i,j,k);
	  dykyy = D2_c4(kyy,i,j,k); 
	  dykyz = D2_c4(kyz,i,j,k); 
	  dykzz = D2_c4(kzz,i,j,k);
	  dzkxx = D3_c4(kxx,i,j,k);
	  dzkxy = D3_c4(kxy,i,j,k);
	  dzkxz = D3_c4(kxz,i,j,k);
	  dzkyy = D3_c4(kyy,i,j,k);
	  dzkyz = D3_c4(kyz,i,j,k);
	  dzkzz = D3_c4(kzz,i,j,k);

	  gxxxL = D1_c4(gxx,i,j,k);
	  gxxyL = D2_c4(gxx,i,j,k);
	  gxxzL = D3_c4(gxx,i,j,k);
	  gxyxL = D1_c4(gxy,i,j,k);
	  gxyyL = D2_c4(gxy,i,j,k);
	  gxyzL = D3_c4(gxy,i,j,k);
	  gxzxL = D1_c4(gxz,i,j,k);
	  gxzyL = D2_c4(gxz,i,j,k);
	  gxzzL = D3_c4(gxz,i,j,k);
	  gyyxL = D1_c4(gyy,i,j,k);
	  gyyyL = D2_c4(gyy,i,j,k);
	  gyyzL = D3_c4(gyy,i,j,k);
	  gyzxL = D1_c4(gyz,i,j,k);
	  gyzyL = D2_c4(gyz,i,j,k);
	  gyzzL = D3_c4(gyz,i,j,k);
	  gzzxL = D1_c4(gzz,i,j,k);
	  gzzyL = D2_c4(gzz,i,j,k);
	  gzzzL = D3_c4(gzz,i,j,k);

	  d2xx_gxx = D11_c4(gxx_phys,i,j,k);
	  d2xx_gxy = D11_c4(gxy_phys,i,j,k);
	  d2xx_gxz = D11_c4(gxz_phys,i,j,k);
	  d2xx_gyy = D11_c4(gyy_phys,i,j,k);
	  d2xx_gyz = D11_c4(gyz_phys,i,j,k);
	  d2xx_gzz = D11_c4(gzz_phys,i,j,k);

	  d2xy_gxx = D21_c4(gxx_phys,i,j,k);
	  d2xy_gxy = D21_c4(gxy_phys,i,j,k);
	  d2xy_gxz = D21_c4(gxz_phys,i,j,k);
	  d2xy_gyy = D21_c4(gyy_phys,i,j,k);
	  d2xy_gyz = D21_c4(gyz_phys,i,j,k);
	  d2xy_gzz = D21_c4(gzz_phys,i,j,k);

	  d2xz_gxx = D31_c4(gxx_phys,i,j,k);
	  d2xz_gxy = D31_c4(gxy_phys,i,j,k);
	  d2xz_gxz = D31_c4(gxz_phys,i,j,k);
	  d2xz_gyy = D31_c4(gyy_phys,i,j,k);
	  d2xz_gyz = D31_c4(gyz_phys,i,j,k);
	  d2xz_gzz = D31_c4(gzz_phys,i,j,k);

	  d2yy_gxx = D22_c4(gxx_phys,i,j,k);
	  d2yy_gxy = D22_c4(gxy_phys,i,j,k);
	  d2yy_gxz = D22_c4(gxz_phys,i,j,k);
	  d2yy_gyy = D22_c4(gyy_phys,i,j,k);
	  d2yy_gyz = D22_c4(gyz_phys,i,j,k);
	  d2yy_gzz = D22_c4(gzz_phys,i,j,k);

	  d2yz_gxx = D32_c4(gxx_phys,i,j,k);
	  d2yz_gxy = D32_c4(gxy_phys,i,j,k);
	  d2yz_gxz = D32_c4(gxz_phys,i,j,k);
	  d2yz_gyy = D32_c4(gyy_phys,i,j,k);
	  d2yz_gyz = D32_c4(gyz_phys,i,j,k);
	  d2yz_gzz = D32_c4(gzz_phys,i,j,k);

	  d2zz_gxx = D33_c4(gxx_phys,i,j,k);
	  d2zz_gxy = D33_c4(gxy_phys,i,j,k);
	  d2zz_gxz = D33_c4(gxz_phys,i,j,k);
	  d2zz_gyy = D33_c4(gyy_phys,i,j,k);
	  d2zz_gyz = D33_c4(gyz_phys,i,j,k);
	  d2zz_gzz = D33_c4(gzz_phys,i,j,k);

	  psixL = D1_c4(psi,i,j,k);
	  psiyL = D2_c4(psi,i,j,k);
	  psizL = D3_c4(psi,i,j,k);

	} else if(local_spatial_order==2) {
	  dxkxx = D1_c2(kxx,i,j,k);
	  dxkxy = D1_c2(kxy,i,j,k);
	  dxkxz = D1_c2(kxz,i,j,k);
	  dxkyy = D1_c2(kyy,i,j,k);
	  dxkyz = D1_c2(kyz,i,j,k);
	  dxkzz = D1_c2(kzz,i,j,k);
	  dykxx = D2_c2(kxx,i,j,k);
	  dykxy = D2_c2(kxy,i,j,k);
	  dykxz = D2_c2(kxz,i,j,k);
	  dykyy = D2_c2(kyy,i,j,k); 
	  dykyz = D2_c2(kyz,i,j,k); 
	  dykzz = D2_c2(kzz,i,j,k);
	  dzkxx = D3_c2(kxx,i,j,k);
	  dzkxy = D3_c2(kxy,i,j,k);
	  dzkxz = D3_c2(kxz,i,j,k);
	  dzkyy = D3_c2(kyy,i,j,k);
	  dzkyz = D3_c2(kyz,i,j,k);
	  dzkzz = D3_c2(kzz,i,j,k);

	  gxxxL = D1_c2(gxx,i,j,k);
	  gxxyL = D2_c2(gxx,i,j,k);
	  gxxzL = D3_c2(gxx,i,j,k);
	  gxyxL = D1_c2(gxy,i,j,k);
	  gxyyL = D2_c2(gxy,i,j,k);
	  gxyzL = D3_c2(gxy,i,j,k);
	  gxzxL = D1_c2(gxz,i,j,k);
	  gxzyL = D2_c2(gxz,i,j,k);
	  gxzzL = D3_c2(gxz,i,j,k);
	  gyyxL = D1_c2(gyy,i,j,k);
	  gyyyL = D2_c2(gyy,i,j,k);
	  gyyzL = D3_c2(gyy,i,j,k);
	  gyzxL = D1_c2(gyz,i,j,k);
	  gyzyL = D2_c2(gyz,i,j,k);
	  gyzzL = D3_c2(gyz,i,j,k);
	  gzzxL = D1_c2(gzz,i,j,k);
	  gzzyL = D2_c2(gzz,i,j,k);
	  gzzzL = D3_c2(gzz,i,j,k);

	  d2xx_gxx = D11_c2(gxx_phys,i,j,k);
	  d2xx_gxy = D11_c2(gxy_phys,i,j,k);
	  d2xx_gxz = D11_c2(gxz_phys,i,j,k);
	  d2xx_gyy = D11_c2(gyy_phys,i,j,k);
	  d2xx_gyz = D11_c2(gyz_phys,i,j,k);
	  d2xx_gzz = D11_c2(gzz_phys,i,j,k);

	  d2xy_gxx = D21_c2(gxx_phys,i,j,k);
	  d2xy_gxy = D21_c2(gxy_phys,i,j,k);
	  d2xy_gxz = D21_c2(gxz_phys,i,j,k);
	  d2xy_gyy = D21_c2(gyy_phys,i,j,k);
	  d2xy_gyz = D21_c2(gyz_phys,i,j,k);
	  d2xy_gzz = D21_c2(gzz_phys,i,j,k);

	  d2xz_gxx = D31_c2(gxx_phys,i,j,k);
	  d2xz_gxy = D31_c2(gxy_phys,i,j,k);
	  d2xz_gxz = D31_c2(gxz_phys,i,j,k);
	  d2xz_gyy = D31_c2(gyy_phys,i,j,k);
	  d2xz_gyz = D31_c2(gyz_phys,i,j,k);
	  d2xz_gzz = D31_c2(gzz_phys,i,j,k);

	  d2yy_gxx = D22_c2(gxx_phys,i,j,k);
	  d2yy_gxy = D22_c2(gxy_phys,i,j,k);
	  d2yy_gxz = D22_c2(gxz_phys,i,j,k);
	  d2yy_gyy = D22_c2(gyy_phys,i,j,k);
	  d2yy_gyz = D22_c2(gyz_phys,i,j,k);
	  d2yy_gzz = D22_c2(gzz_phys,i,j,k);

	  d2yz_gxx = D32_c2(gxx_phys,i,j,k);
	  d2yz_gxy = D32_c2(gxy_phys,i,j,k);
	  d2yz_gxz = D32_c2(gxz_phys,i,j,k);
	  d2yz_gyy = D32_c2(gyy_phys,i,j,k);
	  d2yz_gyz = D32_c2(gyz_phys,i,j,k);
	  d2yz_gzz = D32_c2(gzz_phys,i,j,k);

	  d2zz_gxx = D33_c2(gxx_phys,i,j,k);
	  d2zz_gxy = D33_c2(gxy_phys,i,j,k);
	  d2zz_gxz = D33_c2(gxz_phys,i,j,k);
	  d2zz_gyy = D33_c2(gyy_phys,i,j,k);
	  d2zz_gyz = D33_c2(gyz_phys,i,j,k);
	  d2zz_gzz = D33_c2(gzz_phys,i,j,k);

	  psixL = D1_c2(psi,i,j,k);
	  psiyL = D2_c2(psi,i,j,k);
	  psizL = D3_c2(psi,i,j,k);

	} else {
	  dxkxx = 100.0;
	  dxkxy = 100.0;
	  dxkxz = 100.0;
	  dxkyy = 100.0;
	  dxkyz = 100.0;
	  dxkzz = 100.0;
	  dykxx = 100.0;
	  dykxy = 100.0;
	  dykxz = 100.0;
	  dykyy = 100.0;
	  dykyz = 100.0;
	  dykzz = 100.0;
	  dzkxx = 100.0;
	  dzkxy = 100.0;
	  dzkxz = 100.0;
	  dzkyy = 100.0;
	  dzkyz = 100.0;
	  dzkzz = 100.0;

	  gxxxL = 100.0;
	  gxxyL = 100.0;
	  gxxzL = 100.0;
	  gxyxL = 100.0;
	  gxyyL = 100.0;
	  gxyzL = 100.0;
	  gxzxL = 100.0;
	  gxzyL = 100.0;
	  gxzzL = 100.0;
	  gyyxL = 100.0;
	  gyyyL = 100.0;
	  gyyzL = 100.0;
	  gyzxL = 100.0;
	  gyzyL = 100.0;
	  gyzzL = 100.0;
	  gzzxL = 100.0;
	  gzzyL = 100.0;
	  gzzzL = 100.0;


	  d2xx_gxx = 100.0;
	  d2xx_gxy = 100.0;
	  d2xx_gxz = 100.0;
	  d2xx_gyy = 100.0;
	  d2xx_gyz = 100.0;
	  d2xx_gzz = 100.0;

	  d2xy_gxx = 100.0;
	  d2xy_gxy = 100.0;
	  d2xy_gxz = 100.0;
	  d2xy_gyy = 100.0;
	  d2xy_gyz = 100.0;
	  d2xy_gzz = 100.0;

	  d2xz_gxx = 100.0;
	  d2xz_gxy = 100.0;
	  d2xz_gxz = 100.0;
	  d2xz_gyy = 100.0;
	  d2xz_gyz = 100.0;
	  d2xz_gzz = 100.0;

	  d2yy_gxx = 100.0;
	  d2yy_gxy = 100.0;
	  d2yy_gxz = 100.0;
	  d2yy_gyy = 100.0;
	  d2yy_gyz = 100.0;
	  d2yy_gzz = 100.0;

	  d2yz_gxx = 100.0;
	  d2yz_gxy = 100.0;
	  d2yz_gxz = 100.0;
	  d2yz_gyy = 100.0;
	  d2yz_gyz = 100.0;
	  d2yz_gzz = 100.0;

	  d2zz_gxx = 100.0;
	  d2zz_gxy = 100.0;
	  d2zz_gxz = 100.0;
	  d2zz_gyy = 100.0;
	  d2zz_gyz = 100.0;
	  d2zz_gzz = 100.0;

	  psixL = 100.0;
	  psiyL = 100.0;
	  psizL = 100.0;
	}

	dxgxx = gxxxL*psi_4thpower +  4.0*psixL*gxx_phys[index];
	dxgxy = gxyxL*psi_4thpower +  4.0*psixL*gxy_phys[index];
	dxgxz = gxzxL*psi_4thpower +  4.0*psixL*gxz_phys[index];
	dxgyy = gyyxL*psi_4thpower +  4.0*psixL*gyy_phys[index];
	dxgyz = gyzxL*psi_4thpower +  4.0*psixL*gyz_phys[index];
	dxgzz = gzzxL*psi_4thpower +  4.0*psixL*gzz_phys[index];

	dygxx = gxxyL*psi_4thpower +  4.0*psiyL*gxx_phys[index];
	dygxy = gxyyL*psi_4thpower +  4.0*psiyL*gxy_phys[index];
	dygxz = gxzyL*psi_4thpower +  4.0*psiyL*gxz_phys[index];
	dygyy = gyyyL*psi_4thpower +  4.0*psiyL*gyy_phys[index];
	dygyz = gyzyL*psi_4thpower +  4.0*psiyL*gyz_phys[index];
	dygzz = gzzyL*psi_4thpower +  4.0*psiyL*gzz_phys[index];

	dzgxx = gxxzL*psi_4thpower +  4.0*psizL*gxx_phys[index];
	dzgxy = gxyzL*psi_4thpower +  4.0*psizL*gxy_phys[index];
	dzgxz = gxzzL*psi_4thpower +  4.0*psizL*gxz_phys[index];
	dzgyy = gyyzL*psi_4thpower +  4.0*psizL*gyy_phys[index];
	dzgyz = gyzzL*psi_4thpower +  4.0*psizL*gyz_phys[index];
	dzgzz = gzzzL*psi_4thpower +  4.0*psizL*gzz_phys[index];

	//Compute the Ricci tensor from scratch:
#include "ricci.C"

	//     And calculate
#include "electric_and_magnetic.C"
		 
#include "psiKadelia.C"
		 
	//     If required scale by radial coordinate
	if (scale_with_radius == 1) {
	  fac = pow(r[index],radius_power);
	} else {
	  fac = 1.0;
	}
		 
	//     Psis we do pairwise
	psi0re[index] = fac*psi0re_loc;
	psi0im[index] = fac*psi0im_loc;
      } else { 
	//     Set values to zero without calculating on the z-axis
	psi0re[index] = 0;
	psi0im[index] = 0;
      }
    } else {
      psi0re[index] = -100000;
      psi0im[index] = -100000;
    }

  }
  //////$omp end parallel do

  printf("outside Psikadelia()\n");
}

extern "C" void CCTK_FCALL psikadelia_lowlevel_
(const cGH **cctkGH,int *cctk_lsh,int *cctk_nghostzones,
 double &dx,double &dy,double &dz,double *r,double *x,double *y,double *z,
 double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
 double *gxx_phys,double *gxy_phys,double *gxz_phys,double *gyy_phys,double *gyz_phys,double *gzz_phys,
 double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,
 double *phi,double *psi,
 double *Kxx,double *Kxy,double *Kxz,double *Kyy,double *Kyz,double *Kzz,
 double *Rxx,double *Rxy,double *Rxz,double *Ryy,double *Ryz,double *Rzz,
 double *Gammaxxx,double *Gammayxx,double *Gammazxx,double *Gammaxyy,double *Gammayyy,double *Gammazyy,
 double *Gammaxzz,double *Gammayzz,double *Gammazzz,double *Gammaxxy,double *Gammayxy,double *Gammazxy,
 double *Gammaxxz,double *Gammayxz,double *Gammazxz,double *Gammaxyz,double *Gammayyz,double *Gammazyz,
 double *psi0re,double *psi0im)
{
  PsiKadelia_lowlevel(*cctkGH,cctk_lsh,cctk_nghostzones,
		      dx,dy,dz,r,x,y,z,
		      gxx,gxy,gxz,gyy,gyz,gzz,
		      gxx_phys,gxy_phys,gxz_phys,gyy_phys,gyz_phys,gzz_phys,
		      gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,
		      phi,psi,
		      Kxx,Kxy,Kxz,Kyy,Kyz,Kzz,
		      Rxx,Rxy,Rxz,Ryy,Ryz,Rzz,
		      Gammaxxx,Gammayxx,Gammazxx,Gammaxyy,Gammayyy,Gammazyy,
		      Gammaxzz,Gammayzz,Gammazzz,Gammaxxy,Gammayxy,Gammazxy,
		      Gammaxxz,Gammayxz,Gammazxz,Gammaxyz,Gammayyz,Gammazyz,
		      psi0re,psi0im);
}
