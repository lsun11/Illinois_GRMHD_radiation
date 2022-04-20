#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "stdio.h"
#include "stdlib.h"

extern "C" void CCTK_FCALL ang_freq_of_mag_field_lines_
  (const cGH **cctkGH,int *cctk_lsh,int *cctk_nghostzones,
   double &dx,double &dy,double &dz,double *r,double *x,double *y,double *z,
   double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
   double *lapm1,double *shiftx,double *shifty,double *shiftz,
   double *phi,double *psi,
   double *Bx,double *By,double *Bz,double *Ex,double *Ey,double *Ez,double *Bfreq1, double *Bfreq2, double &bhx, double &bhy, double &bhz);

extern "C" void ang_freq_of_mag_field_lines(const cGH *cctkGH,int *cctk_lsh,int *cctk_nghostzones,
				       double &dx,double &dy,double &dz,double *r,double *x,double *y,double *z,
				       double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
				       double *lapm1,double *shiftx,double *shifty,double *shiftz,
				       double *phi,double *psi,
				       double *Bx,double *By,double *Bz,double *Ex,double *Ey,double *Ez, double *Bfreq1, double *Bfreq2,
				       double &bhx, double &bhy, double &bhz){
  
  DECLARE_CCTK_PARAMETERS

  int istart,iend,jstart,jend,kstart,kend;
  int istart2,iend2,jstart2,jend2,kstart2,kend2;


  printf("em_extraction begins...gxx[11] = %f\n",gxx[11]);
  printf("n_ghostx = %d, n_ghosty = %d, n_ghostz = %d\n",cctk_nghostzones[0], cctk_nghostzones[1], cctk_nghostzones[2]);

    //     Set up variables used in the grid loop for the physical grid points
  istart = cctk_nghostzones[0]+1;
  jstart = cctk_nghostzones[1]+1;
  kstart = cctk_nghostzones[2]+1;
  iend = cctk_lsh[0] - cctk_nghostzones[0];
  jend = cctk_lsh[1] - cctk_nghostzones[1];
  kend = cctk_lsh[2] - cctk_nghostzones[2];

  double PI = 3.14159265358979323846;

#pragma omp parallel for
  //  for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
  for(int k=1;k<cctk_lsh[2]-1;k++) for(int j=1;j<cctk_lsh[1]-1;j++) for(int i=1;i<cctk_lsh[0]-1;i++) {
  //  for(k=kstart;k<kend;k++) for(j=jstart;j<jend;j++) for(i=istart;i<iend;i++){ 
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);


	  double fac;
	  double phixL,phiyL,phizL,func;
	  double psi_4thpower;
	  double gxxL,gxyL,gxzL,gyyL,gyzL,gzzL;
	  double psixL,psiyL,psizL;
	  double shiftxL,shiftyL,shiftzL,lapse;
	  double g4tt, g4tx, g4ty, g4tz, g4xx, g4xy, g4xz, g4yy, g4yz, g4zz;
	  double Ftx, Fty, Ftz, Fxy, Fxz, Fyz, F_tx, F_ty, F_tz, F_xy, F_xz, F_yz;
	  double F_tr, F_tth, F_rph, F_thph;
	  double BxL, ByL, BzL, B_xL, B_yL, B_zL, ExL, EyL, EzL, xL, yL, zL;
	  double r2, dthx, dthy,dthz, dphx, dphy, dphz, rhoL;
	  double drdotdth, drx, dry, drz, thx,thy, thz, th2, phx, phy,phz, ph2,drdotdph, thdotdph;
	  double hatt_t, hatt_x, hatt_y, hatt_z;
	  double r, varpi;
	  double hatph_x, hatph_y, hatph_z;
	  double hatth_x, hatth_y, hatth_z;
	  double norm_hat_r, norm_hat_th, norm_hat_ph;
	  double hatr_dot_hatth, hatr_dot_hatph, hatth_dot_hatph;
	  double detgam, shift_xL, shift_yL, shift_zL;
	  double F_hatt_hatth,F_hatr_hatth,F_hatt_hatph,F_hatr_hatph, Rephi0,Imphi0 ,Rephi2,Imphi2;
	  
	  psi_4thpower = exp(4.0*phi[index]);
	  /*
	  // Vasilis says: We shouldn't need the following 6 lines
	  local_spatial_order=4
	  if(i==1.or.j==1.or.k==1.or.i==cctk_lsh(1).or.j==cctk_lsh(2).or.k==cctk_lsh(3)) then
	  local_spatial_order=0
	  else if(i==2.or.j==2.or.k==2.or.i==cctk_lsh(1)-1.or.j==cctk_lsh(2)-1.or.k==cctk_lsh(3)-1) then
	  local_spatial_order=2
	  end if
	  */
	  
	  //    Define local variables to calculate the orthonormal tetrad
	  

	  gxxL=gxx[index]*psi_4thpower;
	  gxyL=gxy[index]*psi_4thpower;
	  gxzL=gxz[index]*psi_4thpower;
	  gyyL=gyy[index]*psi_4thpower;
	  gyzL=gyz[index]*psi_4thpower;
	  gzzL=gzz[index]*psi_4thpower;
	  // The shift vector
	  shiftxL = shiftx[index];
	  shiftyL = shifty[index];
	  shiftzL = shiftz[index];
	  // The lapse function
	  lapse = lapm1[index] + 1.0;
	
	  // The coordinates
	  xL=x[index] - bhx;
	  yL=y[index] - bhy;
	  zL=z[index] - bhz;
	  
	  // Spherical Coordinates
	  double rL = sqrt(xL*xL + yL*yL + zL*zL);
	  double varpiL = sqrt(xL*xL + yL*yL);
	  


	  // normal timelike vector
	  
	  // \hat t^a = n^a
	  hatt_t = 1.0/lapse;
	  hatt_x = -shiftxL/lapse;
	  hatt_y = -shiftyL/lapse;
	  hatt_z = -shiftzL/lapse;
	  
	  // hatt = (hatt_t, hatt_x hatt_y, hatt_z) = n^a ,i.e., the time like normal vector	
	    
	    // Now, we turn to the calculation of the Faraday tensor F_\mu\nu
	    
	  BxL = Bx[index];   // These are the B^i
	  ByL = By[index];
	  BzL = Bz[index];
	  
	  ExL = Ex[index];   // These are the E^i
	  EyL = Ey[index];
	  EzL = Ez[index];
	
	  // First we need the covariant components of the B-field
	  // B_i = \gamma_ij B^j
	    
	  B_xL = gxxL*BxL + gxyL*ByL + gxzL*BzL;
	  B_yL = gxyL*BxL + gyyL*ByL + gyzL*BzL;
	  B_zL = gxzL*BxL + gyzL*ByL + gzzL*BzL;
	    
	  // We will also need the determinant of the 3-metric
	  
	  detgam = -gxzL*gxzL*gyyL + 2.0*gxyL*gxzL*gyzL 
	    -gxxL*gyzL*gyzL - gxyL*gxyL*gzzL + gxxL*gyyL*gzzL;
	  
	  // Now compute the Maxwell tensor, F^{munu} (contravariant components) 
	  
	  Ftx = hatt_t*ExL;
	  Fty = hatt_t*EyL;
	  Ftz = hatt_t*EzL;
	  Fxy = hatt_x*EyL - hatt_y*ExL + B_zL/sqrt(detgam);
	  Fxz = hatt_x*EzL - hatt_z*ExL - B_yL/sqrt(detgam);
	  Fyz = hatt_y*EzL - hatt_z*EyL + B_xL/sqrt(detgam);
	  
	  
	  // Compute 4-metric
	  shift_xL = gxxL*shiftxL + gxyL*shiftyL + gxzL*shiftzL;
	  shift_yL = gxyL*shiftxL + gyyL*shiftyL + gyzL*shiftzL;
	  shift_zL = gxzL*shiftxL + gyzL*shiftyL + gzzL*shiftzL;
	  
	  
	  g4tt = -lapse*lapse + shift_xL*shiftxL + shift_yL*shiftyL + shift_zL*shiftzL;
	  g4tx = shift_xL;
	  g4ty = shift_yL;
	  g4tz = shift_zL;
	  g4xx = gxxL;
	  g4xy = gxyL;
	  g4xz = gxzL;
	  g4yy = gyyL;
	  g4yz = gyzL;
	  g4zz = gzzL;
	  
	    
	  // Compute F_munu (covariant components)
	  
	  F_tx = -Ftx*g4tx*g4tx - Fty*g4tx*g4ty - Ftz*g4tx*g4tz + Ftx*g4tt*g4xx 
	    -Fxy*g4ty*g4xx - Fxz*g4tz*g4xx + Fty*g4tt*g4xy + Fxy*g4tx*g4xy 
	    -Fyz*g4tz*g4xy + Ftz*g4tt*g4xz + Fxz*g4tx*g4xz + Fyz*g4ty*g4xz;

	  F_ty = -Ftx*g4tx*g4ty - Fty*g4ty*g4ty - Ftz*g4ty*g4tz + Ftx*g4tt*g4xy 
	    -Fxy*g4ty*g4xy - Fxz*g4tz*g4xy + Fty*g4tt*g4yy + Fxy*g4tx*g4yy 
	    -Fyz*g4tz*g4yy + Ftz*g4tt*g4yz + Fxz*g4tx*g4yz + Fyz*g4ty*g4yz;

	  F_tz = -Ftx*g4tx*g4tz - Fty*g4ty*g4tz - Ftz*g4tz*g4tz + Ftx*g4tt*g4xz 
	    -Fxy*g4ty*g4xz - Fxz*g4tz*g4xz + Fty*g4tt*g4yz + Fxy*g4tx*g4yz 
	    -Fyz*g4tz*g4yz + Ftz*g4tt*g4zz + Fxz*g4tx*g4zz + Fyz*g4ty*g4zz;

	  F_xy = -Ftx*g4ty*g4xx + Ftx*g4tx*g4xy - Fty*g4ty*g4xy - Fxy*g4xy*g4xy 
	    -Ftz*g4ty*g4xz - Fxz*g4xy*g4xz + Fty*g4tx*g4yy + Fxy*g4xx*g4yy 
	    -Fyz*g4xz*g4yy + Ftz*g4tx*g4yz + Fxz*g4xx*g4yz + Fyz*g4xy*g4yz;

	  F_xz = -Ftx*g4tz*g4xx - Fty*g4tz*g4xy + Ftx*g4tx*g4xz - Ftz*g4tz*g4xz 
	    -Fxy*g4xy*g4xz - Fxz*g4xz*g4xz + Fty*g4tx*g4yz + Fxy*g4xx*g4yz 
	    -Fyz*g4xz*g4yz + Ftz*g4tx*g4zz + Fxz*g4xx*g4zz + Fyz*g4xy*g4zz;

	  F_yz = -Ftx*g4tz*g4xy + Ftx*g4ty*g4xz - Fty*g4tz*g4yy - Fxy*g4xz*g4yy 
	    +Fty*g4ty*g4yz - Ftz*g4tz*g4yz + Fxy*g4xy*g4yz - Fxz*g4xz*g4yz 
	    -Fyz*g4yz*g4yz + Ftz*g4ty*g4zz + Fxz*g4xy*g4zz + Fyz*g4yy*g4zz;


	  // Compute covariant  Ftr, Ftth, Frph, Fthph using the standard coordinate transformation law

	  F_tr = xL*F_tx/rL + yL*F_ty/rL + zL*F_tz/rL;
	  F_tth = xL*zL*F_tx/varpiL + yL*zL*F_ty/varpiL - varpiL*F_tz;
	  F_rph = -varpiL*varpiL*F_xy/rL - yL*zL*F_xz/rL + xL*zL*F_yz/rL;
	  F_thph = -varpiL*zL*F_xy + yL*varpiL*F_xz - xL*varpiL*F_yz;
	  
	  Bfreq1[index] = - F_tr/F_rph;
	  Bfreq2[index] = - F_tth/F_thph;


      }
  
  printf("outside ang_freq_of_mag_field_lines()\n");
}


extern "C" void CCTK_FCALL ang_freq_of_mag_field_lines_
(const cGH **cctkGH,int *cctk_lsh,int *cctk_nghostzones,
 double &dx,double &dy,double &dz,double *r,double *x,double *y,double *z,
 double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
 double *lapm1,double *shiftx,double *shifty,double *shiftz,
 double *phi,double *psi,
 double *Bx,double *By,double *Bz,double *Ex,double *Ey,double *Ez,double *Bfreq1,double *Bfreq2,
 double &bhx, double &bhy, double &bhz)
{
  ang_freq_of_mag_field_lines(*cctkGH,cctk_lsh,cctk_nghostzones,
			 dx,dy,dz,r,x,y,z,gxx,gxy,gxz,gyy,gyz,gzz,lapm1,shiftx,shifty,shiftz,
			      phi,psi,Bx,By,Bz,Ex,Ey,Ez,Bfreq1,Bfreq2,bhx,bhy,bhz);
}
