#include "cctk.h"
//#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "stdio.h"
#include "stdlib.h"

extern "C" void CCTK_FCALL CCTK_FNAME(em_extraction_lowlevel)
(const cGH **cctkGH,int *cctk_lsh,int *cctk_nghostzones,
 double &dx,double &dy,double &dz,double *r,double *x,double *y,double *z,
 double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
 double *lapm1,double *shiftx,double *shifty,double *shiftz,
 double *phi,double *psi,
 double *Bx,double *By,double *Bz,double *Ex,double *Ey,double *Ez,
 double *phi0re,double *phi0im,double *phi2re,double *phi2im,int *set_up_spherical_EM_wave_in_flat_spacetime,double *compute_phi2_min_radius,double *compute_phi2_max_radius, int *scale_with_radius_em,int *radius_power_em);


extern "C" void em_extraction_lowlevel(const cGH *cctkGH,int *cctk_lsh,int *cctk_nghostzones,
				       double &dx,double &dy,double &dz,double *r,double *x,double *y,double *z,
				       double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
				       double *lapm1,double *shiftx,double *shifty,double *shiftz,
				       double *phi,double *psi,
				       double *Bx,double *By,double *Bz,double *Ex,double *Ey,double *Ez,
				       double *phi0re,double *phi0im,double *phi2re,double *phi2im,int set_up_spherical_EM_wave_in_flat_spacetime, double compute_phi2_min_radius,double compute_phi2_max_radius,int scale_with_radius_em,int radius_power_em){


  /*
subroutine em_extraction_lowlevel(cctkGH,cctk_lsh,cctk_nghostzones, &
     dX,dY,dZ,r,x,y,z, &
     gxx,gxy,gxz,gyy,gyz,gzz, &
     lapm1,shiftx,shifty,shiftz, &
     phi,psi, &
     Bx,By,Bz,Ex,Ey,Ez, &
     phi0re,phi0im,phi2re,phi2im)
  implicit none
  */
  //  DECLARE_CCTK_PARAMETERS;

  int nx,ny,nz;
  CCTK_REAL one,two,fac;
  int i,j,k,local_spatial_order;
  int istart,iend,jstart,jend,kstart,kend;

  double phixL,phiyL,phizL,func;
  double psi_4thpower;
  double gxxL,gxyL,gxzL,gyyL,gyzL,gzzL;
  double psixL,psiyL,psizL;
  double shiftxL,shiftyL,shiftzL,lapse;
  double g4tt, g4tx, g4ty, g4tz, g4xx, g4xy, g4xz, g4yy, g4yz, g4zz;
  double Ftx, Fty, Ftz, Fxy, Fxz, Fyz, F_tx, F_ty, F_tz, F_xy, F_xz, F_yz;
  double BxL, ByL, BzL, B_xL, B_yL, B_zL, ExL, EyL, EzL, xL, yL, zL;
  double r2, dthx, dthy,dthz, dphx, dphy, dphz, rhoL;
  double drdotdth, drx, dry, drz, thx,thy, thz, th2, phx, phy,phz, ph2,drdotdph, thdotdph;
  double hatt_t, hatt_x, hatt_y, hatt_z;
  double hatr_x, hatr_y, hatr_z;
  double hatph_x, hatph_y, hatph_z;
  double hatth_x, hatth_y, hatth_z;
  double norm_hat_r, norm_hat_th, norm_hat_ph;
  double hatr_dot_hatth, hatr_dot_hatph, hatth_dot_hatph;
  double detgam, shift_xL, shift_yL, shift_zL;
  double F_hatt_hatth,F_hatr_hatth,F_hatt_hatph,F_hatr_hatph, Rephi0,Imphi0 ,Rephi2,Imphi2;

  printf("em_extraction begins...gxx[11] = %f\n",gxx[11]);

    //     Set up variables used in the grid loop for the physical grid points
  istart = cctk_nghostzones[0]+1;
  jstart = cctk_nghostzones[1]+1;
  kstart = cctk_nghostzones[2]+1;
  iend = cctk_lsh[0] - cctk_nghostzones[0];
  jend = cctk_lsh[1] - cctk_nghostzones[1];
  kend = cctk_lsh[2] - cctk_nghostzones[2];

  two = 2.0;
  one = 1.0;

#pragma omp parallel for
  for(k=0;k<cctk_lsh[2];k++) for(j=0;j<cctk_lsh[1];j++) for(i=0;i<cctk_lsh[0];i++) {
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	if(r[index] > compute_phi2_min_radius && r[index] < compute_phi2_max_radius) {
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
	  
	  // The 3 metric
	  if (set_up_spherical_EM_wave_in_flat_spacetime==1){
	    gxxL=1.0;
	    gxyL=0.0;
	    gxzL=0.0;
	    gyyL=1.0;
	    gyzL=0.0;
	    gzzL=1.0;
	    // The shift vector
	    shiftxL = 0.0;
	    shiftyL = 0.0;
	    shiftzL = 0.0;
	    // The lapse function
	    lapse = 1.0;

	  }
	  else{
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
	  }
	    // The coordinates
	    xL=x[index];
	    yL=y[index];
	    zL=z[index];
	  
	  //===============================================================================
	  // FIRST: Set up the tetrad vectors assuming spherical coordinates
	  // So that in flat spacetime, \hat r = e_{\hat r},
	  // \hat \phi = e_{\hat \phi}, \hat \theta = e_{\hat \theta}  
	  // centered at the origin of the coordinate system
	  //===============================================================================
	  
	  // Components of the radial vector
	  drx = xL;
	  dry = yL;
	  drz = zL;
	  
	  // Magnitude of (radial) position vector squared
	  r2 = gxxL*xL*xL + gyyL*yL*yL + gzzL*zL*zL + 2.0*(gxyL*xL*yL + gxzL*xL*zL + gyzL*yL*zL)+1.e-300;
	  
	  // Flat-space distance from Z-axis;
	  rhoL = sqrt(xL*xL+yL*yL)+1.e-300;
	  // cartesian components of partial_theta
	  dthx = xL*zL/rhoL;
	  dthy = yL*zL/rhoL;
	  dthz = -rhoL;
	  
	  // Start Gram-Schmidt orhonormalization process
	  // The timelike unit vector normal to t=const hypersurfaces
	  // will serve as the orthonormal timelike \hat t^a vector
	  // The radial vector \partial_r = x^i\partial_i is chosen next
	  // and since it is purely spatial it is already normal to \hat t^a
	  // so we only need to normalize it. 
	  // Next we orthonormalize \partial_\theta
	  
	  // \partial_r \dot \partial_\theta
	  
	  drdotdth = drx*dthx*gxxL + (dry*dthx + drx*dthy)*gxyL + (drz*dthx + drx*dthz)*gxzL 
	    + dry*dthy*gyyL + (drz*dthy + dry*dthz)*gyzL + drz*dthz*gzzL;
	  
	  // Cartesian Components of theta vector normal to \hat r, which is also normal to r
	  thx = dthx - (drdotdth/r2)*drx;
	  thy = dthy - (drdotdth/r2)*dry;
	  thz = dthz - (drdotdth/r2)*drz;
	  
	  // (Magnitude of theta)^2 
	  th2 = gxxL*thx*thx + gyyL*thy*thy + gzzL*thz*thz 
	    + 2.0*(gxyL*thx*thy + gxzL*thx*thz + gyzL*thy*thz);
	  
	  // Cartesian Components of vector \partial_\phi
	  dphx = -yL;
	  dphy = xL;
	  dphz = 0.0;
	  
	  
	  // \partial_r \dot \partial_\phi
	  drdotdph = dphx*drx*gxxL + (dphy*drx + dphx*dry)*gxyL + (dphz*drx + dphx*drz)*gxzL 
	    + dphy*dry*gyyL + (dphz*dry + dphy*drz)*gyzL + dphz*drz*gzzL;
	  
	  // \partial_\phi \dot \theta
	  thdotdph = dphx*thx*gxxL + (dphy*thx + dphx*thy)*gxyL + (dphz*thx + dphx*thz)*gxzL 
	    + dphy*thy*gyyL + (dphz*thy + dphy*thz)*gyzL + dphz*thz*gzzL;
	  
	  // Components of phi vector which is normal to both \hat\theta and \hat r
	  phx = dphx - (drdotdph/r2)*drx - (thdotdph/th2)*thx;
	  phy = dphy - (drdotdph/r2)*dry - (thdotdph/th2)*thy;
	  phz = dphz - (drdotdph/r2)*drz - (thdotdph/th2)*thz;
	  
	  // (magnitude of \phi)^2
	  ph2 = gxxL*phx*phx + gyyL*phy*phy + gzzL*phz*phz 
	    + 2.0*(gxyL*phx*phy + gxzL*phx*phz + gyzL*phy*phz);
	  
	  
	  // Orhonormal tetrad
	  
	  // \hat t^a = n^a
	  hatt_t = 1.0/lapse;
	  hatt_x = -shiftxL/lapse;
	  hatt_y = -shiftyL/lapse;
	  hatt_z = -shiftzL/lapse;
	  
	  // \hat r
	  hatr_x = drx/sqrt(r2);
	  hatr_y = dry/sqrt(r2);
	  hatr_z = drz/sqrt(r2);
	  
	  // \hat \theta
	  hatth_x = thx/sqrt(th2);
	  hatth_y = thy/sqrt(th2);
	  hatth_z = thz/sqrt(th2);
	  
	  // \hat \phi
	  hatph_x = phx/sqrt(ph2);
	  hatph_y = phy/sqrt(ph2);
	  hatph_z = phz/sqrt(ph2);
	  
	  
	  // Check orthonormality of the tetrad
	  
	  norm_hat_r = gxxL*hatr_x*hatr_x + 2.0*gxyL*hatr_x*hatr_y + gyyL*hatr_y*hatr_y 
	    + 2.0*gxzL*hatr_x*hatr_z + 2.0*gyzL*hatr_y*hatr_z + gzzL*hatr_z*hatr_z;
	  
	  norm_hat_th = gxxL*hatth_x*hatth_x + 2.0*gxyL*hatth_x*hatth_y + gyyL*hatth_y*hatth_y 
	    + 2.0*gxzL*hatth_x*hatth_z + 2.0*gyzL*hatth_y*hatth_z + gzzL*hatth_z*hatth_z;
	  
	  norm_hat_ph = gxxL*hatph_x*hatph_x + 2.0*gxyL*hatph_x*hatph_y + gyyL*hatph_y*hatph_y 
	    + 2.0*gxzL*hatph_x*hatph_z + 2.0*gyzL*hatph_y*hatph_z + gzzL*hatph_z*hatph_z;
	  
	  if ( fabs(norm_hat_r-1.0) > 1.e-13 || fabs(norm_hat_th-1.0) > 1.e-13 || fabs(norm_hat_ph-1.0) > 1.e-13) {
	    printf("Stopping: tetrad is not normalized properly in em_extraction_lowlevel\n");
	    printf("|hat r|=%16.14f, |hat th|= %16.14f, |hat ph|= %16.14f",norm_hat_r,norm_hat_th,norm_hat_ph);
	    printf("at grid point x, y, z = %f, %f, %f\n", xL, yL, zL);
	    printf("hat r dot hat th= %16.14f, hat r dot hat ph= %16.14f, hat th dot hat ph= %16.14f\n", hatr_dot_hatth, hatr_dot_hatph,hatth_dot_hatph);
	    printf("hat thx= %16.14f,  hat thy = %16.14f ,   hat thz=  %16.14f\n", hatth_x,hatth_y,hatth_z );
	    printf("hat phx= %16.14f,  hat phy = %16.14f ,   hat phz=  %16.14f\n", hatph_x,hatph_y,hatph_z );
	    printf("hat rx= %16.14f,  hat ry = %16.14f ,   hat rz=  %16.14f\n", hatr_x,hatr_y,hatr_z );
	    printf("dthx=  %16.14f, dthy= %16.14f , dthz= %16.14f\n", dthx,dthy,dthz );
	    printf("thx= %16.14f, thy=%16.14f, thz=%16.14f\n", thx,thy,thz );
	    printf("gxx= %16.14f,  gxy=%16.14f, gxz=%16.14f\n", gxxL,gxyL,gxzL);
	    printf("gyy= %16.14f,  gyz=%16.14f, gzz=%16.14f\n", gyyL,gyzL,gzzL);
	    printf("psi4=%16.14f, psi4_b=%16.14f\n",psi_4thpower,exp(4.0*phi[index]));		  

	    exit(1);
	  }
	  
	  hatr_dot_hatth = gxxL*hatr_x*hatth_x + gxyL*hatr_y*hatth_x + gxzL*hatr_z*hatth_x 
	    + gxyL*hatr_x*hatth_y + gyyL*hatr_y*hatth_y + gyzL*hatr_z*hatth_y 
	    + gxzL*hatr_x*hatth_z + gyzL*hatr_y*hatth_z + gzzL*hatr_z*hatth_z;
	  
	  hatr_dot_hatph = gxxL*hatr_x*hatph_x + gxyL*(hatr_y*hatph_x + hatr_x*hatph_y)
	    + gxzL*(hatr_z*hatph_x + hatr_x*hatph_z) + gyyL*hatr_y*hatph_y 
	    + gyzL*(hatr_z*hatph_y + hatr_y*hatph_z) + gzzL*hatr_z*hatph_z;
	  
	  hatth_dot_hatph = gxxL*hatth_x*hatph_x + gxyL*(hatth_y*hatph_x + hatth_x*hatph_y)
	    + gxzL*(hatth_z*hatph_x + hatth_x*hatph_z) + gyyL*hatth_y*hatph_y 
	    + gyzL*(hatth_z*hatph_y + hatth_y*hatph_z) + gzzL*hatth_z*hatph_z;
	  
	  if (fabs(hatr_dot_hatth) > 1.e-13 || fabs(hatr_dot_hatph) > 1.e-13 || fabs(hatth_dot_hatph) > 1.e-13){
	    printf("Stopping: tetrad is not orthogonal in em_extraction_lowlevel\n");
	    printf("at grid point x, y, z = %f, %f, %f\n", xL, yL, zL);
	    printf("hat r dot hat th= %16.14f, hat r dot hat ph= %16.14f, hat th dot hat ph= %16.14f\n", hatr_dot_hatth, hatr_dot_hatph,hatth_dot_hatph);
	    printf("hat thx= %16.14f,  hat thy = %16.14f ,   hat thz=  %16.14f\n", hatth_x,hatth_y,hatth_z );
	    printf("hat phx= %16.14f,  hat phy = %16.14f ,   hat phz=  %16.14f\n", hatph_x,hatph_y,hatph_z );
	    printf("hat rx= %16.14f,  hat ry = %16.14f ,   hat rz=  %16.14f\n", hatr_x,hatr_y,hatr_z );
	    printf("dthx=  %16.14f, dthy= %16.14f , dthz= %16.14f\n", dthx,dthy,dthz );
	    printf("thx= %16.14f, thy=%16.14f, thz=%16.14f\n", thx,thy,thz );
	    printf("gxx= %16.14f,  gxy=%16.14f, gxz=%16.14f\n", gxxL,gxyL,gxzL);
	    printf("gyy= %16.14f,  gyz=%16.14f, gzz=%16.14f\n", gyyL,gyzL,gzzL);
	    printf("psi4=%16.14f, psi4_b=%16.14f\n",psi_4thpower,exp(4.0*phi[index]));		  
	    exit(1);
	  }
	  // =======================================================================================
	  //          End of orthonormal tetrad calculation
	  // =======================================================================================
	  
	  // We have the components of the orthonormal tetrad in the cartesian coordinate basis
	  // The following are (contravariant) vector components
	  // hatt = (hatt_t, hatt_x hatt_y, hatt_z) = n^a ,i.e., the time like normal vector	
	  // hatr = (0, hatr_x, hatr_y, hatr_z)
	  // hatth = (0, hatth_x, hatth_y, hatth_z)
	  // hatph = (0, hatph_x, hatph_y, hatph_z)
	  
	  if( xL*yL != 0.){
	    //     Calculate only if not on the z-axis for the "radial" tetrad case
	    
	    /*	    printf("r= %16.14f,%16.14f, %16.14f, x y z= %16.14f,%16.14f, %16.14f\n", hatr_x,hatr_y,hatr_z,xL,yL,zL);
	    printf("th= %16.14f,%16.14f, %16.14f, x y z= %16.14f,%16.14f, %16.14f\n", hatth_x,hatth_y,hatth_z,xL,yL,zL);
	    printf("ph= %16.14f,%16.14f, %16.14f, x y z= %16.14f,%16.14f, %16.14f\n", hatph_x,hatph_y,hatph_z,xL,yL,zL);

	    exit(1);*/
	    
	    // Now, we turn to the calculation of the Faraday tensor F_\mu\nu
	    
	    if (set_up_spherical_EM_wave_in_flat_spacetime==1){
	      double norm_k = 0.1155; //Wavenumber. Yields a wavelength of 54.4
	      double kixi = norm_k*sqrt(r2);

	      double sinth = rhoL/sqrt(r2);

	      BxL = -hatph_x*sinth*cos(kixi)/sqrt(r2);   // These are the B^i
	      ByL = -hatph_y*sinth*cos(kixi)/sqrt(r2);
	      BzL = -hatph_z*sinth*cos(kixi)/sqrt(r2);
	      
	      ExL = -hatth_x*sinth*cos(kixi)/sqrt(r2);   // These are the E^i
	      EyL = -hatth_y*sinth*cos(kixi)/sqrt(r2);
	      EzL = -hatth_z*sinth*cos(kixi)/sqrt(r2);

	    }
	    else{
	      BxL = Bx[index];   // These are the B^i
	      ByL = By[index];
	      BzL = Bz[index];
	      
	      ExL = Ex[index];   // These are the E^i
	      EyL = Ey[index];
	      EzL = Ez[index];
	    }
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
	    
	    
	    // Calculate the Faraday tensor F dotted into the orhonormal tetrad 
	    // In slot notation of two-forms F(\hat t,\hat \theta)
	    F_hatt_hatth = F_tx*hatth_x*hatt_t + F_ty*hatth_y*hatt_t + F_tz*hatth_z*hatt_t 
	      + F_xy*hatth_y*hatt_x + F_xz*hatth_z*hatt_x - F_xy*hatth_x*hatt_y 
	      + F_yz*hatth_z*hatt_y - F_xz*hatth_x*hatt_z - F_yz*hatth_y*hatt_z;
	    
	    // In slot notation of two-forms F(\hat r,\hat \theta)
	    F_hatr_hatth = -F_xy*hatr_y*hatth_x - F_xz*hatr_z*hatth_x + F_xy*hatr_x*hatth_y 
	      - F_yz*hatr_z*hatth_y + F_xz*hatr_x*hatth_z + F_yz*hatr_y*hatth_z;
	    
	    // In slot notation of two-forms F(\hat t,\hat \phi)
	    F_hatt_hatph = F_tx*hatph_x*hatt_t + F_ty*hatph_y*hatt_t + F_tz*hatph_z*hatt_t 
	      + F_xy*hatph_y*hatt_x + F_xz*hatph_z*hatt_x - F_xy*hatph_x*hatt_y 
	      + F_yz*hatph_z*hatt_y - F_xz*hatph_x*hatt_z - F_yz*hatph_y*hatt_z;
	    
	    // In slot notation of two-forms F(\hat r,\hat \phi)
	    F_hatr_hatph = -F_xy*hatr_y*hatph_x - F_xz*hatr_z*hatph_x + F_xy*hatr_x*hatph_y 
	      - F_yz*hatr_z*hatph_y + F_xz*hatr_x*hatph_z + F_yz*hatr_y*hatph_z;
	    
	    
	    Rephi0 = 0.5*(F_hatt_hatth + F_hatr_hatth);
	    Imphi0 = 0.5*(F_hatt_hatph + F_hatr_hatph);
	    
	    Rephi2 = 0.5*(-F_hatt_hatth + F_hatr_hatth);
	    Imphi2 = 0.5*(F_hatt_hatph - F_hatr_hatph);
	    	    
	    //     If required scale by radial coordinate
	    if (scale_with_radius_em == 1){
	      fac = pow(r[index],radius_power_em);
	    }
	    else{
	      fac = 1.0;
	    }
	    
	    //     Phis we do pairwise
	    phi0re[index] = fac*Rephi0;
	    phi0im[index] = fac*Imphi0;
	    phi2re[index] = fac*Rephi2;
	    phi2im[index] = fac*Imphi2;


	    /*	    if (sqrt(r2) > 95. && sqrt(r2) < 100.){ 
	    printf("Stopping: Rephi2 = %15.14f , 1/r = %15.14f, Imphi2 = %15.14f, at r = %15.14f \n", phi2re[index],1./sqrt(r2),phi2im[index],sqrt(r2));
	    }
	    */

	  }
	  else{
	    //     Set values to zero without calculating on the z-axis
	    phi0re[index] = 0.0;
	    phi0im[index] = 0.0;
	    phi2re[index] = 0.0;
	    phi2im[index] = 0.0;
	  }
	}
	else{ // (r[index].lt.compute_phi2_min_radius .and. r[index].gt.compute_phi2_max_radius)
	  phi0re[index] = 0.0;
	  phi0im[index] = 0.0;
	  phi2re[index] = 0.0;
	  phi2im[index] = 0.0;
	}
	
	
      }
  
  printf("outside em_extraction_lowlevel()\n");
}


extern "C" void CCTK_FCALL CCTK_FNAME(em_extraction_lowlevel)
(const cGH **cctkGH,int *cctk_lsh,int *cctk_nghostzones,
 double &dx,double &dy,double &dz,double *r,double *x,double *y,double *z,
 double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
 double *lapm1,double *shiftx,double *shifty,double *shiftz,
 double *phi,double *psi,
 double *Bx,double *By,double *Bz,double *Ex,double *Ey,double *Ez,
 double *phi0re,double *phi0im,double *phi2re,double *phi2im, int *set_up_spherical_EM_wave_in_flat_spacetime,double *compute_phi2_min_radius,double *compute_phi2_max_radius,int *scale_with_radius_em,int *radius_power_em)
{
  em_extraction_lowlevel(*cctkGH,cctk_lsh,cctk_nghostzones,
			 dx,dy,dz,r,x,y,z,gxx,gxy,gxz,gyy,gyz,gzz,lapm1,shiftx,shifty,shiftz,
			 phi,psi,Bx,By,Bz,Ex,Ey,Ez,phi0re,phi0im,phi2re,phi2im,*set_up_spherical_EM_wave_in_flat_spacetime,*compute_phi2_min_radius,*compute_phi2_max_radius,*scale_with_radius_em,*radius_power_em);
}

