//-----------------------------------------------------------------------------
// Add third-order accurate curvature terms to S_i_rhs
//-----------------------------------------------------------------------------
#include "math.h"
#include "cctk.h"
#include "stdio.h"
#include "primitives_solver_header.h"
//#include "compute_T_fluid.C"

#define SQR(x) ((x) * (x))
#define D1_o3_cpp(gf,i,j,k) \
  ((   gf[CCTK_GFINDEX3D(cctkGH,i+1,j,k)] \
       - gf[CCTK_GFINDEX3D(cctkGH,i ,j,k)] ) * dxi )
#define D2_o3_cpp(gf,i,j,k) \
  ((   gf[CCTK_GFINDEX3D(cctkGH,i,j+1,k)] \
       - gf[CCTK_GFINDEX3D(cctkGH,i ,j,k)] ) * dyi )
#define D3_o3_cpp(gf,i,j,k) \
  ((   gf[CCTK_GFINDEX3D(cctkGH,i,j,k+1)] \
       - gf[CCTK_GFINDEX3D(cctkGH,i ,j,k)] ) * dzi )

double Fermi4p(double &eta);
double Fermi4m(double &eta);
double Fermi4(double &eta);

double Fermi5p(double &eta);
double Fermi5m(double &eta);
double Fermi5(double &eta);

//void compute_opacity_emissivity(double &rho_star, double &al, double &Psi6, double &u0, double &kappa_a, double &kappa_s, double &eta_gfL, double &T_fluid, double &chi_rad);

extern "C" void CCTK_FCALL mhd_source_z_tau_cpp_
  (int *flux_direction, const cGH **cctkGH,int *cctk_lsh, int *nghostzones, int *Symmetry,
   int *enable_HARM_energyvariable,
   double *dX, double *dY,double *dZ, double *Z,
   double *mhd_st_x_rhs, double *mhd_st_y_rhs, double *mhd_st_z_rhs, 
   double *tau_rhs, double *rho_star, 
   double *P, double *h, double *u0, 
   double *vx, double *vy, double *vz,  
   double *sbt, double *sbx, double *sby, double *sbz, 
   double *alpha, double *betax, double *betay, double *betaz, 
   double *phi, 
   double *gxx, double *gxy, double *gxz, double *gyy, double *gyz, double *gzz, 
   double *gupxx, double *gupxy, double *gupxz, double *gupyy, double *gupyz, double *gupzz, 
   double *alpha_f, double *betax_f, double *betay_f, double *betaz_f, 
   double *phi_f, 				 
   double *gxx_f, double *gxy_f, double *gxz_f, double *gyy_f, double *gyz_f, double *gzz_f,
   double *g4tt,double *g4tx,double *g4ty,double *g4tz,
   double *g4tt_f,double *g4tx_f,double *g4ty_f,double *g4tz_f,
   double *E_rad, double *F_radx, double *F_rady, double *F_radz, double *Y_e, double *optd, double *eta_nue, 
   double *E_rad_nue, double *F_radx_nue, double *F_rady_nue, double *F_radz_nue, 
   double *E_rad_nux, double *F_radx_nux, double *F_rady_nux, double *F_radz_nux,
   int &rad_fourforce_enable, int *rad_closure_scheme, double &rad_opacity_abs, double &rad_opacity_sct, double &rad_const, 
   double *chi_rad, double *chi_rad_nue,
   double *T_fluid,  int &enable_OS_collapse, int &compute_microphysics, int &microphysics_scheme, double &T_fluid_cgs_atm, int &rad_fix);

extern "C" void mhd_source_z_tau_cpp(int flux_direction, const cGH *cctkGH,int *cctk_lsh, int *nghostzones, int Symmetry,
				     int enable_HARM_energyvariable,
				     double dX, double dY,double dZ, double *Z,
				     double *mhd_st_x_rhs, double *mhd_st_y_rhs, double *mhd_st_z_rhs, 
				     double *tau_rhs, double *rho_star, 
				     double *P, double *h, double *u0, 
				     double *vx, double *vy, double *vz,  
				     double *sbt, double *sbx, double *sby, double *sbz, 
				     double *alpha, double *betax, double *betay, double *betaz, 
				     double *phi, 
				     double *gxx, double *gxy, double *gxz, double *gyy, double *gyz, double *gzz, 
				     double *gupxx, double *gupxy, double *gupxz, double *gupyy, double *gupyz, double *gupzz, 
				     double *alpha_f, double *betax_f, double *betay_f, double *betaz_f, 
				     double *phi_f, 				 
				     double *gxx_f, double *gxy_f, double *gxz_f, double *gyy_f, double *gyz_f, double *gzz_f,
				     double *g4tt,double *g4tx,double *g4ty,double *g4tz,
				     double *g4tt_f,double *g4tx_f,double *g4ty_f,double *g4tz_f,
				     double *E_rad, double *F_radx, double *F_rady, double *F_radz, double *Y_e, double *optd, double *eta_nue, 
				     double *E_rad_nue, double *F_radx_nue, double *F_rady_nue, double *F_radz_nue,
				     double *E_rad_nux, double *F_radx_nux, double *F_rady_nux, double *F_radz_nux,
				     int &rad_fourforce_enable, int *rad_closure_scheme, double &rad_opacity_abs, double &rad_opacity_sct, double &rad_const, 
				     double *chi_rad, double *chi_rad_nue,
				     double *T_fluid,  int &enable_OS_collapse, int &compute_microphysics, int &microphysics_scheme, double &T_fluid_cgs_atm, int &rad_fix) {

  //  printf("START mhd_source_z_tau_cpp!!!!!!!!!!!!!!! \n");

  double sfpi = sqrt(4.0*M_PI);

  int AXISYM = 4;

  double f1os4pi = 1.0/sqrt(4.0*M_PI);

  /* Set up variables used in the grid loop for the physical grid points */
  //NOTE: _f variables are defined everywhere except ijkmin, depending on flux_direction
  //     I.e., flux_direction==1 -> _f imin's are not defined
  //   Further, the D123_o3_cpp derivatives cannot read beyond ijk=cctk_lsh[012], depending on flux_direction again.
  int istart = 1;
  int jstart = 1;
  int kstart = 1;
  int iend = cctk_lsh[0]-1;
  int jend = cctk_lsh[1]-1;
  int kend = cctk_lsh[2]-1;

  double dxi = 1.0/dX;
  double dyi = 1.0/dY;
  double dzi = 1.0/dZ;

  if(Symmetry==4) {
    jstart = 0;
    jend = cctk_lsh[1];
    jstart++;
    jend--;
  }
  
#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    double Psi4 = exp(4.0*phi[index]);
    double al = 1.0 + alpha[index];

    double Psi4_f = exp(4.0*phi_f[index]);
    double al_f = 1.0 + alpha_f[index];

    double gxxL = gxx[index];
    double gxyL = gxy[index];
    double gxzL = gxz[index];
    double gyyL = gyy[index];
    double gyzL = gyz[index];
    double gzzL = gzz[index];

    double gxx_fL = gxx_f[index];
    double gxy_fL = gxy_f[index];
    double gxz_fL = gxz_f[index];
    double gyy_fL = gyy_f[index];
    double gyz_fL = gyz_f[index];
    double gzz_fL = gzz_f[index];

    double betaxL = betax[index];
    double betayL = betay[index];
    double betazL = betaz[index];

    double betax_fL = betax_f[index];
    double betay_fL = betay_f[index];
    double betaz_fL = betaz_f[index];
    //----------------------------------------------------------------------------
    // Compute the 4-metric, converting gij's from conformal to physical 
    //----------------------------------------------------------------------------
    g4tt[index] = -SQR(al) + Psi4*(gxxL*SQR(betaxL) + 2.0*gxyL*betaxL*betayL + 2.0*gxzL*betaxL*betazL + 
				   gyyL*SQR(betayL) + 2.0*gyzL*betayL*betazL + gzzL*SQR(betazL));
    g4tx[index] = Psi4*(gxxL*betaxL + gxyL*betayL + gxzL*betazL);
    g4ty[index] = Psi4*(gxyL*betaxL + gyyL*betayL + gyzL*betazL);
    g4tz[index] = Psi4*(gxzL*betaxL + gyzL*betayL + gzzL*betazL);
    gxx[index] = Psi4*gxxL;
    gxy[index] = Psi4*gxyL;
    gxz[index] = Psi4*gxzL;
    gyy[index] = Psi4*gyyL;
    gyz[index] = Psi4*gyzL;
    gzz[index] = Psi4*gzzL;

    g4tt_f[index] = -SQR(al_f) + Psi4_f*(gxx_fL*SQR(betax_fL) + 2.0*gxy_fL*betax_fL*betay_fL + 2.0*gxz_fL*betax_fL*betaz_fL + 
					 gyy_fL*SQR(betay_fL) + 2.0*gyz_fL*betay_fL*betaz_fL + gzz_fL*SQR(betaz_fL));
    g4tx_f[index] = Psi4_f*(gxx_fL*betax_fL + gxy_fL*betay_fL + gxz_fL*betaz_fL);
    g4ty_f[index] = Psi4_f*(gxy_fL*betax_fL + gyy_fL*betay_fL + gyz_fL*betaz_fL);
    g4tz_f[index] = Psi4_f*(gxz_fL*betax_fL + gyz_fL*betay_fL + gzz_fL*betaz_fL);
    gxx_f[index] = Psi4_f*gxx_fL;
    gxy_f[index] = Psi4_f*gxy_fL;
    gxz_f[index] = Psi4_f*gxz_fL;
    gyy_f[index] = Psi4_f*gyy_fL;
    gyz_f[index] = Psi4_f*gyz_fL;
    gzz_f[index] = Psi4_f*gzz_fL;


  }

  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!           
  // To modifiers: The function above have already convert gij to physical. NO ADDTIONAL PSI4 NEEDED!!!!!!           
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

  if(flux_direction==1) {
    //-----------------------------------------------------------------------------
    // Compute the addition to mhd_st_x_rhs and tau_rhs
    //-----------------------------------------------------------------------------
#pragma omp parallel for
    for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
      int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
      double al = 1.0 + alpha[index];
      double al_f = 1.0 + alpha_f[index];
      double betaxL = betax[index];
      double betayL = betay[index];
      double betazL = betaz[index];
      double Psi = exp(phi[index]);
      double Psi4 = Psi*Psi*Psi*Psi;
      double Psim4 = 1.0/Psi4;
      double Psi6 = Psi4*Psi*Psi;
      double gxxL = gxx[index];
      double gxyL = gxy[index];
      double gxzL = gxz[index];
      double gyyL = gyy[index];
      double gyzL = gyz[index];
      double gzzL = gzz[index];

      double gupxxL = gupxx[index];
      double gupxyL = gupxy[index];
      double gupxzL = gupxz[index];
      double gupyyL = gupyy[index];
      double gupyzL = gupyz[index];
      double gupzzL = gupzz[index];
      // Divide b^{\mu} by alpha sqrt(4 pi) 
      double sbtL = sbt[index]/al/sfpi;
      double sbxL = sbx[index]/al/sfpi;
      double sbyL = sby[index]/al/sfpi;
      double sbzL = sbz[index]/al/sfpi;
 
      double rho_starL = rho_star[index];
      double u0L = u0[index];
      double vxL = vx[index];
      double vyL = vy[index];
      double vzL = vz[index];
      double hL = h[index];
      double PL = P[index];

      // We don't need chi_rad_nux since there is no heavy-lepton neutrinos in beta-process.
      double chi_radL = chi_rad[index];
      double chi_rad_nueL = chi_rad_nue[index];
      //F_rad_\alpha = g_{0\alpha} F_rad0 + g_{z\alpha} F_radz + g_{y\alpha} F_rady + g_{z\alpha} F_radz                                                                                                  

      double beta_xL = betaxL*gxxL + betayL*gxyL +betazL*gxzL;
      double beta_yL = betaxL*gxyL + betayL*gyyL +betazL*gyzL;
      double beta_zL = betaxL*gxzL + betayL*gyzL +betazL*gzzL;

      double u_xL = (gxxL*(betaxL+vxL) +
	           gxyL*(betayL+vyL) +
                   gxzL*(betazL+vzL))*u0L;
      double u_yL = (gxyL*(betaxL+vxL) +
                   gyyL*(betayL+vyL) +
                   gyzL*(betazL+vzL))*u0L;
      double u_zL = (gxzL*(betaxL+vxL) +
                   gyzL*(betayL+vyL) +
                   gzzL*(betazL+vzL))*u0L;
      double uxL = u0L*vxL;
      double uyL = u0L*vyL;
      double uzL = u0L*vzL;

      double u_0L = -(1.0 + uxL*u_xL + uyL*u_yL + uzL*u_zL)/u0L;
  
      double T_fluidL = T_fluid[index];
      //double T_fluidL = PL*al*Psi6*u0L/(rho_starL);

      // Compute b^2.  Note that gij is PHYSICAL gij here, so NO Psi4 factor out front!
      double b2 = -SQR(al*sbtL) + ( gxxL*SQR(sbxL+betaxL*sbtL) + 
				    2.0*gxyL*(sbxL+betaxL*sbtL)*(sbyL+betayL*sbtL) + 
				    2.0*gxzL*(sbxL+betaxL*sbtL)*(sbzL+betazL*sbtL) + 
				    gyyL*SQR(sbyL+betayL*sbtL) + 
				    2.0*gyzL*(sbyL+betayL*sbtL)*(sbzL+betazL*sbtL) + 
				    gzzL*SQR(sbzL+betazL*sbtL) );
      
      // take derivatives of the four-metric:
      double g4txm = D1_o3_cpp(g4tx_f,i,j,k);
      double g4xym = D1_o3_cpp(gxy_f,i,j,k);
      double g4xzm = D1_o3_cpp(gxz_f,i,j,k);
      double g4tym = D1_o3_cpp(g4ty_f,i,j,k);
      double g4yzm = D1_o3_cpp(gyz_f,i,j,k);
      
      double g4ttm = D1_o3_cpp(g4tt_f,i,j,k);
      double g4tzm = D1_o3_cpp(g4tz_f,i,j,k);
      double g4xxm = D1_o3_cpp(gxx_f,i,j,k);
      double g4yym = D1_o3_cpp(gyy_f,i,j,k);
      double g4zzm = D1_o3_cpp(gzz_f,i,j,k);


      double mhd_st_x_rhsL = mhd_st_x_rhs[index];


      // add \sqrt{-g}/2 * T^{00}\partial_x g_{00}
      mhd_st_x_rhsL += 0.5*( rho_starL*hL*u0L + al*Psi6*b2*SQR(u0L) 
			     - Psi6/al*(PL+0.5*b2) -al*Psi6*SQR(sbtL))*g4ttm;

      //if(i==1 && j==1 && k==1) printf("gggg00 %.15e\t%.15e\t%.15e\t%.15e\t%.15e\t%.15e\n",gxx_f[index],gxy_f[index],gxz_f[index],gyz_f[index],gyz_f[index],gzz_f[index]);
      
      // add \sqrt{-g} * T^{0i}\partial_x g_{0i} 
      mhd_st_x_rhsL += ( rho_starL*hL*u0L*vxL + al*Psi6*b2*vxL*SQR(u0L) +  
			 Psi6*betaxL/al*(PL+0.5*b2) - al*Psi6*sbtL*sbxL)*g4txm + 
        (rho_starL*hL*u0L*vyL + al*Psi6*b2*vyL*SQR(u0L) + 
         Psi6*betayL/al*(PL+0.5*b2) - al*Psi6*sbtL*sbyL)*g4tym + 
        (rho_starL*hL*u0L*vzL + al*Psi6*b2*vzL*SQR(u0L) + 
	 Psi6*betazL/al*(PL+0.5*b2) - al*Psi6*sbtL*sbzL)*g4tzm;

      //if(i==1 && j==1 && k==1) printf("gggg1 %.15e\n",mhd_st_x_rhsL);
      
      // add \sqrt{-g}/2 * T^{ij}\partial_x g_{ij}
      mhd_st_x_rhsL += 0.5*(rho_starL*hL*u0L*vxL*vxL + al*Psi6*b2*SQR(vxL*u0L) + 
			    (PL+0.5*b2)*(al*Psi6*gupxxL*Psim4 - SQR(betaxL)*Psi6/al) 
			    - al*Psi6*SQR(sbxL))*g4xxm;
      mhd_st_x_rhsL += (rho_starL*hL*u0L*vxL*vyL + al*Psi6*b2*vxL*vyL*SQR(u0L) + 
			(PL+0.5*b2)*(al*Psi6*gupxyL*Psim4 - betaxL*betayL*Psi6/al) 
			- al*Psi6*sbxL*sbyL)*g4xym;
      mhd_st_x_rhsL += (rho_starL*hL*u0L*vxL*vzL + al*Psi6*b2*vxL*vzL*SQR(u0L) + 
			(PL+0.5*b2)*(al*Psi6*gupxzL*Psim4 - betaxL*betazL*Psi6/al) 
			- al*Psi6*sbxL*sbzL)*g4xzm;
      mhd_st_x_rhsL += 0.5*(rho_starL*hL*u0L*vyL*vyL + al*Psi6*b2*SQR(vyL*u0L) + 
			    (PL+0.5*b2)*(al*Psi6*gupyyL*Psim4 - SQR(betayL)*Psi6/al) 
			    - al*Psi6*SQR(sbyL))*g4yym;
      mhd_st_x_rhsL += (rho_starL*hL*u0L*vyL*vzL + al*Psi6*b2*vyL*vzL*SQR(u0L) + 
			(PL+0.5*b2)*(al*Psi6*gupyzL*Psim4 - betayL*betazL*Psi6/al)  
			- al*Psi6*sbyL*sbzL)*g4yzm;
      mhd_st_x_rhsL += 0.5*(rho_starL*hL*u0L*vzL*vzL + al*Psi6*b2*SQR(vzL*u0L) + 
			    (PL+0.5*b2)*(al*Psi6*gupzzL*Psim4 - SQR(betazL)*Psi6/al)  
			    - al*Psi6*SQR(sbzL))*g4zzm;

      // Add radiation part of RHS (\sqrt{-g}*G_i = \sqrt{-g}*rho_star/(u^0*\sqrt{-g})*(opacity_abs*(E_rad - 4\pi*B)*u_x + (opacity_abs+opacity_sct)*F_rad_x)    ;

      if(rad_fourforce_enable==1){
	double kappa_a, kappa_s, eta_gf, kappa_a_nue, kappa_s_nue, eta_gf_nue, kappa_a_nux, kappa_s_nux, eta_gf_nux;
	double Y_eL = Y_e[index];                                                                                                                        
	//	double Y_eL = 0.0;
	double optd_L = optd[index];
	double eta_nueL = eta_nue[index];
	double fourforce_diag;

        if (compute_microphysics==1){
	  compute_opacity_emissivity(rho_starL, al, Psi6, u0L, kappa_a, kappa_s, eta_gf, kappa_a_nue, kappa_s_nue, eta_gf_nue, kappa_a_nux, kappa_s_nux, eta_gf_nux, T_fluidL, chi_radL, chi_rad_nueL, Y_eL, optd_L, eta_nueL, T_fluid_cgs_atm, microphysics_scheme, rad_fix);
	}
        else{ // Using constant opacities                                           
          kappa_a = rad_opacity_abs;
          kappa_s = rad_opacity_sct;
        }


	double E_radL = E_rad[index];
	double F_radxL = F_radx[index];
	double F_radyL = F_rady[index];
	double F_radzL = F_radz[index];
	double F_rad0L = - (F_radxL*u_xL + F_radyL*u_yL + F_radzL*u_zL)/u_0L;
	double F_rad_xL = gxxL * F_radxL + gxyL * F_radyL + gxzL * F_radzL + beta_xL*F_rad0L;

        if(enable_OS_collapse==1){
          mhd_st_x_rhsL += rho_starL/u0L*(kappa_a+kappa_s)*F_rad_xL;
        }
        else {
	  // E_rad_cgs/c has the unit of energy density (~L^-2), which is the same as nkbT
	  if (microphysics_scheme == 0){	    	   	    
	    mhd_st_x_rhsL += al * Psi6 * ((kappa_a*E_radL-eta_gf)*u_xL + (kappa_a+kappa_s)*F_rad_xL);
	    fourforce_diag = al * Psi6 * ((kappa_a*E_radL-eta_gf)*u_xL + (kappa_a+kappa_s)*F_rad_xL);
	  }
	  else{
	    // For full microphysics scheme, we need nue and nux set of radiation variables
	    double E_rad_nueL = E_rad_nue[index];
	    double F_radx_nueL = F_radx_nue[index];
	    double F_rady_nueL = F_rady_nue[index];
	    double F_radz_nueL = F_radz_nue[index];
	    double F_rad0_nueL = - (F_radx_nueL*u_xL + F_rady_nueL*u_yL + F_radz_nueL*u_zL)/u_0L;
	    double F_rad_x_nueL = gxxL * F_radx_nueL + gxyL * F_rady_nueL + gxzL * F_radz_nueL + beta_xL*F_rad0_nueL;
            double E_rad_nuxL = E_rad_nux[index];
            double F_radx_nuxL = F_radx_nux[index];
            double F_rady_nuxL = F_rady_nux[index];
            double F_radz_nuxL = F_radz_nux[index];
            double F_rad0_nuxL = - (F_radx_nuxL*u_xL + F_rady_nuxL*u_yL + F_radz_nuxL*u_zL)/u_0L;
            double F_rad_x_nuxL = gxxL * F_radx_nuxL + gxyL * F_rady_nuxL + gxzL * F_radz_nuxL + beta_xL*F_rad0_nuxL;
	    mhd_st_x_rhsL += al * Psi6 * ((kappa_a*E_radL-eta_gf)*u_xL + (kappa_a+kappa_s)*F_rad_xL);
	    mhd_st_x_rhsL += al * Psi6 * ((kappa_a_nue*E_rad_nueL-eta_gf_nue)*u_xL + (kappa_a_nue+kappa_s_nue)*F_rad_x_nueL);
	    mhd_st_x_rhsL += al * Psi6 * ((kappa_a_nux*E_rad_nuxL-eta_gf_nux)*u_xL + (kappa_a_nux+kappa_s_nux)*F_rad_x_nuxL);

	    if(mhd_st_x_rhsL>1.0)
	      {
                printf("inside source_z_tau.C, mhd_st_x_rhsL is NAN!!! \n");
                printf("kappa_a, kappa_a_nue, kappa_a_nux, eta_gf, eta_gf_nue, eta_gf_nux = %e,%e,%e,%e,%e,%e \n", kappa_a, kappa_a_nue, kappa_a_nux, eta_gf, eta_gf_nue, eta_gf_nux);
		printf("kappa_s, kappa_s_nue, kappa_s_nux = %e,%e,%e\n", kappa_s, kappa_s_nue, kappa_s_nux);
		printf("E_radL, E_rad_nueL, E_rad_nuxL, F_rad_xL, F_rad_x_nueL, F_rad_x_nuxL = %e,%e,%e,%e,%e,%e \n", E_radL, E_rad_nueL, E_rad_nuxL, F_rad_xL, F_rad_x_nueL, F_rad_x_nuxL);
              }
	  }
        }
	

      if (isnan(fourforce_diag))
	{
	  printf("Inside source_z_tau, al * Psi6 * ((kappa_a*E_radL-eta_gf)*u_xL + (kappa_a+kappa_s)*F_rad_xL) is NAN, kappa_a, E_radL, eta_gf, kappa_s, F_rad_xL, T_fluidL = %e,%e,%e,%e,%e,%e \n", kappa_a, E_radL, eta_gf, kappa_s, F_rad_xL, T_fluidL);
	}

      }


      mhd_st_x_rhs[index] = mhd_st_x_rhsL;
      
      // add -\sqrt{-g} * (T^{00} \betaxL + T^{0x})\partial_x alpha to tau_rhs
      // Note that in general, we would need to add time derivative terms of the metric to tau_rhs with the HARM energy variable.
      //  However, the only time we set enable_HARM_energyvariable=1 is with a Cowling (usually disk) run.
      if(enable_HARM_energyvariable==0) {
	double alphax_f = D1_o3_cpp(alpha_f,i,j,k);
  
	tau_rhs[index] += - ( (rho_starL*hL*u0L + al*Psi6*b2*SQR(u0L)
			       - al*Psi6*SQR(sbtL))*betaxL 
			      + rho_starL*hL*u0L*vxL + al*Psi6*b2*vxL*SQR(u0L) - al*Psi6*sbtL*sbxL )*alphax_f;
      }

    }
  } else if(flux_direction==2 && Symmetry != AXISYM) {
    //-----------------------------------------------------------------------------
    // Compute the addition to mhd_st_y_rhs and tau_rhs
    //-----------------------------------------------------------------------------
#pragma omp parallel for
    for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
      int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
      double al = 1.0 + alpha[index];
      double al_f = 1.0 + alpha_f[index];
      double betaxL = betax[index];
      double betayL = betay[index];
      double betazL = betaz[index];
      double Psi = exp(phi[index]);
      double Psi4 = Psi*Psi*Psi*Psi;
      double Psim4 = 1.0/Psi4;
      double Psi6 = Psi4*Psi*Psi;
      double gxxL = gxx[index];
      double gxyL = gxy[index];
      double gxzL = gxz[index];
      double gyyL = gyy[index];
      double gyzL = gyz[index];
      double gzzL = gzz[index];

      double gupxxL = gupxx[index];
      double gupxyL = gupxy[index];
      double gupxzL = gupxz[index];
      double gupyyL = gupyy[index];
      double gupyzL = gupyz[index];
      double gupzzL = gupzz[index];
      // Divide b^{\mu} by alpha sqrt(4 pi) 
      double sbtL = sbt[index]/al/sfpi;
      double sbxL = sbx[index]/al/sfpi;
      double sbyL = sby[index]/al/sfpi;
      double sbzL = sbz[index]/al/sfpi;

      double rho_starL = rho_star[index];
      double u0L = u0[index];
      double vxL = vx[index];
      double vyL = vy[index];
      double vzL = vz[index];
      double hL = h[index];
      double PL = P[index];

      double chi_radL = chi_rad[index];
      double chi_rad_nueL = chi_rad_nue[index];

      double beta_xL = betaxL*gxxL + betayL*gxyL +betazL*gxzL;
      double beta_yL = betaxL*gxyL + betayL*gyyL +betazL*gyzL;
      double beta_zL = betaxL*gxzL + betayL*gyzL +betazL*gzzL;

      double u_xL = (gxxL*(betaxL+vxL) +
                   gxyL*(betayL+vyL) +
                   gxzL*(betazL+vzL))*u0L;
      double u_yL = (gxyL*(betaxL+vxL) +
                   gyyL*(betayL+vyL) +
                   gyzL*(betazL+vzL))*u0L;
      double u_zL = (gxzL*(betaxL+vxL) +
                   gyzL*(betayL+vyL) +
                   gzzL*(betazL+vzL))*u0L;
      double uxL = u0L*vxL;
      double uyL = u0L*vyL;
      double uzL = u0L*vzL;

      double u_0L = -(1.0 + uxL*u_xL + uyL*u_yL + uzL*u_zL)/u0L;
      double T_fluidL = T_fluid[index];
      //double T_fluidL = PL*al*Psi6*u0L/(rho_starL);

      // Compute b^2.  Note that gij is PHYSICAL gij here, so NO Psi4 factor out front!
      double b2 = -SQR(al*sbtL) + ( gxxL*SQR(sbxL+betaxL*sbtL) + 
				    2.0*gxyL*(sbxL+betaxL*sbtL)*(sbyL+betayL*sbtL) + 
				    2.0*gxzL*(sbxL+betaxL*sbtL)*(sbzL+betazL*sbtL) + 
				    gyyL*SQR(sbyL+betayL*sbtL) + 
				    2.0*gyzL*(sbyL+betayL*sbtL)*(sbzL+betazL*sbtL) + 
				    gzzL*SQR(sbzL+betazL*sbtL) );
      //-----------------------------------------------------------------------------
      // Compute the addition to mhd_st_x_rhs and tau_rhs
      //-----------------------------------------------------------------------------
      
      // take derivatives of the four-metric:
      double g4txm = D2_o3_cpp(g4tx_f,i,j,k);
      double g4xym = D2_o3_cpp(gxy_f,i,j,k);
      double g4xzm = D2_o3_cpp(gxz_f,i,j,k);
      double g4tym = D2_o3_cpp(g4ty_f,i,j,k);
      double g4yzm = D2_o3_cpp(gyz_f,i,j,k);
      
      double g4ttm = D2_o3_cpp(g4tt_f,i,j,k);
      double g4tzm = D2_o3_cpp(g4tz_f,i,j,k);
      double g4xxm = D2_o3_cpp(gxx_f,i,j,k);
      double g4yym = D2_o3_cpp(gyy_f,i,j,k);
      double g4zzm = D2_o3_cpp(gzz_f,i,j,k);

     
      double mhd_st_y_rhsL = mhd_st_y_rhs[index];



      // add \sqrt{-g}/2 * T^{00}\partial_y g_{00}
      mhd_st_y_rhsL += 0.5*( rho_starL*hL*u0L + al*Psi6*b2*SQR(u0L) 
			     - Psi6/al*(PL+0.5*b2) -al*Psi6*SQR(sbtL))*g4ttm;

      // add \sqrt{-g} * T^{0i}\partial_y g_{0i} 
      mhd_st_y_rhsL += ( rho_starL*hL*u0L*vxL + al*Psi6*b2*vxL*SQR(u0L) + 
			 Psi6*betaxL/al*(PL+0.5*b2) - al*Psi6*sbtL*sbxL)*g4txm + 
	(rho_starL*hL*u0L*vyL + al*Psi6*b2*vyL*SQR(u0L) + 
	 Psi6*betayL/al*(PL+0.5*b2) - al*Psi6*sbtL*sbyL)*g4tym + 
	(rho_starL*hL*u0L*vzL + al*Psi6*b2*vzL*SQR(u0L) + 
	 Psi6*betazL/al*(PL+0.5*b2) - al*Psi6*sbtL*sbzL)*g4tzm;
      
      // add \sqrt{-g}/2 * T^{ij}\partial_y g_{ij}
      mhd_st_y_rhsL += 0.5*(rho_starL*hL*u0L*vxL*vxL + al*Psi6*b2*SQR(vxL*u0L) + 
			    (PL+0.5*b2)*(al*Psi6*gupxxL*Psim4 - SQR(betaxL)*Psi6/al) 
			    - al*Psi6*SQR(sbxL))*g4xxm;
      mhd_st_y_rhsL += (rho_starL*hL*u0L*vxL*vyL + al*Psi6*b2*vxL*vyL*SQR(u0L) + 
			(PL+0.5*b2)*(al*Psi6*gupxyL*Psim4 - betaxL*betayL*Psi6/al) 
			- al*Psi6*sbxL*sbyL)*g4xym;
      mhd_st_y_rhsL += (rho_starL*hL*u0L*vxL*vzL + al*Psi6*b2*vxL*vzL*SQR(u0L) + 
			(PL+0.5*b2)*(al*Psi6*gupxzL*Psim4 - betaxL*betazL*Psi6/al) 
			- al*Psi6*sbxL*sbzL)*g4xzm;
      mhd_st_y_rhsL += 0.5*(rho_starL*hL*u0L*vyL*vyL + al*Psi6*b2*SQR(vyL*u0L) + 
			    (PL+0.5*b2)*(al*Psi6*gupyyL*Psim4 - SQR(betayL)*Psi6/al) 
			    - al*Psi6*SQR(sbyL))*g4yym;
      mhd_st_y_rhsL += (rho_starL*hL*u0L*vyL*vzL + al*Psi6*b2*vyL*vzL*SQR(u0L) + 
			(PL+0.5*b2)*(al*Psi6*gupyzL*Psim4 - betayL*betazL*Psi6/al) 
			- al*Psi6*sbyL*sbzL)*g4yzm;
      mhd_st_y_rhsL += 0.5*(rho_starL*hL*u0L*vzL*vzL + al*Psi6*b2*SQR(vzL*u0L) + 
			    (PL+0.5*b2)*(al*Psi6*gupzzL*Psim4 - SQR(betazL)*Psi6/al) 
			    - al*Psi6*SQR(sbzL))*g4zzm;


      // Add radiation part of RHS;        
      if(rad_fourforce_enable==1){
	double kappa_a, kappa_s, eta_gf, kappa_a_nue, kappa_s_nue, eta_gf_nue, kappa_a_nux, kappa_s_nux, eta_gf_nux;
	double Y_eL = Y_e[index];                                                                                                                    
	//double Y_eL = 0.0;
	double optd_L = optd[index];
	double eta_nueL = eta_nue[index];
        if (compute_microphysics==1){
          compute_opacity_emissivity(rho_starL, al, Psi6, u0L, kappa_a, kappa_s, eta_gf, kappa_a_nue, kappa_s_nue, eta_gf_nue, kappa_a_nux, kappa_s_nux, eta_gf_nux, T_fluidL, chi_radL, chi_rad_nueL, Y_eL, optd_L, eta_nueL, T_fluid_cgs_atm, microphysics_scheme, rad_fix);
        }
        else{ // Using constant opacities                                                                                            
          kappa_a = rad_opacity_abs;
          kappa_s = rad_opacity_sct;
        }


        double E_radL = E_rad[index];
        double F_radxL = F_radx[index];
        double F_radyL = F_rady[index];
        double F_radzL = F_radz[index];
        double F_rad0L = - (F_radxL*u_xL + F_radyL*u_yL + F_radzL*u_zL)/u_0L;
	double F_rad_yL = gxyL * F_radxL + gyyL * F_radyL + gyzL * F_radzL + beta_yL*F_rad0L;

        if(enable_OS_collapse==1){
          mhd_st_y_rhsL += rho_starL/u0L*(kappa_a+kappa_s)*F_rad_yL;
        }
        else {
	  //mhd_st_y_rhsL += rho_starL/u0L*(kappa_a+kappa_s)*F_rad_yL; 
	  if (microphysics_scheme == 0){
	    mhd_st_y_rhsL += al * Psi6 * ((kappa_a*E_radL-eta_gf)*u_yL + (kappa_a+kappa_s)*F_rad_yL);
	  }
	  else{
	    // For full microphysics scheme, we need nue and nux set of radiation variables                                                                                                                    
            double E_rad_nueL = E_rad_nue[index];
            double F_radx_nueL = F_radx_nue[index];
            double F_rady_nueL = F_rady_nue[index];
            double F_radz_nueL = F_radz_nue[index];
            double F_rad0_nueL = - (F_radx_nueL*u_xL + F_rady_nueL*u_yL + F_radz_nueL*u_zL)/u_0L;
	    double F_rad_y_nueL = gxyL * F_radx_nueL + gyyL * F_rady_nueL + gyzL * F_radz_nueL + beta_yL*F_rad0_nueL;
            double E_rad_nuxL = E_rad_nux[index];
            double F_radx_nuxL = F_radx_nux[index];
            double F_rady_nuxL = F_rady_nux[index];
            double F_radz_nuxL = F_radz_nux[index];
            double F_rad0_nuxL = - (F_radx_nuxL*u_xL + F_rady_nuxL*u_yL + F_radz_nuxL*u_zL)/u_0L;
	    double F_rad_y_nuxL = gxyL * F_radx_nuxL + gyyL * F_rady_nuxL + gyzL * F_radz_nuxL + beta_yL*F_rad0_nuxL;
            mhd_st_y_rhsL += al * Psi6 * ((kappa_a*E_radL-eta_gf)*u_yL + (kappa_a+kappa_s)*F_rad_yL);
            mhd_st_y_rhsL += al * Psi6 * ((kappa_a_nue*E_rad_nueL-eta_gf_nue)*u_yL + (kappa_a_nue+kappa_s_nue)*F_rad_y_nueL);
            mhd_st_y_rhsL += al * Psi6 * ((kappa_a_nux*E_rad_nuxL-eta_gf_nux)*u_yL + (kappa_a_nux+kappa_s_nux)*F_rad_y_nuxL);


	    if(mhd_st_y_rhsL>1.0)
	      {
                printf("inside source_z_tau.C, mhd_st_y_rhsL is NAN!!! \n");
                printf("kappa_a, kappa_a_nue, kappa_a_nux, eta_gf, eta_gf_nue, eta_gf_nux = %e,%e,%e,%e,%e,%e \n", kappa_a, kappa_a_nue, kappa_a_nux, eta_gf, eta_gf_nue, eta_gf_nux);
		printf("kappa_s, kappa_s_nue, kappa_s_nux = %e,%e,%e\n", kappa_s, kappa_s_nue, kappa_s_nux);
		printf("E_radL, E_rad_nueL, E_rad_nuxL, F_rad_yL, F_rad_y_nueL, F_rad_y_nuxL = %e,%e,%e,%e,%e,%e \n", E_radL, E_rad_nueL, E_rad_nuxL, F_rad_yL, F_rad_y_nueL, F_rad_y_nuxL);
              }
	  }
        }
      }
  
      mhd_st_y_rhs[index] = mhd_st_y_rhsL;

      // add -\sqrt{-g} * (T^{00} \betayL + T^{0y})\partial_y alpha to tau_rhs
      // Note that in general, we would need to add time derivative terms of the metric to tau_rhs with the HARM energy variable.
      //  However, the only time we set enable_HARM_energyvariable=1 is with a Cowling (usually disk) run.
      if(enable_HARM_energyvariable==0) {
	double alphay_f = D2_o3_cpp(alpha_f,i,j,k);
        tau_rhs[index] += - ( (rho_starL*hL*u0L + al*Psi6*b2*SQR(u0L)
			       - al*Psi6*SQR(sbtL))*betayL +
                              rho_starL*hL*u0L*vyL + al*Psi6*b2*vyL*SQR(u0L) - al*Psi6*sbtL*sbyL )*alphay_f;
	/*
	tau_rhs[index] += - ( (rho_starL*hL*u0L + al*Psi6*b2*SQR(u0L) 
			       - Psi6/al*(PL+0.5*b2) - al*Psi6*SQR(sbtL))*betayL + 
			      rho_starL*hL*u0L*vyL + al*Psi6*b2*vyL*SQR(u0L) - al*Psi6*sbtL*sbyL )*alphay_f;
	*/      
      }
    }
  } else if(flux_direction==3) {
    //-----------------------------------------------------------------------------
    // Compute the addition to mhd_st_z_rhs and tau_rhs
    //-----------------------------------------------------------------------------
#pragma omp parallel for
    for(int k=kstart;k<kend;k++) for(int j=jstart;j<jend;j++) for(int i=istart;i<iend;i++) {
      int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
      double al = 1.0 + alpha[index];
      double al_f = 1.0 + alpha_f[index];
      double betaxL = betax[index];
      double betayL = betay[index];
      double betazL = betaz[index];
      double Psi = exp(phi[index]);
      double Psi4 = Psi*Psi*Psi*Psi;
      double Psim4 = 1.0/Psi4;
      double Psi6 = Psi4*Psi*Psi;
      double gxxL = gxx[index];
      double gxyL = gxy[index];
      double gxzL = gxz[index];
      double gyyL = gyy[index];
      double gyzL = gyz[index];
      double gzzL = gzz[index];

      double gupxxL = gupxx[index];
      double gupxyL = gupxy[index];
      double gupxzL = gupxz[index];
      double gupyyL = gupyy[index];
      double gupyzL = gupyz[index];
      double gupzzL = gupzz[index];
      // Divide b^{\mu} by alpha sqrt(4 pi) 
      double sbtL = sbt[index]/al/sfpi;
      double sbxL = sbx[index]/al/sfpi;
      double sbyL = sby[index]/al/sfpi;
      double sbzL = sbz[index]/al/sfpi;

      double rho_starL = rho_star[index];
      double u0L = u0[index];
      double vxL = vx[index];
      double vyL = vy[index];
      double vzL = vz[index];
      double hL = h[index];
      double PL = P[index];

      double chi_radL = chi_rad[index];
      double chi_rad_nueL = chi_rad_nue[index];

      double beta_xL = betaxL*gxxL + betayL*gxyL +betazL*gxzL;
      double beta_yL = betaxL*gxyL + betayL*gyyL +betazL*gyzL;
      double beta_zL = betaxL*gxzL + betayL*gyzL +betazL*gzzL;

      double u_xL = (gxxL*(betaxL+vxL) +
                   gxyL*(betayL+vyL) +
                   gxzL*(betazL+vzL))*u0L;
      double u_yL = (gxyL*(betaxL+vxL) +
                   gyyL*(betayL+vyL) +
                   gyzL*(betazL+vzL))*u0L;
      double u_zL = (gxzL*(betaxL+vxL) +
                   gyzL*(betayL+vyL) +
                   gzzL*(betazL+vzL))*u0L;
      double uxL = u0L*vxL;
      double uyL = u0L*vyL;
      double uzL = u0L*vzL;

      double u_0L = -(1.0 + uxL*u_xL + uyL*u_yL + uzL*u_zL)/u0L;
      double T_fluidL = T_fluid[index];
      //double T_fluidL = PL*al*Psi6*u0L/(rho_starL);
      // Compute b^2.  Note that gij is PHYSICAL gij here, so NO Psi4 factor out front!
      double b2 = -SQR(al*sbtL) + ( gxxL*SQR(sbxL+betaxL*sbtL) + 
				    2.0*gxyL*(sbxL+betaxL*sbtL)*(sbyL+betayL*sbtL) + 
				    2.0*gxzL*(sbxL+betaxL*sbtL)*(sbzL+betazL*sbtL) + 
				    gyyL*SQR(sbyL+betayL*sbtL) + 
				    2.0*gyzL*(sbyL+betayL*sbtL)*(sbzL+betazL*sbtL) + 
				    gzzL*SQR(sbzL+betazL*sbtL) );
      
      // take derivatives of the four-metric:
      double g4txm = D3_o3_cpp(g4tx_f,i,j,k);
      double g4xym = D3_o3_cpp(gxy_f,i,j,k);
      double g4xzm = D3_o3_cpp(gxz_f,i,j,k);
      double g4tym = D3_o3_cpp(g4ty_f,i,j,k);
      double g4yzm = D3_o3_cpp(gyz_f,i,j,k);
      
      double g4ttm = D3_o3_cpp(g4tt_f,i,j,k);
      double g4tzm = D3_o3_cpp(g4tz_f,i,j,k);
      double g4xxm = D3_o3_cpp(gxx_f,i,j,k);
      double g4yym = D3_o3_cpp(gyy_f,i,j,k);
      double g4zzm = D3_o3_cpp(gzz_f,i,j,k);

      double mhd_st_z_rhsL = mhd_st_z_rhs[index];
      
      // add \sqrt{-g}/2 * T^{00}\partial_z g_{00}
      mhd_st_z_rhsL += 0.5*( rho_starL*hL*u0L + al*Psi6*b2*SQR(u0L) 
			     - Psi6/al*(PL+0.5*b2) -al*Psi6*SQR(sbtL))*g4ttm;

      // add \sqrt{-g} * T^{0i}\partial_z g_{0i} 
      mhd_st_z_rhsL += ( rho_starL*hL*u0L*vxL + al*Psi6*b2*vxL*SQR(u0L) + 
			 Psi6*betaxL/al*(PL+0.5*b2) - al*Psi6*sbtL*sbxL)*g4txm + 
	(rho_starL*hL*u0L*vyL + al*Psi6*b2*vyL*SQR(u0L) + 
	 Psi6*betayL/al*(PL+0.5*b2) - al*Psi6*sbtL*sbyL)*g4tym + 
	(rho_starL*hL*u0L*vzL + al*Psi6*b2*vzL*SQR(u0L) + 
	 Psi6*betazL/al*(PL+0.5*b2) - al*Psi6*sbtL*sbzL)*g4tzm;

      // add \sqrt{-g}/2 * T^{ij}\partial_z g_{ij}
      mhd_st_z_rhsL += 0.5*(rho_starL*hL*u0L*vxL*vxL + al*Psi6*b2*SQR(vxL*u0L) + 
			    (PL+0.5*b2)*(al*Psi6*gupxxL*Psim4 - SQR(betaxL)*Psi6/al) 
			    - al*Psi6*SQR(sbxL))*g4xxm;
      mhd_st_z_rhsL += (rho_starL*hL*u0L*vxL*vyL + al*Psi6*b2*vxL*vyL*SQR(u0L) + 
			(PL+0.5*b2)*(al*Psi6*gupxyL*Psim4 - betaxL*betayL*Psi6/al) 
			- al*Psi6*sbxL*sbyL)*g4xym;
      mhd_st_z_rhsL += (rho_starL*hL*u0L*vxL*vzL + al*Psi6*b2*vxL*vzL*SQR(u0L) + 
			(PL+0.5*b2)*(al*Psi6*gupxzL*Psim4 - betaxL*betazL*Psi6/al) 
			- al*Psi6*sbxL*sbzL)*g4xzm;
      mhd_st_z_rhsL += 0.5*(rho_starL*hL*u0L*vyL*vyL + al*Psi6*b2*SQR(vyL*u0L) + 
			    (PL+0.5*b2)*(al*Psi6*gupyyL*Psim4 - SQR(betayL)*Psi6/al) 
			    - al*Psi6*SQR(sbyL))*g4yym;
      mhd_st_z_rhsL += (rho_starL*hL*u0L*vyL*vzL + al*Psi6*b2*vyL*vzL*SQR(u0L) + 
			(PL+0.5*b2)*(al*Psi6*gupyzL*Psim4 - betayL*betazL*Psi6/al) 
			- al*Psi6*sbyL*sbzL)*g4yzm;
      mhd_st_z_rhsL += 0.5*(rho_starL*hL*u0L*vzL*vzL + al*Psi6*b2*SQR(vzL*u0L) + 
			    (PL+0.5*b2)*(al*Psi6*gupzzL*Psim4 - SQR(betazL)*Psi6/al) 
			    - al*Psi6*SQR(sbzL))*g4zzm;

      // Add radiation part of RHS;      
      if(rad_fourforce_enable==1){
	double kappa_a, kappa_s, eta_gf, kappa_a_nue, kappa_s_nue, eta_gf_nue, kappa_a_nux, kappa_s_nux, eta_gf_nux;
	double Y_eL = Y_e[index];                                                                                                                 
	//double Y_eL = 0.0;
	double optd_L = optd[index];  
	double eta_nueL = eta_nue[index];
        if (compute_microphysics==1){
	  compute_opacity_emissivity(rho_starL, al, Psi6, u0L, kappa_a, kappa_s, eta_gf, kappa_a_nue, kappa_s_nue, eta_gf_nue, kappa_a_nux, kappa_s_nux, eta_gf_nux, T_fluidL, chi_radL, chi_rad_nueL, Y_eL, optd_L, eta_nueL, T_fluid_cgs_atm, microphysics_scheme, rad_fix);
        }
        else{ // Using constant opacities                                                                                             
          kappa_a = rad_opacity_abs;
          kappa_s = rad_opacity_sct;
        }

        double E_radL = E_rad[index];
        double F_radxL = F_radx[index];
        double F_radyL = F_rady[index];
        double F_radzL = F_radz[index];
        double F_rad0L = - (F_radxL*u_xL + F_radyL*u_yL + F_radzL*u_zL)/u_0L;
	double F_rad_zL = gxzL * F_radxL + gyzL * F_radyL + gzzL * F_radzL + beta_zL*F_rad0L;

        if(enable_OS_collapse==1){
          mhd_st_z_rhsL += rho_starL/u0L*(kappa_a+kappa_s)*F_rad_zL;
        }
        else {
	  if (microphysics_scheme == 0){
	    mhd_st_z_rhsL += al * Psi6 * ((kappa_a*E_radL-eta_gf)*u_zL + (kappa_a+kappa_s)*F_rad_zL);
	  }
	  else{
	    // For full microphysics scheme, we need nue and nux set of radiation variables                                                                                                                    
            double E_rad_nueL = E_rad_nue[index];
            double F_radx_nueL = F_radx_nue[index];
            double F_rady_nueL = F_rady_nue[index];
            double F_radz_nueL = F_radz_nue[index];
            double F_rad0_nueL = - (F_radx_nueL*u_xL + F_rady_nueL*u_yL + F_radz_nueL*u_zL)/u_0L;
	    double F_rad_z_nueL = gxzL * F_radx_nueL + gyzL * F_rady_nueL + gzzL * F_radz_nueL + beta_zL*F_rad0_nueL;
            double E_rad_nuxL = E_rad_nux[index];
            double F_radx_nuxL = F_radx_nux[index];
            double F_rady_nuxL = F_rady_nux[index];
            double F_radz_nuxL = F_radz_nux[index];
            double F_rad0_nuxL = - (F_radx_nuxL*u_xL + F_rady_nuxL*u_yL + F_radz_nuxL*u_zL)/u_0L;
	    double F_rad_z_nuxL = gxzL * F_radx_nuxL + gyzL * F_rady_nuxL + gzzL * F_radz_nuxL + beta_zL*F_rad0_nuxL;
            mhd_st_z_rhsL += al * Psi6 * ((kappa_a*E_radL-eta_gf)*u_zL + (kappa_a+kappa_s)*F_rad_zL);
            mhd_st_z_rhsL += al * Psi6 * ((kappa_a_nue*E_rad_nueL-eta_gf_nue)*u_zL + (kappa_a_nue+kappa_s_nue)*F_rad_z_nueL);
            mhd_st_z_rhsL += al * Psi6 * ((kappa_a_nux*E_rad_nuxL-eta_gf_nux)*u_zL + (kappa_a_nux+kappa_s_nux)*F_rad_z_nuxL);

	    if(mhd_st_z_rhsL>1.0)
	      {
		printf("inside source_z_tau.C, mhd_st_z_rhsL is NAN!!! \n");
		printf("kappa_a, kappa_a_nue, kappa_a_nux, eta_gf, eta_gf_nue, eta_gf_nux = %e,%e,%e,%e,%e,%e \n", kappa_a, kappa_a_nue, kappa_a_nux, eta_gf, eta_gf_nue, eta_gf_nux);
		printf("kappa_s, kappa_s_nue, kappa_s_nux = %e,%e,%e\n", kappa_s, kappa_s_nue, kappa_s_nux);
		printf("E_radL, E_rad_nueL, E_rad_nuxL, F_rad_zL, F_rad_z_nueL, F_rad_z_nuxL = %e,%e,%e,%e,%e,%e \n", E_radL, E_rad_nueL, E_rad_nuxL, F_rad_zL, F_rad_z_nueL, F_rad_z_nuxL);
	      }
	  }
        }
      }	

      mhd_st_z_rhs[index] = mhd_st_z_rhsL;

      // add -\sqrt{-g} * (T^{00} \betazL + T^{0z})\partial_z alpha to tau_rhs
      // Note that in general, we would need to add time derivative terms of the metric to tau_rhs with the HARM energy variable.
      //  However, the only time we set enable_HARM_energyvariable=1 is with a Cowling (usually disk) run.
      if(enable_HARM_energyvariable==0) {
	double alphaz_f = D3_o3_cpp(alpha_f,i,j,k);
	tau_rhs[index] += - ( (rho_starL*hL*u0L + al*Psi6*b2*SQR(u0L) 
			       - al*Psi6*SQR(sbtL))*betazL + 
			      rho_starL*hL*u0L*vzL + al*Psi6*b2*vzL*SQR(u0L) - al*Psi6*sbtL*sbzL )*alphaz_f;
            }
	}
  }


  //----------------------------------------------------------------------------
  // Convert 3-metric back to conformal from physical
  //----------------------------------------------------------------------------
     
#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    double Psim4 = exp(-4.0*phi[index]);
    double Psim4_f = exp(-4.0*phi_f[index]);

    gxx[index] *= Psim4;
    gxy[index] *= Psim4;
    gxz[index] *= Psim4;
    gyy[index] *= Psim4;
    gyz[index] *= Psim4;
    gzz[index] *= Psim4;

    gxx_f[index] *= Psim4_f;
    gxy_f[index] *= Psim4_f;
    gxz_f[index] *= Psim4_f;
    gyy_f[index] *= Psim4_f;
    gyz_f[index] *= Psim4_f;
    gzz_f[index] *= Psim4_f;
   }


  //  printf("END mhd_source_z_tau_cpp!!!!!!!!!!!!!!! \n");

}

extern "C" void CCTK_FCALL mhd_source_z_tau_cpp_
  (int *flux_direction, const cGH **cctkGH,int *cctk_lsh, int *nghostzones, int *Symmetry,
   int *enable_HARM_energyvariable,
   double *dX, double *dY,double *dZ, double *Z,
   double *mhd_st_x_rhs, double *mhd_st_y_rhs, double *mhd_st_z_rhs, 
   double *tau_rhs, double *rho_star, 
   double *P, double *h, double *u0, 
   double *vx, double *vy, double *vz,  
   double *sbt, double *sbx, double *sby, double *sbz, 
   double *alpha, double *betax, double *betay, double *betaz, 
   double *phi, 
   double *gxx, double *gxy, double *gxz, double *gyy, double *gyz, double *gzz, 
   double *gupxx, double *gupxy, double *gupxz, double *gupyy, double *gupyz, double *gupzz, 
   double *alpha_f, double *betax_f, double *betay_f, double *betaz_f, 
   double *phi_f, 				 
   double *gxx_f, double *gxy_f, double *gxz_f, double *gyy_f, double *gyz_f, double *gzz_f,
   double *g4tt,double *g4tx,double *g4ty,double *g4tz,
   double *g4tt_f,double *g4tx_f,double *g4ty_f,double *g4tz_f,
   double *E_rad, double *F_radx, double *F_rady, double *F_radz, double *Y_e, double *optd, double *eta_nue, 
   double *E_rad_nue, double *F_radx_nue, double *F_rady_nue, double *F_radz_nue,
   double *E_rad_nux, double *F_radx_nux, double *F_rady_nux, double *F_radz_nux,
   int &rad_fourforce_enable, int *rad_closure_scheme, double &rad_opacity_abs, double &rad_opacity_sct, double &rad_const, 
   double *chi_rad, double *chi_rad_nue,
   double *T_fluid,  int &enable_OS_collapse, int &compute_microphysics, int &microphysics_scheme, double &T_fluid_cgs_atm, int &rad_fix) 
{
  mhd_source_z_tau_cpp(*flux_direction, *cctkGH,cctk_lsh,nghostzones,*Symmetry,
		       *enable_HARM_energyvariable,
		       *dX, *dY,*dZ, Z,
		       mhd_st_x_rhs, mhd_st_y_rhs, mhd_st_z_rhs, 
		       tau_rhs, rho_star, 
		       P, h, u0, 
		       vx, vy, vz,  
		       sbt, sbx, sby, sbz, 
		       alpha, betax, betay, betaz, 
		       phi, 
		       gxx, gxy, gxz, gyy, gyz, gzz, 
		       gupxx, gupxy, gupxz, gupyy, gupyz, gupzz, 
		       alpha_f, betax_f, betay_f, betaz_f, 
		       phi_f, 				 
		       gxx_f, gxy_f, gxz_f, gyy_f, gyz_f, gzz_f,
		       g4tt,g4tx,g4ty,g4tz,
		       g4tt_f,g4tx_f,g4ty_f,g4tz_f,
		       E_rad,F_radx,F_rady,F_radz, Y_e, optd, eta_nue,
		       E_rad_nue,F_radx_nue,F_rady_nue,F_radz_nue,
		       E_rad_nux,F_radx_nux,F_rady_nux,F_radz_nux,
		       rad_fourforce_enable,rad_closure_scheme,rad_opacity_abs,rad_opacity_sct, rad_const, 
		       chi_rad, chi_rad_nue,
		       T_fluid, enable_OS_collapse, compute_microphysics, microphysics_scheme, T_fluid_cgs_atm, rad_fix);
}





extern "C" void CCTK_FCALL mhd_tau_scalar_rad_cpp_
  (const cGH **cctkGH,int *cctk_lsh,  
   double *tau_rhs, double *rhoYe_rhs, double *rho_star, double *P,
   double *u0, double *vx, double *vy, double *vz,
   double *alpha, double *betax, double *betay, double *betaz,
   double *phi, double *gxx, double *gxy, double *gxz, double *gyy, double *gyz, double *gzz,
   double *E_rad, double *F_radx, double *F_rady, double *F_radz, double *Y_e,double *optd, double *eta_nue, 
   double *E_rad_nue, double *F_radx_nue, double *F_rady_nue, double *F_radz_nue,
   double *E_rad_nux, double *F_radx_nux, double *F_rady_nux, double *F_radz_nux,
   int *rad_closure_scheme, double &rad_opacity_abs, double &rad_opacity_sct, double &rad_const, 
   double *chi_rad, double *chi_rad_nue,
   double *T_fluid,  int &enable_OS_collapse, int &compute_microphysics, int &microphysics_scheme, double &T_fluid_cgs_atm, int &rad_fix);

extern "C" void mhd_tau_scalar_rad_cpp  (const cGH *cctkGH,int *cctk_lsh,
					 double *tau_rhs, double *rhoYe_rhs, double *rho_star, double *P,
					 double *u0, double *vx, double *vy, double *vz,
					 double *alpha, double *betax, double *betay, double *betaz,
					 double *phi, double *gxx, double *gxy, double *gxz, double *gyy, double *gyz, double *gzz,
					 double *E_rad, double *F_radx, double *F_rady, double *F_radz, double *Y_e,double *optd, double *eta_nue, 
					 double *E_rad_nue, double *F_radx_nue, double *F_rady_nue, double *F_radz_nue,
					 double *E_rad_nux, double *F_radx_nux, double *F_rady_nux, double *F_radz_nux,
					 int *rad_closure_scheme, double &rad_opacity_abs, double &rad_opacity_sct, double &rad_const, 
					 double *chi_rad, double *chi_rad_nue,
					 double *T_fluid,  int &enable_OS_collapse, int &compute_microphysics, int &microphysics_scheme, double &T_fluid_cgs_atm, int &rad_fix){
//Finally, add the scalar radiation term on tau_rhs (alpha * (rho_star/u^0) * (opacity_abs*(E-4pi*B)*u^0 +  (opacity_abs + opacity_sct)*F^0)   );

//  printf("START mhd_tau_scalar_rad_cpp!!!!!!!!!!!!!!! \n");

#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
        int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
        double al = 1.0 + alpha[index];
        double Psi = exp(phi[index]);
        double Psi4 = Psi*Psi*Psi*Psi;
	double Psi6 = Psi4*Psi*Psi;

	double rho_starL = rho_star[index];
        double u0L = u0[index];
        double vxL = vx[index];
        double vyL = vy[index];
        double vzL = vz[index];
	double PL = P[index];

        double betaxL = betax[index];
        double betayL = betay[index];
        double betazL = betaz[index];

        double gxxL = gxx[index];
        double gxyL = gxy[index];
        double gxzL = gxz[index];
        double gyyL = gyy[index];
        double gyzL = gyz[index];
        double gzzL = gzz[index];

        double beta_xL = Psi4*(betaxL*gxxL + betayL*gxyL +betazL*gxzL);
        double beta_yL = Psi4*(betaxL*gxyL + betayL*gyyL +betazL*gyzL);
        double beta_zL = Psi4*(betaxL*gxzL + betayL*gyzL +betazL*gzzL);

	double chi_radL = chi_rad[index];
	double chi_rad_nueL = chi_rad_nue[index];

        double u_xL = Psi4*(gxxL*(betaxL+vxL) +
			    gxyL*(betayL+vyL) +
			    gxzL*(betazL+vzL))*u0L;
        double u_yL = Psi4*(gxyL*(betaxL+vxL) +
			    gyyL*(betayL+vyL) +
			    gyzL*(betazL+vzL))*u0L;
        double u_zL = Psi4*(gxzL*(betaxL+vxL) +
			    gyzL*(betayL+vyL) +
			    gzzL*(betazL+vzL))*u0L;
        double uxL = u0L*vxL;
        double uyL = u0L*vyL;
        double uzL = u0L*vzL;
        double u_0L = -(1.0 + uxL*u_xL + uyL*u_yL + uzL*u_zL)/u0L;
	double T_fluidL = T_fluid[index];

	double kappa_a, kappa_s, eta_gf, kappa_a_nue, kappa_s_nue, eta_gf_nue, kappa_a_nux, kappa_s_nux, eta_gf_nux;	
	double Y_eL = Y_e[index];                                                                                                        
	//double Y_eL = 0.0;
	double optd_L = optd[index];
	double eta_nueL = eta_nue[index];
        if (compute_microphysics==1){
	  compute_opacity_emissivity(rho_starL, al, Psi6, u0L, kappa_a, kappa_s, eta_gf, kappa_a_nue, kappa_s_nue, eta_gf_nue, kappa_a_nux, kappa_s_nux, eta_gf_nux, T_fluidL, chi_radL, chi_rad_nueL, Y_eL, optd_L, eta_nueL, T_fluid_cgs_atm, microphysics_scheme, rad_fix);
        }
        else{ // Using constant opacities                                     
          kappa_a = rad_opacity_abs;
          kappa_s = rad_opacity_sct;
        }

	double E_radL = E_rad[index];
        double F_radxL = F_radx[index];
        double F_radyL = F_rady[index];
        double F_radzL = F_radz[index];
        double F_rad0L = - (F_radxL*u_xL + F_radyL*u_yL + F_radzL*u_zL)/u_0L;
	double F_rad_0L = - F_rad0L/al/al + F_radxL*beta_xL + F_radyL*beta_yL + F_radzL*beta_zL;
	
        if(enable_OS_collapse==1){
          tau_rhs[index] += al*rho_starL/u0L*(kappa_a+kappa_s)*F_rad0L;
	}
        else {
	  //tau_rhs[index] += al*rho_starL/u0L*(kappa_a+kappa_s)*F_rad0L;
	  if (microphysics_scheme == 0){
	    tau_rhs[index] += al* al * Psi6 *( (kappa_a*E_radL - eta_gf)*u0L + (kappa_a+kappa_s)*F_rad0L);
	  }
	  else{
	    double E_rad_nueL = E_rad_nue[index];
	    double F_radx_nueL = F_radx_nue[index];
	    double F_rady_nueL = F_rady_nue[index];
	    double F_radz_nueL = F_radz_nue[index];
	    double F_rad0_nueL = - (F_radx_nueL*u_xL + F_rady_nueL*u_yL + F_radz_nueL*u_zL)/u_0L;
	    double F_rad_0_nueL = - F_rad0_nueL/al/al + F_radxL*beta_xL + F_rady_nueL*beta_yL + F_radz_nueL*beta_zL;
            double E_rad_nuxL = E_rad_nux[index];
            double F_radx_nuxL = F_radx_nux[index];
            double F_rady_nuxL = F_rady_nux[index];
            double F_radz_nuxL = F_radz_nux[index];
            double F_rad0_nuxL = - (F_radx_nuxL*u_xL + F_rady_nuxL*u_yL + F_radz_nuxL*u_zL)/u_0L;
            double F_rad_0_nuxL = - F_rad0_nuxL/al/al + F_radxL*beta_xL + F_rady_nuxL*beta_yL + F_radz_nuxL*beta_zL;
	    tau_rhs[index] += al* al * Psi6 *( (kappa_a*E_radL - eta_gf)*u0L + (kappa_a+kappa_s)*F_rad0L);
	    tau_rhs[index] += al* al * Psi6 *( (kappa_a_nue*E_rad_nueL - eta_gf_nue)*u0L + (kappa_a_nue+kappa_s_nue)*F_rad0_nueL);
	    tau_rhs[index] += al* al * Psi6 *( (kappa_a_nux*E_rad_nuxL - eta_gf_nux)*u0L + (kappa_a_nux+kappa_s_nux)*F_rad0_nuxL);

	    if(isnan(tau_rhs[index]))
	      {
                printf("inside source_z_tau.C, tau_rhs[index] is NAN!!! \n");
                printf("kappa_a, kappa_a_nue, kappa_a_nux, eta_gf, eta_gf_nue, eta_gf_nux = %e,%e,%e,%e,%e,%e \n", kappa_a, kappa_a_nue, kappa_a_nux, eta_gf, eta_gf_nue, eta_gf_nux);
		printf("E_radL, E_rad_nueL, E_rad_nuxL, F_rad0L, F_rad0_nueL, F_rad0_nuxL = %e,%e,%e,%e,%e,%e \n", E_radL, E_rad_nueL, E_rad_nuxL, F_rad0L, F_rad0_nueL, F_rad0_nuxL);
              }


	  }
        }
	
        //We have set chemical potential to zero, so the eta in the Fermi intergral is zero.                                                                
	double T_fluid_cgs = T_fluidL/(kb_cgs * G_cgs) * pow(c_cgs,4.0)*pow(k2c, 1.0);
	double rhoYe_rhsL,rhoYe_rhs_nueL,rhoYe_rhs_nuxL;

	double E_rad_nueL,E_rad_nuxL,ave_eng_cgs,ave_eng;
	double ave_eng_nue_cgs,ave_eng_nue,ave_eng_nux_cgs, ave_eng_nux;

	if (T_fluid_cgs > 1.0e8){
	  double kBT_cgs = T_fluid_cgs * kb_cgs; 
	  double eta_nuae = 0.0;
	  double F5o4 = Fermi5(eta_nuae)/Fermi4(eta_nuae);
	  double rho_b = rho_starL/u0L/al/Psi6;
	  double rho_b_cgs = fasterpow_prim(c_cgs,2.0) / G_cgs * rho_b * fasterpow_prim(k2c, -2.0);
	  double n_nucleon = rho_b_cgs/m_n_cgs;
	  double ave_eng_cgs = F5o4*kBT_cgs*n_nucleon;
	  double ave_eng = ave_eng_cgs * fasterpow_prim(c_cgs,-4.0) * G_cgs * fasterpow_prim(c2k, -2.0);
	  rhoYe_rhsL = al*Psi6*(kappa_a*E_radL - eta_gf)* rho_b/ave_eng;
	
	  if (microphysics_scheme == 1){
	    E_rad_nueL = E_rad_nue[index];
	    E_rad_nuxL = E_rad_nux[index];
	    
	    ave_eng_nue_cgs = F5o4*kBT_cgs*n_nucleon;
	    ave_eng_nue = ave_eng_cgs * fasterpow_prim(c_cgs,-4.0) * G_cgs * fasterpow_prim(c2k, -2.0);	   
	    rhoYe_rhs_nueL = al*Psi6*(kappa_a_nue*E_rad_nueL - eta_gf_nue)* rho_b/ave_eng_nue;
	  
	    ave_eng_nux_cgs = F5o4*kBT_cgs*n_nucleon;
            ave_eng_nux = ave_eng_cgs * fasterpow_prim(c_cgs,-4.0) * G_cgs * fasterpow_prim(c2k, -2.0);
            rhoYe_rhs_nuxL = al*Psi6*(kappa_a_nux*E_rad_nuxL - eta_gf_nux)* rho_b/ave_eng_nux;	    
	  }
	  else{
	    rhoYe_rhs_nueL = 0.0;
	    rhoYe_rhs_nuxL = 0.0;
	  }
	}else{
	  rhoYe_rhsL = 0.0;
	  rhoYe_rhs_nueL = 0.0;
	  rhoYe_rhs_nuxL = 0.0;
	}
	rhoYe_rhs[index] +=  rhoYe_rhsL + rhoYe_rhs_nueL + rhoYe_rhs_nuxL;

	if(isnan(rhoYe_rhs[index]))
	  {
	    printf("inside source_z_tau.C, rhoYe_rhs[idnex] is NAN!!! \n");
	    printf("rhoYe_rhsL, rhoYe_rhs_nueL, rhoYe_rhs_nuxL = %e,%e,%e \n", rhoYe_rhsL, rhoYe_rhs_nueL, rhoYe_rhs_nuxL);
	    printf("ave_eng, ave_eng_nue, ave_eng_nux = %e,%e,%e \n", ave_eng, ave_eng_nue, ave_eng_nux);
	  }

      }      
}

extern "C" void CCTK_FCALL mhd_tau_scalar_rad_cpp_
  (const cGH **cctkGH,int *cctk_lsh,
   double *tau_rhs, double *rhoYe_rhs, double *rho_star, double *P,
   double *u0, double *vx, double *vy, double *vz,
   double *alpha, double *betax, double *betay, double *betaz,
   double *phi, double *gxx, double *gxy, double *gxz, double *gyy, double *gyz, double *gzz,
   double *E_rad, double *F_radx, double *F_rady, double *F_radz, double *Y_e,double *optd, double *eta_nue, 
   double *E_rad_nue, double *F_radx_nue, double *F_rady_nue, double *F_radz_nue,
   double *E_rad_nux, double *F_radx_nux, double *F_rady_nux, double *F_radz_nux,
   int *rad_closure_scheme, double &rad_opacity_abs, double &rad_opacity_sct, double &rad_const, 
   double *chi_rad, double *chi_rad_nue,
   double *T_fluid,  int &enable_OS_collapse, int &compute_microphysics, int &microphysics_scheme, double &T_fluid_cgs_atm, int &rad_fix)
{mhd_tau_scalar_rad_cpp(*cctkGH, cctk_lsh,
			tau_rhs, rhoYe_rhs, rho_star, P,
			u0, vx, vy, vz,
			alpha, betax, betay, betaz,
			phi, gxx, gxy, gxz, gyy, gyz, gzz,
			E_rad, F_radx, F_rady, F_radz, Y_e, optd, eta_nue, 
			E_rad_nue,F_radx_nue,F_rady_nue,F_radz_nue,
			E_rad_nux,F_radx_nux,F_rady_nux,F_radz_nux,
			rad_closure_scheme, rad_opacity_abs, rad_opacity_sct, rad_const, 
			chi_rad, chi_rad_nue,
			T_fluid,  enable_OS_collapse, compute_microphysics, microphysics_scheme, T_fluid_cgs_atm, rad_fix);
}
