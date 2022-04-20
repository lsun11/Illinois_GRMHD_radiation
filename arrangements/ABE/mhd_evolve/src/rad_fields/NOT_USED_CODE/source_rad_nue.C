//-----------------------------------------------------------------------------
// Add third-order accurate curvature terms to S_i_rhs
//-----------------------------------------------------------------------------

#include "math.h"
#include "cctk.h"
#include "stdio.h"
#include "primitives_solver_header.h"

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

extern "C" void CCTK_FCALL CCTK_FNAME(rad_source_nue_cpp)
  (int *flux_direction, const cGH **cctkGH,int *cctk_lsh, int *nghostzones, int *Symmetry,
   int *enable_HARM_energyvariable,
   double *dX, double *dY,double *dZ,
   double *S_rad_x_nue_rhs, double *S_rad_y_nue_rhs, double *S_rad_z_nue_rhs, 
   double *tau_rad_nue_rhs, double *rho_star, 
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
   double *E_rad_nue, double *F_radx_nue, double *F_rady_nue, double *F_radz_nue, double *Y_e, double *optd, 
   int &rad_closure_scheme, double &rad_opacity_abs, double &rad_opacity_sct, double &rad_const, 
   double *ka_gf, double *ks_gf, double *emission_gf, double *chi_rad,
   double *T_fluid, double &Erad_atm_cut, int &enable_OS_collapse, int &compute_microphysics, int &rad_fourforce_enable, double &T_fluid_cgs_atm);

extern "C" void rad_source_nue_cpp(int flux_direction, const cGH *cctkGH,int *cctk_lsh, int *nghostzones, int Symmetry,
				     int enable_HARM_energyvariable,
				     double dX, double dY,double dZ,
				     double *S_rad_x_nue_rhs, double *S_rad_y_nue_rhs, double *S_rad_z_nue_rhs, 
				     double *tau_rad_nue_rhs, double *rho_star, 
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
			             double *E_rad_nue, double *F_radx_nue, double *F_rady_nue, double *F_radz_nue, double *Y_e,double *optd,
			             int &rad_closure_scheme, double &rad_opacity_abs, double &rad_opacity_sct, double &rad_const,
			             double *ka_gf, double *ks_gf, double *emission_gf, double *chi_rad,
			       double *T_fluid, double &Erad_atm_cut, int &enable_OS_collapse, int &compute_microphysics, int &rad_fourforce_enable, double &T_fluid_cgs_atm) {

  //  printf("START rad_source_nue_cpp!!!!!!!!!!!!!!!!!\n");

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

  if(flux_direction==1) {
    //-----------------------------------------------------------------------------
    // Compute the addition to S_rad_x_nue_rhs and tau_rad_nue_rhs
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

      double E_rad_nueL = E_rad_nue[index];
      double F_radx_nueL = F_radx_nue[index];
      double F_rady_nueL = F_rady_nue[index];
      double F_radz_nueL = F_radz_nue[index];
      double chi_radL = chi_rad[index];

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
      double F_rad0L = - (F_radx_nueL*u_xL + F_rady_nueL*u_yL + F_radz_nueL*u_zL)/u_0L;

      
      double beta_xL = betaxL*gxxL + betayL*gxyL +betazL*gxzL;
      double beta_yL = betaxL*gxyL + betayL*gyyL +betazL*gyzL;
      double beta_zL = betaxL*gxzL + betayL*gyzL +betazL*gzzL;

      double F_rad_xL = gxxL * F_radx_nueL + gxyL * F_rady_nueL + gxzL * F_radz_nueL + beta_xL*F_rad0L;
      double F_rad_yL = gxyL * F_radx_nueL + gyyL * F_rady_nueL + gyzL * F_radz_nueL + beta_yL*F_rad0L;
      double F_rad_zL = gxzL * F_radx_nueL + gyzL * F_rady_nueL + gzzL * F_radz_nueL + beta_zL*F_rad0L;
      double F_rad_0L = - (F_rad_xL*uxL + F_rad_yL*uyL + F_rad_zL*uzL)/u0L;


      double T_fluidL = T_fluid[index];
      //      double T_fluidL = pow((E_rad_nueL/rad_opacity_abs),0.25);
      //      double T_fluidL = PL*al*Psi6*u0L/(rho_starL);
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

   
      double S_rad_x_nue_rhsL = S_rad_x_nue_rhs[index];
      
      if (rad_closure_scheme==0){
	//	printf("Inside souce_rad.C, rad_closure_scheme  is 0!! \n");
      double  P_radL = E_rad_nueL/3.0;

      // add \sqrt{-g}/2 * R^{00}\partial_x g_{00}                                                                             
      S_rad_x_nue_rhsL += 0.5*al*Psi6*((E_rad_nueL+P_radL)*SQR(u0L)
                                   + 2.0*F_rad0L*u0L-P_radL/al/al)*g4ttm;

  
      // add \sqrt{-g} * R^{0i}\partial_x g_{0i}                                                                              
      S_rad_x_nue_rhsL += al*Psi6*
        (   ((E_rad_nueL+P_radL)*SQR(u0L)*vxL + F_rad0L*u0L*vxL+F_radx_nueL*u0L+P_radL*betaxL/al/al)*g4txm +
            ((E_rad_nueL+P_radL)*SQR(u0L)*vyL + F_rad0L*u0L*vyL+F_rady_nueL*u0L+P_radL*betayL/al/al)*g4tym +
            ((E_rad_nueL+P_radL)*SQR(u0L)*vzL + F_rad0L*u0L*vzL+F_radz_nueL*u0L+P_radL*betazL/al/al)*g4tzm );

      // add \sqrt{-g}/2 * R^{ij}\partial_x g_{ij}                                                                                                                                                         
      S_rad_x_nue_rhsL += 0.5*al*Psi6*((E_rad_nueL+P_radL)*SQR(u0L*vxL) +
                           2.0*F_radx_nueL*u0L*vxL + P_radL*(gupxxL*Psim4 -SQR(betaxL)/al/al))*g4xxm;
      S_rad_x_nue_rhsL += al*Psi6*((E_rad_nueL+P_radL)*SQR(u0L)*vxL*vyL +
                       u0L*(F_radx_nueL*vyL+F_rady_nueL*vxL) + P_radL*(gupxyL*Psim4-betaxL*betayL/al/al))*g4xym;
      S_rad_x_nue_rhsL += al*Psi6*((E_rad_nueL+P_radL)*SQR(u0L)*vxL*vzL +
                       u0L*(F_radx_nueL*vzL+F_radz_nueL*vxL) + P_radL*(gupxzL*Psim4-betaxL*betazL/al/al))*g4xzm;
      S_rad_x_nue_rhsL += 0.5*al*Psi6*((E_rad_nueL+P_radL)*SQR(u0L*vyL) +
                           2.0*F_rady_nueL*u0L*vyL + P_radL*(gupyyL*Psim4 -SQR(betayL)/al/al))*g4yym;
      S_rad_x_nue_rhsL += al*Psi6*((E_rad_nueL+P_radL)*SQR(u0L)*vyL*vzL +
                       u0L*(F_rady_nueL*vzL+F_radz_nueL*vyL) + P_radL*(gupyzL*Psim4-betayL*betazL/al/al))*g4yzm;
      S_rad_x_nue_rhsL += 0.5*al*Psi6*((E_rad_nueL+P_radL)*SQR(u0L*vzL) +
                           2.0*F_radz_nueL*u0L*vzL + P_radL*(gupzzL*Psim4 -SQR(betazL)/al/al))*g4zzm;

      // add -\sqrt{-g} * (R^{00} \betaxL + R^{0x})\partial_x alpha to tau_rad_nue_rhs
      // Note that in general, we would need to add time derivative terms of the metric to tau_rad_nue_rhs with the HARM energy variable.
      //  However, the only time we set enable_HARM_energyvariable=1 is with a Cowling (usually disk) run.
     
      if(enable_HARM_energyvariable==0) {
        double alphax_f = D1_o3_cpp(alpha_f,i,j,k);
        tau_rad_nue_rhs[index] +=  -al*Psi6*u0L*( (betaxL+vxL)*((E_rad_nueL+P_radL)*u0L + F_rad0L) +
                                             (F_rad0L*betaxL+F_radx_nueL) )*alphax_f;
      }      
      }
      else{
	//printf("Inside souce_rad.C, rad_closure_scheme  is 1!! \n");
	double P_radxxL,P_radyyL,P_radzzL,P_radxyL,P_radxzL,P_radyzL;
	double zetaL;
	double zeta_temp = sqrt(fabs(F_rad_0L*F_rad0L + F_rad_xL*F_radx_nueL +  F_rad_yL*F_rady_nueL +  F_rad_zL*F_radz_nueL )/SQR(E_rad_nueL));
        double zeta_cut = Erad_atm_cut*1.5;

        if (E_rad_nueL <= zeta_cut){
          zetaL = 1.0;
        }
        else{
          zetaL = zeta_temp;
        }

        if (zetaL > 1.0){
          zetaL = 1.0;
        }

	
	double chiL = 1/3.0 + SQR(zetaL)*(6.0-2.0*zetaL+6*SQR(zetaL))/15.0;
	double Fksqr =  F_rad_xL*F_radx_nueL +  F_rad_yL*F_rady_nueL +  F_rad_zL*F_radz_nueL;
	double Fasq =   F_rad_0L*F_rad0L + F_rad_xL*F_radx_nueL +  F_rad_yL*F_rady_nueL +  F_rad_zL*F_radz_nueL;
	
	compute_M1(P_radxxL, F_radx_nueL, F_radx_nueL, Fasq, E_rad_nueL, gupxxL, betaxL, betaxL, alpha[index], uxL, uxL, chiL, Psim4, Erad_atm_cut);
	compute_M1(P_radyyL, F_rady_nueL, F_rady_nueL, Fasq, E_rad_nueL, gupyyL, betayL, betayL, alpha[index], uyL, uyL, chiL, Psim4, Erad_atm_cut);
        compute_M1(P_radzzL, F_radz_nueL, F_radz_nueL, Fasq, E_rad_nueL, gupzzL, betazL, betazL, alpha[index], uzL, uzL, chiL, Psim4, Erad_atm_cut);
        compute_M1(P_radxyL, F_radx_nueL, F_rady_nueL, Fasq, E_rad_nueL, gupxyL, betaxL, betayL, alpha[index], uxL, uyL, chiL, Psim4, Erad_atm_cut);
        compute_M1(P_radxzL, F_radx_nueL, F_radz_nueL, Fasq, E_rad_nueL, gupxzL, betaxL, betazL, alpha[index], uxL, uzL, chiL, Psim4, Erad_atm_cut);
        compute_M1(P_radyzL, F_rady_nueL, F_radz_nueL, Fasq, E_rad_nueL, gupyzL, betayL, betazL, alpha[index], uyL, uzL, chiL, Psim4, Erad_atm_cut);	

	double P_rad0xL = - (P_radxxL * u_xL + P_radxyL * u_yL + P_radxzL * u_zL)/u_0L;
	double P_rad0yL = - (P_radxyL * u_xL + P_radyyL * u_yL + P_radyzL * u_zL)/u_0L;
	double P_rad0zL = - (P_radxzL * u_xL + P_radyzL * u_yL + P_radzzL * u_zL)/u_0L;
        double P_rad00L = - (P_rad0xL * u_xL + P_rad0yL * u_yL + P_rad0zL * u_zL)/u_0L;

      // add \sqrt{-g}/2 * R^{00}\partial_x g_{00}                                                                               
     	S_rad_x_nue_rhsL += 0.5*al*Psi6*(E_rad_nueL*SQR(u0L) + 2.0*F_rad0L*u0L + P_rad00L)*g4ttm;

      // add \sqrt{-g} * R^{0i}\partial_x g_{0i}                                                                                      
	S_rad_x_nue_rhsL += al*Psi6*
          ( (E_rad_nueL*SQR(u0L)*vxL + F_rad0L*u0L*vxL+F_radx_nueL*u0L+P_rad0xL)*g4txm +
            (E_rad_nueL*SQR(u0L)*vyL + F_rad0L*u0L*vyL+F_rady_nueL*u0L+P_rad0yL)*g4tym +
            (E_rad_nueL*SQR(u0L)*vzL + F_rad0L*u0L*vzL+F_radz_nueL*u0L+P_rad0zL)*g4tzm );

      // add \sqrt{-g}/2 * R^{ij}\partial_x g_{ij}                                                                                          
	S_rad_x_nue_rhsL += 0.5*al*Psi6*(E_rad_nueL*SQR(u0L*vxL) +
                           2.0*F_radx_nueL*u0L*vxL + P_radxxL)*g4xxm;
	S_rad_x_nue_rhsL += al*Psi6*(E_rad_nueL*SQR(u0L)*vxL*vyL +
                       u0L*(F_radx_nueL*vyL+F_rady_nueL*vxL) + P_radxyL)*g4xym;
	S_rad_x_nue_rhsL += al*Psi6*(E_rad_nueL*SQR(u0L)*vxL*vzL +
                       u0L*(F_radx_nueL*vzL+F_radz_nueL*vxL) + P_radxzL)*g4xzm;
	S_rad_x_nue_rhsL += 0.5*al*Psi6*(E_rad_nueL*SQR(u0L*vyL) +
                           2.0*F_rady_nueL*u0L*vyL + P_radyyL)*g4yym;
	S_rad_x_nue_rhsL += al*Psi6*(E_rad_nueL*SQR(u0L)*vyL*vzL +
                       u0L*(F_rady_nueL*vzL+F_radz_nueL*vyL) + P_radyzL)*g4yzm;
	S_rad_x_nue_rhsL += 0.5*al*Psi6*(E_rad_nueL*SQR(u0L*vzL) +
                           2.0*F_radz_nueL*u0L*vzL + P_radzzL)*g4zzm;
	
      // add -\sqrt{-g} * (R^{00} \betaxL + R^{0x})\partial_x alpha to tau_rad_nue_rhs                                     
      // Note that in general, we would need to add time derivative terms of the metric to tau_rad_nue_rhs with the HARM energy variable.  
      //  However, the only time we set enable_HARM_energyvariable=1 is with a Cowling (usually disk) run.

      if(enable_HARM_energyvariable==0) {
        double alphax_f = D1_o3_cpp(alpha_f,i,j,k);
        tau_rad_nue_rhs[index] +=  -al*Psi6*( (betaxL+vxL)*(E_rad_nueL)*SQR(u0L) +u0L* (F_rad0L*(2.0*betaxL + vxL) + F_radx_nueL) +
					  P_rad00L*betaxL + P_rad0xL)*alphax_f;
      }
      }
      // Add radiation 4-force term in RHS, note that is has the same form as in mhd_st_i but with a minus sign;
      // This part is the same for both closure scheme.
      //      S_rad_x_nue_rhsL += -rho_starL/u0L* (rad_opacity_abs*(E_rad_nueL-rad_const*SQR(T_fluidL)*SQR(T_fluidL))*u_xL +                
 				       //				       (rad_opacity_abs+rad_opacity_sct)*F_rad_xL);

      double kappa_a, kappa_s, eta_gf;
      //      double Y_eL = Y_e[index];
      double Y_eL = 0.0;
      double optd_L = optd[index];  
      if (compute_microphysics==1){
	compute_opacity_emissivity(rho_starL, al, Psi6, u0L, kappa_a, kappa_s, eta_gf, T_fluidL, chi_radL, Y_eL, optd_L, T_fluid_cgs_atm);
      }
      else{ // Using constant opacities                                                 
	kappa_a = rad_opacity_abs;
	kappa_s = rad_opacity_sct;
      }
      
      if(enable_OS_collapse==1){
	S_rad_x_nue_rhsL += -rho_starL/u0L*(kappa_a+kappa_s)*F_rad_xL;
      }
      else if(rad_fourforce_enable==1){
	// E_rad_nue_cgs/c has the unit of energy density (~L^-2), which is the same as nkbT                                                                                        
	//S_rad_x_nue_rhsL += -rho_starL/u0L*(kappa_a+kappa_s)*F_rad_xL;
	S_rad_x_nue_rhsL += -al * Psi6 * ((kappa_a*E_rad_nueL-eta_gf)*u_xL + (kappa_a+kappa_s)*F_rad_xL);
      }
    
      S_rad_x_nue_rhs[index] = S_rad_x_nue_rhsL;
	}
  } else if(flux_direction==2) {
    //-----------------------------------------------------------------------------
    // Compute the addition to S_rad_y_nue_rhs and tau_rad_nue_rhs
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

      double E_rad_nueL = E_rad_nue[index];
      double F_radx_nueL = F_radx_nue[index];
      double F_rady_nueL = F_rady_nue[index];
      double F_radz_nueL = F_radz_nue[index];
      double chi_radL =chi_rad[index];

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
      double F_rad0L = - (F_radx_nueL*u_xL + F_rady_nueL*u_yL + F_radz_nueL*u_zL)/u_0L;

      double beta_xL = betaxL*gxxL + betayL*gxyL +betazL*gxzL; 
      double beta_yL = betaxL*gxyL + betayL*gyyL +betazL*gyzL;
      double beta_zL = betaxL*gxzL + betayL*gyzL +betazL*gzzL;
      
      double F_rad_xL = gxxL * F_radx_nueL + gxyL * F_rady_nueL + gxzL * F_radz_nueL + beta_xL*F_rad0L;
      double F_rad_yL = gxyL * F_radx_nueL + gyyL * F_rady_nueL + gyzL * F_radz_nueL + beta_yL*F_rad0L;
      double F_rad_zL = gxzL * F_radx_nueL + gyzL * F_rady_nueL + gzzL * F_radz_nueL + beta_zL*F_rad0L;
      double F_rad_0L = - (F_rad_xL*uxL + F_rad_yL*uyL + F_rad_zL*uzL)/u0L;

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

      
      double S_rad_y_nue_rhsL = S_rad_y_nue_rhs[index];


      if (rad_closure_scheme==0){
      double  P_radL = E_rad_nueL/3.0;

      
      // add \sqrt{-g}/2 * R^{00}\partial_y g_{00}               
      S_rad_y_nue_rhsL += 0.5*al*Psi6*((E_rad_nueL+P_radL)*SQR(u0L)
                      + 2.0*F_rad0L*u0L-P_radL/al/al)*g4ttm;

      
      // add \sqrt{-g} * R^{0i}\partial_y g_{0i}                                             
      S_rad_y_nue_rhsL += al*Psi6*
        (   ((E_rad_nueL+P_radL)*SQR(u0L)*vxL + F_rad0L*u0L*vxL+F_radx_nueL*u0L+P_radL*betaxL/al/al)*g4txm +
            ((E_rad_nueL+P_radL)*SQR(u0L)*vyL + F_rad0L*u0L*vyL+F_rady_nueL*u0L+P_radL*betayL/al/al)*g4tym +
            ((E_rad_nueL+P_radL)*SQR(u0L)*vzL + F_rad0L*u0L*vzL+F_radz_nueL*u0L+P_radL*betazL/al/al)*g4tzm );


      // add \sqrt{-g}/2 * R^{ij}\partial_y g_{ij}
      S_rad_y_nue_rhsL += 0.5*al*Psi6*((E_rad_nueL+P_radL)*SQR(u0L*vxL) +
                           2.0*F_radx_nueL*u0L*vxL + P_radL*(gupxxL*Psim4 -SQR(betaxL)/al/al))*g4xxm;
      S_rad_y_nue_rhsL += al*Psi6*((E_rad_nueL+P_radL)*SQR(u0L)*vxL*vyL +
                       u0L*(F_radx_nueL*vyL+F_rady_nueL*vxL) + P_radL*(gupxyL*Psim4-betaxL*betayL/al/al))*g4xym;
      S_rad_y_nue_rhsL += al*Psi6*((E_rad_nueL+P_radL)*SQR(u0L)*vxL*vzL +
                       u0L*(F_radx_nueL*vzL+F_radz_nueL*vxL) + P_radL*(gupxzL*Psim4-betaxL*betazL/al/al))*g4xzm;
      S_rad_y_nue_rhsL += 0.5*al*Psi6*((E_rad_nueL+P_radL)*SQR(u0L*vyL) +
                           2.0*F_rady_nueL*u0L*vyL + P_radL*(gupyyL*Psim4 -SQR(betayL)/al/al))*g4yym;
      S_rad_y_nue_rhsL += al*Psi6*((E_rad_nueL+P_radL)*SQR(u0L)*vyL*vzL +
                       u0L*(F_rady_nueL*vzL+F_radz_nueL*vyL) + P_radL*(gupyzL*Psim4-betayL*betazL/al/al))*g4yzm;
      S_rad_y_nue_rhsL += 0.5*al*Psi6*((E_rad_nueL+P_radL)*SQR(u0L*vzL) +
                           2.0*F_radz_nueL*u0L*vzL + P_radL*(gupzzL*Psim4 -SQR(betazL)/al/al))*g4zzm;


       // add -\sqrt{-g} * (R^{00} \betayL + R^{0y})\partial_y alpha to tau_rad_nue_rhs
      // Note that in general, we would need to add time derivative terms of the metric to tau_rad_nue_rhs with the HARM energy variable.
      //  However, the only time we set enable_HARM_energyvariable=1 is with a Cowling (usually disk) run.
      
      if(enable_HARM_energyvariable==0) {
        double alphay_f = D2_o3_cpp(alpha_f,i,j,k);
        tau_rad_nue_rhs[index] +=  -al*Psi6*u0L*( (betayL+vyL)*((E_rad_nueL+P_radL)*u0L + F_rad0L) +
                                             (F_rad0L*betayL+F_rady_nueL) )*alphay_f;
      }
      }
       else{
        double P_radxxL,P_radyyL,P_radzzL,P_radxyL,P_radxzL,P_radyzL;
        double zetaL;
        double zeta_temp = sqrt(fabs(F_rad_0L*F_rad0L + F_rad_xL*F_radx_nueL +  F_rad_yL*F_rady_nueL +  F_rad_zL*F_radz_nueL )/SQR(E_rad_nueL));
        double zeta_cut = Erad_atm_cut*1.5;

        if (E_rad_nueL <= zeta_cut){
          zetaL = 1.0;
        }
        else{
          zetaL = zeta_temp;
        }

        if (zetaL > 1.0){
          zetaL = 1.0;
        }

       	double chiL = 1/3.0 + SQR(zetaL)*(6.0-2.0*zetaL+6*SQR(zetaL))/15.0;        
	double Fksqr =  F_rad_xL*F_radx_nueL +  F_rad_yL*F_rady_nueL +  F_rad_zL*F_radz_nueL;
	double Fasq = F_rad_0L*F_rad0L + F_rad_xL*F_radx_nueL +  F_rad_yL*F_rady_nueL +  F_rad_zL*F_radz_nueL;
   
        compute_M1(P_radxxL, F_radx_nueL, F_radx_nueL, Fasq, E_rad_nueL, gupxxL, betaxL, betaxL, alpha[index], uxL, uxL, chiL, Psim4, Erad_atm_cut);
        compute_M1(P_radyyL, F_rady_nueL, F_rady_nueL, Fasq, E_rad_nueL, gupyyL, betayL, betayL, alpha[index], uyL, uyL, chiL, Psim4, Erad_atm_cut);
        compute_M1(P_radzzL, F_radz_nueL, F_radz_nueL, Fasq, E_rad_nueL, gupzzL, betazL, betazL, alpha[index], uzL, uzL, chiL, Psim4, Erad_atm_cut);
        compute_M1(P_radxyL, F_radx_nueL, F_rady_nueL, Fasq, E_rad_nueL, gupxyL, betaxL, betayL, alpha[index], uxL, uyL, chiL, Psim4, Erad_atm_cut);
        compute_M1(P_radxzL, F_radx_nueL, F_radz_nueL, Fasq, E_rad_nueL, gupxzL, betaxL, betazL, alpha[index], uxL, uzL, chiL, Psim4, Erad_atm_cut);
        compute_M1(P_radyzL, F_rady_nueL, F_radz_nueL, Fasq, E_rad_nueL, gupyzL, betayL, betazL, alpha[index], uyL, uzL, chiL, Psim4, Erad_atm_cut);

        double P_rad0xL = - (P_radxxL * u_xL + P_radxyL * u_yL + P_radxzL * u_zL)/u_0L;
        double P_rad0yL = - (P_radxyL * u_xL + P_radyyL * u_yL + P_radyzL * u_zL)/u_0L;
        double P_rad0zL = - (P_radxzL * u_xL + P_radyzL * u_yL + P_radzzL * u_zL)/u_0L;
        double P_rad00L = - (P_rad0xL * u_xL + P_rad0yL * u_yL + P_rad0zL * u_zL)/u_0L;


	// add \sqrt{-g}/2 * R^{00}\partial_y g_{00}                                                                                                     
        S_rad_y_nue_rhsL += 0.5*al*Psi6*(E_rad_nueL*SQR(u0L) + 2.0*F_rad0L*u0L + P_rad00L)*g4ttm;

      // add \sqrt{-g} * R^{0i}\partial_y g_{0i}
	S_rad_y_nue_rhsL += al*Psi6*
          ( (E_rad_nueL*SQR(u0L)*vxL + F_rad0L*u0L*vxL+F_radx_nueL*u0L+P_rad0xL)*g4txm +
            (E_rad_nueL*SQR(u0L)*vyL + F_rad0L*u0L*vyL+F_rady_nueL*u0L+P_rad0yL)*g4tym +
            (E_rad_nueL*SQR(u0L)*vzL + F_rad0L*u0L*vzL+F_radz_nueL*u0L+P_rad0zL)*g4tzm );

      // add \sqrt{-g}/2 * R^{ij}\partial_y g_{ij}                                                                                                                                                                                               
        S_rad_y_nue_rhsL += 0.5*al*Psi6*(E_rad_nueL*SQR(u0L*vxL) +
                           2.0*F_radx_nueL*u0L*vxL + P_radxxL)*g4xxm;
        S_rad_y_nue_rhsL += al*Psi6*(E_rad_nueL*SQR(u0L)*vxL*vyL +
                       u0L*(F_radx_nueL*vyL+F_rady_nueL*vxL) + P_radxyL)*g4xym;
        S_rad_y_nue_rhsL += al*Psi6*(E_rad_nueL*SQR(u0L)*vxL*vzL +
                       u0L*(F_radx_nueL*vzL+F_radz_nueL*vxL) + P_radxzL)*g4xzm;
        S_rad_y_nue_rhsL += 0.5*al*Psi6*(E_rad_nueL*SQR(u0L*vyL) +
                           2.0*F_rady_nueL*u0L*vyL + P_radyyL)*g4yym;
        S_rad_y_nue_rhsL += al*Psi6*(E_rad_nueL*SQR(u0L)*vyL*vzL +
                       u0L*(F_rady_nueL*vzL+F_radz_nueL*vyL) + P_radyzL)*g4yzm;
        S_rad_y_nue_rhsL += 0.5*al*Psi6*(E_rad_nueL*SQR(u0L*vzL) +
                           2.0*F_radz_nueL*u0L*vzL + P_radzzL)*g4zzm;

      // add -\sqrt{-g} * (R^{00} \betayL + R^{0y})\partial_y alpha to tau_rad_nue_rhs                                                                                
      // Note that in general, we would need to add time derivative terms of the metric to tau_rad_nue_rhs with the HARM energy variable.
      //  However, the only time we set enable_HARM_energyvariable=1 is with a Cowling (usually disk) run.                                                        
       
      if(enable_HARM_energyvariable==0) {
        double alphay_f = D2_o3_cpp(alpha_f,i,j,k);
        tau_rad_nue_rhs[index] +=  -al*Psi6*( (betayL+vyL)*(E_rad_nueL)*SQR(u0L) + u0L*(F_rad0L*(2.0*betayL + vyL) + F_rady_nueL) +
                                          P_rad00L*betayL + P_rad0yL)*alphay_f;
      }
       }
      // Add radiation 4-force term in RHS, note that is has the same form as in mhd_st_i but with a minus sign;  

   
      double kappa_a, kappa_s, eta_gf;
      //      double Y_eL = Y_e[index];                                                                                                                  
      double Y_eL = 0.0;
      double optd_L = optd[index];  
      if (compute_microphysics==1){
	compute_opacity_emissivity(rho_starL, al, Psi6, u0L, kappa_a, kappa_s, eta_gf, T_fluidL, chi_radL, Y_eL, optd_L, T_fluid_cgs_atm);
      }
      else{ // Using constant opacities                                                                                                                                                   
        kappa_a = rad_opacity_abs;
        kappa_s = rad_opacity_sct;
      }

      if(enable_OS_collapse==1){
        S_rad_y_nue_rhsL += -rho_starL/u0L*(kappa_a+kappa_s)*F_rad_yL;
      }
      else if (rad_fourforce_enable==1){
	//S_rad_y_nue_rhsL += -rho_starL/u0L*(kappa_a+kappa_s)*F_rad_yL;
	S_rad_y_nue_rhsL += -al * Psi6 * ((kappa_a*E_rad_nueL-eta_gf)*u_yL + (kappa_a+kappa_s)*F_rad_yL);
      }

      S_rad_y_nue_rhs[index]=S_rad_y_nue_rhsL;        
	}
 } else if(flux_direction==3) {
    //-----------------------------------------------------------------------------
    // Compute the addition to S_rad_z_nue_rhs and tau_rad_nue_rhs
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

      double E_rad_nueL = E_rad_nue[index];
      double F_radx_nueL = F_radx_nue[index];
      double F_rady_nueL = F_rady_nue[index];
      double F_radz_nueL = F_radz_nue[index];
      double chi_radL =chi_rad[index];

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
      double F_rad0L = - (F_radx_nueL*u_xL + F_rady_nueL*u_yL + F_radz_nueL*u_zL)/u_0L;

     
      double beta_xL = betaxL*gxxL + betayL*gxyL +betazL*gxzL;
      double beta_yL = betaxL*gxyL + betayL*gyyL +betazL*gyzL;
      double beta_zL = betaxL*gxzL + betayL*gyzL +betazL*gzzL;

      double F_rad_xL = gxxL * F_radx_nueL + gxyL * F_rady_nueL + gxzL * F_radz_nueL + beta_xL*F_rad0L;
      double F_rad_yL = gxyL * F_radx_nueL + gyyL * F_rady_nueL + gyzL * F_radz_nueL + beta_yL*F_rad0L;
      double F_rad_zL = gxzL * F_radx_nueL + gyzL * F_rady_nueL + gzzL * F_radz_nueL + beta_zL*F_rad0L;
      double F_rad_0L = - (F_rad_xL*uxL + F_rad_yL*uyL + F_rad_zL*uzL)/u0L;
      
      double T_fluidL = T_fluid[index];
      //      double T_fluidL = pow((E_rad_nueL/rad_opacity_abs),0.25);
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


      double S_rad_z_nue_rhsL = S_rad_z_nue_rhs[index];


      if (rad_closure_scheme==0){
      double  P_radL = E_rad_nueL/3.0;

      // add \sqrt{-g}/2 * R^{00}\partial_z g_{00}
      S_rad_z_nue_rhsL +=0.5*al*Psi6*((E_rad_nueL+P_radL)*SQR(u0L)
				  + 2.0*F_rad0L*u0L-P_radL/al/al)*g4ttm;

      // add \sqrt{-g} * R^{0i}\partial_z g_{0i} 
      S_rad_z_nue_rhsL += al*Psi6*
        (   ((E_rad_nueL+P_radL)*SQR(u0L)*vxL + F_rad0L*u0L*vxL+F_radx_nueL*u0L+P_radL*betaxL/al/al)*g4txm +
            ((E_rad_nueL+P_radL)*SQR(u0L)*vyL + F_rad0L*u0L*vyL+F_rady_nueL*u0L+P_radL*betayL/al/al)*g4tym +
            ((E_rad_nueL+P_radL)*SQR(u0L)*vzL + F_rad0L*u0L*vzL+F_radz_nueL*u0L+P_radL*betazL/al/al)*g4tzm );

      // add \sqrt{-g}/2 * R^{ij}\partial_z g_{ij}
      S_rad_z_nue_rhsL += 0.5*al*Psi6*((E_rad_nueL+P_radL)*SQR(u0L*vxL) +
                           2.0*F_radx_nueL*u0L*vxL + P_radL*(gupxxL*Psim4 -SQR(betaxL)/al/al))*g4xxm;
      S_rad_z_nue_rhsL += al*Psi6*((E_rad_nueL+P_radL)*SQR(u0L)*vxL*vyL +
                       u0L*(F_radx_nueL*vyL+F_rady_nueL*vxL) + P_radL*(gupxyL*Psim4-betaxL*betayL/al/al))*g4xym;
      S_rad_z_nue_rhsL += al*Psi6*((E_rad_nueL+P_radL)*SQR(u0L)*vxL*vzL +
                       u0L*(F_radx_nueL*vzL+F_radz_nueL*vxL) + P_radL*(gupxzL*Psim4-betaxL*betazL/al/al))*g4xzm;
      S_rad_z_nue_rhsL += 0.5*al*Psi6*((E_rad_nueL+P_radL)*SQR(u0L*vyL) +
                           2.0*F_rady_nueL*u0L*vyL + P_radL*(gupyyL*Psim4 -SQR(betayL)/al/al))*g4yym;
      S_rad_z_nue_rhsL += al*Psi6*((E_rad_nueL+P_radL)*SQR(u0L)*vyL*vzL +
		       u0L*(F_rady_nueL*vzL+F_radz_nueL*vyL) + P_radL*(gupyzL*Psim4-betayL*betazL/al/al))*g4yzm;
      S_rad_z_nue_rhsL += 0.5*al*Psi6*((E_rad_nueL+P_radL)*SQR(u0L*vzL) +
                           2.0*F_radz_nueL*u0L*vzL + P_radL*(gupzzL*Psim4 -SQR(betazL)/al/al))*g4zzm;
      

    // add -\sqrt{-g} * (R^{00} \betazL + R^{0z})\partial_z alpha to tau_rad_nue_rhs
    // Note that in general, we would need to add time derivative terms of the metric to tau_rad_nue_rhs with the HARM energy variable.
    //  However, the only time we set enable_HARM_energyvariable=1 is with a Cowling (usually disk) run.
      
      if(enable_HARM_energyvariable==0) {
        double alphaz_f = D3_o3_cpp(alpha_f,i,j,k);
	tau_rad_nue_rhs[index] += -al*Psi6*u0L*( (betazL+vzL)*((E_rad_nueL+P_radL)*u0L + F_rad0L) +
					     (F_rad0L*betazL+F_radz_nueL) )*alphaz_f;
      }
      }
      else{
        double P_radxxL,P_radyyL,P_radzzL,P_radxyL,P_radxzL,P_radyzL;
	double zetaL;
	double zeta_temp = sqrt(fabs(F_rad_0L*F_rad0L + F_rad_xL*F_radx_nueL +  F_rad_yL*F_rady_nueL +  F_rad_zL*F_radz_nueL )/SQR(E_rad_nueL));
        double zeta_cut = Erad_atm_cut*1.5;

        if (E_rad_nueL <= zeta_cut){
          zetaL = 1.0;
        }
        else{
          zetaL = zeta_temp;
        }

        if (zetaL > 1.0){
          zetaL = 1.0;
        }

	double chiL = 1/3.0 + SQR(zetaL)*(6.0-2.0*zetaL+6*SQR(zetaL))/15.0;        
	double Fksqr =  F_rad_xL*F_radx_nueL +  F_rad_yL*F_rady_nueL +  F_rad_zL*F_radz_nueL;
	double Fasq = F_rad_0L*F_rad0L + F_rad_xL*F_radx_nueL +  F_rad_yL*F_rady_nueL +  F_rad_zL*F_radz_nueL;

	compute_M1(P_radxxL, F_radx_nueL, F_radx_nueL, Fasq, E_rad_nueL, gupxxL, betaxL, betaxL, alpha[index], uxL, uxL, chiL, Psim4, Erad_atm_cut);
        compute_M1(P_radyyL, F_rady_nueL, F_rady_nueL, Fasq, E_rad_nueL, gupyyL, betayL, betayL, alpha[index], uyL, uyL, chiL, Psim4, Erad_atm_cut);
        compute_M1(P_radzzL, F_radz_nueL, F_radz_nueL, Fasq, E_rad_nueL, gupzzL, betazL, betazL, alpha[index], uzL, uzL, chiL, Psim4, Erad_atm_cut);
        compute_M1(P_radxyL, F_radx_nueL, F_rady_nueL, Fasq, E_rad_nueL, gupxyL, betaxL, betayL, alpha[index], uxL, uyL, chiL, Psim4, Erad_atm_cut);
        compute_M1(P_radxzL, F_radx_nueL, F_radz_nueL, Fasq, E_rad_nueL, gupxzL, betaxL, betazL, alpha[index], uxL, uzL, chiL, Psim4, Erad_atm_cut);
        compute_M1(P_radyzL, F_rady_nueL, F_radz_nueL, Fasq, E_rad_nueL, gupyzL, betayL, betazL, alpha[index], uyL, uzL, chiL, Psim4, Erad_atm_cut);
       
        double P_rad0xL = - (P_radxxL * u_xL + P_radxyL * u_yL + P_radxzL * u_zL)/u_0L;
        double P_rad0yL = - (P_radxyL * u_xL + P_radyyL * u_yL + P_radyzL * u_zL)/u_0L;
        double P_rad0zL = - (P_radxzL * u_xL + P_radyzL * u_yL + P_radzzL * u_zL)/u_0L;
        double P_rad00L = - (P_rad0xL * u_xL + P_rad0yL * u_yL + P_rad0zL * u_zL)/u_0L;
	
	// add \sqrt{-g}/2 * R^{00}\partial_z g_{00}                        
        S_rad_z_nue_rhsL += 0.5*al*Psi6*(E_rad_nueL*SQR(u0L) + 2.0*F_rad0L*u0L + P_rad00L)*g4ttm;
	
	// add \sqrt{-g} * R^{0i}\partial_z g_{0i}                                                                                                                           
        S_rad_z_nue_rhsL += al*Psi6*
          ( (E_rad_nueL*SQR(u0L)*vxL + F_rad0L*u0L*vxL+F_radx_nueL*u0L+P_rad0xL)*g4txm +
            (E_rad_nueL*SQR(u0L)*vyL + F_rad0L*u0L*vyL+F_rady_nueL*u0L+P_rad0yL)*g4tym +
            (E_rad_nueL*SQR(u0L)*vzL + F_rad0L*u0L*vzL+F_radz_nueL*u0L+P_rad0zL)*g4tzm );
	
	// add \sqrt{-g}/2 * R^{ij}\partial_z g_{ij}                                                            
        S_rad_z_nue_rhsL += 0.5*al*Psi6*(E_rad_nueL*SQR(u0L*vxL) +
				     2.0*F_radx_nueL*u0L*vxL + P_radxxL)*g4xxm;
        S_rad_z_nue_rhsL += al*Psi6*(E_rad_nueL*SQR(u0L)*vxL*vyL +
				 u0L*(F_radx_nueL*vyL+F_rady_nueL*vxL) + P_radxyL)*g4xym;
        S_rad_z_nue_rhsL += al*Psi6*(E_rad_nueL*SQR(u0L)*vxL*vzL +
				 u0L*(F_radx_nueL*vzL+F_radz_nueL*vxL) + P_radxzL)*g4xzm;
        S_rad_z_nue_rhsL += 0.5*al*Psi6*(E_rad_nueL*SQR(u0L*vyL) +
				     2.0*F_rady_nueL*u0L*vyL + P_radyyL)*g4yym;
        S_rad_z_nue_rhsL += al*Psi6*(E_rad_nueL*SQR(u0L)*vyL*vzL +
				 u0L*(F_rady_nueL*vzL+F_radz_nueL*vyL) + P_radyzL)*g4yzm;
        S_rad_z_nue_rhsL += 0.5*al*Psi6*(E_rad_nueL*SQR(u0L*vzL) +
				     2.0*F_radz_nueL*u0L*vzL + P_radzzL)*g4zzm;
	
      // add -\sqrt{-g} * (R^{00} \betazL + R^{0z})\partial_z alpha to tau_rad_nue_rhs                                               
      // Note that in general, we would need to add time derivative terms of the metric to tau_rad_nue_rhs with the HARM energy variable.     //  However, the only time we set enable_HARM_energyvariable=1 is with a Cowling (usually disk) run.


	if(enable_HARM_energyvariable==0) {
        double alphaz_f = D3_o3_cpp(alpha_f,i,j,k);
        tau_rad_nue_rhs[index] +=  -al*Psi6*( (betazL+vzL)*(E_rad_nueL)*SQR(u0L) + u0L*(F_rad0L*(2.0*betazL + vzL) + F_radz_nueL) +
                                          P_rad00L*betazL + P_rad0zL)*alphaz_f;
	}
      }
      // Add radiation 4-force term in RHS, note that is has the same form as in mhd_st_i but with a minus sign; 
      // This part is the same for both closure scheme.      

      double kappa_a, kappa_s, eta_gf;
      //      double Y_eL = Y_e[index];                                                                                                                         
      double Y_eL = 0.0;
      double optd_L = optd[index];  
      if (compute_microphysics==1){
	compute_opacity_emissivity(rho_starL, al, Psi6, u0L, kappa_a, kappa_s, eta_gf, T_fluidL, chi_radL, Y_eL, optd_L, T_fluid_cgs_atm);
      }
      else{ // Using constant opacities                                                                        
        kappa_a = rad_opacity_abs;
        kappa_s = rad_opacity_sct;
      }
      if(enable_OS_collapse==1){
        S_rad_z_nue_rhsL += -rho_starL/u0L*(kappa_a+kappa_s)*F_rad_zL;
      }
      else if (rad_fourforce_enable==1){
	//	S_rad_z_nue_rhsL += -rho_starL/u0L*(kappa_a+kappa_s)*F_rad_zL;
	S_rad_z_nue_rhsL += -al * Psi6 * ((kappa_a*E_rad_nueL-eta_gf)*u_zL + (kappa_a+kappa_s)*F_rad_zL);
      }
      S_rad_z_nue_rhs[index] = S_rad_z_nue_rhsL;     
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

  //  printf("END rad_source_nue_cpp!!!!!!!!!!!!!!!!!\n");
}

extern "C" void CCTK_FCALL CCTK_FNAME(rad_source_nue_cpp)
  (int *flux_direction, const cGH **cctkGH,int *cctk_lsh, int *nghostzones, int *Symmetry,
   int *enable_HARM_energyvariable,
   double *dX, double *dY,double *dZ,
   double *S_rad_x_nue_rhs, double *S_rad_y_nue_rhs, double *S_rad_z_nue_rhs, 
   double *tau_rad_nue_rhs, double *rho_star, 
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
   double *E_rad_nue, double *F_radx_nue, double *F_rady_nue, double *F_radz_nue, double *Y_e, double *optd,
   int &rad_closure_scheme, double &rad_opacity_abs, double &rad_opacity_sct, double &rad_const, 
   double *ka_gf, double *ks_gf, double *emission_gf, double *chi_rad,
   double *T_fluid, double &Erad_atm_cut, int &enable_OS_collapse, int &compute_microphysics, int &rad_fourforce_enable, double &T_fluid_cgs_atm) 
{
  rad_source_nue_cpp(*flux_direction, *cctkGH,cctk_lsh,nghostzones,*Symmetry,
		       *enable_HARM_energyvariable,
		       *dX, *dY,*dZ,
		       S_rad_x_nue_rhs, S_rad_y_nue_rhs, S_rad_z_nue_rhs, 
		       tau_rad_nue_rhs, rho_star, 
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
		       E_rad_nue,F_radx_nue,F_rady_nue,F_radz_nue,Y_e, optd,
		       rad_closure_scheme,rad_opacity_abs,rad_opacity_sct, rad_const, 
		       ka_gf,  ks_gf, emission_gf, chi_rad,
		 T_fluid, Erad_atm_cut, enable_OS_collapse, compute_microphysics, rad_fourforce_enable, T_fluid_cgs_atm);
}



extern "C" void CCTK_FCALL CCTK_FNAME(rad_tau_scalar_rad_cpp)
  (const cGH **cctkGH,int *cctk_lsh,
   double *tau_rad_nue_rhs, double *rho_star, double *P,
   double *u0, double *vx, double *vy, double *vz,
   double *alpha, double *betax, double *betay, double *betaz,
   double *phi, double *gxx, double *gxy, double *gxz, double *gyy, double *gyz, double *gzz,
   double *E_rad_nue, double *F_radx_nue, double *F_rady_nue, double *F_radz_nue, double *Y_e,
   int &rad_closure_scheme, double &rad_opacity_abs, double &rad_opacity_sct, double &rad_const, 
   double *ka_gf, double *ks_gf, double*chi_rad,  double *optd,
   double *T_fluid, double *tau_rad_scalar_diag, int &enable_OS_collapse, int &compute_microphysics, double &T_fluid_cgs_atm);

extern "C" void rad_tau_scalar_rad_cpp (const cGH *cctkGH,int *cctk_lsh,
   double *tau_rad_nue_rhs, double *rho_star, double *P,
   double *u0, double *vx, double *vy, double *vz,
   double *alpha, double *betax, double *betay, double *betaz,
   double *phi, double *gxx, double *gxy, double *gxz, double *gyy, double *gyz, double *gzz,
   double *E_rad_nue, double *F_radx_nue, double *F_rady_nue, double *F_radz_nue, double *Y_e,
   int &rad_closure_scheme, double &rad_opacity_abs, double &rad_opacity_sct, double &rad_const, 
   double *ka_gf,  double *ks_gf, double*chi_rad,double *optd,
   double *T_fluid,  double *tau_rad_scalar_diag, int &enable_OS_collapse, int &compute_microphysics, double &T_fluid_cgs_atm){

  //  printf("START rad_tau_scalar_rad_cpp!!!!!!!!!!!!!!!!!\n");

#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) for(int j=0;j<cctk_lsh[1];j++) for(int i=0;i<cctk_lsh[0];i++) {
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	double al = 1.0 + alpha[index];
	double Psi = exp(phi[index]);
	double Psi4 = Psi*Psi*Psi*Psi;
	double Psi6 = Psi4*Psi*Psi;

        double rho_starL = rho_star[index];
        double E_rad_nueL = E_rad_nue[index];
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

	
        double F_radx_nueL = F_radx_nue[index];
        double F_rady_nueL = F_rady_nue[index];
        double F_radz_nueL = F_radz_nue[index];
	double chi_radL =chi_rad[index];

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
        double F_rad0L =  - (F_radx_nueL*u_xL + F_rady_nueL*u_yL + F_radz_nueL*u_zL)/u_0L;
       	double T_fluidL = T_fluid[index];
	//	double T_fluidL = PL*al*Psi6*u0L/(rho_starL);
	//	tau_rad_nue_rhs[index] += -al*rho_starL/u0L*(rad_opacity_abs*(E_rad_nueL-rad_const*SQR(T_fluidL)*SQR(T_fluidL))*u0L +
	//					 (rad_opacity_abs+rad_opacity_sct)*F_rad0L);     
	double kappa_a, kappa_s, eta_gf;
	//	double Y_eL = Y_e[index];                                                                                                                     
	double Y_eL = 0.0;
	double optd_L = optd[index];
	if (compute_microphysics==1){
	  compute_opacity_emissivity(rho_starL, al, Psi6, u0L, kappa_a, kappa_s, eta_gf, T_fluidL, chi_radL, Y_eL, optd_L, T_fluid_cgs_atm);
	}
	else{ // Using constant opacities                                                                                                                                        
	  kappa_a = rad_opacity_abs;
	  kappa_s = rad_opacity_sct;
	}
	
	if(enable_OS_collapse==1){
	  tau_rad_nue_rhs[index] += -al*rho_starL/u0L*(kappa_a+kappa_s)*F_rad0L;
	}
	else {
	  //	  tau_rad_nue_rhs[index] += -al*rho_starL/u0L*(kappa_a+kappa_s)*F_rad0L;
	  tau_rad_nue_rhs[index] += -al* al * Psi6 *( (kappa_a*E_rad_nueL - eta_gf)*u0L + (kappa_a+kappa_s)*F_rad0L);
	}

	tau_rad_scalar_diag[index] = -al* al * Psi6 *( (kappa_a*E_rad_nueL - eta_gf)*u0L + (kappa_a+kappa_s)*F_rad0L);
	//tau_rad_scalar_diag[index] = -al*rho_starL/u0L*(kappa_a+kappa_s)*F_rad0L;
      }

  //  printf("END rad_tau_scalar_rad_cpp!!!!!!!!!!!!!!!!!\n");
}
 
extern "C" void CCTK_FCALL CCTK_FNAME(rad_tau_scalar_rad_cpp)
  (const cGH **cctkGH,int *cctk_lsh,
   double *tau_rad_nue_rhs, double *rho_star, double *P,
   double *u0, double *vx, double *vy, double *vz,
   double *alpha, double *betax, double *betay, double *betaz,
   double *phi, double *gxx, double *gxy, double *gxz, double *gyy, double *gyz, double *gzz,
   double *E_rad_nue, double *F_radx_nue, double *F_rady_nue, double *F_radz_nue,double *Y_e,
   int &rad_closure_scheme, double &rad_opacity_abs, double &rad_opacity_sct, double &rad_const, 
   double *ka_gf, double *ks_gf,  double *chi_rad, double *optd,
   double *T_fluid, double *tau_rad_scalar_diag, int &enable_OS_collapse, int &compute_microphysics, double &T_fluid_cgs_atm)
{rad_tau_scalar_rad_cpp(*cctkGH, cctk_lsh,
                        tau_rad_nue_rhs, rho_star, P,
                        u0, vx, vy, vz,
                        alpha, betax, betay, betaz,
                        phi, gxx, gxy, gxz, gyy, gyz, gzz,
                        E_rad_nue, F_radx_nue, F_rady_nue, F_radz_nue, Y_e,
                        rad_closure_scheme, rad_opacity_abs, rad_opacity_sct, rad_const, 
			ka_gf, ks_gf, chi_rad,optd,
                        T_fluid, tau_rad_scalar_diag, enable_OS_collapse, compute_microphysics, T_fluid_cgs_atm);
}
