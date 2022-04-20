#include "harm_primitives_recompute_conservs_tau_floor.C"
#include "harm_u2p_util.c"
#include "harm_utoprim_2d.c"
#include "eigen.C"
#include "font_fix.C"
#include "font_fix_gamma_law.C"
#include "primitives_solver_header.h"

int harm_primitives_gammalaw_lowlevel(int &index,double *X,double *Y,double *Z,
				      int &enable_OS_collapse,
				      double *phi,double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
				      double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,
				      double *lapm1,double *shiftx,double *shifty,double *shiftz,
				      double *Bx,double *By,double *Bz,
				      double *T_fluid,
				      double *mhd_st_x,double *mhd_st_y,double *mhd_st_z,double *tau,double *rho_star,
				      double *vx,double *vy,double *vz,double *P,double *rho_b,double *h,double *u0,
				      double *ka_gf, double *ks_gf, double *emission_gf, double *chi_rad, double *chi_rad_nue, double *Y_e, double *optd, double *eta_nue,  
				      double *ka_gf_nue, double *ks_gf_nue, double *emission_gf_nue, double *ka_gf_nux, double *ks_gf_nux, double *emission_gf_nux,
				      double *eps_tot, double *eps_thermal, double *eps_cld, double *P_cld,
				      double &rhobatm,double &tau_atm,
				      int &neos,int &ergo_star, double &ergo_sigma,
				      double *rho_tab,double *P_tab,double *eps_tab,double *k_tab,double *gamma_tab,
				      double &rho_max, 
				      struct output_primitives &new_primitives,
				      struct output_stats &stats,
				      double &Psi6threshold, int &horizon_enforce_rho_profile, 
				      int &rad_evolve_enable, int &compute_microphysics, int &microphysics_scheme, double &T_fluid_cgs_atm, int &rad_fix) {
 
  //  printf("inside harm_primitives_lowlevel.C \n");

  // declare the HARM-style arrays
  FTYPE U[NPR]; 
  FTYPE gcov[NDIM][NDIM]; 
  FTYPE gcon[NDIM][NDIM]; 
  FTYPE detg;
  FTYPE prim[NPR];

  double kpoly = k_tab[0];
  double gamma = gamma_tab[0];
  int gamma_equals2 = 1;
  if (fabs(gamma-2.0) > 1.e-10) gamma_equals2 = 0;

  // These are the tentative best guess for our primitives.
  struct output_primitives best_guess_primitives;

  double phiL     = phi[index];
  double Psi      = exp(phiL);
  double Psi4     = Psi*Psi*Psi*Psi;
  double Psi6     = Psi4*Psi*Psi;
  double Psim4    = 1.0/Psi4;
  double alpha    = lapm1[index]+1.0;
  double alpha_inv= 1.0/alpha;
  detg     = alpha*Psi6; // detg already declared above, used for HARM
  double shiftxL  = shiftx[index];
  double shiftyL  = shifty[index];
  double shiftzL  = shiftz[index];

  //Define physical metric & its inverse:
  double gxx_physL     = gxx[index]*Psi4;
  double gxy_physL     = gxy[index]*Psi4;
  double gxz_physL     = gxz[index]*Psi4;
  double gyy_physL     = gyy[index]*Psi4;
  double gyz_physL     = gyz[index]*Psi4;
  double gzz_physL     = gzz[index]*Psi4;
  double gupxx_physL   = gupxx[index]*Psim4;
  double gupxy_physL   = gupxy[index]*Psim4;
  double gupxz_physL   = gupxz[index]*Psim4;
  double gupyy_physL   = gupyy[index]*Psim4;
  double gupyz_physL   = gupyz[index]*Psim4;
  double gupzz_physL   = gupzz[index]*Psim4;

  // Check to see if the metric is positive-defitive
  double lam1,lam2,lam3;
  double M11 = gxx[index], M12=gxy[index], M13=gxz[index], M22=gyy[index], M23=gyz[index], M33=gzz[index];
  eigenvalues_3by3_real_sym_matrix(lam1, lam2, lam3,M11, M12, M13, M22, M23, M33);
  if (lam1 < 0.0 || lam2 < 0.0 || lam3 < 0) {
     // Metric is not positive-defitive, reset the physical metric to be conformally-flat.
     gxx_physL = Psi4;
     gxy_physL = 0.0;
     gxz_physL = 0.0;
     gyy_physL = Psi4;
     gyz_physL = 0.0;
     gzz_physL = Psi4;
     gupxx_physL = Psim4;
     gupxy_physL = 0.0;
     gupxz_physL = 0.0;
     gupyy_physL = Psim4;
     gupyz_physL = 0.0;
     gupzz_physL = Psim4;
  }

  double shift_xL = gxx_physL*shiftxL + gxy_physL*shiftyL + gxz_physL*shiftzL;
  double shift_yL = gxy_physL*shiftxL + gyy_physL*shiftyL + gyz_physL*shiftzL;
  double shift_zL = gxz_physL*shiftxL + gyz_physL*shiftyL + gzz_physL*shiftzL;
  double beta2L   = shift_xL*shiftxL + shift_yL*shiftyL + shift_zL*shiftzL;


  // Note that the factor here gets us to the object
  // referred to as B^i in the Noble et al paper (and
  // apparently also in the comments to their code).
  // This is NOT the \mathcal{B}^i, which differs by 
  // a factor of the lapse.
  double sqrtfourpi_inv      = 1.0/sqrt(4.0*M_PI);
  double BxL       = Bx[index];
  double ByL       = By[index];
  double BzL       = Bz[index];

  double BxL_over_alpha_sqrt_fourpi = BxL*alpha_inv*sqrtfourpi_inv;
  double ByL_over_alpha_sqrt_fourpi = ByL*alpha_inv*sqrtfourpi_inv;
  double BzL_over_alpha_sqrt_fourpi = BzL*alpha_inv*sqrtfourpi_inv;

  // translate 
  double rho_starL  = rho_star[index];
  double tauL       = tau[index];
  double mhd_st_xL  = mhd_st_x[index];
  double mhd_st_yL  = mhd_st_y[index];
  double mhd_st_zL  = mhd_st_z[index];

  double rho_b_oldL = rho_b[index];
  double P_oldL     = P[index];
  double vxL        = vx[index];
  double vyL        = vy[index];
  double vzL        = vz[index];
  double u0L        = u0[index];
  

 
  /*
    -- Driver for new prim. var. solver.  The driver just translates
    between the two sets of definitions for U and P.  The user may 
    wish to alter the translation as they see fit.  


    //         /  rho u^t           \                             //
    //    U =  |  T^t_t   + rho u^t |  sqrt(-det(g_{\mu\nu}))     //
    //         |  T^t_i             |                             //
    //         \   B^i              /                             //
    //                                                            //
    //        /    rho        \                                   //
    //    P = |    uu         |                                   //
    //        | \tilde{u}^i   |                                   //
    //        \   B^i         /                                   //

    (above equations have been fixed by Yuk Tung & Zach)
  */
  
  // U[NPR]    = conserved variables (current values on input/output);
  // gcov[NDIM][NDIM] = covariant form of the metric ;
  // gcon[NDIM][NDIM] = contravariant form of the metric ;
  // gdet             = sqrt( - determinant of the metric) ;
  // prim[NPR] = primitive variables (guess on input, calculated values on 
  //	 			    output if there are no problems);

  // U[1]   =  
  // U[2-4] =  stildei + rhostar
		    
  gcov[0][0] = -alpha*alpha + beta2L;
  gcov[0][1] = shift_xL;
  gcov[0][2] = shift_yL;
  gcov[0][3] = shift_zL;
  gcov[1][1] = gxx_physL;
  gcov[1][2] = gxy_physL;
  gcov[1][3] = gxz_physL;
  gcov[2][2] = gyy_physL;
  gcov[2][3] = gyz_physL;
  gcov[3][3] = gzz_physL;
  gcov[1][0] = gcov[0][1];
  gcov[2][0] = gcov[0][2];
  gcov[3][0] = gcov[0][3];
  gcov[2][1] = gcov[1][2];
  gcov[3][1] = gcov[1][3];
  gcov[3][2] = gcov[2][3];
		
  if (isnan(gcov[1][1])||isinf(gcov[1][1])||gcov[1][1]>1.0e100){ printf("gcov[1][1] is nan or inf, gcov[1][1], Psi4, gxx[index] are %e,%e,%e \n",  gcov[1][1], Psi4, gxx[index] );} 

    
  gcon[0][0] = -1.0/(alpha*alpha);
  gcon[0][1] = shiftxL/(alpha*alpha);
  gcon[0][2] = shiftyL/(alpha*alpha);
  gcon[0][3] = shiftzL/(alpha*alpha);
  gcon[1][1] = gupxx_physL - shiftxL*shiftxL/(alpha*alpha);
  gcon[1][2] = gupxy_physL - shiftxL*shiftyL/(alpha*alpha);
  gcon[1][3] = gupxz_physL - shiftxL*shiftzL/(alpha*alpha);
  gcon[2][2] = gupyy_physL - shiftyL*shiftyL/(alpha*alpha);
  gcon[2][3] = gupyz_physL - shiftyL*shiftzL/(alpha*alpha);
  gcon[3][3] = gupzz_physL - shiftzL*shiftzL/(alpha*alpha);
  gcon[1][0] = gcon[0][1];
  gcon[2][0] = gcon[0][2];
  gcon[3][0] = gcon[0][3];
  gcon[2][1] = gcon[1][2];
  gcon[3][1] = gcon[1][3];
  gcon[3][2] = gcon[2][3];

  //printf("rhobatm= %e\n",rhobatm); 

  double tau_orig = tauL;
  double rho_star_orig = rho_starL;
  double mhd_st_x_orig = mhd_st_xL;
  double mhd_st_y_orig = mhd_st_yL;
  double mhd_st_z_orig = mhd_st_zL;

  double GAMMA_SPEED_LIMIT = 100.0;

  //if(Psi6>Psi6threshold) GAMMA_SPEED_LIMIT=500.0;
  if(Psi6>Psi6threshold) GAMMA_SPEED_LIMIT=100.0;

  for(int which_guess=0;which_guess<3;which_guess++) {
	    
    int check;
    double eps_oldL;
    
    if(which_guess==1) {
      //Use a different initial guess:
      rho_b_oldL = rho_starL/Psi6;

      compute_pcold_cpp(rho_b_oldL, P_oldL, neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab);
      compute_epscold_cpp(rho_b_oldL, P_oldL, eps_oldL, neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab, enable_OS_collapse);
      u0L = alpha_inv;
      vxL = -shiftxL;
      vyL = -shiftyL;
      vzL = -shiftzL;
    }

    if(which_guess==2) {
      //Use atmosphere as initial guess:
      rho_b_oldL = 100.0*rhobatm;

      compute_pcold_cpp(rho_b_oldL, P_oldL, neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab);
      compute_epscold_cpp(rho_b_oldL, P_oldL, eps_oldL, neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab, enable_OS_collapse);
      u0L = alpha_inv;
      vxL = -shiftxL;
      vyL = -shiftyL;
      vzL = -shiftzL;
    }

    // Fill the array of conserved variables according to the wishes of Utoprim_2d.
    U[RHO]    = rho_starL;
    U[UU]     = -tauL*alpha - lapm1[index]*rho_starL + shiftxL*mhd_st_xL + shiftyL*mhd_st_yL  + shiftzL*mhd_st_zL ; // note the minus sign on tauL!
    U[UTCON1] = mhd_st_xL;
    U[UTCON2] = mhd_st_yL;
    U[UTCON3] = mhd_st_zL;
    U[BCON1]  = detg*BxL_over_alpha_sqrt_fourpi;
    U[BCON2]  = detg*ByL_over_alpha_sqrt_fourpi;
    U[BCON3]  = detg*BzL_over_alpha_sqrt_fourpi;


    double P_coldL, eps_coldL;
    //compute P_cold using rho_b_oldL and EOS:
    compute_pcold_cpp(rho_b_oldL, P_coldL, neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab);
    //compute eps_cold using P_cold, rho_b_oldL and EOS:
    compute_epscold_cpp(rho_b_oldL, P_coldL, eps_coldL, neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab, enable_OS_collapse); 
    
    // P = P_cold + rhob(Gamma_th-1)(eps-eps_cold) = P_cold + (Gamma_th-1)(u-u_cold)                                                    
    // --> u = (P-P_cold)/(Gamma_th-1) + u_cold = (P-P_cold)/(Gamma_th-1) + eps_cold * rhob

    double uL;
    if (compute_microphysics==1){
    eps_oldL = h[index] - 1.0 - P_oldL/rho_b_oldL;
    uL = rho_b_oldL*eps_oldL;
    }else{
    uL = (P_oldL-P_coldL)/(GAMMA_th - 1.0) + eps_coldL*rho_b_oldL;
    }
    double utxL = u0L*(vxL + shiftxL);
    double utyL = u0L*(vyL + shiftyL);
    double utzL = u0L*(vzL + shiftzL);

    prim[RHO]    = rho_b_oldL;
    prim[UU]     = uL;
    prim[UTCON1] = utxL;
    prim[UTCON2] = utyL;
    prim[UTCON3] = utzL;
    prim[BCON1]  = BxL_over_alpha_sqrt_fourpi;
    prim[BCON2]  = ByL_over_alpha_sqrt_fourpi;
    prim[BCON3]  = BzL_over_alpha_sqrt_fourpi;

    double vsq = SQR(vxL) + SQR(vyL) + SQR(vzL);
    double shiftsq = SQR(shiftxL) + SQR(shiftyL) + SQR(shiftzL);
    
    /*
    if(vsq==0.0){printf("before Utoprim_2d, vx=vy=vz= 0 \n");}
    
    if(shiftsq==0.0){printf("before Utoprim_2d, shiftx=shifty=shiftz= 0 \n");}

    if (u0L==0.0){printf("before Utoprim_2d, u0= 0 \n");}

    if(prim[UTCON1]*prim[UTCON2]*prim[UTCON3] == 0.0) {

      if(prim[UTCON1]==0.0){printf("before Utoprim_2d, prim[UTCON1] is 0 \n");}
      if(prim[UTCON2]==0.0){printf("before Utoprim_2d, prim[UTCON2] is 0 \n");}
      if(prim[UTCON3]==0.0){printf("before Utoprim_2d, prim[UTCON3] is 0 \n");}
    
      printf("u0L=%.17g, vxL=%.17g, shiftxL=%.17g, vyL=%.17g, shiftyL=%.17g, vzL=%.17g, shiftzL=%.17g \n", u0L, vxL, shiftxL,vyL, shiftyL,vzL, shiftzL);
      printf("which_guess = %d \n", which_guess);
    }
    */
    
    if(detg < 1.0e-40) { printf("BEFORE Utoprim_2d, detg is too small to be correct.  gdet, alpha, Psi6 are %e, %e, %e, \n", detg, alpha, Psi6);}

    double T_fluidL;
    /*************************************************************/
    // CALL HARM PRIMITIVES SOLVER:
    check = Utoprim_2d(U, gcov, gcon, detg, prim, neos, ergo_star, ergo_sigma, T_fluidL, rho_tab, P_tab, eps_tab, k_tab, gamma_tab, enable_OS_collapse, compute_microphysics);
    // Note that we have modified this solver, so that nearly 100% 
    // of the time it yields either a good root, or a root with 
    // negative epsilon (i.e., pressure).
    /*************************************************************/
    /*************************************************************/
    // CALL FONT FIX, IF HARM PRIMITIVES SOLVER DIDN'T FIND A ROOT
    // Use the new Font fix subroutine 
   
    //    if (check != 0) { printf ( "primitive_solver returns: %d\n", check);}



    int font_fix_applied=0;
    if(check!=0) {
      font_fix_applied=1;
      //printf("Font fix applied!\n");
      double u_xl, u_yl, u_zl, rhob;
      if (gamma_equals2==1) {
         check = font_fix_gamma_equals2(u_xl,u_yl,u_zl,rhob,rho_starL,mhd_st_xL,mhd_st_yL,mhd_st_zL,
             BxL,ByL,BzL,gxx_physL,gxy_physL,gxz_physL,gyy_physL,gyz_physL,gzz_physL,
	     gupxx_physL,gupxy_physL,gupxz_physL,gupyy_physL,gupyz_physL,gupzz_physL,Psi6,kpoly,
	     neos, ergo_star, ergo_sigma, rho_tab, P_tab, eps_tab, k_tab, gamma_tab,enable_OS_collapse);
      } else {
         check = font_fix_general_gamma(u_xl,u_yl,u_zl,rhob,
             rho_starL,mhd_st_xL,mhd_st_yL,mhd_st_zL,BxL,ByL,BzL,
             gxx_physL,gxy_physL,gxz_physL,gyy_physL,gyz_physL,gzz_physL,
             gupxx_physL,gupxy_physL,gupxz_physL,gupyy_physL,gupyz_physL,gupzz_physL,
	     Psi6, gamma, kpoly, neos, ergo_star, ergo_sigma, rho_tab, P_tab, eps_tab, k_tab, gamma_tab,enable_OS_collapse);
      }

      //Translate to HARM primitive now:
      prim[UTCON1] = gupxx_physL*u_xl + gupxy_physL*u_yl + gupxz_physL*u_zl;
      prim[UTCON2] = gupxy_physL*u_xl + gupyy_physL*u_yl + gupyz_physL*u_zl;
      prim[UTCON3] = gupxz_physL*u_xl + gupyz_physL*u_yl + gupzz_physL*u_zl;

      //      if(isnan(prim[UTCON1])) { printf("prim[UTCON1] is nan, u_xl, u_yl, u_zl are %e,%e,%e \n",  u_xl, u_yl, u_zl);}

      if (check==1) {
	printf("Font fix failed!\n");
	printf("x,y,z = %e %e %e, st_i = %e %e %e, rhostar = %e, Bi = %e %e %e, gij = %e %e %e %e %e %e, Psi6 = %e\n",X[index],Y[index],Z[index],mhd_st_x_orig,mhd_st_y_orig,mhd_st_z_orig,rho_star_orig,BxL,ByL,BzL,gxx_physL,gxy_physL,gxz_physL,gyy_physL,gyz_physL,gzz_physL,Psi6);
      }
      //printf("after font_fix, check is %d \n", check);
    }
    stats.font_fixed=font_fix_applied;
    /*************************************************************/

    if(check==0) {
      //Now that we have found some solution, we first set the best_guess_primitives struct & limit velocity:

      double u_x_new,u_y_new,u_z_new,au0m1; // <- needed below for recomputing conservatives
    
      double utx_new = prim[UTCON1];
      double uty_new = prim[UTCON2];
      double utz_new = prim[UTCON3];
     
      //Velocity limiter:
      double gijuiuj = gxx_physL*SQR(utx_new ) +
	2.0*gxy_physL*utx_new*uty_new + 2.0*gxz_physL*utx_new*utz_new +
	gyy_physL*SQR(uty_new) + 2.0*gyz_physL*uty_new*utz_new +
	gzz_physL*SQR(utz_new);
      au0m1 = gijuiuj/( 1.0+sqrt(1.0+gijuiuj) );
      
      if(isnan(au0m1)) { printf("CHECK==0 au0m1 is nan, gijuiuj, utx_new, uty_new, utz_new are %e,%e,%e,%e \n",  gijuiuj, utx_new, uty_new, utz_new);}
      
      best_guess_primitives.u0_new = (au0m1+1.0)*alpha_inv;
      // *** Limit velocity
      if (au0m1 > GAMMA_SPEED_LIMIT-1.0) {
	double fac = sqrt((SQR(GAMMA_SPEED_LIMIT)-1.0)/(SQR(1.0+au0m1) - 1.0));
	utx_new *= fac;
	uty_new *= fac;
	utz_new *= fac;
	gijuiuj = gijuiuj * SQR(fac);
	au0m1 = gijuiuj/( 1.0+sqrt(1.0+gijuiuj) );
	// Reset rho_b and u0
	best_guess_primitives.u0_new = (au0m1+1.0)*alpha_inv;
	prim[RHO] =  rho_star_orig/(alpha*best_guess_primitives.u0_new*Psi6);
      }//Finished limiting velocity

      //The Font fix only sets the velocities.  Here we set the pressure & density HARM primitives.

           
      if(font_fix_applied==1) {
	prim[RHO] = rho_star_orig/(alpha*best_guess_primitives.u0_new*Psi6);
	if(isnan(best_guess_primitives.u0_new)) { printf("font_fix_applied==1 & best_guess_primitives.u0_new \n");}
	//Next set P = P_cold:
	double P_cold, eps_cold;
	compute_pcold_cpp(prim[RHO], P_cold,neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab);
	compute_epscold_cpp(prim[RHO], P_cold, eps_cold, neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab, enable_OS_collapse);
	prim[UU] = eps_cold * prim[RHO];
      } //Finished setting remaining primitives if there was a Font fix.
      
      double P_cold, eps_cold;
      best_guess_primitives.rho_b_new = prim[RHO];     
      compute_pcold_cpp(best_guess_primitives.rho_b_new, P_cold, neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab);
      compute_epscold_cpp(best_guess_primitives.rho_b_new, P_cold, eps_cold, neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab, enable_OS_collapse);
     
      //double T_fluidL = T_fluid[index];
      best_guess_primitives.vx_new  = utx_new/best_guess_primitives.u0_new - shiftxL;
      best_guess_primitives.vy_new  = uty_new/best_guess_primitives.u0_new - shiftyL;
      best_guess_primitives.vz_new  = utz_new/best_guess_primitives.u0_new - shiftzL;
      u_x_new = gxx_physL*utx_new + gxy_physL*uty_new + gxz_physL*utz_new;
      u_y_new = gxy_physL*utx_new + gyy_physL*uty_new + gyz_physL*utz_new;
      u_z_new = gxz_physL*utx_new + gyz_physL*uty_new + gzz_physL*utz_new;

      //we floor rho_b if necessary:
      if(best_guess_primitives.rho_b_new<rhobatm) best_guess_primitives.rho_b_new=rhobatm;

      //Also set horizon cap density 
      if(Psi6 > Psi6threshold) {
	double rho_horiz_cap = 1.0e5*rhobatm;
	if(best_guess_primitives.rho_b_new >= rho_horiz_cap){
	  best_guess_primitives.rho_b_new = rho_horiz_cap;
	}
      }

      double eps_new = prim[UU]/prim[RHO];
      double T_fluidL;

      if (compute_microphysics == 1){
        compute_P_T_microphys(best_guess_primitives.P_new, T_fluidL, P_cold, eps_new, eps_cold, best_guess_primitives.rho_b_new);
        //!!!!The line below is just to reproduced the old way of getting P. i.e. compute_P_T_microphys is used only for Temperature.!!!!!                                       
	//        best_guess_primitives.P_new = P_cold + (GAMMA_th-1.0)*(prim[UU]-best_guess_primitives.rho_b_new*eps_cold);		      
      }
      else {
	//The line below is for diagnostic Temperature only!
	//	compute_P_T_microphys(best_guess_primitives.P_new, T_fluidL, P_cold, eps_new, eps_cold, best_guess_primitives.rho_b_new);
        best_guess_primitives.P_new = P_cold + (GAMMA_th-1.0)*(prim[UU]-best_guess_primitives.rho_b_new*eps_cold);
      }

      // We may have reset the density but we need to reset the pressure anyway
      // in general so that we place a cap and floor on the ratio P/Pcold      
	
      // compute K 
      double knew  = best_guess_primitives.P_new/P_cold; 
      double klow=0.9;
      if(knew< klow){knew = klow;}      
      
      // Set potentially a different cap on K inside the horizon
      //      double khigh_BH=100.;
      double khigh_BH=100.;
      double khigh=100.;
      if(Psi6 > Psi6threshold) {
	if(knew> khigh_BH){knew = khigh_BH;}
      }else{
	if(knew> khigh){knew = khigh_BH;}
      }


      //Next change P according to the rho_b value and the new K:
      //      compute_pcold_cpp(best_guess_primitives.rho_b_new, best_guess_primitives.P_new, neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab);

      best_guess_primitives.P_new = P_cold*knew;
 
      double eps_totL, eps_cldL, P_cldL;
      //Now that P and possibly rho_b have been changed, recompute conservatives:
      recompute_conservatives_fast(best_guess_primitives,
				   enable_OS_collapse, compute_microphysics,
				   Psi6,alpha,
				   gxx_physL,gxy_physL,gxz_physL,gyy_physL,gyz_physL,gzz_physL,
				   gupxx_physL,gupxy_physL,gupxz_physL,gupyy_physL,gupyz_physL,gupzz_physL,
				   BxL,ByL,BzL,T_fluidL,eps_totL, eps_cldL, P_cldL, 
				   rho_starL,tauL,mhd_st_xL,mhd_st_yL,mhd_st_zL,
				   u_x_new,u_y_new,u_z_new,au0m1,
				   neos,ergo_star,ergo_sigma,rho_tab,P_tab,eps_tab,
				   k_tab,gamma_tab,rho_max,rhobatm);
    
      
      eps_tot[index] = eps_totL;
      eps_thermal[index] = eps_totL - eps_cldL;
      eps_cld[index] = eps_cldL;
      P_cld[index] = P_cldL;
      T_fluid[index] = T_fluidL;      

      double kappa_a, kappa_s, eta_gf, kappa_a_nue, kappa_s_nue, eta_gf_nue, kappa_a_nux, kappa_s_nux, eta_gf_nux;
      double chi_radL = chi_rad[index];
      double chi_rad_nueL = chi_rad_nue[index];
      double Y_eL = Y_e[index];
      //double Y_eL = 0.0;
      double optd_L = optd[index];
      double eta_nueL = eta_nue[index];
      if (rad_evolve_enable==1){
	compute_opacity_emissivity(rho_starL, alpha, Psi6, best_guess_primitives.u0_new, kappa_a, kappa_s, eta_gf, kappa_a_nue, kappa_s_nue, eta_gf_nue, kappa_a_nux, kappa_s_nux, eta_gf_nux, T_fluidL, chi_radL, chi_rad_nueL, Y_eL, optd_L, eta_nueL, T_fluid_cgs_atm, microphysics_scheme, rad_fix);
	ka_gf[index] = kappa_a;
	ks_gf[index] = kappa_s;
	emission_gf[index] = eta_gf;
	ka_gf_nue[index] = kappa_a_nue;
	ks_gf_nue[index] = kappa_s_nue;
	emission_gf_nue[index] = eta_gf_nue;
	ka_gf_nux[index] = kappa_a_nux;
	ks_gf_nux[index] = kappa_s_nux;
	emission_gf_nux[index] = eta_gf_nux;
      }
      else{
	ka_gf[index] = 0.0;
        ks_gf[index] = 0.0;
        emission_gf[index] = 0.0;
        ka_gf_nue[index] = 0.0;
        ks_gf_nue[index] = 0.0;
        emission_gf_nue[index] = 0.0;
        ka_gf_nux[index] = 0.0;
        ks_gf_nux[index] = 0.0;
        emission_gf_nux[index] = 0.0;
      }

      // We save the latest conservatives & primitives here.
      rho_star[index] = rho_starL;
      tau[index] = tauL;
      mhd_st_x[index] = mhd_st_xL;
      mhd_st_y[index] = mhd_st_yL;
      mhd_st_z[index] = mhd_st_zL;

      //if(best_guess_primitives.rho_b_new < rhobatm*0.99) { printf("ERROR: rhob<rhobatm ->  %e < %e\n",best_guess_primitives.rho_b_new,rhobatm); }
      new_primitives.rho_b_new  = best_guess_primitives.rho_b_new;
      new_primitives.P_new      = best_guess_primitives.P_new;
      new_primitives.u0_new     = best_guess_primitives.u0_new;
      new_primitives.vx_new     = best_guess_primitives.vx_new;
      new_primitives.vy_new     = best_guess_primitives.vy_new;
      new_primitives.vz_new     = best_guess_primitives.vz_new;

      return 0;
    } else {
      //If we didn't find a root, then try again with a different guess.
    }
  }
  printf("Couldn't find root from: %e %e %e %e %e, rhob approx=%e, rhobatm=%e, Bx=%e, By=%e, Bz=%e, gij_phys=%e %e %e %e %e %e, alpha=%e\n",
	 tau_orig,rho_star_orig,mhd_st_x_orig,mhd_st_y_orig,mhd_st_z_orig,rho_star_orig/Psi6,rhobatm,BxL,ByL,BzL,gxx_physL,gxy_physL,gxz_physL,gyy_physL,gyz_physL,gzz_physL,alpha);
  return 1;
}

