void quartic_solver(double &rho_s,double &st_xL,double &st_yL,double &st_zL,double &betaxL,double &betayL,double &betazL,
		    double &gxxL,double &gxyL,double &gxzL,double &gyyL,double &gyzL,double &gzzL,
		    double &gupxx_pL,double &gupxy_pL,double &gupxz_pL,double &gupyy_pL,double &gupyz_pL,double &gupzz_pL,
		    double &Psi2,double &Psi4,double &Psi6,double &Psim6,double &Psim4,
		    double &rho_b_atm,double &tau_atm,int &enable_shocktest_mode,int &tau_stilde_fix_enable,int &Symmetry,
		    int &neos, double *rho_tab, double *P_tab, double *eps_tab, double *k_tab, double *gamma_tab, double &gamma_th,
		    int &index,double *tau,double *st_x,double *st_y,double *st_z,double &w_oldL,double *h,double *alpha, double *X,double *Y,double *Z,
		    double *w,double *rho_b,double *P,double *rho,
		    double *Sx,double *Sy,double *Sz,double *Sxx,double *Sxy,double *Sxz,double *Syy,double *Syz,double *Szz,
		    double *u0,double *vx,double *vy,double *vz,double *rho_star) {
  double sdots_o_rhot = 0.98;
  double k_poly = 1.0;
  double n_poly = 1.0;
  int i=1,j=1,k=1;
  double w_l = w[index];
  if(rho_s>=rho_b_atm*Psi6*0.99 || enable_shocktest_mode==1) {
    //           u = (h_old[index]-1.0)*rho_star_old[index]**2 
    //               / (Gamma * w_old[index])
    bool tau_stilde_fix_applied=false;
    double tau_orig = tau[index];
    double st_x_orig = st_x[index];
    double st_y_orig = st_y[index];
    double st_z_orig = st_z[index];

    double s_dot_s=gupxx_pL*st_xL*st_xL+gupyy_pL*st_yL*st_yL+gupzz_pL*st_zL*st_zL+
      2.0*(gupxy_pL*st_xL*st_yL+gupxz_pL*st_xL*st_zL+gupyz_pL*st_yL*st_zL);

    if(isnan(s_dot_s)) printf("before tau_stilde_fix: BAD SDOTS: %e\n",s_dot_s);

    // The following is the tau_stilde fix (see arXiv:0708.2436v3, pg 4, after Eq 12)
    if(tau_stilde_fix_enable==1) {
      if(tau[index]<0.0) {
	tau[index] = tau_atm;
	tau_stilde_fix_applied=true;
      }
      double rho_t=tau[index]*(tau[index]+2.0*rho_s);
      if(s_dot_s > sdots_o_rhot*rho_t) {
	double r_fact=sqrt(s_dot_s/(sdots_o_rhot*rho_t));
	st_xL /= r_fact;
	st_yL /= r_fact;
	st_zL /= r_fact;
	st_x[index] = st_xL;
	st_y[index] = st_yL;
	st_z[index] = st_zL;
	//Recompute s_dot_s!
	s_dot_s = gupxx_pL*st_xL*st_xL+gupyy_pL*st_yL*st_yL+gupzz_pL*st_zL*st_zL+
	  2.0*(gupxy_pL*st_xL*st_yL+gupxz_pL*st_xL*st_zL+gupyz_pL*st_yL*st_zL);
      }
    }

    if(isnan(s_dot_s)) printf("BAD SDOTS: %e %e. gup's: %e %e %e %e %e %e. st_xyz: %e %e %e\n",s_dot_s,sdots_o_rhot,
			      gupxx_pL,gupxy_pL,gupxz_pL,gupyy_pL,gupyz_pL,gupzz_pL,
			      st_xL,st_yL,st_zL);

    double tau_l = tau[index];

    // Find w {
    double alpha_l = 1.0+alpha[index];
    bool fixed = false;
    double h_l = h[index];

    double P_l,eps;
    double first_roots[5];
    double poly_coeffs[4];
    int num_first_roots;
    same_w_solver_tau_n1(w_l,s_dot_s,tau_l,rho_s,h_l,P_l,eps,Psi6,alpha_l,
			 rho_b_atm,k_poly,n_poly,i,j,k,enable_shocktest_mode,fixed,
			 num_first_roots,first_roots,poly_coeffs,tau_stilde_fix_applied);

    //If Font fix was imposed in w_solver_tau_n1(), set tau[index], etc.:
    if (fixed==true) {
      tau_l = w_l*h_l - rho_s - Psi6*P_l;
      tau[index] = tau_l;
    }

    // Next, we compute the 3-velocity and apply a gamma cap if necessary (gamma = alpha u0)!
    double gammaL = w_l/rho_s;
    double fac;
    if (gammaL > MAX_GAMMAM1 && enable_shocktest_mode==0) { 
      w_l = MAX_GAMMAM1 * rho_s;
      fac = rho_s*h_l * sqrt( (MAX_GAMMAM1*MAX_GAMMAM1-1.0)/s_dot_s );
      st_xL = st_xL *fac; 
      st_yL = st_yL *fac;
      st_zL = st_zL *fac;
      st_x[index] = st_xL;
      st_y[index] = st_yL;
      st_z[index] = st_zL;
      tau_l = w_l*h_l - rho_s - Psi6*P_l;
      tau[index]  = tau_l;
      // The following hides the fact that we've done a Font fix, since the Font fix caused
      //    the velocities to go crazy.  That way the ERROR message below is skipped.
      fixed=false;
    }
    double u0L = w_l/(alpha_l*rho_s);
    fac = alpha_l/(w_l*h_l);
    double vxL = ( gupxx_pL*st_xL + gupxy_pL*st_yL + gupxz_pL*st_zL )*fac - betaxL;
    double vyL = ( gupxy_pL*st_xL + gupyy_pL*st_yL + gupyz_pL*st_zL )*fac - betayL;
    double vzL = ( gupxz_pL*st_xL + gupyz_pL*st_yL + gupzz_pL*st_zL )*fac - betazL;

    u0[index] = u0L;
    vx[index] = vxL;
    vy[index] = vyL;
    vz[index] = vzL;

    double rho_bl = rho_s*rho_s/w_l*Psim6;

    // Limit P by imposing both a floor and ceiling:
    double Pmax;
    // We apply the floor/ceiling on P (really K) when we are in the atmosphere, or where the
    //    grav. fields are REALLY strong (e.g., near a puncture).
    //if (rho_bl < 100.0*rho_b_atm || Psi6>1e5) {
    //if (rho_bl < 100.0*rho_b_atm || Psi6>10) {
    if (rho_bl < 100.0*rho_b_atm || Psi6>15) {
      //if (rho_bl < 100.0*rho_b_atm || Psi6>5e1) {
      //if (rho_bl < 100.0*rho_b_atm || Psi6>1e4) {
      double Pcold = k_poly*pow(rho_bl,gamma_th);
      Pmax = 10.0*Pcold;
      //double Pmax = 300.0*Pcold;
      bool hitpmax=false;
      if (P_l > Pmax && enable_shocktest_mode==0) {
	hitpmax=true;
	P_l = Pmax;
	eps = Pmax/((gamma_th-1.0)*rho_bl);
	double h_m = 1.0 + eps + Pmax/rho_bl;
	st_xL = st_xL*h_m/h_l;
	st_yL = st_yL*h_m/h_l;
	st_zL = st_zL*h_m/h_l;
	h_l = h_m;
	tau_l = w_l*h_m - rho_s - Psi6*Pmax ;
	st_x[index] = st_xL;
	st_y[index] = st_yL;
	st_z[index] = st_zL;
	tau[index]  = tau_l;
      }
      double Pmin = 0.5*Pcold;
      //double Pmin = Pmax/10.0;
      bool hitpmin=false;
      if (P_l < Pmin && enable_shocktest_mode==0) {
	hitpmin=true;
	P_l = Pmin;
	eps = Pmin/((gamma_th-1.0)*rho_bl);
	double h_m = 1.0 + eps + Pmin/rho_bl;
	st_xL = st_xL*h_m/h_l;
	st_yL = st_yL*h_m/h_l;
	st_zL = st_zL*h_m/h_l;
	h_l = h_m;
	tau_l = w_l*h_m - rho_s - Psi6*Pmin;
	st_x[index] = st_xL;
	st_y[index] = st_yL;
	st_z[index] = st_zL;
	tau[index]  = tau_l;
      }
    }
    // Set hydro sources:
    w[index] = w_l;
    rho_b[index] = rho_bl;
    h[index] = h_l;
    fac  = 1.0 / ( Psi6 * w_l * h_l ) ;
    P[index]   = P_l;
    rho[index] = h_l * w_l * Psim6 - P_l;
    Sx[index]  = st_xL * Psim6 ;
    Sy[index]  = st_yL * Psim6 ;
    Sz[index]  = st_zL * Psim6 ;
    Sxx[index] = fac * st_xL*st_xL + Psi4 * gxxL * P_l;
    Sxy[index] = fac * st_xL*st_yL + Psi4 * gxyL * P_l;
    Sxz[index] = fac * st_xL*st_zL + Psi4 * gxzL * P_l;
    Syy[index] = fac * st_yL*st_yL + Psi4 * gyyL * P_l;
    Syz[index] = fac * st_yL*st_zL + Psi4 * gyzL * P_l;
    Szz[index] = fac * st_zL*st_zL + Psi4 * gzzL * P_l;

    if(tau_stilde_fix_applied==true && fixed==true) {
      /*
	printf("ERROR: bug detected at (%d,%d,%d): tau_stilde_fix has already been applied.  Font fix should never be necessary!\n",i,j,k);
	printf("dxyz %e %e %e\n",dX,dY,dZ);
	printf("xyz %e %e %e\n",X[index],Y[index],Z[index]);
	printf("polynomial before font fix: %.15e %.15e %.15e %.15e\n",poly_coeffs[0],poly_coeffs[1],poly_coeffs[2],poly_coeffs[3]);
	printf("ROOTS BEFORE FONT FIX: # of roots = %d... roots= %.15e %.15e %.15e %.15e\n",
	num_first_roots,first_roots[1],first_roots[2],first_roots[3],first_roots[4]);

	printf("gij %e %e %e %e %e %e\n",gxx[index]*Psi4,gxy[index]*Psi4,gxz[index]*Psi4,gyy[index]*Psi4,gyz[index]*Psi4,gzz[index]*Psi4);
	//printf("lapbeta %e %e %e %e\n",alpha[index],betax[index],betay[index],betaz[index]);
	printf("rhoPv %e %e %e %e %e\n",rho_b[index],P[index],vx[index],vy[index],vz[index]);
	printf("w, w_old, h %e %e %e\n",w[index],w_old[index],h[index]);
	printf("orig conserv: rhostar_tau_sti %e %e %e %e %e\n",rho_star[index],tau_orig,st_x_orig,st_y_orig,st_z_orig);
	printf("conserv: rhostar_tau_sti %e %e %e %e %e\n",rho_star[index],tau[index],st_x[index],st_y[index],st_z[index]);
      */
      printf("ERROR: bug detected at xyz %e %e %e: tau_stilde_fix has already been applied.  Font fix should never be necessary!\n",X[index],Y[index],Z[index]);

      exit(1);
    }
  } else {
    // Set rho_b to the atmosphere density and u_i=0
    double rho_bl = rho_b_atm;
    //FIXME maybe: Might want to get this from k_tab:
    //FIXME maybe: Might want to get gamma factor from gamma_tab:
    double P_cold = k_poly*pow(rho_bl,gamma_th);
    double eps_cold;
    if (rho_bl != 0.0) {
      eps_cold = P_cold/rho_bl/(gamma_th-1.0);
    } else {
      eps_cold = 0.0;
    }

    double P_l = P_cold;
    double eps = eps_cold;
    double h_l = 1.0 + P_l/rho_bl + eps;
    double alpha_l = 1.0 + alpha[index];
    u0[index] = 1.0/alpha_l;
    vx[index] = -betaxL;
    vy[index] = -betayL;
    vz[index] = -betazL;
    rho_star[index] = Psi6*rho_bl;
    tau[index]   = Psi6*rho_bl*eps;
    st_x[index]  = 0.0;
    st_y[index]  = 0.0;
    st_z[index]  = 0.0;
    h[index]     = h_l;
    w[index]     = rho_star[index];
    rho_b[index] = rho_bl;
    P[index]     = P_l;
    rho[index]   = rho_bl*(1.0+eps);
    Sx[index]    = 0.0;
    Sy[index]    = 0.0;
    Sz[index]    = 0.0;
    Sxx[index]   = Psi4 * gxxL * P_l;
    Sxy[index]   = Psi4 * gxyL * P_l;
    Sxz[index]   = Psi4 * gxzL * P_l;
    Syy[index]   = Psi4 * gyyL * P_l;
    Syz[index]   = Psi4 * gyzL * P_l;
    Szz[index]   = Psi4 * gzzL * P_l;
  }
}
inline int same_gsl_poly_solve_quartic (double a, double b, double c, double d,
					double *x0, double *x1, double *x2, double *x3);

//------------------------------------------------------------------
//  Solve a quartic equation for w.
//  Warning:  This routine is specific to \Gamma = 2.
//------------------------------------------------------------------
void same_w_solver_tau_n1(double &w,double &s_dot_s,double &tau,double &rho_s,double &h,double &P,double &eps,
		     double Psi6,double alpha,double rho_atm, double k_poly,double n_poly,int ii,int j,int k,int enable_shocktest_mode,bool &fixed,
		     int &num_first_roots,double *first_roots,double *poly_coeffs,bool &tau_stilde_fix_applied) {
  int polish=1;

  int num_valid_roots = 0;
  double ws,wr;
  double dif0 = 1.e300;
  double n_p1 = n_poly+1.0;
  double tau2 = tau*tau;
  double rho2 = rho_s*rho_s;
  double rho4 = rho2*rho2;

  double ar[5];

  ar[4] = n_p1*n_p1*((tau+rho_s)*(tau+rho_s) - s_dot_s);
  ar[3] = 2.0*n_p1*rho_s*( (2.0*n_poly+1.0)*rho2+(4.0*n_poly+3.0)*rho_s*tau + 
			   2.0*n_p1*(tau2-s_dot_s) );
  ar[2] = rho2* ( n_poly*(5.0*n_poly+4.0)*rho2+2.0*n_p1*(5.0*n_poly+2.0)*rho_s*tau + 
		  n_p1*( 5.0*n_p1*tau2-2.0*(3.0*n_poly+2.0)*s_dot_s ) );
  ar[1] = 2.0*rho_s*rho2*( n_poly*n_poly*rho2+2.0*n_poly*n_p1*rho_s*tau + 
			   n_p1*(n_p1*tau2-2.0*n_poly*s_dot_s) );
  ar[0] = -n_poly*n_poly*s_dot_s*rho4;

  double ar4inv = 1.0/ar[4];
  double a=ar[3]*ar4inv;
  double b=ar[2]*ar4inv;
  double c=ar[1]*ar4inv;
  double d=ar[0]*ar4inv;

  double roots[5];
  int numroots=same_gsl_poly_solve_quartic (a,b,c,d,&roots[1],&roots[2],&roots[3],&roots[4]);


  //The following chunk of code is used for debugging the Font+tau_stilde_fix failure:
  num_first_roots=numroots;
  poly_coeffs[0] = a;
  poly_coeffs[1] = b;
  poly_coeffs[2] = c;
  poly_coeffs[3] = d;
  for(int i=1;i<=numroots;i++) first_roots[i]=roots[i];

  for(int i=1;i<=numroots;i++) {
    double current_root = roots[i];
    wr = current_root+rho_s;
    if (wr > 0.0) {
      double h_m1 = n_p1*wr*(tau-current_root)/( n_poly*wr*wr + current_root*(current_root+2.0*rho_s) );
      if (h_m1 >= 0.0) {
	num_valid_roots = num_valid_roots + 1;
	double dif = fabs(w-wr);
	if (dif < dif0) {
	  ws = wr;
	  h = h_m1 + 1.0;
	  P = rho2*h_m1/Psi6/n_p1/wr;
	  eps = h_m1/(1.0+1.0/n_poly);
	  dif0 = dif;
	}
      }
    }
  }

  //********************************************************************************************************************//
  // Note that the quartic solver solves for x:
  //        x^4 + a x^3 + b x^2 + c x + d = 0
  //    When |x^4| + |a x^3| + |b x^2| < 1e-12 * (|c x| + |d|), (i.e., x very small) we find that the quartic solver gives
  //    very inaccurate results.  In this case, we just set x=-d/c, which is very close to the actual root, since x is small.
  //    We only worry about this when we've applied the tau_stilde_fix and tau has been reset to tau_atm, since
  //    the Font fix should _NOT_ be needed in this case.  In fact, the code will die if tau has been reset
  //    to tau_atm and the Font fix has been applied.  The following lines of code are designed to avoid that fate.
  // Another way of seeing this is by looking at h_m1.  Here, tau-current_root could be negative if current_root is
  //    inaccurate.
  //********************************************************************************************************************//
  if(tau_stilde_fix_applied==true && num_valid_roots==0) {
    for(int i=1;i<=numroots;i++) {
      double current_root = roots[i];
      wr = current_root+rho_s;
      if (wr > 0.0) {
	double rootval=current_root;
	double rootval2=rootval*rootval;
	double first_three_terms_fabs=rootval2*rootval2 + fabs(a*rootval2*rootval) + fabs(b*rootval2);
	double last_two_terms_fabs=fabs(c*rootval) + fabs(d);
	if(first_three_terms_fabs<1e-12*last_two_terms_fabs) {
	  rootval=-d/c;
	  printf("Non-Fatal Warning: Fixing a root at (i,j,k)=(%d,%d,%d): old root=%.17e, fixed root=%.17e, tau=%.17e\n",
		 ii,j,k,roots[i],rootval,tau);
	  printf("a=%.17e b=%.17e c=%.17e d=%.17e, first_three_terms_fabs=%.17e,last_two_terms_fabs=%.17e \n sdots=%e,rho_s=%e,n_poly=%e,ar4inv=%e\n",a,b,c,d,first_three_terms_fabs,last_two_terms_fabs,s_dot_s,rho_s,n_poly,ar4inv);
	  printf("w,s_dot_s,tau,rho_s,h,P: %e %e %e %e %e %e\n",w,s_dot_s,tau,rho_s,h,P);
	  printf("Psi6,alpha,rho_atm,n_poly,fixed: %e %e %e %e %d\n",Psi6,alpha,rho_atm,n_poly,fixed);
	  printf("If you don't want to see the above warning message as often, increase the value of tau_fact in your .par file!\n");
	  roots[i]=rootval;
	  // Now reset current_root, wr:
	  current_root=roots[i];
	  wr = current_root+rho_s;
	}
      }
      double h_m1 = n_p1*wr*(tau-current_root)/( n_poly*wr*wr + current_root*(current_root+2.0*rho_s) );
      if (h_m1 >= 0.0) {
	num_valid_roots = num_valid_roots + 1;
	double dif = fabs(w-wr);
	if (dif < dif0) {
	  ws = wr;
	  h = h_m1 + 1.0;
	  P = rho2*h_m1/Psi6/n_p1/wr;
	  eps = h_m1/(1.0+1.0/n_poly);
	  dif0 = dif;
	}
      }
    }
  }

  //==================================================================
  // FONT FIX IF num_valid_roots==0 (no physically valide roots found)
  // NOTE: This fix only applies for n=1 polytropic EOS
  //==================================================================
  if (num_valid_roots==0) {
    if(enable_shocktest_mode==1) {
      printf("Had to fontfix, NO ROOTS FOUND: i,j,k = %d %d %d\n", ii,j,k);
      printf("numroots: %d\n",numroots);
      printf("abcd: %e %e %e %e\n",a,b,c,d);
      if(numroots>0) printf("ROOTS: %e %e %e %e\n",roots[1],roots[2],roots[3],roots[4]);
      printf("w,s_dot_s,tau,rho_s,h,P: %e %e %e %e %e %e\n",w,s_dot_s,tau,rho_s,h,P);
      printf("Psi6,alpha,rho_atm,n_poly,fixed: %e %e %e %e %d\n",Psi6,alpha,rho_atm,n_poly,fixed);
      exit(1);
    }
    fixed = true;
    dif0 = 1.e300;
    ar[4] = 1.0;
    ar[3] = 4.0 * rho2 / Psi6;
    ar[2] = 4.0 * rho4 /(Psi6*Psi6) - rho2 - s_dot_s;
    ar[1] = -4.0* rho4 / Psi6;
    ar[0] = -4.0* rho4*rho2 / (Psi6 * Psi6);

    double ar4inv = 1.0/ar[4];
    double a=ar[3]*ar4inv;
    double b=ar[2]*ar4inv;
    double c=ar[1]*ar4inv;
    double d=ar[0]*ar4inv;

    int numroots=same_gsl_poly_solve_quartic (a,b,c,d,&roots[1],&roots[2],&roots[3],&roots[4]);

    // Look for the physical solution
    for(int i=1;i<=numroots;i++) {
      double current_root = roots[i];
      wr = current_root;
      if(wr > 0.0) {
	num_valid_roots = num_valid_roots + 1;
	double dif = fabs(w-wr);
	if (dif < dif0) {
	  ws = wr;
	  dif0 = dif;
	}
      }
    }

    double rho_b = rho2 / (Psi6 * ws);
    P = k_poly*pow(rho_b,(1.0+1.0/n_poly));
    eps = n_poly*P/rho_b;
    h = 1.0 + n_p1*k_poly*pow(rho_b,(1.0/n_poly));
    if (num_valid_roots==0) {
      //w,s_dot_s,tau,rho_s,h,P: nan nan 7.618771e-12 nan nan nan
      //Psi6,alpha,rho_atm,n_poly,fixed: 1.154094e+00 9.527063e-01 1.252265e-11 1.000000e+00 1

      printf("Secondary w_solver also failed at i,j,k = %d %d %d\n", ii,j,k);
      printf("w,s_dot_s,tau,rho_s,h,P: %e %e %e %e %e %e\n",w,s_dot_s,tau,rho_s,h,P);
      printf("Psi6,alpha,rho_atm,n_poly,fixed: %e %e %e %e %d\n",Psi6,alpha,rho_atm,n_poly,fixed);

      printf("\n I'm setting everything to NaN, and I hope this gets overwritten!\n");

      //Note: We set the following variables to NaN's so that the run will continue.  
      //      In some cases (e.g., near a puncture) the primitives solver cannot find 
      //      the root, but only on lower resolution refinement levels.  These low 
      //      resolution regions should get overwritten during the RESRICTION operation.

      P = 1.0/0.0;
      eps = 1.0/0.0;
      h = 1.0/0.0;

      //exit(1);
    }
  }
  // Font fix ends
  w=ws;
}

/* poly/solve_quartic.c
 * 
 * Copyright (C) 2003 CERN and K.S. K\"{o}lbig
 *
 * Converted to C and implemented into the GSL Library 
 * by Andrew W. Steiner and Andy Buckley
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/* solve_quartic.c - finds the real roots of 
 *  x^4 + a x^3 + b x^2 + c x + d = 0
 */

#define SWAPD(a,b) do { double tmp = b ; b = a ; a = tmp ; } while(0)
#define GSL_MAX(a,b) ((a) > (b) ? (a) : (b))

inline int same_gsl_poly_solve_quartic (double a, double b, double c, double d,
				   double *x0, double *x1, double *x2, double *x3)
{
  /* 
   * This code is based on a simplification of
   * the algorithm from zsolve_quartic.c for real roots
   */
  double u[3];
  double aa, pp, qq, rr, rc, sc, tc, mt;
  double w1r, w1i, w2r, w2i, w3r;
  double v[3], v1, v2, arg, theta;
  double disc, h;
  int k1, k2;
  double zarr[4];

  /* For non-degenerate solutions, proceed by constructing and
   * solving the resolvent cubic */
  aa = a * a;
  pp = b - (3.0/8.0) * aa;
  qq = c - (1.0/2.0) * a * (b - (1.0/4.0) * aa);
  rr = d - (1.0/4.0) * (a * c - (1.0/4.0) * aa * (b - (3.0/16.0) * aa));
  rc = (1.0/2.0) * pp;
  sc = (1.0/4.0) * ((1.0/4.0) * pp * pp - rr);
  tc = -((1.0/8.0) * qq * (1.0/8.0) * qq);

  /* This code solves the resolvent cubic in a convenient fashion
   * for this implementation of the quartic. If there are three real
   * roots, then they are placed directly into u[].  If two are
   * complex, then the real root is put into u[0] and the real
   * and imaginary part of the complex roots are placed into
   * u[1] and u[2], respectively. Additionally, this
   * calculates the discriminant of the cubic and puts it into the
   * variable disc. */
  {
    double qcub = (rc * rc - 3 * sc);
    double rcub = (2 * rc * rc * rc - 9 * rc * sc + 27 * tc);

    double Q = qcub / 9;
    double R = rcub / 54;

    double Q3 = Q * Q * Q;
    double R2 = R * R;

    double CR2 = 729 * rcub * rcub;
    double CQ3 = 2916 * qcub * qcub * qcub;

    disc = (CR2 - CQ3) / 2125764.0;

    if (0 == R && 0 == Q)
      {
	u[0] = -rc / 3;
	u[1] = -rc / 3;
	u[2] = -rc / 3;
      }
    else if (CR2 == CQ3)
      {
	double sqrtQ = sqrt (Q);
	if (R > 0)
	  {
	    u[0] = -2 * sqrtQ - rc / 3;
	    u[1] = sqrtQ - rc / 3;
	    u[2] = sqrtQ - rc / 3;
	  }
	else
	  {
	    u[0] = -sqrtQ - rc / 3;
	    u[1] = -sqrtQ - rc / 3;
	    u[2] = 2 * sqrtQ - rc / 3;
	  }
      }
    else if (CR2 < CQ3)
      {
	double sqrtQ = sqrt (Q);
	double sqrtQ3 = sqrtQ * sqrtQ * sqrtQ;
	double theta = acos (R / sqrtQ3);
	if (R / sqrtQ3 >= 1.0) theta = 0.0;
	{
	  double norm = -2 * sqrtQ;
	  
	  u[0] = norm * cos (theta / 3) - rc / 3;
	  u[1] = norm * cos ((theta + 2.0 * M_PI) / 3) - rc / 3;
	  u[2] = norm * cos ((theta - 2.0 * M_PI) / 3) - rc / 3;
	}
      }
    else
      {
	double sgnR = (R >= 0 ? 1 : -1);
	double modR = fabs (R);
	//Zach says: Occasionally, R2-Q3 < 0 when R2, Q3 ~ 1e-100, so we add fabs()
	double sqrt_disc = sqrt(fabs(R2 - Q3));
	//Original line of code:
	//double sqrt_disc = sqrt (R2 - Q3);
	double A = -sgnR * pow (modR + sqrt_disc, 1.0 / 3.0);
	double B = Q / A;
	double mod_diffAB = fabs (A - B);

	u[0] = A + B - rc / 3;
	u[1] = -0.5 * (A + B) - rc / 3;
	u[2] = -(sqrt (3.0) / 2.0) * mod_diffAB;
      }
  }
  /* End of solution to resolvent cubic */

  /* Combine the square roots of the roots of the cubic 
   * resolvent appropriately. Also, calculate 'mt' which 
   * designates the nature of the roots:
   * mt=1 : 4 real roots (disc == 0)
   * mt=2 : 0 real roots (disc < 0)
   * mt=3 : 2 real roots (disc > 0)
   */

  if (0.0 == disc) 
    u[2] = u[1];

  if (0 >= disc)
    {
      mt = 2; 

      /* One would think that we could return 0 here and exit,
       * since mt=2. However, this assignment is temporary and
       * changes to mt=1 under certain conditions below.  
       */
	  
      v[0] = fabs (u[0]);
      v[1] = fabs (u[1]);
      v[2] = fabs (u[2]);

      v1 = GSL_MAX (GSL_MAX (v[0], v[1]), v[2]);
      /* Work out which two roots have the largest moduli */
      k1 = 0, k2 = 0;
      if (v1 == v[0])
	{
	  k1 = 0;
	  v2 = GSL_MAX (v[1], v[2]);
	}
      else if (v1 == v[1])
	{
	  k1 = 1;
	  v2 = GSL_MAX (v[0], v[2]);
	}
      else
	{
	  k1 = 2;
	  v2 = GSL_MAX (v[0], v[1]);
	}

      if (v2 == v[0])
	{
	  k2 = 0;
	}
      else if (v2 == v[1])
	{
	  k2 = 1;
	}
      else
	{
	  k2 = 2;
	}
	  
      if (0.0 <= u[k1]) 
	{
	  w1r=sqrt(u[k1]);
	  w1i=0.0;
	} 
      else 
	{
	  w1r=0.0;
	  w1i=sqrt(-u[k1]);
	}
      if (0.0 <= u[k2]) 
	{
	  w2r=sqrt(u[k2]);
	  w2i=0.0;
	} 
      else 
	{
	  w2r=0.0;
	  w2i=sqrt(-u[k2]);
	}
    }
  else
    {
      mt = 3;

      if (0.0 == u[1] && 0.0 == u[2]) 
	{
	  arg = 0.0;
	} 
      else 
	{
	  arg = sqrt(sqrt(u[1] * u[1] + u[2] * u[2]));
	}
      theta = atan2(u[2], u[1]);
	  
      w1r = arg * cos(theta / 2.0);
      w1i = arg * sin(theta / 2.0);
      w2r = w1r;
      w2i = -w1i;
    }
  
  /* Solve the quadratic to obtain the roots to the quartic */
  w3r = qq / 8.0 * (w1i * w2i - w1r * w2r) / 
    (w1i * w1i + w1r * w1r) / (w2i * w2i + w2r * w2r);
  h = a / 4.0;

  zarr[0] = w1r + w2r + w3r - h;
  zarr[1] = -w1r - w2r + w3r - h;
  zarr[2] = -w1r + w2r - w3r - h;
  zarr[3] = w1r - w2r - w3r - h;
      
  /* Arrange the roots into the variables z0, z1, z2, z3 */
  if (2 == mt)
    {
      if (u[k1] >= 0 && u[k2] >= 0)
	{
	  mt = 1;
	  *x0 = zarr[0];
	  *x1 = zarr[1];
	  *x2 = zarr[2];
	  *x3 = zarr[3];
	}
      else
	{
	  return 0;
	}
    }
  else 
    {
      *x0 = zarr[0];
      *x1 = zarr[1];
    }
  
  /* Sort the roots as usual */
  if (1 == mt)
    {
      /* Roots are all real, sort them by the real part */
      if (*x0 > *x1)
	SWAPD (*x0, *x1);
      if (*x0 > *x2)
	SWAPD (*x0, *x2);
      if (*x0 > *x3)
	SWAPD (*x0, *x3);

      if (*x1 > *x2)
	SWAPD (*x1, *x2);
      if (*x2 > *x3)
	{
	  SWAPD (*x2, *x3);
	  if (*x1 > *x2)
	    SWAPD (*x1, *x2);
	}
      return 4;
    }
  else
    {
      /* 2 real roots */
      if (*x0 > *x1)
	SWAPD (*x0, *x1);
    }

  return 2;
}

