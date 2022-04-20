#define SQR(x) ((x) * (x))


void compute_pcold_cpp(double &rhob, double &P_cold, 
                       int &neos, int &ergo_star, double &ergo_sigma, double *rho_tab,double *P_tab,double *eps_tab,
		       double *k_tab,double *gamma_tab);

void compute_epscold_cpp(double &rhob, double &P_cold, double &eps_cold,
                       int &neos, int &ergo_star, double &ergo_sigma, double *rho_tab,double *P_tab,double *eps_tab,
			 double *k_tab,double *gamma_tab,int &enable_OS_collapse);

void compute_eps_T_microphys (double &eps, double &T_fluid, double &eps_cold,  double &P, double &P_cold, double &rho_b);

void eigenvalues_3by3_real_sym_matrix(double & lam1, double & lam2, double & lam3,
         double M11, double M12, double M13, double M22, double M23, double M33);

int apply_tau_floor(int &index,double &tau_atm,double &rhobatm,
		    double *Bx,double *By,double *Bz,
		    double *tau,double *rho_star,double *mhd_st_x,double *mhd_st_y,double *mhd_st_z,
		    double *phi,double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
		    double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,
		    struct output_primitives &new_primitives,
		    double *lapm1,double *shiftx,double *shifty,double *shiftz,
		    int &neos, int &ergo_star, double &ergo_sigma, double *rho_tab,double *P_tab,double *eps_tab,double *k_tab,double *gamma_tab,
		    double &Psi6threshold) {
  
  double Psi2 = exp(2.0*phi[index]);
  double Psi4 = Psi2*Psi2;
  double Psi6 = Psi4*Psi2;
  double Psim4 = 1.0/Psi4;

  double gamma = gamma_tab[0];
  double kpoly = k_tab[0];
  int gamma_equals2 = 1;
  if (fabs(gamma-2.0) > 1.e-10) gamma_equals2 = 0;

  //First apply the rho_star floor:

  //rho_star = alpha u0 Psi6 rho_b, alpha u0 > 1, so if rho_star < Psi6 rhobatm, then we are GUARANTEED that we can reset to atmosphere.
  //if(rho_star[index] < 1e4*Psi6*rhobatm) {
  //if(rho_star[index] < 2*Psi6*rhobatm) {
  double gxxL = gxx[index];
  double gxyL = gxy[index];
  double gxzL = gxz[index];
  double gyyL = gyy[index];
  double gyzL = gyz[index];
  double gzzL = gzz[index];

  double lam1,lam2,lam3;
  eigenvalues_3by3_real_sym_matrix(lam1, lam2, lam3,gxxL, gxyL, gxzL, gyyL, gyzL, gzzL);

  double gxxi=gupxx[index]*Psim4,gyyi=gupyy[index]*Psim4,gzzi=gupzz[index]*Psim4,gxyi=gupxy[index]*Psim4,gxzi=gupxz[index]*Psim4,gyzi=gupyz[index]*Psim4;

  if (lam1 < 0.0 || lam2 < 0.0 || lam3 < 0.0) {
     // Metric is not positive-defitive, reset the metric to be conformally-flat.
     gxxL = 1.0;
     gxyL = 0.0;
     gxzL = 0.0;
     gyyL = 1.0;
     gyzL = 0.0;
     gzzL = 1.0;
     gxxi = Psim4;
     gxyi = 0.0;
     gxzi = 0.0;
     gyyi = Psim4;
     gyzi = 0.0;
     gzzi = Psim4;
  }

  //Next, prepare for the tau and stilde fixes:

  double f1os4p = 1.0/sqrt(4.0*M_PI);
  double Bxbar = Bx[index]*f1os4p,Bybar = By[index]*f1os4p,Bzbar = Bz[index]*f1os4p;

  double gxx_physL=gxxL*Psi4,gyy_physL=gyyL*Psi4,gzz_physL=gzzL*Psi4,gxy_physL=gxyL*Psi4,gxz_physL=gxzL*Psi4,gyz_physL=gyzL*Psi4;
  double Bbar_x = gxx_physL*Bxbar + gxy_physL*Bybar + gxz_physL*Bzbar;
  double Bbar_y = gxy_physL*Bxbar + gyy_physL*Bybar + gyz_physL*Bzbar;
  double Bbar_z = gxz_physL*Bxbar + gyz_physL*Bybar + gzz_physL*Bzbar;
  double Bbar2 = Bxbar*Bbar_x + Bybar*Bbar_y + Bzbar*Bbar_z;
  double Bbar = sqrt(Bbar2);
  double check_B_small = fabs(Bxbar)+fabs(Bybar)+fabs(Bzbar);
  if (check_B_small>0 && check_B_small<1.e-150) {
     // need to compute Bbar specially to prevent floating-point underflow
     double Bmax = fabs(Bxbar);
     if (Bmax < fabs(Bybar)) Bmax=fabs(Bybar);
     if (Bmax < fabs(Bzbar)) Bmax=fabs(Bzbar);
     double Bxtmp=Bxbar/Bmax, Bytemp=Bybar/Bmax, Bztemp=Bzbar/Bmax;
     double B_xtemp=Bbar_x/Bmax, B_ytemp=Bbar_y/Bmax, B_ztemp=Bbar_z/Bmax;
     Bbar = sqrt(Bxtmp*B_xtemp + Bytemp*B_ytemp + Bztemp*B_ztemp)*Bmax;
  }
  double stxi=mhd_st_x[index],styi=mhd_st_y[index],stzi=mhd_st_z[index];
  double BbardotS = Bxbar*stxi + Bybar*styi + Bzbar*stzi;
  double hatBbardotS = BbardotS/Bbar;
  if (Bbar<1.e-300) hatBbardotS = 0.0;

  double rho_s = rho_star[index];

  double sdots= gxxi*SQR(stxi)+gyyi*SQR(styi)+gzzi*SQR(stzi)+2.0*
      (gxyi*stxi*styi+gxzi*stxi*stzi+gyzi*styi*stzi);

  double Wm = sqrt(SQR(hatBbardotS)+ SQR(rho_s))/Psi6;
  double Sm2 = (SQR(Wm)*sdots + SQR(BbardotS)*(Bbar2+2.0*Wm))/SQR(Wm+Bbar2);
  double Wmin = sqrt(Sm2 + SQR(rho_s))/Psi6;
  double sdots_fluid_max = sdots;

  //tau fix, applicable when B==0 and B!=0:
  if(tau[index] < 0.5*Psi6*Bbar2) {
    tau[index] = tau_atm+0.5*Psi6*Bbar2;
  }

  double tau_fluid_min = tau[index] - 0.5*Psi6*Bbar2 - (Bbar2*sdots - SQR(BbardotS))*0.5/(Psi6*SQR(Wmin+Bbar2));
  
  //Apply Stilde fix when B==0.
  //if(Bx[index]==0 && By[index]==0 && Bz[index]==0 && (Psi6>30.0 || rho_star[index]/Psi6<100*rhobatm)) {
  //if(check_B_small < 1.e-300) {
  double Patm;
  compute_pcold_cpp(rhobatm, Patm, neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab);
 
  
  if(check_B_small*check_B_small < Patm*1e-32) {
    double rhot=tau[index]*(tau[index]+2.0*rho_s);
    double safetyfactor = 0.999999;
    //if(Psi6>Psi6threshold) safetyfactor=0.99;

    if(sdots > safetyfactor*rhot) {
      double rfactm1 = sqrt((safetyfactor*rhot)/sdots);
      mhd_st_x[index]=mhd_st_x[index]*rfactm1;
      mhd_st_y[index]=mhd_st_y[index]*rfactm1;
      mhd_st_z[index]=mhd_st_z[index]*rfactm1;
    }
   } else if(Psi6>Psi6threshold) {
     //Apply new Stilde fix.
     if (tau_fluid_min < tau_atm*1.001) { 
        tau_fluid_min = tau_atm*1.001;
        tau[index] = tau_fluid_min + 0.5*Psi6*Bbar2 + (Bbar2*sdots - SQR(BbardotS))*0.5/(Psi6*SQR(Wmin+Bbar2));
     }
     double tauL = tau[index];
     double rho_starL = rho_star[index];

     double LHS = tau_fluid_min*(tau_fluid_min+2.0*rho_starL);
     double RHS = sdots_fluid_max;

     double safetyfactor = 0.999999;
     if(safetyfactor*LHS < RHS) {
       double rfactm1 = sqrt((safetyfactor*LHS)/RHS);
       mhd_st_x[index]=mhd_st_x[index]*rfactm1;
       mhd_st_y[index]=mhd_st_y[index]*rfactm1;
       mhd_st_z[index]=mhd_st_z[index]*rfactm1;
     }
   }
  


  return 0;
}

void recompute_conservatives_fast(struct output_primitives &new_primitives,
				  int &enable_OS_collapse, int &compute_microphysics,
				  double &Psi6,double &alpha,
				  double &gxx_phys,double &gxy_phys,double &gxz_phys,double &gyy_phys,double &gyz_phys,double &gzz_phys,
				  double &gupxx_phys,double &gupxy_phys,double &gupxz_phys,double &gupyy_phys,double &gupyz_phys,double &gupzz_phys,
				  double &Bx,double &By,double &Bz, double &T_fluid, double &eps_tot, double &eps_cld, double &P_cld,
				  double &rho_star,double &tau,double &mhd_st_x,double &mhd_st_y,double &mhd_st_z,
				  double &u_xl,double &u_yl, double &u_zl,double &au0m1, 
				  int &neos, int &ergo_star, double &ergo_sigma, double *rho_tab,double *P_tab,double *eps_tab,
				  double *k_tab,double *gamma_tab, double &rho_max, double &rho_b_atm) {
  
  double rho_b_new, P_new;
  double P_cold, eps_cold;
  double eps,T_fluidl;

  rho_b_new = new_primitives.rho_b_new;  
  compute_pcold_cpp(rho_b_new, P_cold, neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab);
  compute_epscold_cpp(rho_b_new, P_cold, eps_cold, neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab, enable_OS_collapse);


  if (rho_b_new < 1.001*rho_b_atm){
    P_new = P_cold;
    eps = eps_cold;
  }
  else{
    P_new =  new_primitives.P_new;
 
    double P_min = 0.9*P_cold;
    if (P_new < P_min) {
      P_new = P_min;
    }
    double P_max = 1.0e3*P_cold;
    if (P_new > P_max) {
      P_new = P_max;
    }


    if (enable_OS_collapse == 1){
      eps   = 0.0;
    }
    else {  
      if (compute_microphysics==1)
	{
	  compute_eps_T_microphys(eps, T_fluidl, eps_cold, P_new, P_cold, rho_b_new);
	}
      else{
	// The line below is for diagnostic Temperature only! 
	//	compute_eps_T_microphys(eps, T_fluidl, eps_cold, P_new, P_cold, rho_b_new);
	eps   = (P_new - P_cold)/(GAMMA_th-1.0)/rho_b_new + eps_cold;
      }
    }
  }
  eps_tot = eps;
  eps_cld = eps_cold;
  P_cld = P_cold;
  T_fluid = T_fluidl;
  // in compute_eps_T_microphys, P is possibly changed. Here reassign the value.
  new_primitives.P_new = P_new;


  double h_new = 1.0 + P_new/rho_b_new + eps;

  double u0_new     = new_primitives.u0_new;

  double alpn1     = alpha;
  double alpn1_inv = 1.0/alpha;
  double f1o4p     = 1.0/(4.0*M_PI);
  double f1o4pa    = sqrt(f1o4p)*alpn1_inv;

  double Bx_f1o4pa = Bx*f1o4pa;
  double By_f1o4pa = By*f1o4pa;
  double Bz_f1o4pa = Bz*f1o4pa;
  
  double B_xl  = (gxx_phys *Bx + gxy_phys * By +
		  gxz_phys * Bz );
  double B_yl  = (gxy_phys *Bx + gyy_phys * By +
		  gyz_phys * Bz);
  double B_zl  = (gxz_phys *Bx + gyz_phys * By +
		  gzz_phys * Bz);
  double B_xl_f1o4pa = B_xl*f1o4pa;
  double B_yl_f1o4pa = B_yl*f1o4pa;
  double B_zl_f1o4pa = B_zl*f1o4pa;

  double B2 = ( gxx_phys*SQR(Bx_f1o4pa) +
		2.0*gxy_phys*Bx_f1o4pa*By_f1o4pa + 2.0*gxz_phys*Bx_f1o4pa*Bz_f1o4pa +
		gyy_phys*SQR(By_f1o4pa) + 2.0*gyz_phys*By_f1o4pa*Bz_f1o4pa +
		gzz_phys*SQR(Bz_f1o4pa) );
	
  double sb0 = u_xl*Bx_f1o4pa + u_yl*By_f1o4pa + u_zl*Bz_f1o4pa;
  double sb2 = (B2 + SQR(sb0))/SQR(u0_new);
  double sb_x = (B_xl_f1o4pa + u_xl*sb0)/u0_new;
  double sb_y = (B_yl_f1o4pa + u_yl*sb0)/u0_new;
  double sb_z = (B_zl_f1o4pa + u_zl*sb0)/u0_new;

  /*
    double tau_old = tau;
    double mhd_st_x_old = mhd_st_x;
  */

  rho_star = alpn1*u0_new*rho_b_new*Psi6;
  mhd_st_x = rho_star*h_new*u_xl +
    alpn1*Psi6*u0_new*sb2*u_xl - alpn1*Psi6*sb0*sb_x;
  mhd_st_y = rho_star*h_new*u_yl +
    alpn1*Psi6*u0_new*sb2*u_yl - alpn1*Psi6*sb0*sb_y;
  mhd_st_z = rho_star*h_new*u_zl +
    alpn1*Psi6*u0_new*sb2*u_zl - alpn1*Psi6*sb0*sb_z;
  tau = (au0m1+(P_new/rho_b_new+eps)*alpn1*u0_new)*rho_star +
    Psi6*sb2*SQR(alpn1*u0_new)
    - Psi6*(P_new+sb2*0.5)-Psi6*SQR(alpn1*sb0);
}




void compute_pcold_cpp(double &rhob, double &P_cold,
                       int &neos, int &ergo_star, double &ergo_sigma, double *rho_tab,double *P_tab,double *eps_tab,
		       double *k_tab,double *gamma_tab) {
  bool exit_do;
  int i = 0;
  exit_do = 0;
  while(exit_do==0) {
    if (rhob <= rho_tab[i]) {
      exit_do = 1;
      P_cold = k_tab[i]*pow(rhob,gamma_tab[i]);
    }
    if (i==neos-1) exit_do=1;
    i++;
  }
  if (rhob > rho_tab[neos-1]) {
    if (ergo_star==0){
      P_cold = k_tab[neos]*pow(rhob,gamma_tab[neos]);
    }
    else {
      P_cold = ((ergo_sigma* (1.0+eps_tab[neos-1]+P_tab[neos-1]/rho_tab[neos-1])/pow(rho_tab[neos-1],ergo_sigma)) * pow(rhob, ergo_sigma+1.0) + P_tab[neos-1] - ergo_sigma*((1.0+eps_tab[neos-1])*rho_tab[neos-1]))/(ergo_sigma+1.0);
    }
  }
}



void compute_epscold_cpp(double &rhob, double &P_cold, double &eps_cold,
                       int &neos, int &ergo_star, double &ergo_sigma, double *rho_tab,double *P_tab,double *eps_tab,
			 double *k_tab,double *gamma_tab, int &enable_OS_collapse){
  int i = 0;
  bool exit_do;

  if(enable_OS_collapse == 1){
    eps_cold = P_cold/rhob/(gamma_tab[i]-1.0);
  }
  else {
  exit_do = 0;
  while(exit_do==0) {
    if (rhob <= rho_tab[i]){
      exit_do = 1;
      if (i==0){
	if (rhob != 0.0){
	  eps_cold = P_cold/rhob/(gamma_tab[i]-1.0);
	} else {
	  eps_cold = 0.0;
	}
      } else {
	eps_cold = eps_tab[i-1] + (P_cold/rhob - P_tab[i-1]/rho_tab[i-1])/(gamma_tab[i]-1.0);
      }
    }
    if (i==neos-1) exit_do=1;
    i++;
  }


  if (rhob > rho_tab[neos-1]) {
    if (ergo_star == 0){
      eps_cold = eps_tab[neos-1] + (P_cold/rhob - P_tab[neos-1]/rho_tab[neos-1])/(gamma_tab[neos]-1.0);
    }
    else {
      eps_cold = ((1.0+eps_tab[neos-1]+P_tab[neos-1]/rho_tab[neos-1])/pow(rho_tab[neos-1],ergo_sigma) * pow(rhob, ergo_sigma+1.0) - P_tab[neos-1] + ergo_sigma*(1.0+eps_tab[neos-1])*rho_tab[neos-1] )/((ergo_sigma+1.0)*rhob)-1.0;
    }
  }
  
  }


}
