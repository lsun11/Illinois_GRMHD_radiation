// Lunan: This is very similar to the original conservative and vars calculator
// except that for the radiation part we take radiation primitive variables as INPUT!!!!
// NOT output like in the original script. This is called only in mhd_postinitialdata.

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sys/time.h>
#include <math.h>
#include "cctk.h"
#include "harm_u2p_defs.h"
#include "harm_primitives_headers.h"
#define SQR(x) ((x) * (x))


void compute_pcold_epscold_cpp(double &rhob, double &P_cold, double &eps_cold,
                               int &neos, int &ergo_star, double &ergo_sigma, double *rho_tab,double *P_tab,double *eps_tab,double *k_tab,double *gamma_tab, int &enable_OS_collapse);

void compute_M1(double &Pij,double &Fi,double &Fj,double &Fasq,double &E,double &gupij, double &shifti, double &shiftj, double &lapse, double &ui, double &uj, double &chi, double &psim4, double &Erad_atm_cut);


extern "C" void CCTK_FCALL CCTK_FNAME(metric_source_terms_and_misc_vars_initial)
  (const cGH **cctkGH,int *ext,
   double *rho,double *Sx,double *Sy,double *Sz,
   double *Sxx,double *Sxy,double *Sxz,double *Syy,double *Syz,double *Szz,
   double *tau_rad, double *S_rad_x, double *S_rad_y, double *S_rad_z,
   double *E_rad,  double *F_radx, double *F_rady,double *F_radz, double *F_rad0, double *F,
   double *h,double *w,
   double *st_x,double *st_y,double *st_z,
   double *Ex,double *Ey,double *Ez,
   double *sbt,double *sbx,double *sby,double *sbz,
   double *lapm1,double *shiftx,double *shifty,double *shiftz,
   double *phi,double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
   double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,
   double *P,double *rho_b,double *u0,double *vx,double *vy,double *vz,
   double *Bx,double *By,double *Bz,
   double *rho_star,double *mhd_st_x,double *mhd_st_y,double *mhd_st_z,
   int &neos,int &ergo_star, double &ergo_sigma,
   double *rho_tab,double *P_tab,double *eps_tab,double *k_tab,double *gamma_tab,
   double &rho_b_atm, int &enable_OS_collapse, int &rad_evolve_enable, int &rad_closure_scheme, double &Psi6threshold, double &Erad_atm_cut);

void metric_source_terms_and_misc_vars_initial(const cGH *cctkGH,int *ext,
					       double *rho,double *Sx,double *Sy,double *Sz,
					       double *Sxx,double *Sxy,double *Sxz,double *Syy,double *Syz,double *Szz,
					       double *tau_rad, double *S_rad_x, double *S_rad_y, double *S_rad_z,
					       double *E_rad,  double *F_radx, double *F_rady,double *F_radz, double *F_rad0, double *F,
					       double *h,double *w,
					       double *st_x,double *st_y,double *st_z,
					       double *Ex,double *Ey,double *Ez,
					       double *sbt,double *sbx,double *sby,double *sbz,
					       double *lapm1,double *shiftx,double *shifty,double *shiftz,
					       double *phi,double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
					       double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,
					       double *P,double *rho_b,double *u0,double *vx,double *vy,double *vz,
					       double *Bx,double *By,double *Bz,
					       double *rho_star,double *mhd_st_x,double *mhd_st_y,double *mhd_st_z,
					       int &neos,int &ergo_star, double &ergo_sigma,
					       double *rho_tab,double *P_tab,double *eps_tab,double *k_tab,double *gamma_tab, 
					       double &rho_b_atm, int &enable_OS_collapse, int &rad_evolve_enable, int &rad_closure_scheme, double &Psi6threshold, double &Erad_atm_cut) {
printf ("start metric_source_terms_and_misc_vars_initial!!!!!! \n");
#pragma omp parallel for
  for(int k=0;k<ext[2];k++)
    for(int j=0;j<ext[1];j++)
      for(int i=0;i<ext[0];i++) {
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

	

        if(isnan(mhd_st_x[index])) { printf("Inside metric source, found nan in msx %d %d %d\n",i,j,k); }
        if(isnan(mhd_st_y[index])) { printf("Inside metric source, found nan in msy %d %d %d\n",i,j,k); }
        if(isnan(mhd_st_z[index])) { printf("Inside metric source, found nan in msz %d %d %d\n",i,j,k); }

	//	if(h[index]==0.0) {printf("Inside metric source, h[index] is 0!!!! %d %d %d\n",i,j,k);}
       
	//	if(i==0&&j==0) { printf("Inside metric source, i=j=0, k=%d\n", k);
	// printf("Inside metric source, i,j, k, h[index] is %d, %d, %d, %e \n",i,j,k,h[index]);}

        double f1o8p = 1.0/(8.0*M_PI);
        double f1o4p = 2.0*f1o8p;

	double rho_starL = rho_star[index];

	double alpn1 = lapm1[index] + 1.0;
        double alpn1_inv = 1.0/alpn1;

	double au0m1 = alpn1*u0[index]-1.0;

	double BxL = Bx[index];
        double ByL = By[index];
        double BzL = Bz[index];


	double gxxL = gxx[index];                                                                                              
	double gxyL = gxy[index];                                                                                           
	double gxzL = gxz[index];                                                     
	double gyyL = gyy[index];                   
	double gyzL = gyz[index];                                                   
        double gzzL = gzz[index];                                                                                                
	
	double shiftxL = shiftx[index];                                          
	double shiftyL = shifty[index];                              
        double shiftzL = shiftz[index];

	double Psi2 = exp(2.0*phi[index]);
	double Psi4 = Psi2*Psi2;
        double Psim4 = 1.0/Psi4;
        double Psi6 = Psi4*Psi2;
        double Psim6 = 1.0/Psi6;

	double u_xl, u_yl, u_zl;
	double u_xll,u_yll,u_zll;
	double u0L, PL, rho_bL;
	double ux_l, uy_l, uz_l;

	if (rho_starL < 0.0){
	  rho_bL = rho_b_atm;
	  printf ("negative rho_star found!!!");
	  double P_cold, eps_cold;
	  compute_pcold_epscold_cpp(rho_bL, P_cold, eps_cold,
				    neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab,enable_OS_collapse);
	  PL = P_cold;
	  double eps= eps_cold;
	  double h_l = 1.0 + PL/rho_bL + eps;
	  
	  u0[index] = 1.0*alpn1_inv;
	  vx[index] = -shiftxL;
	  vy[index] = -shiftyL;
	  vz[index] = -shiftzL;
	  

	  ux_l = vx[index]*u0[index];
	  uy_l = vy[index]*u0[index];
	  uz_l = vz[index]*u0[index];

	  double B2 = Psi4*( gxxL*SQR(BxL) +
                      2.0*gxyL*BxL*ByL +
                      2.0*gxzL*BxL*BzL +
                      gyyL*SQR(ByL) +
                      2.0*gyzL*ByL*BzL +
                      gzzL*SQR(BzL) );
	  double xsb2 = B2*f1o4p;
	  rho_star[index] = Psi6*rho_bL;
	  st_x[index]  = 0.0;
	  st_y[index]  = 0.0;
	  st_z[index]  = 0.0;
	  mhd_st_x[index] = 0.0;
	  mhd_st_y[index] = 0.0;
	  mhd_st_z[index] = 0.0;
	  u_xll   = 0.0;
	  u_yll   = 0.0;
	  u_zll   = 0.0;
	  h[index]     = h_l;
	  w[index]     = rho_star[index];
	  rho_b[index] = rho_bL;
	  P[index]     = PL;
	  rho[index]   = rho_bL*(1.0+eps);
	  Sx[index]    = 0.0;
	  Sy[index]    = 0.0;
	  Sz[index]    = 0.0;
	  Sxx[index]   = Psi4 * gxxL * PL;
	  Sxy[index]   = Psi4 * gxyL * PL;
	  Sxz[index]   = Psi4 * gxzL * PL;
	  Syy[index]   = Psi4 * gyyL * PL;
	  Syz[index]   = Psi4 * gyzL * PL;
	  Szz[index]   = Psi4 * gzzL * PL;	  
	}
	else{
	u0L = u0[index];
	PL = P[index];
	rho_bL = rho_b[index];
	//	double au0m1 = alpn1*u0L-1.0;

	u_xll = (gxxL*(shiftxL+vx[index]) + 
			gxyL*(shiftyL+vy[index]) + 
			gxzL*(shiftzL+vz[index]))*Psi4*u0L;
	u_yll = (gxyL*(shiftxL+vx[index]) + 
			gyyL*(shiftyL+vy[index]) + 
			gyzL*(shiftzL+vz[index]))*Psi4*u0L;
        u_zll = (gxzL*(shiftxL+vx[index]) + 
			gyzL*(shiftyL+vy[index]) + 
			gzzL*(shiftzL+vz[index]))*Psi4*u0L;

	u_xl = u_xll;
	u_yl = u_yll;
	u_zl = u_zll;

	w[index] = (1.0+lapm1[index])*u0L*rho_starL;
	double st_x_l = rho_starL*h[index]*u_xl;
	double st_y_l = rho_starL*h[index]*u_yl;
	double st_z_l = rho_starL*h[index]*u_zl;

	st_x[index] = st_x_l;
	st_y[index] = st_y_l;
	st_z[index] = st_z_l;

	ux_l = vx[index]*u0L;
	uy_l = vy[index]*u0L;
	uz_l = vz[index]*u0L;

	double fac;

	if(w[index]==0.0 || h[index] == 0.0){
	  fac = 0.0; 
	}
	else{
	  fac = 1.0 / ( Psi6 * w[index]* h[index] ) ;
	}

	rho[index] = h[index] * w[index] * Psim6 - PL;
	Sx[index]  = st_x_l * Psim6;
	Sy[index]  = st_y_l * Psim6; 
	Sz[index]  = st_z_l * Psim6;
	Sxx[index] = fac * st_x_l*st_x_l + Psi4 * gxxL * PL;
	Sxy[index] = fac * st_x_l*st_y_l + Psi4 * gxyL * PL;
	Sxz[index] = fac * st_x_l*st_z_l + Psi4 * gxzL * PL;
	Syy[index] = fac * st_y_l*st_y_l + Psi4 * gyyL * PL;
	Syz[index] = fac * st_y_l*st_z_l + Psi4 * gyzL * PL;
	Szz[index] = fac * st_z_l*st_z_l + Psi4 * gzzL * PL;

	
	if(isnan(Sxx[index])) { printf("BAD sxx BEFORE MHD: %d %d %d %e.  %e %e %e %e %e %e %e %e %e %e\n",
				       i,j,k,Sxx[index],Psi4,w[index],gxxL,gxyL,gyyL,rho_starL,h[index],u_xl,PL,rho_bL);
          printf("%e,%e,%e,%e,%e\n", st_x_l, st_y_l, st_z_l, rho_star[index], u_xl);
	}
	

	}


	// MHD metric sources

	double E_xl = Psi6 * ( ByL*(vz[index]+shiftzL) 
			       - BzL*(vy[index]+shiftyL) )*alpn1_inv;
	double E_yl = Psi6 * ( BzL*(vx[index]+shiftxL) 
			       - BxL*(vz[index]+shiftzL) )*alpn1_inv;
	double E_zl = Psi6 * ( BxL*(vy[index]+shiftyL) 
			       - ByL*(vx[index]+shiftxL) )*alpn1_inv;
	Ex[index] = (gupxx[index]*E_xl + gupxy[index]*E_yl + 
		     gupxz[index]*E_zl)*Psim4;
	Ey[index] = (gupxy[index]*E_xl + gupyy[index]*E_yl + 
		     gupyz[index]*E_zl)*Psim4;
	Ez[index] = (gupxz[index]*E_xl + gupyz[index]*E_yl + 
		     gupzz[index]*E_zl)*Psim4;
	double B_xl  = Psi4 * (gxxL * BxL + gxyL * ByL 
			       + gxzL * BzL);
	double B_yl  = Psi4 * (gxyL * BxL + gyyL * ByL 
			       + gyzL * BzL);
	double B_zl  = Psi4 * (gxzL * BxL + gyzL * ByL 
			       + gzzL * BzL);
        double temp = f1o8p*(Ex[index]*E_xl + Ey[index]*E_yl + Ez[index]*E_zl 
			     + BxL*B_xl + ByL*B_yl + BzL*B_zl);
	rho[index]   = rho[index] + temp;
	Sxx[index]   = Sxx[index] + temp*Psi4*gxxL 
	  - f1o4p*(E_xl*E_xl + B_xl*B_xl);
	Sxy[index]   = Sxy[index] + temp*Psi4*gxyL 
	  - f1o4p*(E_xl*E_yl + B_xl*B_yl);
	Sxz[index]   = Sxz[index] + temp*Psi4*gxzL 
	  - f1o4p*(E_xl*E_zl + B_xl*B_zl);
	Syy[index]   = Syy[index] + temp*Psi4*gyyL 
	  - f1o4p*(E_yl*E_yl + B_yl*B_yl);
	Syz[index]   = Syz[index] + temp*Psi4*gyzL 
	  - f1o4p*(E_yl*E_zl + B_yl*B_zl);
	Szz[index]   = Szz[index] + temp*Psi4*gzzL 
	  - f1o4p*(E_zl*E_zl + B_zl*B_zl);
	Sx[index]    = Sx[index] + f1o4p*Psi6*(Ey[index]*BzL - Ez[index]*ByL);
	Sy[index]    = Sy[index] + f1o4p*Psi6*(Ez[index]*BxL - Ex[index]*BzL);
	Sz[index]    = Sz[index] + f1o4p*Psi6*(Ex[index]*ByL - Ey[index]*BxL);
	sbt[index] = u_xll*BxL + u_yll*ByL + u_zll*BzL;
	sbx[index] = BxL/u0L + vx[index]*sbt[index];
	sby[index] = ByL/u0L + vy[index]*sbt[index];
	sbz[index] = BzL/u0L + vz[index]*sbt[index];


	if(isnan(Sxx[index])) { printf("BAD sxx BEFORE rad: %d %d %d %e.  %e %e %e %e %e %e %e %e %e %e\n",
				      i,j,k,Sxx[index],Psi6,B_xl,BzL,Ez[index],ByL,rho_starL,h[index],u_xl,PL,rho_bL); 
	  printf("%e,%e,%e,%e,%e\n", E_xl, vx[index], temp, alpn1_inv, rho_b_atm);
	}
      
	//Add radiation terms to the source
	if (rad_evolve_enable == 1){
	  
	  double E_radl = E_rad[index];
	  double F_radxl = F_radx[index];
	  double F_radyl = F_rady[index];
	  double F_radzl = F_radz[index];
	  
	  double shift_xl = Psi4* (shiftx[index]*gxxL + shifty[index]*gxyL +shiftz[index]*gxzL);
	  double shift_yl = Psi4* (shiftx[index]*gxyL + shifty[index]*gyyL +shiftz[index]*gyzL);
	  double shift_zl = Psi4* (shiftx[index]*gxzL + shifty[index]*gyzL +shiftz[index]*gzzL);
	  
	  double v_xl = Psi4* (vx[index]*gxxL + vy[index]*gxyL +vz[index]*gxzL);
	  double v_yl = Psi4* (vx[index]*gxyL + vy[index]*gyyL +vz[index]*gyzL);
	  double v_zl = Psi4* (vx[index]*gxzL + vy[index]*gyzL +vz[index]*gzzL);
	  
	  double beta2 = shiftx[index]*shift_xl + shifty[index]*shift_yl + shiftz[index]*shift_zl;
	  double udotbeta = u0L*(vx[index]*shift_xl + vy[index]*shift_yl + vz[index]*shift_zl);
	  double g_00L =beta2-alpn1*alpn1;
	  double u_0L = g_00L*u0L + udotbeta;
	  
	  //Apply atm floor of rad primitives:                                                                                                                                      
	  if (E_radl < Erad_atm_cut)
	    {
	      E_radl = Erad_atm_cut;
	      F_radxl = 0.0;
	      F_radyl = 0.0;
	      F_radzl = 0.0;
	    }
	  
	  //Apply Psi6threshold fix on E_radl, and F_radl:                                                                                                                          
	  double F_rad0l;
	  if(Psi6 > Psi6threshold) {
	    double Erad_horiz_cap = 1.0e2*Erad_atm_cut;
	    if(E_radl < Erad_horiz_cap){
	      E_radl = Erad_horiz_cap;
	      F_radxl = 0.0;
	      F_radyl = 0.0;
	      F_radzl = 0.0;
	      F_rad0l = 0.0;
	    }
	  }
	  
	  F_rad0l = - (F_radxl*u_xl + F_radyl*u_yl + F_radzl*u_zl)/u_0L;
	  
	  double F_rad_xl = Psi4 * (gxxL * F_radxl + gxyL * F_radyl + gxzL * F_radzl) + shift_xl*F_rad0l;
	  double F_rad_yl = Psi4 * (gxyL * F_radxl + gyyL * F_radyl + gyzL * F_radzl) + shift_yl*F_rad0l;
	  double F_rad_zl = Psi4 * (gxzL * F_radxl + gyzL * F_radyl + gzzL * F_radzl) + shift_zl*F_rad0l;
	  double F_rad_0l = - (F_rad_xl*ux_l + F_rad_yl*uy_l + F_rad_zl*uz_l)/u0L;
	  
	  double Fasq = F_rad_0l*F_rad0l + F_rad_xl*F_radxl +  F_rad_yl*F_radyl +  F_rad_zl*F_radzl;
	  
	  F[index] = sqrt(Psi4*(gxxL*SQR(F_radxl)+gyyL*SQR(F_radyl)+gzzL*SQR(F_radzl) + 2.0*( gxyL*F_radxl*F_radyl + gxzL*F_radxl*F_radzl + gyzL*F_radyl*F_radzl))  );
	  if (rad_closure_scheme==0){
	    double P_radl = E_radl/3.0;
	    
	    double temp_rad = alpn1*u0[index];
	    double temp_rad1 = temp_rad*temp_rad*(E_radl + P_radl) - P_radl + 2.0*alpn1*u0[index]*F_rad0l;
	    
	    //F_rad_\alpha = g_{0\alpha} F_rad0 + g_{z\alpha} F_radz + g_{y\alpha} F_rady + g_{z\alpha} F_radz  
	    
	    rho[index]   = rho[index] + temp_rad1;
	    
	    Sx[index] = Sx[index] + temp_rad*( ( (E_radl+P_radl)*u0L + F_rad0l) * (shift_xl + v_xl) + F_rad_xl);
	    Sy[index] = Sy[index] + temp_rad*( ( (E_radl+P_radl)*u0L + F_rad0l) * (shift_yl + v_yl) + F_rad_yl);
	    Sz[index] = Sz[index] + temp_rad*( ( (E_radl+P_radl)*u0L + F_rad0l) * (shift_zl + v_zl) + F_rad_zl);
	    
	    Sxx[index] = Sxx[index] + (E_radl+P_radl)*SQR(u0L*(shift_xl+v_xl)) + 2.0*F_rad_xl*u0L*(shift_xl+v_xl) + Psi4 * P_radl * gxxL;
	    Syy[index] = Syy[index] + (E_radl+P_radl)*SQR(u0L*(shift_yl+v_yl)) + 2.0*F_rad_yl*u0L*(shift_yl+v_yl) + Psi4 * P_radl * gyyL;
	    Szz[index] = Szz[index] + (E_radl+P_radl)*SQR(u0L*(shift_zl+v_zl)) + 2.0*F_rad_yl*u0L*(shift_zl+v_zl) + Psi4 * P_radl * gzzL;
	    Sxy[index] = Sxy[index] + (E_radl+P_radl)*SQR(u0L)*(shift_xl+v_xl)*(shift_yl+v_yl) + u0L*(F_rad_xl*(shift_yl+v_yl)+F_rad_yl*(shift_xl+v_xl)) + Psi4 * P_radl * gxyL;
	    Sxz[index] = Sxz[index] + (E_radl+P_radl)*SQR(u0L)*(shift_xl+v_xl)*(shift_zl+v_zl) + u0L*(F_rad_xl*(shift_zl+v_zl)+F_rad_zl*(shift_xl+v_xl)) + Psi4 * P_radl * gxzL;
	    Syz[index] = Syz[index] + (E_radl+P_radl)*SQR(u0L)*(shift_yl+v_yl)*(shift_zl+v_zl) + u0L*(F_rad_yl*(shift_zl+v_zl)+F_rad_zl*(shift_yl+v_yl)) + Psi4 * P_radl * gyzL;
	    
	    //assign diag stuff:                       
	    F_rad0[index] = F_rad0l;
	    
	    // recompute_conservatives:                                                                                                                                               
	    tau_rad[index] = SQR(alpn1)*Psi6*(E_radl*SQR(u0[index]) + 2.0*F_rad0l*u0[index] + P_radl*SQR(u0[index])) - Psi6*P_radl;
	    S_rad_x[index] = alpn1*Psi6*((E_radl+P_radl)*u0[index]*u_xll + F_rad0l*u_xll + F_rad_xl*u0[index]);
	    S_rad_x[index] = alpn1*Psi6*((E_radl+P_radl)*u0[index]*u_yll + F_rad0l*u_yll + F_rad_yl*u0[index]);
	    S_rad_z[index] = alpn1*Psi6*((E_radl+P_radl)*u0[index]*u_zll + F_rad0l*u_zll + F_rad_zl*u0[index]);
	  }
	  else { // M1 closure!!!
	    double Fksq = F_rad_xl*F_radxl +  F_rad_yl*F_radyl +  F_rad_zl*F_radzl;
	    double zeta;
	    double zeta_temp = sqrt(fabs(Fasq)/SQR(E_radl));
	    //	    double zeta_cut = Erad_atm_cut*1.5;
	    double zeta_cut = 1.0e-40;

	    if (E_radl<=Erad_atm_cut){
	      zeta = 1.0;
	    }else{
	      zeta = zeta_temp;
	    }	
	    if (zeta > 1.0){
	      zeta = 1.0;
	    }
	    
	    double chi = 1/3.0 + SQR(zeta)*(6.0-2.0*zeta+6.0*SQR(zeta))/15.0;
	    double P_radxxl,P_radyyl,P_radzzl,P_radxyl,P_radxzl,P_radyzl;
	    compute_M1(P_radxxl, F_radxl, F_radxl, Fasq, E_radl, gupxx[index], shiftx[index], shiftx[index], lapm1[index], ux_l, ux_l, chi, Psim4, Erad_atm_cut);
	    compute_M1(P_radyyl, F_radyl, F_radyl, Fasq, E_radl, gupyy[index], shifty[index], shifty[index], lapm1[index], uy_l, uy_l, chi, Psim4, Erad_atm_cut);
	    compute_M1(P_radzzl, F_radzl, F_radzl, Fasq, E_radl, gupzz[index], shiftz[index], shiftz[index], lapm1[index], uz_l, uz_l, chi, Psim4, Erad_atm_cut);
	    compute_M1(P_radxyl, F_radxl, F_radyl, Fasq, E_radl, gupxy[index], shiftx[index], shifty[index], lapm1[index], ux_l, uy_l, chi, Psim4, Erad_atm_cut);
	    compute_M1(P_radxzl, F_radxl, F_radzl, Fasq, E_radl, gupxz[index], shiftx[index], shiftz[index], lapm1[index], ux_l, uz_l, chi, Psim4, Erad_atm_cut);
	    compute_M1(P_radyzl, F_radyl, F_radzl, Fasq, E_radl, gupyz[index], shifty[index], shiftz[index], lapm1[index], uy_l, uz_l, chi, Psim4, Erad_atm_cut);
	    
	    double P_rad0xl = - (P_radxxl * u_xll + P_radxyl * u_yll + P_radxzl * u_zll)/u_0L;
	    double P_rad0yl = - (P_radxyl * u_xll + P_radyyl * u_yll + P_radyzl * u_zll)/u_0L;
	    double P_rad0zl = - (P_radxzl * u_xll + P_radyzl * u_yll + P_radzzl * u_zll)/u_0L;
	    double P_rad00l = - (P_rad0xl * u_xll + P_rad0yl * u_yll + P_rad0zl * u_zll)/u_0L;
	    
	    //Asign radiation variables.
	    rho[index] = rho[index] + SQR(alpn1)*(E_radl*SQR(u0L) + 2.0*F_rad0l*u0L + P_rad00l);
	    
	    Sx[index] = Sx[index] + alpn1*(u0L*E_radl*u_xll + F_rad0l*u_xll + u0L*F_rad_xl +
					   P_rad00l*shift_xl + Psi4*(P_rad0xl*gxxL + P_rad0yl*gxyL + P_rad0zl*gxzL));
	    
	    Sy[index] = Sy[index] + alpn1*(u0L*E_radl*u_yll + F_rad0l*u_yll + u0L*F_rad_yl +
					   P_rad00l*shift_yl + Psi4*(P_rad0xl*gxyL + P_rad0yl*gyyL + P_rad0zl*gyzL));
	    
	    Sz[index] = Sz[index] + alpn1*(u0L*E_radl*u_zll + F_rad0l*u_zll + u0L*F_rad_zl +
					   P_rad00l*shift_zl + Psi4*(P_rad0xl*gxzL + P_rad0yl*gyzL + P_rad0zl*gzzL));
	    
	    
	    
	    Sxx[index] = Sxx[index] + E_radl*u_xll*u_xll + 2*F_rad_xl*u_xll +
	      SQR(shift_xl)*P_rad00l + shift_xl*2.0*(gxxL*P_rad0xl+gxyL*P_rad0yl+gxzL*P_rad0zl) +
	      SQR(Psi4)*( SQR(gxxL)*P_radxxl + SQR(gxyL)*P_radyyl + SQR(gxzL)*P_radzzl +
			  2.0*(gxxL*gxyL*P_radxyl+gxxL*gxzL*P_radxzl+gxyL*gxzL*P_radyzl) );
	    
	    Syy[index] = Syy[index] + E_radl*u_yll*u_yll + 2*F_rad_yl*u_yll +
	      SQR(shift_yl)*P_rad00l + shift_yl*2.0*(gxyL*P_rad0xl+gyyL*P_rad0yl+gyzL*P_rad0zl) +
	      SQR(Psi4)*( SQR(gxyL)*P_radxxl + SQR(gyyL)*P_radyyl + SQR(gyzL)*P_radzzl +
			  2.0*(gxyL*gyyL*P_radxyl+gxyL*gyzL*P_radxzl+gyyL*gyzL*P_radyzl) );
	    
	    Szz[index] = Szz[index] + E_radl*u_zll*u_zll + 2*F_rad_zl*u_zll +
	      SQR(shift_zl)*P_rad00l + shift_zl*2.0*(gxzL*P_rad0xl+gyzL*P_rad0yl+gzzL*P_rad0zl) +
	      SQR(Psi4)*( SQR(gxzL)*P_radxxl + SQR(gyzL)*P_radyyl + SQR(gzzL)*P_radzzl +
			  2.0*(gxzL*gyzL*P_radxyl+gxzL*gzzL*P_radxzl+gyzL*gzzL*P_radyzl) );
	    
	    Sxy[index] = Sxy[index] + E_radl*u_xll*u_yll + F_rad_xl*u_yll + F_rad_yl*u_xll +
	      shift_xl*shift_yl*P_rad00l +
	      Psi4*(shift_xl*(gxyL*P_rad0xl + gyyL*P_rad0yl + gyzL*P_rad0zl) + (shift_yl*(gxxL*P_rad0xl + gxyL*P_rad0yl + gxzL*P_rad0zl)) )+
	      SQR(Psi4)*(gxxL*gxyL*P_radxxl + gxyL*gyyL*P_radyyl + gxzL*gyzL*P_radzzl +
			 (gxxL*gyyL + gxyL*gxyL)*P_radxyl + (gxxL*gyzL + gxzL*gxyL)*P_radxzl + (gxyL*gyzL + gxzL*gyyL)*P_radyzl);
	    
	    Sxz[index] = Sxz[index] + E_radl*u_xll*u_zll + F_rad_xl*u_zll + F_rad_zl*u_xll +
	      shift_xl*shift_zl*P_rad00l +
	      Psi4*(shift_xl*(gxzL*P_rad0xl + gyzL*P_rad0yl + gzzL*P_rad0zl) + (shift_zl*(gxxL*P_rad0xl + gxyL*P_rad0yl + gxzL*P_rad0zl)) )+
	      SQR(Psi4)*(gxxL*gxzL*P_radxxl + gxyL*gyzL*P_radyyl + gxzL*gzzL*P_radzzl +
			 (gxxL*gyzL + gxyL*gxzL)*P_radxyl + (gxxL*gzzL + gxzL*gxzL)*P_radxzl + (gxyL*gzzL + gxzL*gyzL)*P_radyzl);
	    
	    Syz[index] = Syz[index] + E_radl*u_yll*u_zll + F_rad_yl*u_zll + F_rad_zl*u_yll +
	      shift_yl*shift_zl*P_rad00l +
	      Psi4*(shift_yl*(gxzL*P_rad0xl + gyzL*P_rad0yl + gzzL*P_rad0zl) + (shift_zl*(gxyL*P_rad0xl + gyyL*P_rad0yl + gyzL*P_rad0zl)) )+
	      SQR(Psi4)*(gxyL*gxzL*P_radxxl + gyyL*gyzL*P_radyyl + gyzL*gzzL*P_radzzl +
			 (gxyL*gyzL + gxzL*gyyL)*P_radxyl + (gxyL*gzzL + gyzL*gxzL)*P_radxzl + (gyyL*gzzL + gyzL*gyzL)*P_radyzl);
	    
	    if(isnan(Sxx[index])) { printf("BAD sxx in M1: %d %d %d %e.  %e %e %e %e %e %e %e %e %e %e\n",
					   i,j,k,Sxx[index],Psi6,E_radl,u_xll,F_rad_xl,P_radxxl,P_radxyl,P_radxzl,u_xl,PL,rho_bL);
	      printf("%e,%e,%e,%e,%e,%e,%e\n", zeta, chi, P_rad0yl, P_rad0zl, u_0L, Fasq, lapm1[index]);
	    }
	    F_rad0[index] = F_rad0l;
	    
	    tau_rad[index] = alpn1*alpn1*Psi6*(E_radl*u0L*u0L+2.0*F_rad0l*u0L+P_rad00l);
	    S_rad_x[index] = alpn1*Psi6*(E_radl*u0L*u_xll + F_rad0l*u_xll + F_rad_xl * u0L + P_rad00l*shift_xl + Psi4*(P_rad0xl*gxxL + P_rad0yl*gxyL + P_rad0zl*gxzL));
	    S_rad_y[index] = alpn1*Psi6*(E_radl*u0L*u_yll + F_rad0l*u_yll + F_rad_yl * u0L + P_rad00l*shift_yl + Psi4*(P_rad0xl*gxyL + P_rad0yl*gyyL + P_rad0zl*gyzL));
	    S_rad_z[index] = alpn1*Psi6*(E_radl*u0L*u_zll + F_rad0l*u_zll + F_rad_zl * u0L + P_rad00l*shift_zl + Psi4*(P_rad0xl*gxzL + P_rad0yl*gyzL + P_rad0zl*gzzL));
	  }
	}
      }
  printf ("End metric_source_terms_and_misc_vars_initial \n");
}


extern "C" void CCTK_FCALL CCTK_FNAME(metric_source_terms_and_misc_vars_initial)
  (const cGH **cctkGH,int *ext,
   double *rho,double *Sx,double *Sy,double *Sz,
   double *Sxx,double *Sxy,double *Sxz,double *Syy,double *Syz,double *Szz,
   double *tau_rad, double *S_rad_x, double *S_rad_y, double *S_rad_z,
   double *E_rad,  double *F_radx, double *F_rady,double *F_radz, double *F_rad0, double *F,
   double *h,double *w,
   double *st_x,double *st_y,double *st_z,
   double *Ex,double *Ey,double *Ez,
   double *sbt,double *sbx,double *sby,double *sbz,
   double *lapm1,double *shiftx,double *shifty,double *shiftz,
   double *phi,double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
   double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,
   double *P,double *rho_b,double *u0,double *vx,double *vy,double *vz,
   double *Bx,double *By,double *Bz,
   double *rho_star,double *mhd_st_x,double *mhd_st_y,double *mhd_st_z,
   int &neos,int &ergo_star, double &ergo_sigma,
   double *rho_tab,double *P_tab,double *eps_tab,double *k_tab,double *gamma_tab,
   double &rho_b_atm, int &enable_OS_collapse, int &rad_evolve_enable,  int &rad_closure_scheme, double &Psi6threshold, double &Erad_atm_cut) {
  metric_source_terms_and_misc_vars_initial (*cctkGH,ext,
				     rho,Sx,Sy,Sz,
				     Sxx,Sxy,Sxz,Syy,Syz,Szz,
				     tau_rad, S_rad_x,S_rad_y, S_rad_z,
				     E_rad, F_radx, F_rady, F_radz, F_rad0, F,
				     h,w,
				     st_x,st_y,st_z,
				     Ex,Ey,Ez,
				     sbt,sbx,sby,sbz,
				     lapm1,shiftx,shifty,shiftz,
				     phi,gxx,gxy,gxz,gyy,gyz,gzz,
				     gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,
				     P,rho_b,u0,vx,vy,vz,
				     Bx,By,Bz,
				     rho_star,mhd_st_x,mhd_st_y,mhd_st_z,
				     neos,ergo_star, ergo_sigma,
				     rho_tab, P_tab,eps_tab,k_tab,gamma_tab,
				     rho_b_atm, enable_OS_collapse, rad_evolve_enable, rad_closure_scheme, Psi6threshold, Erad_atm_cut);

}
