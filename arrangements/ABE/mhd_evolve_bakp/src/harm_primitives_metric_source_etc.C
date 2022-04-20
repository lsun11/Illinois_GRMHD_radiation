// Here is a wrapper for the 2d solver described in Noble et al.  
// It is meant to be an interface with their code.
// Note that it assumes a simple gamma law for the moment.  No hybrid.
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
                               int &neos, int &ergo_star, double &ergo_sigma, double *rho_tab,double *P_tab,double *eps_tab,double *k_tab,double *gamma_tab);

extern "C" void CCTK_FCALL CCTK_FNAME(metric_source_terms_and_misc_vars)
  (const cGH **cctkGH,int *ext,
   double *rho,double *Sx,double *Sy,double *Sz,
   double *Sxx,double *Sxy,double *Sxz,double *Syy,double *Syz,double *Szz,
   double *tau_rad, double *S_rad_x, double *S_rad_y, double *S_rad_z,
   double *E_rad, double *F_rad0, double *F_radx, double *F_rady,double *F_radz, double *P_rad,
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
   double &rho_b_atm);

void metric_source_terms_and_misc_vars(const cGH *cctkGH,int *ext,double *rho,double *Sx,double *Sy,double *Sz,
				       double *Sxx,double *Sxy,double *Sxz,double *Syy,double *Syz,double *Szz,
				       double *tau_rad, double *S_rad_x, double *S_rad_y, double *S_rad_z,
				       double *E_rad, double *F_rad0, double *F_radx, double *F_rady,double *F_radz, double *P_rad,
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
				       double &rho_b_atm) {
  //printf ("start metric_source_terms_and_misc_vars \n");
#pragma omp parallel for
  for(int k=0;k<ext[2];k++)
    for(int j=0;j<ext[1];j++)
      for(int i=0;i<ext[0];i++) {
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);


        if(isnan(mhd_st_x[index])) { printf("Inside metric source, found nan in msx %d %d %d\n",i,j,k); }
        if(isnan(mhd_st_y[index])) { printf("Inside metric source, found nan in msy %d %d %d\n",i,j,k); }
        if(isnan(mhd_st_z[index])) { printf("Inside metric source, found nan in msz %d %d %d\n",i,j,k); }

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
				    neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab);
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

	double fac   = 1.0 / ( Psi6 * w[index]* h[index] ) ;

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

	//	if (j==22 && k==1) {printf("%d %d %d %e %e %e\n", i,j,k, vx[index],vy[index],vz[index]); }

	if(isnan(Sx[index])) { printf("BAD sx: %d %d %d %e.  %e %e %e %e %e %e %e %e %e %e\n",
				      i,j,k,Sx[index],Psi6,Ey[index],BzL,Ez[index],ByL,rho_starL,h[index],u_xl,PL,rho_bL); 
	  printf("%e,%e,%e,%e,%e\n", E_yl, vx[index], shiftxL, alpn1_inv, rho_b_atm);
	  printf("%e,%e,%e,%e,%e\n", vy[index], vz[index], shiftx[index], shifty[index],shiftz[index]);
	} // exit(1); }
      


	//Add radiation terms to the source

	/*
	double tau_radl = tau_rad[index];
	double S_rad_xl = S_rad_x[index];
	double S_rad_yl = S_rad_y[index];
	double S_rad_zl = S_rad_z[index];

	double E_rad_denom = -(au0m1+1.0)*Psi6+ SQR(Psi6)*(2.0*(au0m1+1.0)/3.0-1.0/(au0m1+1.0)/6.0);
	double E_rad_num = (Psi6/(2.0*(au0m1+1.0)) - (au0m1+1.0)) * tau_radl
	  +(u0[index]*shiftxL+ux_l)*S_rad_xl
	  +(u0[index]*shiftyL+uy_l)*S_rad_yl
	  +(u0[index]*shiftzL+uz_l)*S_rad_zl;


	double E_radl = E_rad_num/E_rad_denom;
	double P_radl = E_radl/3.0;

	double F_rad0l = (tau_radl/Psi6- (4*SQR(au0m1+1.0) - 1.0)/3.0 * E_radl)/(2.0*(au0m1+1.0)*alpn1);

	double F_radxl = ((gupxx[index]*S_rad_xl +
			   gupxy[index]*S_rad_yl +
			   gupxz[index]*S_rad_zl)*Psim4)/(alpn1*Psi6*u0[index])-
	                  (E_radl+P_radl)*u0[index]*(vx[index]+shiftx[index]) -
	                   2.0*(F_rad0l * shiftx[index]) - F_rad0l*vx[index];

	double F_radyl = ((gupxy[index]*S_rad_xl +
			   gupyy[index]*S_rad_yl +
			   gupyz[index]*S_rad_zl)*Psim4)/(alpn1*Psi6*u0[index])-
	                 (E_radl+P_radl)*u0[index]*(vy[index]+shifty[index]) -
	                  2.0*(F_rad0l * shifty[index]) - F_rad0l*vy[index];

	double F_radzl = ((gupxz[index]*S_rad_xl +
			   gupyz[index]*S_rad_yl +
			   gupzz[index]*S_rad_zl)*Psim4)/(alpn1*Psi6*u0[index])-
	                  (E_radl+P_radl)*u0[index]*(vz[index]+shiftz[index]) -
	                   2.0*(F_rad0l * shiftz[index]) - F_rad0l*vz[index];
	E_rad[index] = E_radl;
	P_rad[index] = P_radl;
	F_rad0[index] = F_rad0l;
	F_radx[index] = F_radxl;
	F_rady[index] = F_radyl;
	F_radz[index] = F_radzl;

	double temp_rad = alpn1*u0[index];
	double temp_rad1 = temp_rad*temp_rad*(E_radl+P_radl) - P_radl + 2*alpn1*u0[index]*F_rad0l;


	double F_rad_xl = Psi4 * (gxxL * F_radxl + gxyL * F_radyl + gxzL * F_radzl);
	double F_rad_yl = Psi4 * (gxyL * F_radxl + gyyL * F_radyl + gyzL * F_radzl);
	double F_rad_zl = Psi4 * (gxzL * F_radxl + gyzL * F_radyl + gzzL * F_radzl);


	rho[index]   = rho[index] + temp_rad1;

	Sx[index] = Sx[index] + temp_rad*(E_radl + P_radl) * u_xll + alpn1*(F_rad0l*u_xll + u0L*F_rad_xl);
	Sy[index] = Sy[index] + temp_rad*(E_radl + P_radl) * u_yll + alpn1*(F_rad0l*u_yll + u0L*F_rad_yl);
	Sz[index] = Sz[index] + temp_rad*(E_radl + P_radl) * u_zll + alpn1*(F_rad0l*u_zll + u0L*F_rad_zl);

	Sxx[index] = Sxx[index] + (E_radl+P_radl)*u_xll*u_xll + Psi4 * P_radl * gxxL + 2*F_rad_xl*u_xll;
	Syy[index] = Syy[index] + (E_radl+P_radl)*u_yll*u_yll + Psi4 * P_radl * gyyL + 2*F_rad_yl*u_yll;
	Szz[index] = Szz[index] + (E_radl+P_radl)*u_zll*u_zll + Psi4 * P_radl * gzzL + 2*F_rad_zl*u_zll;
	Sxy[index] = Sxy[index] + (E_radl+P_radl)*u_xll*u_yll + Psi4 * P_radl * gxyL + F_rad_xl*u_yll + F_rad_yl*u_xll;
	Sxz[index] = Sxz[index] + (E_radl+P_radl)*u_xll*u_zll + Psi4 * P_radl * gxzL + F_rad_xl*u_zll + F_rad_zl*u_xll;
	Syz[index] = Syz[index] + (E_radl+P_radl)*u_yll*u_zll + Psi4 * P_radl * gyzL + F_rad_yl*u_zll + F_rad_zl*u_yll;
	*/
      }

  // printf ("End metric_source_terms_and_misc_vars \n");
}


extern "C" void CCTK_FCALL CCTK_FNAME(metric_source_terms_and_misc_vars)
  (const cGH **cctkGH,int *ext,
   double *rho,double *Sx,double *Sy,double *Sz,
   double *Sxx,double *Sxy,double *Sxz,double *Syy,double *Syz,double *Szz,
   double *tau_rad, double *S_rad_x, double *S_rad_y, double *S_rad_z,
   double *E_rad, double *F_rad0, double *F_radx, double *F_rady,double *F_radz, double *P_rad,
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
   double &rho_b_atm) {
  metric_source_terms_and_misc_vars (*cctkGH,ext,
				     rho,Sx,Sy,Sz,
				     Sxx,Sxy,Sxz,Syy,Syz,Szz,
				     tau_rad, S_rad_x,S_rad_y, S_rad_z,
				     E_rad, F_rad0, F_radx, F_rady, F_radz, P_rad,
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
				     rho_b_atm);

}
