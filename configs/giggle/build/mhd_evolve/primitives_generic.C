// Here is a wrapper for the 2d solver described in Noble et al.  
// It is meant to be an interface with their code.
// Note that it assumes a simple gamma law for the moment.  No hybrid.
#include <iostream>
#include <iomanip>
#include <fstream>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include "cctk.h"

#include "harm_primitives_headers.h"
#include "primitives_solver_header.h"
#include "harm_primitives_lowlevel.C"


#define SQR(x) ((x) * (x))

double min_prim(double val1,double val2);
double min_prim(double val1, double val2) {
  if(val1<val2) return val1;
  return val2;
}

extern "C" void CCTK_FCALL primitives_generic_
  (const cGH **cctkGH,int *ext,int *nghostzones,double *X,double *Y,double *Z,
   double *phi,double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
   double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,
   double *lapm1,double *shiftx,double *shifty,double *shiftz,
   double *Bx,double *By,double *Bz,
   double *eps_tot, double *eps_thermal, double *eps_cld, double *P_cld, 
   double *ka_gf, double *ks_gf, double *emission_gf, double *chi_rad, double *chi_rad_nue, double *Y_e, double *optd,  double *eta_nue,  
   double *ka_gf_nue, double *ks_gf_nue, double *emission_gf_nue, double *ka_gf_nux, double *ks_gf_nux, double *emission_gf_nux,
   double *mhd_st_x,double *mhd_st_y,double *mhd_st_z,double *tau,double *rho_star, double *rhoYe,
   double *vx,double *vy,double *vz,double *P,double *rho_b,double *h,double *u0,
   double &rhobatm,double &tau_atm, double &rho_max,
   double &xNS1,double &yNS1, double &xNS2,double &yNS2,  double &M_B, double &rad_T_fac, double &rad_T_cutoff, double &rad_T_pow, double &rad_T_floor, double *T_fluid,
   int &neos,int &ergo_star, double &ergo_sigma,double *rho_tab,double *P_tab,double *eps_tab,double *k_tab,double *gamma_tab,
   double *mhd_st_x_old,double *mhd_st_y_old,double *mhd_st_z_old,double *tau_old,double *rho_star_old,
   int &primitives_debug,double &Psi6threshold, double &rad_const, int &horizon_enforce_rho_profile, int &enable_OS_collapse, 
   int &rad_evolve_enable, int &compute_microphysics, int &microphysics_scheme, double &T_fluid_cgs_atm, int &rad_fix);

void primitives_generic(const cGH *cctkGH,int *ext,int *nghostzones,double *X,double *Y,double *Z,
			double *phi,double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
			double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,
			double *lapm1,double *shiftx,double *shifty,double *shiftz,
			double *Bx,double *By,double *Bz,
			double *eps_tot, double *eps_thermal, double *eps_cld, double *P_cld,
			double *ka_gf, double *ks_gf, double *emission_gf, double *chi_rad, double *chi_rad_nue, double *Y_e, double *optd,  double *eta_nue,  
			double *ka_gf_nue, double *ks_gf_nue, double *emission_gf_nue, double *ka_gf_nux, double *ks_gf_nux, double *emission_gf_nux,
			double *mhd_st_x,double *mhd_st_y,double *mhd_st_z,double *tau,double *rho_star, double *rhoYe,
			double *vx,double *vy,double *vz,double *P,double *rho_b,double *h,double *u0,
			double &rhobatm,double &tau_atm, double &rho_max, 
			double &xNS1,double &yNS1, double &xNS2,double &yNS2,  double &M_B, double &rad_T_fac, double &rad_T_cutoff, double &rad_T_pow, double &rad_T_floor, double *T_fluid,
			int &neos,int &ergo_star, double &ergo_sigma, double *rho_tab,double *P_tab,double *eps_tab,double *k_tab,double *gamma_tab,
			double *mhd_st_x_old,double *mhd_st_y_old,double *mhd_st_z_old,double *tau_old,double *rho_star_old,
			int &primitives_debug,double &Psi6threshold, double &rad_const, int &horizon_enforce_rho_profile,int &enable_OS_collapse, 
			int &rad_evolve_enable, int &compute_microphysics, int &microphysics_scheme, double &T_fluid_cgs_atm, int &rad_fix) {
  
  using namespace std;

  //Start the timer, so we can benchmark the primitives solver during evolution.
  //  Slower solver -> harder to find roots -> things may be going crazy!
  struct timeval start, end;
  long mtime, seconds, useconds;
  gettimeofday(&start, NULL);

  int failures=0,font_fixes=0;
  int pointcount=0;
  int failures_inhoriz=0;
  int pointcount_inhoriz=0;
  
  //  printf("START primitives_generic!!!!!!!!!!!!!!! compute_microphysics is %d \n", compute_microphysics);

  if (fabs(gamma_tab[0]-GAMMA) > 1.e-10) {
    printf("GAMMA=%e must be equal to gamma_tab[0]=%e",GAMMA,gamma_tab[0]);
     exit(1);
  }

  int pressure_cap_hit=0;

  double error_int_numer=0,error_int_denom=0;
  
  int rho_star_fix_applied=0;
#pragma omp parallel for reduction(+:failures,font_fixes,pointcount,failures_inhoriz,pointcount_inhoriz,error_int_numer,error_int_denom,pressure_cap_hit,rho_star_fix_applied) schedule(static)
  for(int k=0;k<ext[2];k++)
    for(int j=0;j<ext[1];j++)
      for(int i=0;i<ext[0];i++) {
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

	struct output_primitives new_primitives;
	struct output_stats stats;

	if(isnan(tau[index])) { printf("found nan in tau %d %d %d\n",i,j,k); }
	if(isnan(mhd_st_x[index])) { printf("found nan in msx %d %d %d\n",i,j,k); }
	if(isnan(mhd_st_y[index])) { printf("found nan in msy %d %d %d\n",i,j,k); }
	if(isnan(mhd_st_z[index])) { printf("found nan in msz %d %d %d\n",i,j,k); }

     	if(isnan(rho_star[index])) { printf("found nan in rho_star %d %d %d\n",i,j,k); }


	/*
	double vsq = SQR(vx[index]) + SQR(vy[index]) + SQR(vz[index]);
	double shiftsq = SQR(shiftx[index]) + SQR(shifty[index]) + SQR(shiftz[index]);
	if(vsq==0.0){printf("At the beginning of prim generic, vx=vy=vz= 0 \n");}
	if(shiftsq==0.0){printf("At the beginnning of prim generic, shiftx=shifty=shiftz= 0 \n");}
	if (u0[index]==0.0){printf("At the beginnning of prim generic, u0= 0 \n");}
	*/



	tau_old[index] = tau[index];
	rho_star_old[index] = rho_star[index];
	mhd_st_x_old[index] = mhd_st_x[index];
	mhd_st_y_old[index] = mhd_st_y[index];
	mhd_st_z_old[index] = mhd_st_z[index];

	double tau_orig = tau[index];
	double rho_star_orig = rho_star[index];
	double mhd_st_x_orig = mhd_st_x[index];
	double mhd_st_y_orig = mhd_st_y[index];
	double mhd_st_z_orig = mhd_st_z[index];

	int check=0;
	
	//double eps_tot; This is used to compute the total eps = eps_cold + eps_th + eps_nucl for microphysics run. 

	// VASILIS SAYS: CHECK THE PART WITH rho_star[index]<=0.0 
	// AJK's moto is that every five lines there's a bug 
	// and I coded more than 5 lines :-)

	if(rho_star[index]<=0.0) {
	  rho_star[index]=rhobatm*exp(6.0*phi[index]);
	  rho_star_fix_applied++;
	  

	  if(rhobatm<0) {
	    double xL=X[index];
	    double yL=Y[index];
	    double zL=Z[index];
	    double rL=sqrt(xL*xL+yL*yL+zL*zL);
	    double r_falloff_radius = 12.0;
	    double r_falloff_dr     = 5.0;
	    
	    double rhobatm_initial_value = 1.0e-8;
	    double rhobatm_final_value = 1.0e-11;
	    
	    rhobatm=-erf((rL-r_falloff_radius-r_falloff_dr*2.0)/r_falloff_dr)*(rhobatm_initial_value-rhobatm_final_value)*0.5 + (rhobatm_initial_value-rhobatm_final_value)*0.5 + rhobatm_final_value;
          }
	  

	  new_primitives.rho_b_new = rhobatm;
	  	  
	  // when rho_star<=0.0, we set P_tot = P_cold, and compute eps entirely based on EOS (no thermal part).
 
	  compute_pcold_cpp(new_primitives.rho_b_new, new_primitives.P_new, neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab);

	  double alpha     = lapm1[index]+1;
	  double alpha_inv = 1.0/alpha;
	  double betax = shiftx[index];
	  double betay = shifty[index];
	  double betaz = shiftz[index];

	  new_primitives.u0_new=1.0/alpha;
	  new_primitives.vx_new=-betax;
	  new_primitives.vy_new=-betay;
	  new_primitives.vz_new=-betaz;

	  rho_b[index] = new_primitives.rho_b_new;
	  P[index] = new_primitives.P_new;
	  u0[index] = new_primitives.u0_new;
	  vx[index] = new_primitives.vx_new;
	  vy[index] = new_primitives.vy_new;
	  vz[index] = new_primitives.vz_new;


	  double f1o4p     = 1.0/(4.0*M_PI);
	  double f1o4pa    = sqrt(f1o4p)*alpha_inv;

	  //double eps = new_primitives.P_new/(GAMMA_th - 1.0)/new_primitives.rho_b_new;
	  double eps, eps_cold;
	  compute_epscold_cpp(new_primitives.rho_b_new, new_primitives.P_new, eps_cold, neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab, enable_OS_collapse);
          if (compute_microphysics==0){
            eps_thermal[index] = new_primitives.P_new/(GAMMA_th - 1.0)/new_primitives.rho_b_new;
            eps = new_primitives.P_new/(GAMMA_th - 1.0)/new_primitives.rho_b_new + eps_cold;
          }
          else{
            double eps_thermalL;
            double T_fluid_floor = 1.0e1* (kb_cgs * G_cgs)/( fasterpow_prim(c_cgs,4.0) * fasterpow_prim(k2c, 1.0));
            compute_eps_th(T_fluid_floor, new_primitives.rho_b_new, eps_thermalL);
            eps_thermal[index] = eps_thermalL;
            eps = eps_thermalL + eps_cold;
            //      printf("T_fluid_floor=%e, eps_thermal=%e, eps_cold=%e, eps=%e \n", T_fluid_floor, eps_thermal, eps_cold, eps); 
            T_fluid[index] = T_fluid_floor;
          }

	  double h   = 1.0 + new_primitives.P_new/new_primitives.rho_b_new + eps;
          eps_tot[index] = eps;
          eps_cld[index] = eps_cold;
          P_cld[index] = new_primitives.P_new;

	  double phiL     = phi[index];
	  double Psi2     = exp(2.0*phiL);
	  double Psi4     = Psi2*Psi2;
	  double Psi6     = Psi4*Psi2;
	  double u_xl = 0.0;
	  double u_yl = 0.0;
	  double u_zl = 0.0;

	  rho_star[index] = alpha*rhobatm*Psi6;

	  double Bx_f1o4pa = Bx[index]*f1o4pa;
	  double By_f1o4pa = By[index]*f1o4pa;
	  double Bz_f1o4pa = Bz[index]*f1o4pa;

	  double gxx_phys = gxx[index]*Psi4;
	  double gxy_phys = gxy[index]*Psi4;
	  double gxz_phys = gxz[index]*Psi4;
	  double gyy_phys = gyy[index]*Psi4;
	  double gyz_phys = gyz[index]*Psi4;
	  double gzz_phys = gzz[index]*Psi4;
  
	  double B_xl  = (gxx_phys *Bx[index] + gxy_phys * By[index] +
			  gxz_phys * Bz[index] );
	  double B_yl  = (gxy_phys *Bx[index] + gyy_phys * By[index] +
			  gyz_phys * Bz[index]);
	  double B_zl  = (gxz_phys *Bx[index] + gyz_phys * By[index] +
			  gzz_phys * Bz[index]);
	  double B_xl_f1o4pa = B_xl*f1o4pa;
	  double B_yl_f1o4pa = B_yl*f1o4pa;
	  double B_zl_f1o4pa = B_zl*f1o4pa;

	  double B2 = ( gxx_phys*SQR(Bx_f1o4pa) +
                2.0*gxy_phys*Bx_f1o4pa*By_f1o4pa + 2.0*gxz_phys*Bx_f1o4pa*Bz_f1o4pa +
			gyy_phys*SQR(By_f1o4pa) + 2.0*gyz_phys*By_f1o4pa*Bz_f1o4pa +
			gzz_phys*SQR(Bz_f1o4pa) );
        
	  double sb0 = u_xl*Bx_f1o4pa + u_yl*By_f1o4pa + u_zl*Bz_f1o4pa;
	  double sb2 = (B2 + SQR(sb0))/SQR(new_primitives.u0_new);
	  double sb_x = (B_xl_f1o4pa + u_xl*sb0)/new_primitives.u0_new;
	  double sb_y = (B_yl_f1o4pa + u_yl*sb0)/new_primitives.u0_new;
	  double sb_z = (B_zl_f1o4pa + u_zl*sb0)/new_primitives.u0_new;


	  mhd_st_x[index] = rho_star[index]*h*u_xl +
	    alpha*Psi6*new_primitives.u0_new*sb2*u_xl - alpha*Psi6*sb0*sb_x;
	  mhd_st_y[index] = rho_star[index]*h*u_yl +
	    alpha*Psi6*new_primitives.u0_new*sb2*u_yl - alpha*Psi6*sb0*sb_y;
	  mhd_st_z[index] = rho_star[index]*h*u_zl +
	    alpha*Psi6*new_primitives.u0_new*sb2*u_zl - alpha*Psi6*sb0*sb_z;
	  tau[index] = ((alpha*new_primitives.u0_new-1.0)+
			(new_primitives.P_new/new_primitives.rho_b_new+eps)*alpha*new_primitives.u0_new)*rho_star[index] +
	                Psi6*sb2*SQR(alpha*new_primitives.u0_new)
	               - Psi6*(new_primitives.P_new+sb2*0.5)-Psi6*SQR(alpha*sb0);
	  

	  if(isnan(tau[index])) { printf("Check1: found nan in tau %d %d %d\n",i,j,k); }
	} else {
	  // !!!!!!!!!!!!!!! rho_star > 0 !!!!!!!!!!!!!!!!!!!!
	  // Lunan: OS collapse ID cannot be applied to tau floor. Nans are generated.
	  if(enable_OS_collapse==0){
	  // Apply the tau floor 
	  apply_tau_floor(index,tau_atm,rhobatm,
			  Bx,By,Bz,
			  tau,rho_star,mhd_st_x,mhd_st_y,mhd_st_z,
			  phi,gxx,gxy,gxz,gyy,gyz,gzz,
			  gupxx,gupxy,gupxz,gupyy,gupyz,gupzz, 
			  new_primitives,
			  lapm1,shiftx,shifty,shiftz,
			  neos,ergo_star,ergo_sigma,rho_tab,P_tab,eps_tab,k_tab,gamma_tab,
			  Psi6threshold);
	  }
	  tau_orig = tau[index];
	  rho_star_orig = rho_star[index];
	  mhd_st_x_orig = mhd_st_x[index];
	  mhd_st_y_orig = mhd_st_y[index];
	  mhd_st_z_orig = mhd_st_z[index];

	  if(isnan(tau[index])) { printf("2found nan in tau %d %d %d\n",i,j,k); }
	  //if(isnan(mhd_st_x[index])) { printf("2found nan in msx %d %d %d\n",i,j,k); exit(1); }
	  //if(isnan(mhd_st_y[index])) { printf("2found nan in msy %d %d %d\n",i,j,k); exit(1); }
	  //if(isnan(mhd_st_z[index])) { printf("2found nan in msz %d %d %d\n",i,j,k); exit(1); }  \
	  //if(isnan(rho_star[index])) { printf("2found nan in rho_star %d %d %d\n",i,j,k); exit(1); }
	  stats.font_fixed=0;

	  if(isnan(vx[index])) { printf("BEFORE prim_gamma_law BAD vx inside primitives generic: %d %d %d %e %e %e %e %e %e %e %e\n",
					i,j,k,vx[index],vy[index],vz[index],rho_b[index],P[index],u0[index],rho_star[index],tau[index]); }

	  
	  double vsq = SQR(vx[index]) + SQR(vy[index]) + SQR(vz[index]);
	  double shiftsq = SQR(shiftx[index]) + SQR(shifty[index]) + SQR(shiftz[index]);
	  	  
	  double Y_eL = min_prim(0.99, rhoYe[index]/rho_star[index]);
	  Y_e[index] = Y_eL;

	  for(int ii=0;ii<3;ii++) {
	    check = harm_primitives_gammalaw_lowlevel(index,X,Y,Z,
						      enable_OS_collapse,
						      phi,gxx,gxy,gxz,gyy,gyz,gzz,
						      gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,
						      lapm1,shiftx,shifty,shiftz,
						      Bx,By,Bz,
						      T_fluid, 
						      mhd_st_x,mhd_st_y,mhd_st_z,tau,rho_star,
						      vx,vy,vz,P,rho_b,h,u0,
						      ka_gf, ks_gf, emission_gf, chi_rad, chi_rad_nue, Y_e, optd, eta_nue,
						      ka_gf_nue, ks_gf_nue, emission_gf_nue,ka_gf_nux, ks_gf_nux, emission_gf_nux,
						      eps_tot, eps_thermal, eps_cld, P_cld,
						      rhobatm,tau_atm,
						      neos,ergo_star, ergo_sigma,
						      rho_tab,P_tab,eps_tab,k_tab,gamma_tab,	     
						      rho_max,
						      new_primitives,stats,
						      Psi6threshold, horizon_enforce_rho_profile, 
						      rad_evolve_enable,compute_microphysics, microphysics_scheme, T_fluid_cgs_atm, rad_fix);
	  
	    // Set primitives, and/or provide a better guess.
	    P[index] = new_primitives.P_new;
	    rho_b[index] = new_primitives.rho_b_new;
	    u0[index] = new_primitives.u0_new;
	    vx[index] = new_primitives.vx_new;
	    vy[index] = new_primitives.vy_new;
	    vz[index] = new_primitives.vz_new;

	    if(check==0) ii=4;
	    if(isnan(vx[index])) { printf("AFTER prim_gamma_law BAD vx inside primitives generic: %d %d %d %e %e %e %e %e %e %e %e\n",
					  i,j,k,vx[index],vy[index],vz[index],rho_b[index],P[index],u0[index],rho_star[index],tau[index]); }

	  }
	}

	//Finally, we set h, the enthalpy:
	h[index] = 1.0 + P[index]/rho_b[index] + eps_tot[index];

	if(isnan(tau[index])) { printf("Check2: found nan in tau %d %d %d\n",i,j,k); }
	/*************************************************************************************************************************/
	// DIAGNOSTICS:
	//Pressure cap hit?
	//double P_cold = rho_b[index]*rho_b[index];
	/*
	if(exp(phi[index]*6.0) > Psi6threshold) {
	  //if(P[index]/P_cold > 0.99*1e6) pressure_cap_hit++;
	} else {
	  if(P[index]/P_cold > 0.99*1e3 && rho_b[index]>100.0*rhobatm) pressure_cap_hit++;
	}
	*/
        //Now we compute the difference between original & new conservatives, for diagnostic purposes:
        error_int_numer += fabs(tau[index] - tau_orig) + fabs(rho_star[index] - rho_star_orig) + 
          fabs(mhd_st_x[index] - mhd_st_x_orig) + fabs(mhd_st_y[index] - mhd_st_y_orig) + fabs(mhd_st_z[index] - mhd_st_z_orig);
        error_int_denom += tau_orig + rho_star_orig + fabs(mhd_st_x_orig) + fabs(mhd_st_y_orig) + fabs(mhd_st_z_orig);

	if(stats.font_fixed==1) font_fixes++;
        if(check!=0) failures++;
        pointcount++;
        if(exp(phi[index]*6.0)>15.0) {
          if(check!=0) failures_inhoriz++;
          pointcount_inhoriz++;
        }
	/***************************************************************************************************************************/
	

      }

  gettimeofday(&end, NULL);

  seconds  = end.tv_sec  - start.tv_sec;
  useconds = end.tv_usec - start.tv_usec;

  mtime = ((seconds) * 1000 + useconds/1000.0) + 0.999;  // We add 0.999 since mtime is a long int; this rounds up the result before setting the value.  Here, rounding down is incorrect.
  printf("Pointcount: %d Font fixes: %d Failures: %d, InHoriz: %d / %d = %.3e\t%f solutions/second, Error: %e, ErrDenom: %e, rho*fixes: %d\n",
         pointcount,font_fixes,
         failures,
         failures_inhoriz,pointcount_inhoriz,failures_inhoriz/((double)pointcount_inhoriz+1e-10), 
         ext[0]*ext[1]*ext[2] / ((double)mtime/1000.0),
         error_int_numer/error_int_denom,error_int_denom,
	 rho_star_fix_applied);

  if(pressure_cap_hit!=0) {
    //printf("PRESSURE CAP HIT %d TIMES!  Outputting debug file!\n",pressure_cap_hit);
  }



  //if(pressure_cap_hit!=0 || (primitives_debug==1 && (failures>0 || error_int_numer/error_int_denom>0.8))) {
  //if((primitives_debug==1 && (failures>0 || error_int_numer/error_int_denom>0.8))) {
  if(primitives_debug==1) {

    ofstream myfile;
    char filename[100];
    srand(time(NULL));
    sprintf(filename,"primitives_debug-%d.dat",rand());
    myfile.open (filename, ios::out | ios::binary);
    //myfile.open ("data.bin", ios::out | ios::binary);
    myfile.write((char*)ext, 3*sizeof(int));

    myfile.write((char*)&rhobatm, 1*sizeof(double));
    myfile.write((char*)&tau_atm, 1*sizeof(double));

    myfile.write((char*)&Psi6threshold, 1*sizeof(double));

    double gamma_th=GAMMA_th;
    myfile.write((char*)&gamma_th, 1*sizeof(double));
    myfile.write((char*)&neos,     1*sizeof(int));

    myfile.write((char*)gamma_tab, (neos+1)*sizeof(double));
    myfile.write((char*)k_tab,     (neos+1)*sizeof(double));

    myfile.write((char*)eps_tab,   neos*sizeof(double));
    myfile.write((char*)rho_tab,   neos*sizeof(double));
    myfile.write((char*)P_tab,     neos*sizeof(double));

    int fullsize=ext[0]*ext[1]*ext[2];
    myfile.write((char*)X,   (fullsize)*sizeof(double));
    myfile.write((char*)Y,   (fullsize)*sizeof(double));
    myfile.write((char*)Z,   (fullsize)*sizeof(double));
    myfile.write((char*)phi, (fullsize)*sizeof(double));
    myfile.write((char*)gxx, (fullsize)*sizeof(double));
    myfile.write((char*)gxy, (fullsize)*sizeof(double));
    myfile.write((char*)gxz, (fullsize)*sizeof(double));
    myfile.write((char*)gyy, (fullsize)*sizeof(double));
    myfile.write((char*)gyz, (fullsize)*sizeof(double));
    myfile.write((char*)gzz, (fullsize)*sizeof(double));

    myfile.write((char*)gupxx, (fullsize)*sizeof(double));
    myfile.write((char*)gupxy, (fullsize)*sizeof(double));
    myfile.write((char*)gupxz, (fullsize)*sizeof(double));
    myfile.write((char*)gupyy, (fullsize)*sizeof(double));
    myfile.write((char*)gupyz, (fullsize)*sizeof(double));
    myfile.write((char*)gupzz, (fullsize)*sizeof(double));

    myfile.write((char*)shiftx, (fullsize)*sizeof(double));
    myfile.write((char*)shifty, (fullsize)*sizeof(double));
    myfile.write((char*)shiftz, (fullsize)*sizeof(double));

    myfile.write((char*)lapm1, (fullsize)*sizeof(double));
 
    myfile.write((char*)tau_old,      (fullsize)*sizeof(double));
    myfile.write((char*)mhd_st_x_old, (fullsize)*sizeof(double));
    myfile.write((char*)mhd_st_y_old, (fullsize)*sizeof(double));
    myfile.write((char*)mhd_st_z_old, (fullsize)*sizeof(double));

    myfile.write((char*)rho_star_old, (fullsize)*sizeof(double));

    myfile.write((char*)Bx,   (fullsize)*sizeof(double));
    myfile.write((char*)By,   (fullsize)*sizeof(double));
    myfile.write((char*)Bz,   (fullsize)*sizeof(double));

    myfile.write((char*)vx,   (fullsize)*sizeof(double));
    myfile.write((char*)vy,   (fullsize)*sizeof(double));
    myfile.write((char*)vz,   (fullsize)*sizeof(double));
    myfile.write((char*)P,    (fullsize)*sizeof(double));
    myfile.write((char*)rho_b,(fullsize)*sizeof(double));
    myfile.write((char*)h,    (fullsize)*sizeof(double));
    myfile.write((char*)u0,   (fullsize)*sizeof(double));

    int checker=1063; myfile.write((char*)&checker,sizeof(int));

    myfile.close();
    printf("Finished writing...\n");
    //sleep(2);
    //exit(1);
  }

  printf("END primitives_generic!!!!!!!!!!!!!!! \n");

}



extern "C" void CCTK_FCALL primitives_generic_
  (const cGH **cctkGH,int *ext,int *nghostzones,double *X,double *Y,double *Z,
   double *phi,double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
   double *gupxx,double *gupxy,double *gupxz,double *gupyy,double *gupyz,double *gupzz,
   double *lapm1,double *shiftx,double *shifty,double *shiftz,
   double *Bx,double *By,double *Bz,
   double *eps_tot, double *eps_thermal,double *eps_cld, double *P_cld,
   double *ka_gf, double *ks_gf, double *emission_gf, double *chi_rad, double *chi_rad_nue, double *Y_e, double *optd, double *eta_nue,  
   double *ka_gf_nue, double *ks_gf_nue, double *emission_gf_nue, double *ka_gf_nux, double *ks_gf_nux, double *emission_gf_nux,
   double *mhd_st_x,double *mhd_st_y,double *mhd_st_z,double *tau,double *rho_star, double *rhoYe,
   double *vx,double *vy,double *vz,double *P,double *rho_b,double *h,double *u0,
   double &rhobatm,double &tau_atm, double &rho_max,
   double &xNS1,double &yNS1, double &xNS2,double &yNS2,  double &M_B, double &rad_T_fac, double &rad_T_cutoff, double &rad_T_pow, double &rad_T_floor, double *T_fluid,
   int &neos,int &ergo_star, double &ergo_sigma, double *rho_tab,double *P_tab,double *eps_tab,double *k_tab,double *gamma_tab,
   double *mhd_st_x_old,double *mhd_st_y_old,double *mhd_st_z_old,double *tau_old,double *rho_star_old,
   int &primitives_debug,double &Psi6threshold, double &rad_const, int &horizon_enforce_rho_profile, int &enable_OS_collapse, 
   int &rad_evolve_enable, int &compute_microphysics, int &microphysics_scheme, double &T_fluid_cgs_atm, int &rad_fix) 
{
  primitives_generic(*cctkGH,ext,nghostzones,X,Y,Z,
		     phi,gxx,gxy,gxz,gyy,gyz,gzz,
		     gupxx,gupxy,gupxz,gupyy,gupyz,gupzz,
		     lapm1,shiftx,shifty,shiftz,
		     Bx,By,Bz,
		     eps_tot, eps_thermal, eps_cld, P_cld,
		     ka_gf, ks_gf, emission_gf, chi_rad, chi_rad_nue, Y_e, optd, eta_nue,
		     ka_gf_nue, ks_gf_nue, emission_gf_nue, ka_gf_nux,ks_gf_nux, emission_gf_nux,
		     mhd_st_x,mhd_st_y,mhd_st_z,tau,rho_star,rhoYe,
		     vx,vy,vz,P,rho_b,h,u0,
		     rhobatm,tau_atm,rho_max,
		     xNS1, yNS1, xNS2, yNS2, M_B, rad_T_fac, rad_T_cutoff, rad_T_pow, rad_T_floor, T_fluid,
		     neos,ergo_star,ergo_sigma,rho_tab,P_tab,eps_tab,k_tab,gamma_tab,
		     mhd_st_x_old,mhd_st_y_old,mhd_st_z_old,tau_old,rho_star_old,
		     primitives_debug, Psi6threshold, rad_const, horizon_enforce_rho_profile, 
		     enable_OS_collapse, rad_evolve_enable, compute_microphysics, 
		     microphysics_scheme, T_fluid_cgs_atm, rad_fix);

}
