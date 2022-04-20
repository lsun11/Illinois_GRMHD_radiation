#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <sys/time.h>
#include "cctk.h"

#include "primitives_solver_header.h"




#define SQR(x) ((x) * (x))
#define AXISYM 4
#define MAX_GAMMA 100.0



struct auxarray_rad {
  double tau_rad,S_rad_x, S_rad_y, S_rad_z;
  double ux, uy, uz, u_x, u_y, u_z;
  double gupxx_phys,gupxy_phys,gupxz_phys,gupyy_phys,gupyz_phys,gupzz_phys;
  double Psi6, u_0, alpn1, g_00;
  double shiftx, shifty, shiftz;
  double rho_s, Erad_atm_cut, zeta_cut;
  double sri_scal_inv, tau_rad_scal_inv;
};


void newt2_cpp_rad(double x[5],struct auxarray_rad &aux,
	       void (*function_rad)(int &n,double *x,double *fvec,struct auxarray_rad &aux),
	       void (*jacobian_rad)
	       (int &n,double *x,struct auxarray_rad &aux, double *fvec,int &np, double fjac[][5]),
	       int &n,bool &check,
	       int *indx,double *g,double *p,double *xold,double *fvec,int &MAXITS,double &STPMX);

void function_rad(int &n,double *x,double *fvec,struct auxarray_rad &aux);
void function_rad_font_fix(int &n,double *x,double *fvec,struct auxarray_rad &aux);
void jacobian_rad
(int &n,double *x,struct auxarray_rad &aux, double *fvec,int &np, double fjac[][5]);
void jacobian_rad_font_fix
(int &n,double *x,struct auxarray_rad &aux, double *fvec,int &np, double fjac[][5]);
double max_val_rad(double val1,double val2);
double **dmatrix_newt(long nrl, long nrh, long ncl, long nch);
void free_dmatrix_newt(double **m, long nrl, long nrh, long ncl, long nch);
void ludcmp_rad(double a[][5], int n, int *indx, double *d);
void lubksb_rad(double a[][5], int n, int *indx, double b[]);
void lnsrch_rad(int n, double xold[], double fold, double g[], double p[], double x[5],
		 double *f, double stpmax, bool &check, 
		 double (*fminrad)
		 (double *x,struct auxarray_rad &aux,double *fvec,
		  void (*nrfuncv_rad)(int &n,double *x,double *fvec, struct auxarray_rad &aux),
		  int &n),
		 void (*nrfuncv_rad)(int &n,double *x,double *fvec,struct auxarray_rad &aux),
		 struct auxarray_rad &aux,double *fvec);
double fmin_rad(double *x,struct auxarray_rad &aux,double *fvec,
		 void (*nrfuncv_rad)(int &n,double *x,double *fvec,struct auxarray_rad &aux),
		 int &n);




extern "C" void CCTK_FCALL CCTK_FNAME(primitive_vars_rad_cpp)
  (int *ext,int *nghostzones, double *X, double *Y, double *Z, 
   double *rho_star,
   double *tau_rad, double *S_rad_x, double *S_rad_y, double *S_rad_z,
   double *Sx, double *Sy, double *Sz, double *rho, 
   double *Sxx, double *Sxy, double *Sxz, double *Syy, double *Syz, double *Szz,
   double *E_rad,  double *F_radx, double *F_rady, double *F_radz, double *F_rad0, double *F,
   double *P_radxx, double *P_radyy, double *P_radzz, double *P_radxy, double *P_radxz, double *P_radyz,
   double *phi, double *alpha, double *shiftx, double *shifty, double *shiftz, double *gxx, double *gxy, double *gxz, double *gyy, double *gyz, double *gzz,
   double *gupxx, double *gupxy, double *gupxz, double *gupyy, double *gupyz, double *gupzz,
   double *vx, double *vy, double *vz, double *u0, double *chi_rad, double *zeta_rad,
   const cGH **cctkGH, double *failure_tracker_rad, int &ignore_ghostzones, int &repairs_rad_needed, double &Psi6threshold, double &Erad_atm_cut);


//-----------------------------------------------------------------------------
//
// reconstruct primitive variables, compute sources for "hybrid" EOS
//
//-----------------------------------------------------------------------------
void primitive_vars_rad_cpp
(int *ext,int *nghostzones, double *X, double *Y, double *Z, 
 double *rho_star,
 double *tau_rad, double *S_rad_x, double *S_rad_y, double *S_rad_z,
 double *Sx, double *Sy, double *Sz, double *rho,
 double *Sxx, double *Sxy, double *Sxz, double *Syy, double *Syz, double *Szz,
 double *E_rad, double *F_radx, double *F_rady, double *F_radz, double *F_rad0, double *F,
 double *P_radxx, double *P_radyy, double *P_radzz, double *P_radxy, double *P_radxz,double *P_radyz,
 double *phi, double *alpha, double *shiftx, double *shifty, double *shiftz, double *gxx, double *gxy, double *gxz, double *gyy, double *gyz, double *gzz,
 double *gupxx, double *gupxy, double *gupxz, double *gupyy, double *gupyz, double *gupzz,
 double *vx, double *vy, double *vz, double *u0, double *chi_rad,double *zeta_rad,
 const cGH *cctkGH, double *failure_tracker_rad, int &ignore_ghostzones, int &repairs_rad_needed, double &Psi6threshold, double &Erad_atm_cut) {

  //int count = 0;
  double f1o4p = 1.0/(4.0*M_PI);
  double f1o8p = f1o4p*0.5;

  // Timer:
  // reset the clock
  struct timeval start, end;
  long mtime, seconds, useconds;    
  gettimeofday(&start, NULL);

  printf("===============Start Primitives_rad.C========================= \n");

  //
  // Input translation
  //
  /* Set up variables used in the grid loop for the physical grid points */
  int istart = 0;
  int jstart = 0;
  int kstart = 0;
  int iend = ext[0];
  int jend = ext[1];
  int kend = ext[2];

  if(ignore_ghostzones==1) {
    istart = nghostzones[0];
    jstart = nghostzones[1];
    kstart = nghostzones[2];
    iend   = ext[0]-nghostzones[0];
    jend   = ext[1]-nghostzones[1];
    kend   = ext[2]-nghostzones[2];
  } else if(ignore_ghostzones==0 || ignore_ghostzones==-1) {
    //do nothing
  } else if(ignore_ghostzones==2) {
    istart = nghostzones[0]-1;
    jstart = nghostzones[1]-1;
    kstart = nghostzones[2]-1;
    iend   = ext[0]-nghostzones[0]+1;
    jend   = ext[1]-nghostzones[1]+1;
    kend   = ext[2]-nghostzones[2]+1;
  } else {
    printf("YOU MUST SET ignore_ghostzones (last variable in primitives solver function call) to zero (compute everywhere on grid) or one (ignore ghostzones) or two (ignore all but one ghostzone, so that averaging works when fontfix fails.) or -1 (only compute inside the ghostzones)\n");
    exit(1);
  }

  double dX = X[CCTK_GFINDEX3D(cctkGH,1,0,0)]-X[CCTK_GFINDEX3D(cctkGH,0,0,0)];
  double dY = Y[CCTK_GFINDEX3D(cctkGH,0,1,0)]-Y[CCTK_GFINDEX3D(cctkGH,0,0,0)];
  double dZ = Z[CCTK_GFINDEX3D(cctkGH,0,0,1)]-Z[CCTK_GFINDEX3D(cctkGH,0,0,0)];

  double dV = dX * dY * dZ;


  //-----------------------------------------------------------------------------
  // Funny parameters...  (See Shibata, Oohara & Nakamura...)
  //-----------------------------------------------------------------------------

  // Initialize the failure_tracker_rad array
#pragma omp parallel for
  for(int i=0;i<ext[0];i++) for(int j=0;j<ext[1];j++) for(int k=0;k<ext[2];k++) {
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	failure_tracker_rad[index] = 0.0;
      }
  repairs_rad_needed = 0;

  //
  // Now go to each gridpoint
  //
#pragma omp parallel for
  for(int k=kstart;k<kend;k++)
    for(int j=jstart;j<jend;j++)
      for(int i=istart;i<iend;i++) {
        int compute_primitives = 1;

        // Exclude regions not in the ghostzones when ignore_ghostzones==-1
        if ( ignore_ghostzones==-1) {
	  if ( (k>nghostzones[2]-1 && k<kend-nghostzones[2]) &&
                  (j>nghostzones[1]-1 && j<jend-nghostzones[1]) &&
                  (i>nghostzones[0]-1 && i<iend-nghostzones[0]) )
                compute_primitives = 0;
        }

	double UU_radguess[5];
	bool recom_rad = false;

        if (compute_primitives==1) { 
	   int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

	   int MAXITS_beforefontfix=50;
	   double STPMX_beforefontfix=20.0; // <== faster than 100.0, when maxits set to 50
	   

	   double au0m1;

	   //Structure that holds auxilliary data for Newton-Raphson solver
	   
	   double Psi2 = exp(2.0*phi[index]);
	   double Psim2= 1.0/Psi2;
	   double Psi4 = Psi2*Psi2;
	   double Psim4= 1.0/(Psi2*Psi2);
	   double Psi6 = Psi2*Psi4;
	   double Psim6= 1.0/(Psi2*Psi4);
	   double gxxL = gxx[index];
	   double gxyL = gxy[index];
	   double gxzL = gxz[index];
	   double gyyL = gyy[index];
	   double gyzL = gyz[index];
	   double gzzL = gzz[index];

	   double shiftxL = shiftx[index];
	   double shiftyL = shifty[index];
	   double shiftzL = shiftz[index];

           double shift_xL = Psi4*(shiftxL*gxxL + shiftyL*gxyL + shiftzL*gxzL);
           double shift_yL = Psi4*(shiftxL*gxyL + shiftyL*gyyL + shiftzL*gyzL);
           double shift_zL = Psi4*(shiftxL*gxzL + shiftyL*gyzL + shiftzL*gzzL);

           double tau_radl = tau_rad[index];
           double S_rad_xl = S_rad_x[index];
           double S_rad_yl = S_rad_y[index];
           double S_rad_zl = S_rad_z[index];


           double gupxx_phys = gupxx[index]*Psim4;
           double gupxy_phys = gupxy[index]*Psim4;
           double gupxz_phys = gupxz[index]*Psim4;
           double gupyy_phys = gupyy[index]*Psim4;
           double gupyz_phys = gupyz[index]*Psim4;
           double gupzz_phys = gupzz[index]*Psim4;


	   double u0L = u0[index];

	   double u_xL = (gxxL*(shiftxL+vx[index]) +
			  gxyL*(shiftyL+vy[index]) +
			  gxzL*(shiftzL+vz[index]))*Psi4*u0L;
	   double u_yL = (gxyL*(shiftxL+vx[index]) +
			  gyyL*(shiftyL+vy[index]) +
			  gyzL*(shiftzL+vz[index]))*Psi4*u0L;
	   double u_zL = (gxzL*(shiftxL+vx[index]) +
			  gyzL*(shiftyL+vy[index]) +
			  gzzL*(shiftzL+vz[index]))*Psi4*u0L;

     	   double alpn1 = alpha[index] + 1.0;
	   double alpn1_inv = 1.0/(alpha[index] + 1.0);
	   double beta2 = shiftxL*shift_xL + shiftyL*shift_yL + shiftzL*shift_zL;
	   double udotbeta = u0L*(vx[index]*shift_xL + vy[index]*shift_yL + vz[index]*shift_zL);
	   double g_00L =beta2-alpn1*alpn1;
	   double u_0L = g_00L*u0L + udotbeta;
	   
	   //double F_rad0L = F_rad0[index];	     //- (F_radx[index]*u_xL + F_rady[index]*u_yL + F_radz[index]*u_zL)/u_0L;	   
	   double F_rad0L = - (F_radx[index]*u_xL + F_rady[index]*u_yL + F_radz[index]*u_zL)/u_0L;
	   double F_rad_x = (gxxL*F_radx[index] +
			     gxyL*F_rady[index] +
			     gxzL*F_radz[index])*Psi4 + shift_xL*F_rad0L;
	   double F_rad_y = (gxyL*F_radx[index] +
			     gyyL*F_rady[index] +
			     gyzL*F_radz[index])*Psi4 + shift_yL*F_rad0L;
	   double F_rad_z = (gxzL*F_radx[index] +
			     gyzL*F_rady[index] +
			     gzzL*F_radz[index])*Psi4 + shift_zL*F_rad0L;
	   
	   double rho_s = rho_star[index];	  
	   
	   int nn = 4;
	   double UU_rad[5];
	   
           double uxL = -shiftxL*u0[index] + gupxx_phys*u_xL +
             gupxy_phys*u_yL + gupxz_phys*u_zL;
           double uyL = -shiftyL*u0[index] + gupxy_phys*u_xL +
             gupyy_phys*u_yL + gupyz_phys*u_zL;
           double uzL = -shiftzL*u0[index] + gupxz_phys*u_xL +
             gupyz_phys*u_yL + gupzz_phys*u_zL;

	   double E_radl = E_rad[index];	     

           // Apply a atmopheric floor before the solver.                                                                                    
	   if (E_radl < Erad_atm_cut)
             {
	       recom_rad = true;
               E_radl = Erad_atm_cut;
               F_rad_x = 0.0;
               F_rad_y = 0.0;
               F_rad_z = 0.0;
             }
	   
	   double F_rad_0l, F_rad0l;
           F_rad_0l = - (F_rad_x * uxL + F_rad_y * uyL + F_rad_z * uzL)/u0L;
	   F_rad0l = (-F_rad_0l+F_rad_x*shiftxL+F_rad_y*shiftyL+F_rad_z*shiftzL)/SQR(alpn1);

           // Limit Flux                                                                                                                                                          
           double FaFa = F_radx[index]*F_rad_x + F_rady[index]*F_rad_y + F_radz[index]*F_rad_z + F_rad0l*F_rad_0l;
	   double fac_rad = sqrt(SQR(E_radl)/FaFa);
           if (FaFa > SQR(E_radl)){
             recom_rad = true;
             F_rad_x *= fac_rad;
             F_rad_y *= fac_rad;
             F_rad_z *= fac_rad;
           }

	   double maxF_i = max_val_rad(fabs(F_rad_x),
				       max_val_rad(fabs(F_rad_y),fabs(F_rad_z) ));
	   
	   double sri_scal =max_val_rad(0.001*Psi6*(E_radl+F_rad0L+maxF_i),
					max_val_rad(fabs(S_rad_x[index]),
						    max_val_rad(fabs(S_rad_y[index]),fabs(S_rad_z[index]) )) );
	   if (sri_scal==0.0) 
	     {
	       printf("Inside primitive rad, sri_scal is zero!!! E_radl, maxF_i, Psi6 = %e,%e,%e \n", E_radl, maxF_i, Psi6);
	     }

	   if(isnan(sri_scal))
	     {
               printf("Inside primitive rad, sri_scal is NAN!!! E_radl, maxF_i, Psi6, F_rad0L = %e,%e,%e,%e \n", E_radl, maxF_i, Psi6, F_rad0L);
	       printf("0.001*Psi6*(E_radl+F_rad0L+maxF_i), fabs(S_rad_x[index]), fabs(S_rad_y[index]), fabs(S_rad_z[index]) = %e,%e,%e,%e \n", 0.001*Psi6*(E_radl+F_rad0L+maxF_i), fabs(S_rad_x[index]), fabs(S_rad_y[index]), fabs(S_rad_z[index])  );
             }

	   double tau_rad_scal = max_val_rad(0.01*Psi6*(E_radl+2.0*F_rad0L), fabs(tau_rad[index]));
	     
	     
	     //	     double sri_scal = max_val_rad(fabs(S_rad_x[index]),
	     //				 max_val_rad(fabs(S_rad_y[index]),fabs(S_rad_z[index]) ));
	     //	   double tau_rad_scal = fabs(tau_rad[index]);
	     
	   double zeta_cut = Erad_atm_cut*1.5;

	   UU_rad[0] = 0.0;
 	   UU_rad[1] = F_rad_x;
	   UU_rad[2] = F_rad_y;
	   UU_rad[3] = F_rad_z;
	   UU_rad[4] = E_radl;

	  
	   // CHECK ABNORMAL VALUES!!!
	   if (abs(UU_rad[1]) > 100.0 || abs(UU_rad[2]) > 100.0 || abs(UU_rad[3]) > 100.0)
	     {
	       printf ("Inside primitives_rad.C, before newt2, UU_rad[0]=%e, UU_rad[1]=%e, UU_rad[2]=%e, UU_rad[3]=%e, UU_rad[5]=%e, \n", UU_rad[0], UU_rad[1], UU_rad[2], UU_rad[3], UU_rad[4]);
	       printf ("E_rad[index], tau_rad[index], S_rad_x[index], S_rad_y[index], S_rad_z[index] = %e,%e,%e,%e,%e \n", E_rad[index], tau_rad[index], S_rad_x[index], S_rad_y[index], S_rad_z[index]);
	       printf ("F_radx[index], F_rady[index], F_radz[index], F_rad0[index], FaFa = %e,%e,%e,%e,%e \n", F_radx[index], F_rady[index], F_radz[index], F_rad0[index], FaFa);
	       printf ("gxxL,gxyL,gxzL,gyyL,gxzL,gyzL,shift_xL,shift_yL,shift_zL,u0[index],uxL,uyL,uzL=%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e \n", gxxL,gxyL,gxzL,gyyL,gxzL,gyzL,shift_xL,shift_yL,shift_zL,u0[index],uxL,uyL,uzL);
	       printf ("rho_star[index], vx[index], vy[index], vz[index] = %e,%e,%e,%e \n", rho_star[index], vx[index], vy[index], vz[index]);
	       printf ("X,Y,Z=%e,%e,%e \n", X[index], Y[index], Z[index]);	
	     }       




 
	   struct auxarray_rad AUX_rad;
	   AUX_rad.rho_s = rho_s;
	   AUX_rad.tau_rad = tau_rad[index];
	   AUX_rad.S_rad_x = S_rad_x[index];
	   AUX_rad.S_rad_y = S_rad_y[index];
	   AUX_rad.S_rad_z = S_rad_z[index];
	   
	   AUX_rad.ux = u0L * vx[index];
	   AUX_rad.uy = u0L * vy[index];
	   AUX_rad.uz = u0L * vz[index];
	   AUX_rad.u_x = u_xL;
	   AUX_rad.u_y = u_yL;
	   AUX_rad.u_z = u_zL;
	   
	   AUX_rad.gupxx_phys = gupxx[index]*Psim4;
	   AUX_rad.gupxy_phys = gupxy[index]*Psim4;
	   AUX_rad.gupxz_phys = gupxz[index]*Psim4;
	   AUX_rad.gupyy_phys = gupyy[index]*Psim4;
	   AUX_rad.gupyz_phys = gupyz[index]*Psim4;
	   AUX_rad.gupzz_phys = gupzz[index]*Psim4;
	   
	   //Note: Psi6 = sqrtg !
	   AUX_rad.Psi6 = Psi6;
	   AUX_rad.u_0 = u_0L;
	   AUX_rad.alpn1 = alpn1;
	   AUX_rad.g_00 = g_00L;
	   
	   AUX_rad.shiftx = shiftxL;
	   AUX_rad.shifty = shiftyL;
	   AUX_rad.shiftz = shiftzL;
	   
	   AUX_rad.Erad_atm_cut = Erad_atm_cut;
	   AUX_rad.zeta_cut = zeta_cut;
	     

	   AUX_rad.sri_scal_inv = 1.0/sri_scal;
	   if(isnan(AUX_rad.sri_scal_inv))
	     {
	       printf("Inside primitive rad, AUX_rad.sri_scal_inv is NAN!!! sri_scal = %e \n", sri_scal);
	     }
	   
	   
	   AUX_rad.tau_rad_scal_inv = 1.0/tau_rad_scal;
	  
	   UU_radguess[0] = 0.0;
	   UU_radguess[1] = UU_rad[1];
	   UU_radguess[2] = UU_rad[2];
	   UU_radguess[3] = UU_rad[3];
	   UU_radguess[4] = UU_rad[4];
	   
	     
	   bool check;
	   int newt_rad_indx[5];
	   double newt_rad_g[5],newt_rad_p[5],newt_rad_xold[5],newt_rad_fvec[5];	    
	   
	   newt2_cpp_rad(UU_rad,AUX_rad,function_rad,jacobian_rad,nn,check, 
			 newt_rad_indx,newt_rad_g,newt_rad_p,newt_rad_xold,newt_rad_fvec,MAXITS_beforefontfix,STPMX_beforefontfix);
	   
	   int MAXITS_taustildefix=MAXITS_beforefontfix*10; //  <-- Try harder when doing a Font fix; we don't want to fail here!!!
	   double STPMX_taustildefix=STPMX_beforefontfix*10; //  <-- Try harder when doing a Font fix; we don't want to fail here!!!
	   
	   
	   while(check || UU_rad[4] < 0.0) {
	     check = false;
	     newt2_cpp_rad(UU_rad,AUX_rad,function_rad,jacobian_rad,nn,check,
			   newt_rad_indx,newt_rad_g,newt_rad_p,newt_rad_xold,newt_rad_fvec,MAXITS_taustildefix,STPMX_taustildefix);
	     MAXITS_taustildefix*=10;  // <-- Try even harder!  Don't give up!
	     //if(MAXITS_taustildefix>=MAXITS_beforefontfix*1e6) break;
	     if(MAXITS_taustildefix>=MAXITS_beforefontfix*1e4) break;
	   }
	   
	   
	     //****************************************************************
	     //                          FONT FIX
	     // Impose Font fix when Newt-Raph fails (check = true) or when it
	     //  gives negative E_rad [UU_rad[4] < 0]
	     //****************************************************************
	     if(check || UU_rad[4] < 0.0) {

	       //	       failure_tracker_rad[index] = 1.0;
	       //	       repairs_rad_needed = 1;   
	     
	       printf("Problem at (x,y,z) = %e %e %e, %d %d %d\n", X[index],Y[index],Z[index],i,j,k);
	       printf("tau_rad = %e\n",tau_rad[index]);
	       printf("S_rad_x = %e\n",S_rad_x[index]);
	       printf("S_rad_y = %e\n",S_rad_y[index]);
	       printf("S_rad_z = %e\n",S_rad_z[index]);
	       printf("gamma^xx = %e\n",AUX_rad.gupxx_phys);
	       printf("gamma^xy = %e\n",AUX_rad.gupxy_phys);
	       printf("gamma^xz = %e\n",AUX_rad.gupxz_phys);
	       printf("gamma^yy = %e\n",AUX_rad.gupyy_phys);
	       printf("gamma^yz = %e\n",AUX_rad.gupyz_phys);
	       printf("gamma^zz = %e\n",AUX_rad.gupzz_phys);
	       printf("exp(6 phi) = %e\n",AUX_rad.Psi6);
	       exit(0);



	       double UU_rad_font_fix[5];
	       UU_rad_font_fix[0] = 0.0;
	       UU_rad_font_fix[1] = UU_radguess[1];
	       UU_rad_font_fix[2] = UU_radguess[2];
	       UU_rad_font_fix[3] = UU_radguess[3];
	       UU_rad_font_fix[4] = 0.0;

	       nn = 3;
	         
	       check = false;
	       int MAXITS_fontfix=MAXITS_beforefontfix*10; //  <-- Try harder when doing a Font fix; we don't want to fail here!!!
	       double STPMX_fontfix=STPMX_beforefontfix*10; //  <-- Try harder when doing a Font fix; we don't want to fail here!!!
	       newt2_cpp_rad(UU_rad_font_fix,AUX_rad,function_rad_font_fix,jacobian_rad_font_fix,nn,check,
			 newt_rad_indx,newt_rad_g,newt_rad_p,newt_rad_xold,newt_rad_fvec,MAXITS_fontfix,STPMX_fontfix);
	       while(check) {
		 check = false;
		 newt2_cpp_rad(UU_rad_font_fix,AUX_rad,function_rad_font_fix,jacobian_rad_font_fix,nn,check,
			   newt_rad_indx,newt_rad_g,newt_rad_p,newt_rad_xold,newt_rad_fvec,MAXITS_fontfix,STPMX_fontfix);
		 MAXITS_fontfix*=10;  // <-- Try even harder!  Don't give up!
		 if(MAXITS_fontfix>=MAXITS_beforefontfix*1e4) break;
		 //if(MAXITS_fontfix>=MAXITS_beforefontfix*1e4) break;
	       }
	       if(check) {
	         printf("ERROR: FONT FIX (secondary solver) JUST FAILED\n");
	         printf("Problem at (x,y,z) = %e %e %e, %d %d %d\n", X[index],Y[index],Z[index],i,j,k);
	         printf("tau_rad = %e\n",tau_rad[index]);
	         printf("S_rad_x = %e\n",S_rad_x[index]);
	         printf("S_rad_y = %e\n",S_rad_y[index]);
	         printf("S_rad_z = %e\n",S_rad_z[index]);
	         printf("gamma^xx = %e\n",AUX_rad.gupxx_phys);
	         printf("gamma^xy = %e\n",AUX_rad.gupxy_phys);
	         printf("gamma^xz = %e\n",AUX_rad.gupxz_phys);
	         printf("gamma^yy = %e\n",AUX_rad.gupyy_phys);
	         printf("gamma^yz = %e\n",AUX_rad.gupyz_phys);
	         printf("gamma^zz = %e\n",AUX_rad.gupzz_phys);
	         printf("exp(6 phi) = %e\n",AUX_rad.Psi6);
		 
		 failure_tracker_rad[index] = 1.0;
		 repairs_rad_needed = 1;
	      
	         // Set everything to 0.0 before calculating 
	         // these quantities from averages.
	         F_rad_x         = 0.0;
	         F_rad_y         = 0.0;
	         F_rad_z         = 0.0;
		 E_radl          = 0.0;
		 recom_rad        = true;
	       } else {
	         //if the Font fix worked, do the following: 
	         F_rad_x = UU_rad_font_fix[1]*(SQR(UU_rad_font_fix[1]) + 1.0);
	         F_rad_y = UU_rad_font_fix[2]*(SQR(UU_rad_font_fix[2]) + 1.0);
	         F_rad_z = UU_rad_font_fix[3]*(SQR(UU_rad_font_fix[3]) + 1.0);
	         
		 recom_rad = true;    
	       }	     
	     } else {
	       // Inversion worked without the Font fix!  Now set the primtives
       	       F_rad_x = UU_rad[1];
	       F_rad_y = UU_rad[2];
	       F_rad_z = UU_rad[3];
	       E_radl  = UU_rad[4];
	     }
	



	     //check drastic change of UU_rad:
	     double ratio_1 = fabs(UU_rad[1]/UU_radguess[1]);
	     double ratio_2 = UU_rad[2]/UU_radguess[2];
	     double ratio_3 = UU_rad[3]/UU_radguess[3];
	     double ratio_4 = UU_rad[4]/UU_radguess[4];

	     double check_1 = fabs(1.0 - ratio_1);
	     double check_2 = fabs(1.0 - ratio_2);
	     double check_3 = fabs(1.0 - ratio_3);
	     double check_4 = fabs(1.0 - ratio_4);

	   
	     double F_radxl, F_radyl, F_radzl;

	     ////////////////////////Add radiation terms to the source.////////////////////////////
	     //Apply atm floor of rad primitives:  	  
	     //	     if (E_radl < Erad_atm_cut)
	     if (E_radl < Erad_atm_cut) 
               {                    
		 recom_rad = true;
		 E_radl = Erad_atm_cut;
		 F_rad_x = 0.0;
		 F_rad_y = 0.0;
		 F_rad_z = 0.0;
	       }
	         
             //Apply Psi6threshold fix on E_radl:                                                                                                         
             if(AUX_rad.Psi6 > Psi6threshold) {
	       recom_rad = true;
               double Erad_horiz_cap = 2.0*Erad_atm_cut;
               double Frad_horiz_cap = 0.0;
               if(E_radl >= Erad_horiz_cap){
                 E_radl = Erad_horiz_cap;
               }
               if(abs(F_rad_x) >= Frad_horiz_cap){
                 F_rad_x = Frad_horiz_cap;
               }
               if(abs(F_rad_y) >= Frad_horiz_cap){
                 F_rad_y = Frad_horiz_cap;
               }
               if(abs(F_rad_z) >= Frad_horiz_cap){
                 F_rad_z = Frad_horiz_cap;
               }
             }


             F_rad_0l = - (F_rad_x * AUX_rad.ux + F_rad_y * AUX_rad.uy + F_rad_z * AUX_rad.uz)/u0L;
             F_rad0l = (-F_rad_0l+F_rad_x*AUX_rad.shiftx+F_rad_y*AUX_rad.shifty+F_rad_z*AUX_rad.shiftz)/SQR(AUX_rad.alpn1);

             double gupxx4 = AUX_rad.gupxx_phys-AUX_rad.shiftx*AUX_rad.shiftx/SQR(AUX_rad.alpn1);
             double gupyy4 = AUX_rad.gupyy_phys-AUX_rad.shifty*AUX_rad.shifty/SQR(AUX_rad.alpn1);
             double gupzz4 = AUX_rad.gupzz_phys-AUX_rad.shiftz*AUX_rad.shiftz/SQR(AUX_rad.alpn1);
             double gupxy4 = AUX_rad.gupxy_phys-AUX_rad.shiftx*AUX_rad.shifty/SQR(AUX_rad.alpn1);
             double gupxz4 = AUX_rad.gupxz_phys-AUX_rad.shiftx*AUX_rad.shiftz/SQR(AUX_rad.alpn1);
             double gupyz4 = AUX_rad.gupyz_phys-AUX_rad.shifty*AUX_rad.shiftz/SQR(AUX_rad.alpn1);

             F_radxl = F_rad_0l*AUX_rad.shiftx/SQR(AUX_rad.alpn1) + F_rad_x*gupxx4 + F_rad_y*gupxy4 +F_rad_z*gupxz4;
             F_radyl = F_rad_0l*AUX_rad.shifty/SQR(AUX_rad.alpn1) + F_rad_x*gupxy4 + F_rad_y*gupyy4 +F_rad_z*gupyz4;
             F_radzl = F_rad_0l*AUX_rad.shiftz/SQR(AUX_rad.alpn1) + F_rad_x*gupxz4 + F_rad_y*gupyz4 +F_rad_z*gupzz4;	     
	     
	     // Limit E_rad & Flux                                                                                                  
	     if (E_radl > rho_s)
	       {
		 recom_rad = true;
		 E_radl = rho_s;
	       }
             FaFa = F_radxl*F_rad_x + F_radyl*F_rad_y + F_radzl*F_rad_z + F_rad0l*F_rad_0l;
             if (FaFa > SQR(E_radl)){
               recom_rad = true;
	       fac_rad = sqrt(SQR(E_radl)/FaFa);
               F_rad_x *= fac_rad;
               F_rad_y *= fac_rad;
               F_rad_z *= fac_rad;
               F_radxl *= fac_rad;
               F_radyl *= fac_rad;
               F_radzl *= fac_rad;
               F_rad_0l = - (F_rad_x * AUX_rad.ux + F_rad_y * AUX_rad.uy + F_rad_z * AUX_rad.uz)/u0L;
               F_rad0l = (-F_rad_0l+F_rad_x*AUX_rad.shiftx+F_rad_y*AUX_rad.shifty+F_rad_z*AUX_rad.shiftz)/SQR(AUX_rad.alpn1);
             }

	     F[index] = sqrt(Psi4*(gxxL*SQR(F_radxl)+gyyL*SQR(F_radyl)+gzzL*SQR(F_radzl) + 2.0*( gxyL*F_radxl*F_radyl + gxzL*F_radxl*F_radzl + gyzL*F_radyl*F_radzl))  );

	     double Fasq = F_rad_0l*F_rad0l + F_rad_x*F_radxl +  F_rad_y*F_radyl +  F_rad_z*F_radzl;

	     double zeta;	   
	     double zeta_temp =  sqrt(fabs(Fasq)/SQR(E_radl));

	     if (E_radl <= zeta_cut){
	       zeta = 1.0;
	     }
	     else{
	       zeta = zeta_temp;
	     }

	     if (zeta > 1.0){
	       zeta = 1.0;
	     }	   
	     double chi = 1/3.0 + SQR(zeta)*(6.0-2.0*zeta+6.0*SQR(zeta))/15.0;
	     
	     zeta_rad[index] = zeta;
	     chi_rad[index] = chi;

	     double P_radxxl,P_radyyl,P_radzzl,P_radxyl,P_radxzl,P_radyzl;
	     // Since we use gupij_phys = gupij * psim4, we set psim4 = 1.0 in this function.
	     double psim4_temp = 1.0;
	     compute_M1(P_radxxl, F_radxl, F_radxl, Fasq, E_radl, AUX_rad.gupxx_phys, AUX_rad.shiftx, AUX_rad.shiftx, alpha[index], AUX_rad.ux, AUX_rad.ux, chi, psim4_temp, Erad_atm_cut);
	     compute_M1(P_radyyl, F_radyl, F_radyl, Fasq, E_radl, AUX_rad.gupyy_phys, AUX_rad.shifty, AUX_rad.shifty, alpha[index], AUX_rad.uy, AUX_rad.uy, chi, psim4_temp, Erad_atm_cut);
	     compute_M1(P_radzzl, F_radzl, F_radzl, Fasq, E_radl, AUX_rad.gupzz_phys, AUX_rad.shiftz, AUX_rad.shiftz, alpha[index], AUX_rad.uz, AUX_rad.uz, chi, psim4_temp, Erad_atm_cut);
	     compute_M1(P_radxyl, F_radxl, F_radyl, Fasq, E_radl, AUX_rad.gupxy_phys, AUX_rad.shiftx, AUX_rad.shifty, alpha[index], AUX_rad.ux, AUX_rad.uy, chi, psim4_temp, Erad_atm_cut);
	     compute_M1(P_radxzl, F_radxl, F_radzl, Fasq, E_radl, AUX_rad.gupxz_phys, AUX_rad.shiftx, AUX_rad.shiftz, alpha[index], AUX_rad.ux, AUX_rad.uz, chi, psim4_temp, Erad_atm_cut);	       
             compute_M1(P_radyzl, F_radyl, F_radzl, Fasq, E_radl, AUX_rad.gupyz_phys, AUX_rad.shifty, AUX_rad.shiftz, alpha[index], AUX_rad.uy, AUX_rad.uz, chi, psim4_temp, Erad_atm_cut);

	 
	     E_rad[index] = E_radl;
	     F_radx[index] = F_radxl;
	     F_rady[index] = F_radyl;
	     F_radz[index] = F_radzl;
	     F_rad0[index] = F_rad0l;
	     
	     P_radxx[index] = P_radxxl;
	     P_radyy[index] = P_radyyl;
	     P_radzz[index] = P_radzzl;
	     P_radxy[index] = P_radxyl;
	     P_radxz[index] = P_radxzl;
	     P_radyz[index] = P_radyzl;
	   
	   
	   double P_rad0xl = - (P_radxxl * AUX_rad.u_x + P_radxyl * AUX_rad.u_y + P_radxzl * AUX_rad.u_z)/AUX_rad.u_0;
	   double P_rad0yl = - (P_radxyl * AUX_rad.u_x + P_radyyl * AUX_rad.u_y + P_radyzl * AUX_rad.u_z)/AUX_rad.u_0;
	   double P_rad0zl = - (P_radxzl * AUX_rad.u_x + P_radyzl * AUX_rad.u_y + P_radzzl * AUX_rad.u_z)/AUX_rad.u_0;
	   
	   double P_rad00l = - (P_rad0xl * AUX_rad.u_x + P_rad0yl * AUX_rad.u_y + P_rad0zl * AUX_rad.u_z)/AUX_rad.u_0;
	   
	   
	   
	   
	   rho[index] = rho[index] + SQR(alpn1)*(E_radl*SQR(u0L) + 2.0*F_rad0l*u0L + P_rad00l);
	   
	   Sx[index] = Sx[index] + alpn1*(u0L*E_radl*AUX_rad.u_x + F_rad0l*AUX_rad.u_x + u0L*F_rad_x + 
					  P_rad00l*shift_xL + Psi4*(P_rad0xl*gxxL + P_rad0yl*gxyL + P_rad0zl*gxzL));
	   
	   Sy[index] = Sy[index] + alpn1*(u0L*E_radl*AUX_rad.u_y + F_rad0l*AUX_rad.u_y + u0L*F_rad_y +
					  P_rad00l*shift_yL + Psi4*(P_rad0xl*gxyL + P_rad0yl*gyyL + P_rad0zl*gyzL));
	   
	   Sz[index] = Sz[index] + alpn1*(u0L*E_radl*AUX_rad.u_z + F_rad0l*AUX_rad.u_z + u0L*F_rad_z +
					  P_rad00l*shift_zL + Psi4*(P_rad0xl*gxzL + P_rad0yl*gyzL + P_rad0zl*gzzL));
	   
	   
	   
	   Sxx[index] = Sxx[index] + E_radl*AUX_rad.u_x*AUX_rad.u_x + 2*F_rad_x*AUX_rad.u_x + 
	     SQR(shift_xL)*P_rad00l + shift_xL*2.0*(gxxL*P_rad0xl+gxyL*P_rad0yl+gxzL*P_rad0zl) + 
	     SQR(Psi4)*( SQR(gxxL)*P_radxxl + SQR(gxyL)*P_radyyl + SQR(gxzL)*P_radzzl +
			 2.0*(gxxL*gxyL*P_radxyl+gxxL*gxzL*P_radxzl+gxyL*gxzL*P_radyzl) );
	   
	   Syy[index] = Syy[index] + E_radl*AUX_rad.u_y*AUX_rad.u_y + 2*F_rad_y*AUX_rad.u_y +
	     SQR(shift_yL)*P_rad00l + shift_yL*2.0*(gxyL*P_rad0xl+gyyL*P_rad0yl+gyzL*P_rad0zl) + 
	     SQR(Psi4)*( SQR(gxyL)*P_radxxl + SQR(gyyL)*P_radyyl + SQR(gyzL)*P_radzzl +
			 2.0*(gxyL*gyyL*P_radxyl+gxyL*gyzL*P_radxzl+gyyL*gyzL*P_radyzl) );
	   
	   Szz[index] = Szz[index] + E_radl*AUX_rad.u_z*AUX_rad.u_z + 2*F_rad_z*AUX_rad.u_z +
	     SQR(shift_zL)*P_rad00l + shift_zL*2.0*(gxzL*P_rad0xl+gyzL*P_rad0yl+gzzL*P_rad0zl) + 
	     SQR(Psi4)*( SQR(gxzL)*P_radxxl + SQR(gyzL)*P_radyyl + SQR(gzzL)*P_radzzl +
			 2.0*(gxzL*gyzL*P_radxyl+gxzL*gzzL*P_radxzl+gyzL*gzzL*P_radyzl) );
	   
	   Sxy[index] = Sxy[index] + E_radl*AUX_rad.u_x*AUX_rad.u_y + F_rad_x*AUX_rad.u_y + F_rad_y*AUX_rad.u_x +
	     shift_xL*shift_yL*P_rad00l + 
	     Psi4*(shift_xL*(gxyL*P_rad0xl + gyyL*P_rad0yl + gyzL*P_rad0zl) + (shift_yL*(gxxL*P_rad0xl + gxyL*P_rad0yl + gxzL*P_rad0zl)) )+
	     SQR(Psi4)*(gxxL*gxyL*P_radxxl + gxyL*gyyL*P_radyyl + gxzL*gyzL*P_radzzl +
			(gxxL*gyyL + gxyL*gxyL)*P_radxyl + (gxxL*gyzL + gxzL*gxyL)*P_radxzl + (gxyL*gyzL + gxzL*gyyL)*P_radyzl);
	   
	   Sxz[index] = Sxz[index] + E_radl*AUX_rad.u_x*AUX_rad.u_z + F_rad_x*AUX_rad.u_z + F_rad_z*AUX_rad.u_x +
	     shift_xL*shift_zL*P_rad00l +
	     Psi4*(shift_xL*(gxzL*P_rad0xl + gyzL*P_rad0yl + gzzL*P_rad0zl) + (shift_zL*(gxxL*P_rad0xl + gxyL*P_rad0yl + gxzL*P_rad0zl)) )+
	     SQR(Psi4)*(gxxL*gxzL*P_radxxl + gxyL*gyzL*P_radyyl + gxzL*gzzL*P_radzzl +
			(gxxL*gyzL + gxyL*gxzL)*P_radxyl + (gxxL*gzzL + gxzL*gxzL)*P_radxzl + (gxyL*gzzL + gxzL*gyzL)*P_radyzl);
	     
	   Syz[index] = Syz[index] + E_radl*AUX_rad.u_y*AUX_rad.u_z + F_rad_y*AUX_rad.u_z + F_rad_z*AUX_rad.u_y +
	     shift_yL*shift_zL*P_rad00l +
	     Psi4*(shift_yL*(gxzL*P_rad0xl + gyzL*P_rad0yl + gzzL*P_rad0zl) + (shift_zL*(gxyL*P_rad0xl + gyyL*P_rad0yl + gyzL*P_rad0zl)) )+
	     SQR(Psi4)*(gxyL*gxzL*P_radxxl + gyyL*gyzL*P_radyyl + gyzL*gzzL*P_radzzl +
			(gxyL*gyzL + gxzL*gyyL)*P_radxyl + (gxyL*gzzL + gyzL*gxzL)*P_radxzl + (gyyL*gzzL + gyzL*gyzL)*P_radyzl);
	        	   
	   
	   if (recom_rad){
	     tau_rad[index] = alpn1*alpn1*Psi6*(E_radl*u0L*u0L+2.0*F_rad0l*u0L+P_rad00l);
	     S_rad_x[index] = alpn1*Psi6*(E_radl*u0L*AUX_rad.u_x + F_rad0l*AUX_rad.u_x + F_rad_x * u0L + P_rad00l*shift_xL + Psi4*(P_rad0xl*gxxL + P_rad0yl*gxyL + P_rad0zl*gxzL));
	     S_rad_y[index] = alpn1*Psi6*(E_radl*u0L*AUX_rad.u_y + F_rad0l*AUX_rad.u_y + F_rad_y * u0L + P_rad00l*shift_yL + Psi4*(P_rad0xl*gxyL + P_rad0yl*gyyL + P_rad0zl*gyzL));
	     S_rad_z[index] = alpn1*Psi6*(E_radl*u0L*AUX_rad.u_z + F_rad0l*AUX_rad.u_z + F_rad_z * u0L + P_rad00l*shift_zL + Psi4*(P_rad0xl*gxzL + P_rad0yl*gyzL + P_rad0zl*gzzL));


	     // CHECK ABNORMAL VALUES!!!                                                                                                                                         
	     if (abs(F_radx[index]) > 100.0 || abs(F_rady[index]) > 100.0 || abs(F_radz[index]) > 100.0)
               {
                 printf ("Inside primitives_rad.C, AFTER newt2 F_rad_x, F_rad_y, F_rad_z = %e,%e,%e \n", F_rad_x, F_rad_y, F_rad_z);
                 printf ("E_rad[index], tau_rad[index], S_rad_x[index], S_rad_y[index], S_rad_z[index] = %e,%e,%e,%e,%e \n", E_rad[index], tau_rad[index], S_rad_x[index], S_rad_y[index], S_rad_z[index]);
                 printf ("F_radx[index], F_rady[index], F_radz[index], F_rad0[index], FaFa = %e,%e,%e,%e,%e \n", F_radx[index], F_rady[index], F_radz[index], F_rad0[index], FaFa);
                 printf ("gxxL,gxyL,gxzL,gyyL,gxzL,gyzL,shift_xL,shift_yL,shift_zL,u0[index],uxL,uyL,uzL=%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e,%e \n", gxxL,gxyL,gxzL,gyyL,gxzL,gyzL,shift_xL,shift_yL,shift_zL,u0[index],uxL,uyL,uzL);
                 printf ("rho_star[index], vx[index], vy[index], vz[index] = %e,%e,%e,%e \n", rho_star[index], vx[index], vy[index], vz[index]);
                 printf ("X,Y,Z=%e,%e,%e \n", X[index], Y[index], Z[index]);
               }


	   }
	}
      }
  
  
  gettimeofday(&end, NULL);
  
  seconds  = end.tv_sec  - start.tv_sec;
  useconds = end.tv_usec - start.tv_usec;

  mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;  // We add 0.5 since mtime is a long int; this rounds up the result before setting the value.  Here, rounding down is incorrect.


  printf("Primitives radiation solver: %f solutions/second\n",(iend-istart)*(jend-jstart)*(kend-kstart) / ((double)mtime/1000.0));

}



double max_val_rad(double val1, double val2) {
  if(val1>val2) return val1;
  return val2;
}

#define NP 5
//#define MAXITS 50
#define TOLF 1.0e-12
#define TOLMIN 1.0e-14
#define TOLX 3.0e-16

//#define STPMX 10000.0
//#define STPMX 100.0
#define FREERETURN {return;}
//#define FREERETURN {free_dmatrix_newt(fjac,1,NP,1,NP);return;}

void newt2_cpp_rad(double x[5],struct auxarray_rad &aux,
	       void (*funcv2)(int &n,double *x,double *fvec,struct auxarray_rad &aux),
	       void (*fdjac2)
	       (int &n,double *x,struct auxarray_rad &aux, double *fvec,int &np,double fjac[][5]),
	       int &n,bool &check,
	       int *indx,double *g,double *p,double *xold,double *fvec,int &MAXITS,double &STPMX)
{
  void lubksb_rad(double a[][5], int n, int *indx, double b[]);
  void ludcmp_rad(double a[][5], int n, int *indx, double *d);
  int i,its,j;
  double d,den,f,fold,stpmax,sum,temp,test,fjac[5][5];


  //printf ("INSIDE newt2_cpp_rad, before fmin_rad, fvec[1]=%e, fvec[2]=%e,fvec[3]=%e, fvec[4]=%e, ux=x[1]=%e \n", fvec[1], fvec[2],fvec[3], fvec[4], x[1]);
 
  f=fmin_rad(x,aux,fvec,funcv2,n);

  // printf ("INSIDE newt2_cpp_rad, after fmin_rad, fvec[1]=%e, fvec[2]=%e,fvec[3]=%e, fvec[4]=%e, ux=x[1]=%e \n", fvec[1], fvec[2],fvec[3], fvec[4], x[1]);

  test=0.0;
  for (i=1;i<=n;i++)
    if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
 
  //printf("test=%.16g, 0.01*TOLF=%.16g !!!!!!\n", test, 0.01*TOLF);
  if (test < 0.01*TOLF) {
    check=0;
    FREERETURN
      }

  //printf("test > 0.01*TOLF JACOBIAN!!!!!!\n");
  for (sum=0.0,i=1;i<=n;i++) sum += SQR(x[i]);
  stpmax=STPMX*max_val_rad(sqrt(sum),(double)n);
  for (its=1;its<=MAXITS;its++) {
    int np=NP;
    fdjac2(n,x,aux,fvec,np,fjac);

    for (i=1;i<=n;i++) {
      for (sum=0.0,j=1;j<=n;j++) sum += fjac[j][i]*fvec[j];
      g[i]=sum;
    }

    for (i=1;i<=n;i++) xold[i]=x[i];
    fold=f;
    for (i=1;i<=n;i++) p[i] = -fvec[i];

    // The two following lines solve the equation fjac * x = p, where p[i] = -fvec[i] = -\delta f(i).
    // After two functions are called, the vector p(i) is replaced by the new solution x(i), used as the Newton directon of new solution. 
    ludcmp_rad(fjac,n,indx,&d);
    lubksb_rad(fjac,n,indx,p);
    
    
    lnsrch_rad(n,xold,fold,g,p,x,&f,stpmax,check,fmin_rad,funcv2,aux,fvec);

    test=0.0;
    for (i=1;i<=n;i++)
      if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
    if (test < TOLF) {
      check=0;
      FREERETURN
	}
    if (check) {
      test=0.0;
      den=max_val_rad(f,0.5*n);
      for (i=1;i<=n;i++) {
	temp=fabs(g[i])*max_val_rad(fabs(x[i]),1.0)/den;
	if (temp > test) test=temp;
      }
      check=(test < TOLMIN ? 1 : 0);
      //if(check!=0) printf("BAD CHECK %d: %e\t%e\t%e\t%e\t%e\t%e\n",check,test,TOLMIN,fabs(x[1]),fabs(x[2]),fabs(x[3]),fabs(x[4]) );
      FREERETURN
	}
    test=0.0;
    for (i=1;i<=n;i++) {
      temp=(fabs(x[i]-xold[i]))/max_val_rad(fabs(x[i]),1.0);
      if (temp > test) test=temp;
    }
    if (test < TOLX) FREERETURN
		       }
  check=1;
  //if(check!=0) printf("BAD STUFF %d: %e\t%e.  ITS=%d\n",check,test,TOLMIN,its);
  //printf("MAXITS exceeded in newtn\n");
  //exit(0);
}

double fmin_rad(double *x,struct auxarray_rad &aux,double *fvec,
		 void (*nrfuncv_rad)(int &n,double *x,double *fvec,struct auxarray_rad &aux),
		 int &n)
{
  int i;
  double sum;

  (*nrfuncv_rad)(n,x,fvec,aux);
  for (sum=0.0,i=1;i<=n;i++) sum += SQR(fvec[i]);
  return 0.5*sum;
}

#define ALF 1.0e-4
//#define TOLX 1.0e-7
#define TOLX2 1.0e-16
//#define TOLX 1.0e-15

void lnsrch_rad(int n, double xold[], double fold, double g[], double p[], double x[5],
		 double *f, double stpmax, bool &check, 
		 double (*fminrad)
		 (double *x,struct auxarray_rad &aux,double *fvec,
		  void (*nrfuncv_rad)(int &n,double *x,double *fvec,struct auxarray_rad &aux),
		  int &n),
		 void (*nrfuncv_rad)(int &n,double *x,double *fvec,struct auxarray_rad &aux),
		 struct auxarray_rad &aux,double *fvec)
{
  int i;
  double a,alam,alam2,alamin,b,disc,f2,rhs1,rhs2,slope,sum,temp,
    test,tmplam;

  check=0;
  for (sum=0.0,i=1;i<=n;i++) sum += p[i]*p[i];
  sum=sqrt(sum);
  if (sum > stpmax)
    for (i=1;i<=n;i++) p[i] *= stpmax/sum;
  for (slope=0.0,i=1;i<=n;i++)
    slope += g[i]*p[i];
  if (slope >= 0.0) {printf("OUCH!  Roundoff problem in lnsrch.  Note that this function has been updated to correct a bug in Numerical Recipes 2.06.  This is version 2.08, with the 2.06->2.08 diff obtained from http://www.numerical-recipes.com/upgrade/upgrade-208.html  x[1-4] = %e,%e,%e,%e \n", xold[1],xold[2],xold[3],xold[4]);
    printf("g[1-4] = %e,%e,%e,%e, p[1-4] = %e,%e,%e,%e \n", g[1],g[2],g[3],g[4],p[1],p[2],p[3],p[4]);
  }
  test=0.0;
  for (i=1;i<=n;i++) {
    temp=fabs(p[i])/max_val_rad(fabs(xold[i]),1.0);
    if (temp > test) test=temp;
  }
  alamin=TOLX2/test;
  alam=1.0;
  for (;;) {
    for (i=1;i<=n;i++) x[i]=xold[i]+alam*p[i];
    *f=(*fminrad)(x,aux,fvec,nrfuncv_rad,n);
    if (alam < alamin) {
      for (i=1;i<=n;i++) x[i]=xold[i];
      check=1;
      return;
    } else if (*f <= fold+ALF*alam*slope) return;
    else {
      if (alam == 1.0)
        tmplam = -slope/(2.0*(*f-fold-slope));
      else {
        rhs1 = *f-fold-alam*slope;
        rhs2=f2-fold-alam2*slope;
        a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
        b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
        if (a == 0.0) tmplam = -slope/(2.0*b);
        else {
          disc=b*b-3.0*a*slope;
          if (disc < 0.0) tmplam=0.5*alam;
          else if (b <= 0.0) tmplam=(-b+sqrt(disc))/(3.0*a);
          else tmplam=-slope/(b+sqrt(disc));
	}
        if (tmplam > 0.5*alam)
          tmplam=0.5*alam;
      }
    }
    alam2=alam;
    f2 = *f;
    alam=max_val_rad(tmplam,0.1*alam);
    //    alam=FMAX(tmplam,0.1*alam);
  }

}
#undef ALF
#undef TOLX2

void lubksb_rad(double a[][5], int n, int *indx, double b[])
{
  int i,ii=0,ip,j;
  double sum;

  for (i=1;i<=n;i++) {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    //printf("at the beginning of lubksb, ip=%d, sum=%e, b[ip]=%e \n", ip, sum, b[ip]);
    if (ii)
      for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
    else if (sum) ii=i;
    b[i]=sum;
  }
  //printf("In the middle of lubksb, b[1]=%e, b[2]=%e, b[3]=%e, b[4]=%e \n", b[1], b[2], b[3], b[4]);
  for (i=n;i>=1;i--) {
    //printf ("inside lubksb_rad n=%d", n);
    sum=b[i];
    for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
  
}

#define TINY 1.0e-20;
//#define TINY 1.0e-20;
void ludcmp_rad(double a[][5], int n, int *indx, double *d)
{
  int i,imax,j,k;
  double big,dum,sum,temp;
  double vv[5];

  *d=1.0;
  for (i=1;i<=n;i++) {
    big=0.0;
    for (j=1;j<=n;j++)
      if ((temp=fabs(a[i][j])) > big) big=temp;
    if (big == 0.0) {
	for(int www=0;www<10;www++) {
        printf("Singular matrix in routine ludcmp_rad\n"); 
        for (j=1;j<=n;j++) printf("i, j, a[i][j]: %d %d %e\n",i,j,a[i][j]);
	}
        exit(1);
    }
    vv[i]=1.0/big;
  }
  for (j=1;j<=n;j++) {
    for (i=1;i<j;i++) {
      sum=a[i][j];
      for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for (i=j;i<=n;i++) {
      sum=a[i][j];
      for (k=1;k<j;k++)
	sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big) {
	big=dum;
	imax=i;
      }
    }
    if (j != imax) {
      for (k=1;k<=n;k++) {
	dum=a[imax][k];
	a[imax][k]=a[j][k];
	a[j][k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j] == 0.0) a[j][j]=TINY;
    if (j != n) {
      dum=1.0/(a[j][j]);
      for (i=j+1;i<=n;i++) a[i][j] *= dum;
    }
  }
}
#undef TINY

void function_rad(int &n,double *x,double *fvec,struct auxarray_rad &aux){
  double F_rad_x    = x[1];
  double F_rad_y    = x[2];
  double F_rad_z    = x[3];
  double E_rad     = x[4];


  double gijuiuj = aux.gupxx_phys*SQR(aux.u_x) + 2.0*aux.gupxy_phys*aux.u_x*aux.u_y +
    2.0*aux.gupxz_phys*aux.u_x*aux.u_z + aux.gupyy_phys*SQR(aux.u_y) + 2.0*aux.gupyz_phys*aux.u_y*aux.u_z +
    aux.gupzz_phys*SQR(aux.u_z);
  double au0m1 = gijuiuj/( 1.0+sqrt(1.0+gijuiuj) );
  if (aux.rho_s < 0.0) au0m1 = gijuiuj/( 1.0-sqrt(1.0+gijuiuj) );
  double u0 = (au0m1+1.0)/aux.alpn1;

  double F_rad_0 = - (F_rad_x * aux.ux + F_rad_y * aux.uy + F_rad_z * aux.uz)/u0; 
  double F_rad0 = (-F_rad_0+F_rad_x*aux.shiftx+F_rad_y*aux.shifty+F_rad_z*aux.shiftz)/SQR(aux.alpn1);

  double gupxx4 = aux.gupxx_phys-aux.shiftx*aux.shiftx/SQR(aux.alpn1);
  double gupyy4 = aux.gupyy_phys-aux.shifty*aux.shifty/SQR(aux.alpn1);
  double gupzz4 = aux.gupzz_phys-aux.shiftz*aux.shiftz/SQR(aux.alpn1);
  double gupxy4 = aux.gupxy_phys-aux.shiftx*aux.shifty/SQR(aux.alpn1);
  double gupxz4 = aux.gupxz_phys-aux.shiftx*aux.shiftz/SQR(aux.alpn1);
  double gupyz4 = aux.gupyz_phys-aux.shifty*aux.shiftz/SQR(aux.alpn1);

  double F_radx = F_rad_0*aux.shiftx/SQR(aux.alpn1) + F_rad_x*gupxx4 + F_rad_y*gupxy4 +F_rad_z*gupxz4;
  double F_rady = F_rad_0*aux.shifty/SQR(aux.alpn1) + F_rad_x*gupxy4 + F_rad_y*gupyy4 +F_rad_z*gupyz4;
  double F_radz = F_rad_0*aux.shiftz/SQR(aux.alpn1) + F_rad_x*gupxz4 + F_rad_y*gupyz4 +F_rad_z*gupzz4;

  double FaFa = F_radx*F_rad_x + F_rady*F_rad_y + F_radz*F_rad_z + F_rad0*F_rad_0; // F_rad^alpha*F_rad_alpha

  double FaFa_o_E2;
  
  if (E_rad <= aux.zeta_cut){
      FaFa_o_E2 = 1.0;
  }
  else{  
    FaFa_o_E2= fabs(FaFa)/SQR(E_rad);
  }
  
  if (FaFa_o_E2 > 1.0){
    FaFa_o_E2 = 1.0;
  }
  
  double A = FaFa_o_E2*(3.0 - sqrt(FaFa_o_E2) + 3.0*FaFa_o_E2)/5.0;

  //  double A = 0.0;

  double FaFa_inv;
  if (FaFa <= SQR(aux.zeta_cut) || E_rad <= aux.zeta_cut){
    FaFa_inv == 0.0;
  }else{
    FaFa_inv = 1.0/FaFa;
  }

  //
  // fvec(1): Eq. for S_rad_x; fvec[2]: Eq. for S_rad_y; fvec[3]: Eq. for S_rad_z; 
  // fvec[4]: Eq. for tau_rad.
  //

  fvec[1] = (aux.alpn1*aux.Psi6*(E_rad*u0*aux.u_x + F_rad0*aux.u_x + F_rad_x*u0 + E_rad*(F_rad0*F_rad_x*FaFa_inv*A + u0*aux.u_x/3.0*(1.0-A))) - aux.S_rad_x)*aux.sri_scal_inv;
  fvec[2] = (aux.alpn1*aux.Psi6*(E_rad*u0*aux.u_y + F_rad0*aux.u_y + F_rad_y*u0 + E_rad*(F_rad0*F_rad_y*FaFa_inv*A + u0*aux.u_y/3.0*(1.0-A))) - aux.S_rad_y)*aux.sri_scal_inv;
  fvec[3] = (aux.alpn1*aux.Psi6*(E_rad*u0*aux.u_z + F_rad0*aux.u_z + F_rad_z*u0 + E_rad*(F_rad0*F_rad_z*FaFa_inv*A + u0*aux.u_z/3.0*(1.0-A))) - aux.S_rad_z)*aux.sri_scal_inv;
  fvec[4] = (SQR(aux.alpn1)*aux.Psi6*(E_rad*SQR(u0) + 2.0*F_rad0*u0 + E_rad*( SQR(F_rad0)*FaFa_inv*A +  (-1.0/SQR(aux.alpn1)+SQR(u0))/3.0*(1.0-A)) )- aux.tau_rad)*aux.tau_rad_scal_inv;
} 
  


void jacobian_rad
(int &n,double *x,struct auxarray_rad &aux, double *fvec,int &np, double fjac[][5]){
  double F_rad_x    = x[1];
  double F_rad_y    = x[2];
  double F_rad_z    = x[3];
  double E_rad     = x[4];

  double gijuiuj = aux.gupxx_phys*SQR(aux.u_x) + 2.0*aux.gupxy_phys*aux.u_x*aux.u_y +
    2.0*aux.gupxz_phys*aux.u_x*aux.u_z + aux.gupyy_phys*SQR(aux.u_y) + 2.0*aux.gupyz_phys*aux.u_y*aux.u_z +
    aux.gupzz_phys*SQR(aux.u_z);
  double au0m1 = gijuiuj/( 1.0+sqrt(1.0+gijuiuj) );
  if (aux.rho_s < 0.0) au0m1 = gijuiuj/( 1.0-sqrt(1.0+gijuiuj) );
  double u0 = (au0m1+1.0)/aux.alpn1;
  double u0_inv = aux.alpn1/(au0m1+1.0);

  double F_rad_0 = - (F_rad_x * aux.ux + F_rad_y * aux.uy + F_rad_z * aux.uz)/u0;
  double F_rad0 = (-F_rad_0+F_rad_x*aux.shiftx+F_rad_y*aux.shifty+F_rad_z*aux.shiftz)/SQR(aux.alpn1);

  double gupxx4 = aux.gupxx_phys-aux.shiftx*aux.shiftx/SQR(aux.alpn1);
  double gupyy4 = aux.gupyy_phys-aux.shifty*aux.shifty/SQR(aux.alpn1);
  double gupzz4 = aux.gupzz_phys-aux.shiftz*aux.shiftz/SQR(aux.alpn1);
  double gupxy4 = aux.gupxy_phys-aux.shiftx*aux.shifty/SQR(aux.alpn1);
  double gupxz4 = aux.gupxz_phys-aux.shiftx*aux.shiftz/SQR(aux.alpn1);
  double gupyz4 = aux.gupyz_phys-aux.shifty*aux.shiftz/SQR(aux.alpn1);

  double F_radx = F_rad_0*aux.shiftx/SQR(aux.alpn1) + F_rad_x*gupxx4 + F_rad_y*gupxy4 +F_rad_z*gupxz4;
  double F_rady = F_rad_0*aux.shifty/SQR(aux.alpn1) + F_rad_x*gupxy4 + F_rad_y*gupyy4 +F_rad_z*gupyz4;
  double F_radz = F_rad_0*aux.shiftz/SQR(aux.alpn1) + F_rad_x*gupxz4 + F_rad_y*gupyz4 +F_rad_z*gupzz4;

  double FaFa = F_radx*F_rad_x + F_rady*F_rad_y + F_radz*F_rad_z + F_rad0*F_rad_0; // F_rad^alpha*F_rad_alpha                                  
  double FkFk = F_radx*F_rad_x + F_rady*F_rad_y + F_radz*F_rad_z;

  double zetasq;
  double zeta_temp = sqrt(fabs(FaFa/SQR(E_rad)));

  /*if ( zeta_temp < 1.0e-20){
    if (E_rad > aux.Erad_atm_cut){
      zetasq=0.0;
    }
    else{
      zetasq = 1.0;
    }
  }
  */
  if (E_rad <= aux.zeta_cut){
    zetasq = 1.0;
  }
  else{
    zetasq = SQR(zeta_temp);
  }
  
  if (zetasq > 1.0){
    zetasq = 1.0;
  }
  double zeta = sqrt(zetasq);
  //double zetasq = 0.0;
  //double zeta =0.0;

  double A = zetasq*(3.0 - zeta + 3.0*zetasq)/5.0;
  double B = zeta*(6.0 - 3.0*zeta + 12.0*zetasq)/5.0;

  double vx = aux.ux/u0;
  double vy = aux.uy/u0;
  double vz = aux.uz/u0;

  double g00 = -1.0/SQR(aux.alpn1);
  double gx0 = aux.shiftx/SQR(aux.alpn1);
  double gy0 = aux.shifty/SQR(aux.alpn1);
  double gz0 = aux.shiftz/SQR(aux.alpn1);
  // f(1) = S_rad_x; f(2) = S_rad_y; f(3) = S_rad_z; f(4) = tau_rad;
  // x[1] = F_rad_x; x[2] = F_rad_y; x[3] = F_rad_z; x[4] = E_rad;
  // fjac(i,j) = partial f(i) / partial x(j) 
  //

  /*
  fjac[4][4] = SQR(aux.alpn1)*aux.Psi6*(4.0/3.0*SQR(u0) - 1.0/(3.0*SQR(aux.alpn1)) );

  fjac[1][4] = aux.alpn1*aux.Psi6*(4.0*u0*aux.u_x/3.0);

  fjac[2][4] = aux.alpn1*aux.Psi6*(4.0*u0*aux.u_y/3.0);

  fjac[3][4] = aux.alpn1*aux.Psi6*(4.0*u0*aux.u_z/3.0);


  fjac[4][1] = SQR(aux.alpn1)*aux.Psi6*(2.0*(gx0 + vx*g00)*u0);

  fjac[4][2] = SQR(aux.alpn1)*aux.Psi6*(2.0*(gy0 + vy*g00)*u0);
  
  fjac[4][3] = SQR(aux.alpn1)*aux.Psi6*(2.0*(gz0 + vz*g00)*u0);


  fjac[1][1] = aux.alpn1*aux.Psi6*(u0 + (gx0 + vx*g00)*(aux.u_x));
						
  fjac[2][2] = aux.alpn1*aux.Psi6*(u0 + (gy0 + vy*g00)*(aux.u_y));
				
  fjac[3][3] = aux.alpn1*aux.Psi6*(u0 + (gz0 + vz*g00)*(aux.u_z));
				
							

  fjac[1][2] = aux.alpn1*aux.Psi6*((gx0 + vx*g00)*(aux.u_y));
				
  fjac[1][3] = aux.alpn1*aux.Psi6*((gx0 + vx*g00)*(aux.u_z));

  fjac[2][1] = aux.alpn1*aux.Psi6*((gy0 + vy*g00)*(aux.u_x));
                                  
  fjac[3][1] = aux.alpn1*aux.Psi6*((gz0 + vz*g00)*(aux.u_x));

  fjac[2][3] = aux.alpn1*aux.Psi6*((gy0 + vy*g00)*(aux.u_z));
			       
  fjac[3][2] = aux.alpn1*aux.Psi6*((gz0 + vz*g00)*(aux.u_y));
  */
  
  
  double FaFa_inv;
  if (FaFa <= SQR(aux.zeta_cut) || E_rad <= aux.zeta_cut){
    FaFa_inv == 0.0;
  }else{
    FaFa_inv = 1.0/FaFa;
  }
  
  fjac[4][4] = (SQR(aux.alpn1)*aux.Psi6*(4.0*SQR(u0)/3.0 + (A-B*zeta)*(SQR(F_rad0)*FaFa_inv - (g00+SQR(u0))/3.0)  + g00/3.0))*aux.tau_rad_scal_inv;

  fjac[1][4] = (aux.alpn1*aux.Psi6*(4.0*u0*aux.u_x/3.0 + (A-B*zeta)*( F_rad0*F_rad_x*FaFa_inv - u0*aux.u_x/3.0)))*aux.sri_scal_inv;

  fjac[2][4] = (aux.alpn1*aux.Psi6*(4.0*u0*aux.u_y/3.0 + (A-B*zeta)*( F_rad0*F_rad_y*FaFa_inv - u0*aux.u_y/3.0)))*aux.sri_scal_inv;

  fjac[3][4] = (aux.alpn1*aux.Psi6*(4.0*u0*aux.u_z/3.0 + (A-B*zeta)*( F_rad0*F_rad_z*FaFa_inv - u0*aux.u_z/3.0)))*aux.sri_scal_inv;

  fjac[4][1] = (SQR(aux.alpn1)*aux.Psi6*(2.0*(gx0 + vx*g00)*(u0 + A*E_rad*F_rad0*FaFa_inv) +
					 (F_radx - vx*F_rad0)*( B*sqrt(FaFa_inv)* (SQR(F_rad0)*FaFa_inv - (g00+SQR(u0))/3.0) - A*E_rad*2.0*SQR(F_rad0)*SQR(FaFa_inv)) ) )*aux.tau_rad_scal_inv;

  fjac[4][2] = (SQR(aux.alpn1)*aux.Psi6*(2.0*(gy0 + vy*g00)*(u0 + A*E_rad*F_rad0*FaFa_inv) +
					 (F_rady - vy*F_rad0)*( B*sqrt(FaFa_inv)* (SQR(F_rad0)*FaFa_inv - (g00+SQR(u0))/3.0) - A*E_rad*2.0*SQR(F_rad0)*SQR(FaFa_inv)) ) )*aux.tau_rad_scal_inv;

  fjac[4][3] = (SQR(aux.alpn1)*aux.Psi6*(2.0*(gz0 + vz*g00)*(u0 + A*E_rad*F_rad0*FaFa_inv) +
					 (F_radz - vz*F_rad0)*( B*sqrt(FaFa_inv)* (SQR(F_rad0)*FaFa_inv - (g00+SQR(u0))/3.0) - A*E_rad*2.0*SQR(F_rad0)*SQR(FaFa_inv)) ) )*aux.tau_rad_scal_inv;
     
  fjac[1][1] = (aux.alpn1*aux.Psi6*((u0+A*E_rad*F_rad0*FaFa_inv) + (gx0 + vx*g00)*(aux.u_x + A*E_rad*F_rad_x*FaFa_inv) +
				    (F_radx - vx*F_rad0)*( B*sqrt(FaFa_inv)* (F_rad0*F_rad_x*FaFa_inv - u0*aux.u_x/3.0) - A*E_rad*2.0*F_rad0*F_rad_x*SQR(FaFa_inv))))*aux.sri_scal_inv;

  fjac[2][2] = (aux.alpn1*aux.Psi6*((u0+A*E_rad*F_rad0*FaFa_inv) + (gy0 + vy*g00)*(aux.u_y + A*E_rad*F_rad_y*FaFa_inv) +
				    (F_rady - vy*F_rad0)*( B*sqrt(FaFa_inv)* (F_rad0*F_rad_y*FaFa_inv - u0*aux.u_y/3.0) - A*E_rad*2.0*F_rad0*F_rad_y*SQR(FaFa_inv))))*aux.sri_scal_inv;

  fjac[3][3] = (aux.alpn1*aux.Psi6*((u0+A*E_rad*F_rad0*FaFa_inv) + (gz0 + vz*g00)*(aux.u_z + A*E_rad*F_rad_z*FaFa_inv) +
				    (F_radz - vz*F_rad0)*( B*sqrt(FaFa_inv)* (F_rad0*F_rad_z*FaFa_inv - u0*aux.u_z/3.0) - A*E_rad*2.0*F_rad0*F_rad_z*SQR(FaFa_inv))))*aux.sri_scal_inv;
  
  fjac[2][1] = (aux.alpn1*aux.Psi6*((gx0 + vx*g00)*(aux.u_y + A*E_rad*F_rad_y*FaFa_inv) +
				    (F_radx - vx*F_rad0)*( B*sqrt(FaFa_inv)* (F_rad0*F_rad_y*FaFa_inv - u0*aux.u_y/3.0) - A*E_rad*2.0*F_rad0*F_rad_y*SQR(FaFa_inv))))*aux.sri_scal_inv;
  
  fjac[3][1] = (aux.alpn1*aux.Psi6*((gx0 + vx*g00)*(aux.u_z + A*E_rad*F_rad_z*FaFa_inv) +
				    (F_radx - vx*F_rad0)*( B*sqrt(FaFa_inv)* (F_rad0*F_rad_z*FaFa_inv - u0*aux.u_z/3.0) - A*E_rad*2.0*F_rad0*F_rad_z*SQR(FaFa_inv))))*aux.sri_scal_inv;
  
  fjac[1][2] = (aux.alpn1*aux.Psi6*((gy0 + vy*g00)*(aux.u_x + A*E_rad*F_rad_x*FaFa_inv) +
				    (F_rady - vy*F_rad0)*( B*sqrt(FaFa_inv)* (F_rad0*F_rad_x*FaFa_inv - u0*aux.u_x/3.0) - A*E_rad*2.0*F_rad0*F_rad_x*SQR(FaFa_inv))))*aux.sri_scal_inv;
  
  fjac[1][3] = (aux.alpn1*aux.Psi6*((gz0 + vz*g00)*(aux.u_x + A*E_rad*F_rad_x*FaFa_inv) +
				    (F_radz - vz*F_rad0)*( B*sqrt(FaFa_inv)* (F_rad0*F_rad_x*FaFa_inv - u0*aux.u_x/3.0) - A*E_rad*2.0*F_rad0*F_rad_x*SQR(FaFa_inv))))*aux.sri_scal_inv;
  
  fjac[3][2] = (aux.alpn1*aux.Psi6*((gy0 + vy*g00)*(aux.u_z + A*E_rad*F_rad_z*FaFa_inv) +
				    (F_rady - vy*F_rad0)*( B*sqrt(FaFa_inv)* (F_rad0*F_rad_z*FaFa_inv - u0*aux.u_z/3.0) - A*E_rad*2.0*F_rad0*F_rad_z*SQR(FaFa_inv))))*aux.sri_scal_inv;
  
  fjac[2][3] = (aux.alpn1*aux.Psi6*((gz0 + vz*g00)*(aux.u_y + A*E_rad*F_rad_y*FaFa_inv) +
				    (F_radz - vz*F_rad0)*( B*sqrt(FaFa_inv)* (F_rad0*F_rad_y*FaFa_inv - u0*aux.u_y/3.0) - A*E_rad*2.0*F_rad0*F_rad_y*SQR(FaFa_inv))))*aux.sri_scal_inv;
  
  if(isnan(fjac[1][1])) printf("Inside jacobian_rad, fjac[1][1] is nan. aux.sri_scal_inv, FaFa_inv, A, F_rad_x, F_rad0, E_rad, F_radx, F_rady, F_radz  = %e, %e, %e, %e,%e, %e, %e, %e, %e\n", aux.sri_scal_inv, FaFa_inv, A, F_rad_x, F_rad0, E_rad, F_radx, F_rady, F_radz);

}



void function_rad_font_fix(int &n,double *x,double *fvec,struct auxarray_rad &aux){

    double F_rad_x = x[1]*(SQR(x[1]) + 1.0);
    double F_rad_y = x[2]*(SQR(x[2]) + 1.0);
    double F_rad_z = x[3]*(SQR(x[3]) + 1.0);
    //    double E_rad = x[4];


    double gijuiuj = aux.gupxx_phys*SQR(aux.u_x) + 2.0*aux.gupxy_phys*aux.u_x*aux.u_y +
      2.0*aux.gupxz_phys*aux.u_x*aux.u_z + aux.gupyy_phys*SQR(aux.u_y) + 2.0*aux.gupyz_phys*aux.u_y*aux.u_z +
      aux.gupzz_phys*SQR(aux.u_z);
    double au0m1 = gijuiuj/( 1.0+sqrt(1.0+gijuiuj) );
    if (aux.rho_s < 0.0) au0m1 = gijuiuj/( 1.0-sqrt(1.0+gijuiuj) );
    double u0 = (au0m1+1.0)/aux.alpn1;

    double rho_b =  fabs(aux.rho_s/(aux.alpn1*aux.Psi6*u0));
    double E_rad = aux.Erad_atm_cut;

    double F_rad_0 = - (F_rad_x * aux.ux + F_rad_y * aux.uy + F_rad_z * aux.uz)/u0;
    double F_rad0 = (-F_rad_0+F_rad_x*aux.shiftx+F_rad_y*aux.shifty+F_rad_z*aux.shiftz)/SQR(aux.alpn1);
    double gupxx4 = aux.gupxx_phys-aux.shiftx*aux.shiftx/SQR(aux.alpn1);
    double gupyy4 = aux.gupyy_phys-aux.shifty*aux.shifty/SQR(aux.alpn1);
    double gupzz4 = aux.gupzz_phys-aux.shiftz*aux.shiftz/SQR(aux.alpn1);
    double gupxy4 = aux.gupxy_phys-aux.shiftx*aux.shifty/SQR(aux.alpn1);
    double gupxz4 = aux.gupxz_phys-aux.shiftx*aux.shiftz/SQR(aux.alpn1);
    double gupyz4 = aux.gupyz_phys-aux.shifty*aux.shiftz/SQR(aux.alpn1);

    double F_radx = F_rad_0*aux.shiftx/SQR(aux.alpn1) + F_rad_x*gupxx4 + F_rad_y*gupxy4 +F_rad_z*gupxz4;
    double F_rady = F_rad_0*aux.shifty/SQR(aux.alpn1) + F_rad_x*gupxy4 + F_rad_y*gupyy4 +F_rad_z*gupyz4;
    double F_radz = F_rad_0*aux.shiftz/SQR(aux.alpn1) + F_rad_x*gupxz4 + F_rad_y*gupyz4 +F_rad_z*gupzz4;

    double FaFa = F_radx*F_rad_x + F_rady*F_rad_y + F_radz*F_rad_z + F_rad0*F_rad_0; // F_rad^alpha*F_rad_alpha                                                             
    double FkFk = F_radx*F_rad_x + F_rady*F_rad_y + F_radz*F_rad_z;

    double FaFa_o_E2;
    if (E_rad <= aux.zeta_cut){
      FaFa_o_E2 = 1.0;
    }  
    else{
      FaFa_o_E2= fabs(FaFa)/SQR(E_rad);
    }   
    
    if (FaFa_o_E2 > 1.0){
      FaFa_o_E2 = 1.0;
    }
    
    double A = FaFa_o_E2*(3.0 - sqrt(FaFa_o_E2) + 3.0*FaFa_o_E2)/5.0;
    //double A = 0.0;

    double FaFa_inv;
    if (FaFa <= SQR(aux.zeta_cut) || E_rad <= aux.zeta_cut){
      FaFa_inv == 0.0;
    }else{
      FaFa_inv = 1.0/FaFa;
    }

    //                                                                                                                  
    // fvec(1): Eq. for S_rad_x; fvec[2]: Eq. for S_rad_y; fvec[3]: Eq. for S_rad_z;  
    // fvec[4]: Eq. for tau_rad.                                                                                                          
    //
    /*  fvec[1] = aux.alpn1*aux.Psi6*(E_rad*u0*aux.u_x + F_rad0*aux.u_x + F_rad_x*u0 + E_rad/3.0*u0*aux.u_x) - aux.S_rad_x;
    fvec[2] = aux.alpn1*aux.Psi6*(E_rad*u0*aux.u_y + F_rad0*aux.u_y + F_rad_y*u0 + E_rad/3.0*u0*aux.u_y) - aux.S_rad_y;
    fvec[3] = aux.alpn1*aux.Psi6*(E_rad*u0*aux.u_z + F_rad0*aux.u_z + F_rad_z*u0 + E_rad/3.0*u0*aux.u_z) - aux.S_rad_z;
    */
    
    fvec[1] = (aux.alpn1*aux.Psi6*(E_rad*u0*aux.u_x + F_rad0*aux.u_x + F_rad_x*u0 + E_rad*(F_rad0*F_rad_x*FaFa_inv*A + u0*aux.u_x/3.0*(1.0-A))) - aux.S_rad_x)*aux.sri_scal_inv;
    fvec[2] = (aux.alpn1*aux.Psi6*(E_rad*u0*aux.u_y + F_rad0*aux.u_y + F_rad_y*u0 + E_rad*(F_rad0*F_rad_y*FaFa_inv*A + u0*aux.u_y/3.0*(1.0-A))) - aux.S_rad_y)*aux.sri_scal_inv;
    fvec[3] = (aux.alpn1*aux.Psi6*(E_rad*u0*aux.u_z + F_rad0*aux.u_z + F_rad_z*u0 + E_rad*(F_rad0*F_rad_z*FaFa_inv*A + u0*aux.u_z/3.0*(1.0-A))) - aux.S_rad_z)*aux.sri_scal_inv;
    
}



void jacobian_rad_font_fix(int &n,double *x,struct auxarray_rad &aux, double *fvec,int &np, double fjac[][5]){
  //
  double F_rad_x  = x[1]*(SQR(x[1]) + 1.0);
  double F_rad_y  = x[2]*(SQR(x[2]) + 1.0);
  double F_rad_z  = x[3]*(SQR(x[3]) + 1.0);
  //  double E_rad = x[4];
  double fac[4];
  fac[1] = (3.0*SQR(x[1]) + 1.0)*aux.sri_scal_inv;;
  fac[2] = (3.0*SQR(x[2]) + 1.0)*aux.sri_scal_inv;;
  fac[3] = (3.0*SQR(x[3]) + 1.0)*aux.sri_scal_inv;;

  double gijuiuj = aux.gupxx_phys*SQR(aux.u_x) + 2.0*aux.gupxy_phys*aux.u_x*aux.u_y +
    2.0*aux.gupxz_phys*aux.u_x*aux.u_z + aux.gupyy_phys*SQR(aux.u_y) + 2.0*aux.gupyz_phys*aux.u_y*aux.u_z +
    aux.gupzz_phys*SQR(aux.u_z);
  double au0m1 = gijuiuj/( 1.0+sqrt(1.0+gijuiuj) );
  if (aux.rho_s < 0.0) au0m1 = gijuiuj/( 1.0-sqrt(1.0+gijuiuj) );
  double u0 = (au0m1+1.0)/aux.alpn1;
  double u0_inv = aux.alpn1/(au0m1+1.0);

  double rho_b =  fabs(aux.rho_s/(aux.alpn1*aux.Psi6*u0));
  double E_rad = aux.Erad_atm_cut;

  double F_rad_0 = - (F_rad_x * aux.ux + F_rad_y * aux.uy + F_rad_z * aux.uz)/u0;
  double F_rad0 = (-F_rad_0+F_rad_x*aux.shiftx+F_rad_y*aux.shifty+F_rad_z*aux.shiftz)/SQR(aux.alpn1);
  double gupxx4 = aux.gupxx_phys-aux.shiftx*aux.shiftx/SQR(aux.alpn1);
  double gupyy4 = aux.gupyy_phys-aux.shifty*aux.shifty/SQR(aux.alpn1);
  double gupzz4 = aux.gupzz_phys-aux.shiftz*aux.shiftz/SQR(aux.alpn1);
  double gupxy4 = aux.gupxy_phys-aux.shiftx*aux.shifty/SQR(aux.alpn1);
  double gupxz4 = aux.gupxz_phys-aux.shiftx*aux.shiftz/SQR(aux.alpn1);
  double gupyz4 = aux.gupyz_phys-aux.shifty*aux.shiftz/SQR(aux.alpn1);

  double F_radx = F_rad_0*aux.shiftx/SQR(aux.alpn1) + F_rad_x*gupxx4 + F_rad_y*gupxy4 +F_rad_z*gupxz4;
  double F_rady = F_rad_0*aux.shifty/SQR(aux.alpn1) + F_rad_x*gupxy4 + F_rad_y*gupyy4 +F_rad_z*gupyz4;
  double F_radz = F_rad_0*aux.shiftz/SQR(aux.alpn1) + F_rad_x*gupxz4 + F_rad_y*gupyz4 +F_rad_z*gupzz4;

  double FaFa = F_radx*F_rad_x + F_rady*F_rad_y + F_radz*F_rad_z + F_rad0*F_rad_0; // F_rad^alpha*F_rad_alpha                                                        
  double FkFk = F_radx*F_rad_x + F_rady*F_rad_y + F_radz*F_rad_z;

  double zetasq;
  double zeta_temp = sqrt(fabs(FaFa/SQR(E_rad)));
  /*
  if ( zeta_temp < 1.0e-20){
    if (E_rad > aux.Erad_atm_cut){
      zetasq=0.0;
    }
    else{
      zetasq = 1.0;
    }
  }
  */

  if (E_rad <= aux.zeta_cut){
    zetasq = 1.0;
  }
  else{
    zetasq = SQR(zeta_temp);
  }

  if (zetasq > 1.0){
    zetasq = 1.0;
  }

  double zeta = sqrt(zetasq);

  //double zetasq =0.0;
  //double zeta =0.0;
  
  double A = zetasq*(3.0 - zeta + 3.0*zetasq)/5.0;
  double B = zeta*(6.0 - 3.0*zeta + 12.0*zetasq)/5.0;

  double vx = aux.ux/u0;
  double vy = aux.uy/u0;
  double vz = aux.uz/u0;

  double g00 = -1.0/SQR(aux.alpn1);
  double gx0 = aux.shiftx/SQR(aux.alpn1);
  double gy0 = aux.shifty/SQR(aux.alpn1);
  double gz0 = aux.shiftz/SQR(aux.alpn1);
  
  double FaFa_inv;
  if (FaFa <= SQR(aux.zeta_cut) || E_rad <= aux.zeta_cut){
    FaFa_inv == 0.0;
  }else{
    FaFa_inv = 1.0/FaFa;
  }
  /*
  fjac[1][1] = aux.alpn1*aux.Psi6*(u0 + (gx0 + vx*g00)*(aux.u_x));

  fjac[2][2] = aux.alpn1*aux.Psi6*(u0 + (gy0 + vy*g00)*(aux.u_y));

  fjac[3][3] = aux.alpn1*aux.Psi6*(u0 + (gz0 + vz*g00)*(aux.u_z));



  fjac[1][2] = aux.alpn1*aux.Psi6*((gx0 + vx*g00)*(aux.u_y));

  fjac[1][3] = aux.alpn1*aux.Psi6*((gx0 + vx*g00)*(aux.u_z));

  fjac[2][1] = aux.alpn1*aux.Psi6*((gy0 + vy*g00)*(aux.u_x));

  fjac[3][1] = aux.alpn1*aux.Psi6*((gz0 + vz*g00)*(aux.u_x));

  fjac[2][3] = aux.alpn1*aux.Psi6*((gy0 + vy*g00)*(aux.u_z));

  fjac[3][2] = aux.alpn1*aux.Psi6*((gz0 + vz*g00)*(aux.u_y));
  */
  
  fjac[1][1] = aux.alpn1*aux.Psi6*((u0+A*E_rad*F_rad0*FaFa_inv) + (gx0 + vx*g00)*(aux.u_x + A*E_rad*F_rad_x*FaFa_inv) +
                                   (F_radx - vx*F_rad0)*( B*sqrt(FaFa_inv)* (F_rad0*F_rad_x*FaFa_inv - u0*aux.u_x/3.0) - A*E_rad*2.0*F_rad0*F_rad_x*SQR(FaFa_inv)))*fac[1];

  fjac[2][2] = aux.alpn1*aux.Psi6*((u0+A*E_rad*F_rad0*FaFa_inv) + (gy0 + vy*g00)*(aux.u_y + A*E_rad*F_rad_y*FaFa_inv) +
                                   (F_rady - vy*F_rad0)*( B*sqrt(FaFa_inv)* (F_rad0*F_rad_y*FaFa_inv - u0*aux.u_y/3.0) - A*E_rad*2.0*F_rad0*F_rad_y*SQR(FaFa_inv)))*fac[2];

  fjac[3][3] = aux.alpn1*aux.Psi6*((u0+A*E_rad*F_rad0*FaFa_inv) + (gz0 + vz*g00)*(aux.u_z + A*E_rad*F_rad_z*FaFa_inv) +
                                   (F_radz - vz*F_rad0)*( B*sqrt(FaFa_inv)* (F_rad0*F_rad_z*FaFa_inv - u0*aux.u_z/3.0) - A*E_rad*2.0*F_rad0*F_rad_z*SQR(FaFa_inv)))*fac[3];

  fjac[2][1] = aux.alpn1*aux.Psi6*((gx0 + vx*g00)*(aux.u_y + A*E_rad*F_rad_y*FaFa_inv) +
                                   (F_radx - vx*F_rad0)*( B*sqrt(FaFa_inv)* (F_rad0*F_rad_y*FaFa_inv - u0*aux.u_y/3.0) - A*E_rad*2.0*F_rad0*F_rad_y*SQR(FaFa_inv)))*fac[2];

  fjac[3][1] = aux.alpn1*aux.Psi6*((gx0 + vx*g00)*(aux.u_z + A*E_rad*F_rad_z*FaFa_inv) +
                                   (F_radx - vx*F_rad0)*( B*sqrt(FaFa_inv)* (F_rad0*F_rad_z*FaFa_inv - u0*aux.u_z/3.0) - A*E_rad*2.0*F_rad0*F_rad_z*SQR(FaFa_inv)))*fac[3];

  fjac[1][2] = aux.alpn1*aux.Psi6*((gy0 + vy*g00)*(aux.u_x + A*E_rad*F_rad_x*FaFa_inv) +
                                   (F_rady - vy*F_rad0)*( B*sqrt(FaFa_inv)* (F_rad0*F_rad_x*FaFa_inv - u0*aux.u_x/3.0) - A*E_rad*2.0*F_rad0*F_rad_x*SQR(FaFa_inv)))*fac[1];

  fjac[1][3] = aux.alpn1*aux.Psi6*((gz0 + vz*g00)*(aux.u_x + A*E_rad*F_rad_x*FaFa_inv) +
                                   (F_radz - vz*F_rad0)*( B*sqrt(FaFa_inv)* (F_rad0*F_rad_x*FaFa_inv - u0*aux.u_x/3.0) - A*E_rad*2.0*F_rad0*F_rad_x*SQR(FaFa_inv)))*fac[1];

  fjac[3][2] = aux.alpn1*aux.Psi6*((gy0 + vy*g00)*(aux.u_z + A*E_rad*F_rad_z*FaFa_inv) +
                                   (F_rady - vy*F_rad0)*( B*sqrt(FaFa_inv)* (F_rad0*F_rad_z*FaFa_inv - u0*aux.u_z/3.0) - A*E_rad*2.0*F_rad0*F_rad_z*SQR(FaFa_inv)))*fac[3];

  fjac[2][3] = aux.alpn1*aux.Psi6*((gz0 + vz*g00)*(aux.u_y + A*E_rad*F_rad_y*FaFa_inv) +
				   (F_radz - vz*F_rad0)*( B*sqrt(FaFa_inv)* (F_rad0*F_rad_y*FaFa_inv - u0*aux.u_y/3.0) - A*E_rad*2.0*F_rad0*F_rad_y*SQR(FaFa_inv)))*fac[2];

  
  if (isnan(fjac[1][1])){ printf ("fjac[1][1] is nan!!! zeta=%e, zetasq=%e, FaFa_inv=%e, FaFa=%e \n", zeta,zetasq,FaFa_inv, FaFa );}
  if (zetasq > 1.0){
    printf ("zetasq too large!!! F_rad0,F_rad_0,F_radx,F_rad_x,F_rady,F_rad_y,F_radz,F_rad_z,E_rad= %e,%e,%e,%e,%e,%e,%e,%e,%e\n", F_rad0,F_rad_0,F_radx,F_rad_x,F_rady,F_rad_y,F_radz,F_rad_z,E_rad);
    printf("x[1],x[2],x[3],fac[1],fac[2],fac[3]=%e,%e,%e,%e,%e,%e \n", x[1],x[2],x[3],fac[1],fac[2],fac[3]);
  }



}





extern "C" void CCTK_FCALL CCTK_FNAME(primitive_vars_rad_cpp)
  (int *ext,int *nghostzones, double *X, double *Y, double *Z, 
   double *rho_star,
   double *tau_rad, double *S_rad_x, double *S_rad_y, double *S_rad_z,
   double *Sx, double *Sy, double *Sz, double *rho,
   double *Sxx, double *Sxy, double *Sxz, double *Syy, double *Syz, double *Szz,
   double *E_rad, double *F_radx, double *F_rady, double *F_radz, double *F_rad0, double *F,
   double *P_radxx, double *P_radyy, double *P_radzz, double *P_radxy, double *P_radxz,double *P_radyz,
   double *phi, double *alpha, double *shiftx, double *shifty, double *shiftz, double *gxx, double *gxy, double *gxz, double *gyy, double *gyz, double *gzz,
   double *gupxx, double *gupxy, double *gupxz, double *gupyy, double *gupyz, double *gupzz,
   double *vx, double *vy, double *vz, double *u0, double *chi_rad,double *zeta_rad,
   const cGH **cctkGH, double *failure_tracker_rad, int &ignore_ghostzones, int &repairs_rad_needed, double &Psi6threshold, double &Erad_atm_cut) {
  primitive_vars_rad_cpp(ext,nghostzones, X, Y, Z,
			     rho_star,
			     tau_rad,S_rad_x,S_rad_y,S_rad_z,
			     Sx, Sy, Sz, rho,
			     Sxx, Sxy, Sxz, Syy, Syz, Szz,
			     E_rad,F_radx,F_rady,F_radz,F_rad0,F,
			     P_radxx, P_radyy, P_radzz, P_radxy, P_radxz, P_radyz,
			     phi, alpha, shiftx, shifty, shiftz, gxx, gxy, gxz, gyy, gyz, gzz,
			     gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,
			     vx, vy, vz, u0, chi_rad,zeta_rad,
			    *cctkGH, failure_tracker_rad, ignore_ghostzones, repairs_rad_needed, Psi6threshold, Erad_atm_cut);
}
