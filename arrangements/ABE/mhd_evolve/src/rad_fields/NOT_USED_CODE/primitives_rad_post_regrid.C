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


double max_val_rad_pr(double val1,double val2);  

double max_val_rad_pr(double val1, double val2) {
  if(val1>val2) return val1;
  return val2;
}



void newt2_cpp_rad_pr(double x[5],struct auxarray_rad &aux,void (*function_rad)(int &n,double *x,double *fvec,struct auxarray_rad &aux),void (*jacobian_rad)(int &n,double *x,struct auxarray_rad &aux, double *fvec,int &np, double fjac[][5]),int &n,bool &check,int *indx,double *g,double *p,double *xold,double *fvec,int &MAXITS,double &STPMX);

void lnsrch_rad_pr(int n, double xold[], double fold, double g[], double p[], double x[5],double *f, double stpmax, bool &check, double (*fmin)(double *x,struct auxarray_rad &aux,double *fvec,void (*nrfuncv)(int &n,double *x,double *fvec, struct auxarray_rad &aux),int &n),void (*nrfuncv)(int &n,double *x,double *fvec,struct auxarray_rad &aux),struct auxarray_rad &aux,double *fvec);

double fmin_rad_pr(double *x,struct auxarray_rad &aux,double *fvec,void (*nrfuncv)(int &n,double *x,double *fvec,struct auxarray_rad &aux),int &n);


void function_rad_pr(int &n,double *x,double *fvec,struct auxarray_rad &aux); 
void function_rad_pr_font_fix(int &n,double *x,double *fvec,struct auxarray_rad &aux);        
void jacobian_rad_pr(int &n,double *x,struct auxarray_rad &aux, double *fvec,int &np, double fjac[][5]);
void jacobian_rad_pr_font_fix(int &n,double *x,struct auxarray_rad &aux, double *fvec,int &np, double fjac[][5]); 
void ludcmp_rad_pr(double a[][5], int n, int *indx, double *d);                        
void lubksb_rad_pr(double a[][5], int n, int *indx, double b[]);    


//double fmin_rad_pr(double *x,struct auxarray_rad &aux,double *fvec, int &n); 




extern "C" void CCTK_FCALL CCTK_FNAME(primitive_vars_rad_pr)
  (int *ext,int *nghostzones, 
   double *tau_rad, double *S_rad_x, double *S_rad_y, double *S_rad_z,
   double *Sx, double *Sy, double *Sz, double *rho, 
   double *Sxx, double *Sxy, double *Sxz, double *Syy, double *Syz, double *Szz,
   double *E_rad,  double *F_radx, double *F_rady, double *F_radz, 
   double *phi, double *alpha, double *shiftx, double *shifty, double *shiftz, double *gxx, double *gxy, double *gxz, double *gyy, double *gyz, double *gzz,
   double *gupxx, double *gupxy, double *gupxz, double *gupyy, double *gupyz, double *gupzz,
   double *vx, double *vy, double *vz, double *u0, 
   const cGH **cctkGH,int &ignore_ghostzones, int &repairs_rad_needed, double &Psi6threshold, double &Erad_atm_cut);


//-----------------------------------------------------------------------------
//
// reconstruct primitive variables, compute sources for "hybrid" EOS
//
//-----------------------------------------------------------------------------
void primitive_vars_rad_pr
(int *ext,int *nghostzones, 
 double *tau_rad, double *S_rad_x, double *S_rad_y, double *S_rad_z,
 double *Sx, double *Sy, double *Sz, double *rho,
 double *Sxx, double *Sxy, double *Sxz, double *Syy, double *Syz, double *Szz,
 double *E_rad, double *F_radx, double *F_rady, double *F_radz, 
 double *phi, double *alpha, double *shiftx, double *shifty, double *shiftz, double *gxx, double *gxy, double *gxz, double *gyy, double *gyz, double *gzz,
 double *gupxx, double *gupxy, double *gupxz, double *gupyy, double *gupyz, double *gupzz,
 double *vx, double *vy, double *vz, double *u0,
 const cGH *cctkGH,int &ignore_ghostzones, int &repairs_rad_needed, double &Psi6threshold, double &Erad_atm_cut) {

  //int count = 0;
  double f1o4p = 1.0/(4.0*M_PI);
  double f1o8p = f1o4p*0.5;

  // Timer:
  // reset the clock
  struct timeval start, end;
  long mtime, seconds, useconds;    
  gettimeofday(&start, NULL);

  printf("===============Start Primitives_rad_PR.C=========================+++ \n");

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

  //-----------------------------------------------------------------------------
  // Funny parameters...  (See Shibata, Oohara & Nakamura...)
  //-----------------------------------------------------------------------------

  // Initialize the failure_tracker_rad array
  /*    #pragma omp parallel for
    for(int i=0;i<ext[0];i++) for(int j=0;j<ext[1];j++) for(int k=0;k<ext[2];k++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    failure_tracker_rad[index] = 0.0;
    }*/
  //  repairs_rad_needed = 0;

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

	double UUguess[5];

	   int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

	   int MAXITS_beforefontfix=50;
	   double STPMX_beforefontfix=20.0; // <== faster than 100.0, when maxits set to 50
	   bool recom_rad = false;

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
           double alpn1 = alpha[index] + 1.0;
           double alpn1_inv = 1.0/alpn1;

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
	        
	   double beta2 = shiftxL*shift_xL + shiftyL*shift_yL + shiftzL*shift_zL;
	   double udotbeta = u0L*(vx[index]*shift_xL + vy[index]*shift_yL + vz[index]*shift_zL);
	   double g_00L =beta2-alpn1*alpn1;
	   double u_0L = g_00L*u0L + udotbeta;
	   
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
           double gijuiuj = gupxx_phys*SQR(u_xL) +
                  2.0*gupxy_phys*u_xL*u_yL + 2.0*gupxz_phys*u_xL*u_zL +
             gupyy_phys*SQR(u_yL) + 2.0*gupyz_phys*u_yL*u_zL +
             gupzz_phys*SQR(u_zL);
           au0m1 = gijuiuj/( 1.0+sqrt(1.0+gijuiuj) );

           // *** Limit velocity                                                                                                                                                 
           if (au0m1 > MAX_GAMMA-1.0) {
             double fac = sqrt((SQR(MAX_GAMMA) - 1.0)/(SQR(1.0+au0m1)-1.0));
             u_xL = fac*u_xL;
             u_yL = fac*u_yL;
             u_zL = fac*u_zL;
             gijuiuj = gijuiuj * SQR(fac);
             au0m1 = gijuiuj/( 1.0+sqrt(1.0+gijuiuj) );
           }
           u0[index] = (au0m1+1.0)*alpn1_inv;
	   u0L = u0[index];

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

	   double maxF_i = max_val_rad_pr(fabs(F_rad_x),
                                          max_val_rad_pr(fabs(F_rad_y),fabs(F_rad_z) ));

           double sri_scal =max_val_rad_pr(0.001*Psi6*(E_radl+F_rad0L+maxF_i),
                                           max_val_rad_pr(fabs(S_rad_x[index]),
                                                          max_val_rad_pr(fabs(S_rad_y[index]),fabs(S_rad_z[index]) )) );

           double tau_rad_scal = max_val_rad_pr(0.01*Psi6*(E_radl+2.0*F_rad0L), fabs(tau_rad[index]));
           double zeta_cut = Erad_atm_cut*1.5;

	   //	   double rho_s = rho_star[index];	  
	   int nn = 4;
	   double UU[5];
	  	  	         
           UU[1] = F_rad_x;
           UU[2] = F_rad_y;
           UU[3] = F_rad_z;
           UU[4] = E_radl;
 

	   struct auxarray_rad AUX;
	   //	   AUX.rho_s = 0.0;
	   AUX.tau_rad = tau_rad[index];
	   AUX.S_rad_x = S_rad_x[index];
	   AUX.S_rad_y = S_rad_y[index];
	   AUX.S_rad_z = S_rad_z[index];
	   
	   AUX.ux = u0L * vx[index];
	   AUX.uy = u0L * vy[index];
	   AUX.uz = u0L * vz[index];
	   AUX.u_x = u_xL;
	   AUX.u_y = u_yL;
	   AUX.u_z = u_zL;
	   
	   AUX.gupxx_phys = gupxx[index]*Psim4;
	   AUX.gupxy_phys = gupxy[index]*Psim4;
	   AUX.gupxz_phys = gupxz[index]*Psim4;
	   AUX.gupyy_phys = gupyy[index]*Psim4;
	   AUX.gupyz_phys = gupyz[index]*Psim4;
	   AUX.gupzz_phys = gupzz[index]*Psim4;
	   
	   //Note: Psi6 = sqrtg !
	   AUX.Psi6 = Psi6;
	   AUX.u_0 = u_0L;
	   AUX.alpn1 = alpn1;
	   AUX.g_00 = g_00L;
	   
	   AUX.shiftx = shiftxL;
	   AUX.shifty = shiftyL;
	   AUX.shiftz = shiftzL;
	   
	   AUX.Erad_atm_cut = Erad_atm_cut;
	   AUX.zeta_cut = zeta_cut;
	     

	   AUX.sri_scal_inv = 1.0/sri_scal;
	   if(isnan(AUX.sri_scal_inv))
	     {
	       printf("Inside primitive rad, AUX.sri_scal_inv is NAN!!! sri_scal = %e \n", sri_scal);
	     }
	   
	   
	   AUX.tau_rad_scal_inv = 1.0/tau_rad_scal;

	   

	   UUguess[1] = UU[1];
	   UUguess[2] = UU[2];
	   UUguess[3] = UU[3];
	   UUguess[4] = UU[4];	   	   

	   bool check;
	   int newt_indx[5];
	   double newt_g[5],newt_p[5],newt_xold[5],newt_fvec[5];	    
	   	  
	   newt2_cpp_rad_pr(UU,AUX,function_rad_pr,jacobian_rad_pr,nn,check, 
			    newt_indx,newt_g,newt_p,newt_xold,newt_fvec,MAXITS_beforefontfix,STPMX_beforefontfix);


	   int MAXITS_taustildefix=MAXITS_beforefontfix*10; //  <-- Try harder when doing a Font fix; we don't want to fail here!!!
	   double STPMX_taustildefix=STPMX_beforefontfix*10; //  <-- Try harder when doing a Font fix; we don't want to fail here!!!
	   	  

	   while(check || UU[4] < 0.0) {
	     check = false;
	     newt2_cpp_rad_pr(UU,AUX,function_rad_pr,jacobian_rad_pr,nn,check,
			   newt_indx,newt_g,newt_p,newt_xold,newt_fvec,MAXITS_taustildefix,STPMX_taustildefix);
	     MAXITS_taustildefix*=10;  // <-- Try even harder!  Don't give up!
	     //if(MAXITS_taustildefix>=MAXITS_beforefontfix*1e6) break;
	     if(MAXITS_taustildefix>=MAXITS_beforefontfix*1e4) break;
	   }
	   
	   
	     //****************************************************************
	     //                          FONT FIX
	     // Impose Font fix when Newt-Raph fails (check = true) or when it
	     //  gives negative eps [UU[4] < 0]
	     //****************************************************************
	     if(check || UU[4] < 0.0) {
	       //  failure_tracker_rad[index] = 1.0;
	       repairs_rad_needed = 1;   
	       
	       
	       double UU_font_fix[5];
	       
	       UU_font_fix[1] = UUguess[1];
	       UU_font_fix[2] = UUguess[2];
	       UU_font_fix[3] = UUguess[3];

	       nn = 3;
	         
	       check = false;
	       int MAXITS_fontfix=MAXITS_beforefontfix*10; //  <-- Try harder when doing a Font fix; we don't want to fail here!!!
	       double STPMX_fontfix=STPMX_beforefontfix*10; //  <-- Try harder when doing a Font fix; we don't want to fail here!!!
	       newt2_cpp_rad_pr(UU_font_fix,AUX,function_rad_pr_font_fix,jacobian_rad_pr_font_fix,nn,check,
			 newt_indx,newt_g,newt_p,newt_xold,newt_fvec,MAXITS_fontfix,STPMX_fontfix);
	       while(check) {
		 check = false;
		 newt2_cpp_rad_pr(UU_font_fix,AUX,function_rad_pr_font_fix,jacobian_rad_pr_font_fix,nn,check,
			   newt_indx,newt_g,newt_p,newt_xold,newt_fvec,MAXITS_fontfix,STPMX_fontfix);
		 MAXITS_fontfix*=10;  // <-- Try even harder!  Don't give up!
		 if(MAXITS_fontfix>=MAXITS_beforefontfix*1e4) break;
		 //if(MAXITS_fontfix>=MAXITS_beforefontfix*1e4) break;
	       }
	       if(check) {
	         printf("ERROR: FONT FIX (secondary solver) JUST FAILED\n");
		 //	         printf("Problem at (x,y,z) = %e %e %e, %d %d %d\n", X[index],Y[index],Z[index],i,j,k);
	         printf("tau_rad = %e\n",tau_rad[index]);
	         printf("S_rad_x = %e\n",S_rad_x[index]);
	         printf("S_rad_y = %e\n",S_rad_y[index]);
	         printf("S_rad_z = %e\n",S_rad_z[index]);
	         printf("gamma^xx = %e\n",AUX.gupxx_phys);
	         printf("gamma^xy = %e\n",AUX.gupxy_phys);
	         printf("gamma^xz = %e\n",AUX.gupxz_phys);
	         printf("gamma^yy = %e\n",AUX.gupyy_phys);
	         printf("gamma^yz = %e\n",AUX.gupyz_phys);
	         printf("gamma^zz = %e\n",AUX.gupzz_phys);
	         printf("exp(6 phi) = %e\n",AUX.Psi6);
		 
		 //		 failure_tracker_rad[index] = 1.0;
		 repairs_rad_needed = 1;
	      
	         // Set everything to 0.0 before calculating 
	         // these quantities from averages.
	         F_rad_x         = 0.0;
	         F_rad_y         = 0.0;
	         F_rad_z         = 0.0;
	         recom_rad        = true;
	       } else {
	         //if the Font fix worked, do the following: 
	         F_rad_x = UU_font_fix[1]*(SQR(UU_font_fix[1]) + 1.0);
	         F_rad_y = UU_font_fix[2]*(SQR(UU_font_fix[2]) + 1.0);
	         F_rad_z = UU_font_fix[3]*(SQR(UU_font_fix[3]) + 1.0);
	         
		 recom_rad = true;    
	       }	     
	     } else {
	       // Inversion worked without the Font fix!  Now set the primtives
       	       F_rad_x = UU[1];
	       F_rad_y = UU[2];
	       F_rad_z = UU[3];
	       E_radl  = UU[4];
	     }	 
	   
               

 
	     ////////////////////////Add radiation terms to the source.////////////////////////////
	     //Apply atm floor of rad primitives:  	  
	     //	     if (E_radl < Erad_atm_cut)	     
	     double F_radxl, F_radyl, F_radzl;
	     if (E_radl < zeta_cut) 
               {                    
		 recom_rad = true;
		 E_radl = Erad_atm_cut;
		 F_rad_x = 0.0;
		 F_rad_y = 0.0;
		 F_rad_z = 0.0;
	       }
	     	    
             //Apply Psi6threshold fix on E_radl:                                                 
             if(Psi6 > Psi6threshold) {
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

	      
	     	     //	     double F_rad0l, F_rad_0l, F_radxl, F_radyl, F_radzl;	   	    
             F_rad_0l = - (F_rad_x * uxL + F_rad_y * uyL + F_rad_z * uzL)/u0L;
             F_rad0l = (-F_rad_0l+F_rad_x*shiftx[index]+F_rad_y*shifty[index]+F_rad_z*shiftz[index])/SQR(alpn1);

             double gupxx4 = gupxx_phys-shiftx[index]*shiftx[index]/SQR(alpn1);
             double gupyy4 = gupyy_phys-shifty[index]*shifty[index]/SQR(alpn1);
             double gupzz4 = gupzz_phys-shiftz[index]*shiftz[index]/SQR(alpn1);
             double gupxy4 = gupxy_phys-shiftx[index]*shifty[index]/SQR(alpn1);
             double gupxz4 = gupxz_phys-shiftx[index]*shiftz[index]/SQR(alpn1);
             double gupyz4 = gupyz_phys-shifty[index]*shiftz[index]/SQR(alpn1);

             F_radxl = F_rad_0l*shiftx[index]/SQR(alpn1) + F_rad_x*gupxx4 + F_rad_y*gupxy4 +F_rad_z*gupxz4;
             F_radyl = F_rad_0l*shifty[index]/SQR(alpn1) + F_rad_x*gupxy4 + F_rad_y*gupyy4 +F_rad_z*gupyz4;
             F_radzl = F_rad_0l*shiftz[index]/SQR(alpn1) + F_rad_x*gupxz4 + F_rad_y*gupyz4 +F_rad_z*gupzz4;	     
	     

	     // Limit Flux                                                                                                  
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
               F_rad_0l = - (F_rad_x * uxL + F_rad_y * uyL + F_rad_z * uzL)/u0L;
               F_rad0l = (-F_rad_0l+F_rad_x*shiftx[index]+F_rad_y*shifty[index]+F_rad_z*shiftz[index])/SQR(alpn1);
             }
	     
	   
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
	    
	    
	     double P_radxxl,P_radyyl,P_radzzl,P_radxyl,P_radxzl,P_radyzl;
	     // Since we use gupij_phys = gupij * psim4, we set psim4 = 1.0 in this function.
	     double psim4_temp = 1.0;
	     double lapse = alpha[index];
	     compute_M1(P_radxxl, F_radxl, F_radxl, Fasq, E_radl, gupxx_phys, shiftx[index], shiftx[index], lapse, uxL, uxL, chi, psim4_temp, Erad_atm_cut);
	     compute_M1(P_radyyl, F_radyl, F_radyl, Fasq, E_radl, gupyy_phys, shifty[index], shifty[index], lapse, uyL, uyL, chi, psim4_temp, Erad_atm_cut);
	     compute_M1(P_radzzl, F_radzl, F_radzl, Fasq, E_radl, gupzz_phys, shiftz[index], shiftz[index], lapse, uzL, uzL, chi, psim4_temp, Erad_atm_cut);
	     compute_M1(P_radxyl, F_radxl, F_radyl, Fasq, E_radl, gupxy_phys, shiftx[index], shifty[index], lapse, uxL, uyL, chi, psim4_temp, Erad_atm_cut);
	     compute_M1(P_radxzl, F_radxl, F_radzl, Fasq, E_radl, gupxz_phys, shiftx[index], shiftz[index], lapse, uxL, uzL, chi, psim4_temp, Erad_atm_cut);	       
             compute_M1(P_radyzl, F_radyl, F_radzl, Fasq, E_radl, gupyz_phys, shifty[index], shiftz[index], lapse, uyL, uzL, chi, psim4_temp, Erad_atm_cut);

	
	     
	     double P_rad0xl = - (P_radxxl * u_xL + P_radxyl * u_yL + P_radxzl * u_zL)/u_0L;
	     double P_rad0yl = - (P_radxyl * u_xL + P_radyyl * u_yL + P_radyzl * u_zL)/u_0L;
	     double P_rad0zl = - (P_radxzl * u_xL + P_radyzl * u_yL + P_radzzl * u_zL)/u_0L;	   
	     double P_rad00l = - (P_rad0xl * u_xL + P_rad0yl * u_yL + P_rad0zl * u_zL)/u_0L;
	     
	    

	     	   
	     E_rad[index] = E_radl;
	     F_radx[index] = F_radxl;
	     F_rady[index] = F_radyl;
	     F_radz[index] = F_radzl;
	   
	   
	     rho[index] = rho[index] + SQR(alpn1)*(E_radl*SQR(u0L) + 2.0*F_rad0l*u0L + P_rad00l);
	     
	     Sx[index] = Sx[index] + alpn1*(u0L*E_radl*AUX.u_x + F_rad0l*AUX.u_x + u0L*F_rad_x + 
					    P_rad00l*shift_xL + Psi4*(P_rad0xl*gxxL + P_rad0yl*gxyL + P_rad0zl*gxzL));
	   
	     Sy[index] = Sy[index] + alpn1*(u0L*E_radl*AUX.u_y + F_rad0l*AUX.u_y + u0L*F_rad_y +
					    P_rad00l*shift_yL + Psi4*(P_rad0xl*gxyL + P_rad0yl*gyyL + P_rad0zl*gyzL));
	   
	     Sz[index] = Sz[index] + alpn1*(u0L*E_radl*AUX.u_z + F_rad0l*AUX.u_z + u0L*F_rad_z +
					    P_rad00l*shift_zL + Psi4*(P_rad0xl*gxzL + P_rad0yl*gyzL + P_rad0zl*gzzL));
	     
	   	   
	     Sxx[index] = Sxx[index] + E_radl*AUX.u_x*AUX.u_x + 2*F_rad_x*AUX.u_x + 
	       SQR(shift_xL)*P_rad00l + shift_xL*2.0*(gxxL*P_rad0xl+gxyL*P_rad0yl+gxzL*P_rad0zl) + 
	       SQR(Psi4)*( SQR(gxxL)*P_radxxl + SQR(gxyL)*P_radyyl + SQR(gxzL)*P_radzzl +
			   2.0*(gxxL*gxyL*P_radxyl+gxxL*gxzL*P_radxzl+gxyL*gxzL*P_radyzl) );
	     
	     Syy[index] = Syy[index] + E_radl*AUX.u_y*AUX.u_y + 2*F_rad_y*AUX.u_y +
	       SQR(shift_yL)*P_rad00l + shift_yL*2.0*(gxyL*P_rad0xl+gyyL*P_rad0yl+gyzL*P_rad0zl) + 
	       SQR(Psi4)*( SQR(gxyL)*P_radxxl + SQR(gyyL)*P_radyyl + SQR(gyzL)*P_radzzl +
			   2.0*(gxyL*gyyL*P_radxyl+gxyL*gyzL*P_radxzl+gyyL*gyzL*P_radyzl) );
	     
	     Szz[index] = Szz[index] + E_radl*AUX.u_z*AUX.u_z + 2*F_rad_z*AUX.u_z +
	       SQR(shift_zL)*P_rad00l + shift_zL*2.0*(gxzL*P_rad0xl+gyzL*P_rad0yl+gzzL*P_rad0zl) + 
	       SQR(Psi4)*( SQR(gxzL)*P_radxxl + SQR(gyzL)*P_radyyl + SQR(gzzL)*P_radzzl +
			   2.0*(gxzL*gyzL*P_radxyl+gxzL*gzzL*P_radxzl+gyzL*gzzL*P_radyzl) );
	     
	     Sxy[index] = Sxy[index] + E_radl*AUX.u_x*AUX.u_y + F_rad_x*AUX.u_y + F_rad_y*AUX.u_x +
	       shift_xL*shift_yL*P_rad00l + 
	       Psi4*(shift_xL*(gxyL*P_rad0xl + gyyL*P_rad0yl + gyzL*P_rad0zl) + (shift_yL*(gxxL*P_rad0xl + gxyL*P_rad0yl + gxzL*P_rad0zl)) )+
	       SQR(Psi4)*(gxxL*gxyL*P_radxxl + gxyL*gyyL*P_radyyl + gxzL*gyzL*P_radzzl +
			  (gxxL*gyyL + gxyL*gxyL)*P_radxyl + (gxxL*gyzL + gxzL*gxyL)*P_radxzl + (gxyL*gyzL + gxzL*gyyL)*P_radyzl);
	     
	     Sxz[index] = Sxz[index] + E_radl*AUX.u_x*AUX.u_z + F_rad_x*AUX.u_z + F_rad_z*AUX.u_x +
	       shift_xL*shift_zL*P_rad00l +
	       Psi4*(shift_xL*(gxzL*P_rad0xl + gyzL*P_rad0yl + gzzL*P_rad0zl) + (shift_zL*(gxxL*P_rad0xl + gxyL*P_rad0yl + gxzL*P_rad0zl)) )+
	       SQR(Psi4)*(gxxL*gxzL*P_radxxl + gxyL*gyzL*P_radyyl + gxzL*gzzL*P_radzzl +
			  (gxxL*gyzL + gxyL*gxzL)*P_radxyl + (gxxL*gzzL + gxzL*gxzL)*P_radxzl + (gxyL*gzzL + gxzL*gyzL)*P_radyzl);
	     
	     Syz[index] = Syz[index] + E_radl*AUX.u_y*AUX.u_z + F_rad_y*AUX.u_z + F_rad_z*AUX.u_y +
	       shift_yL*shift_zL*P_rad00l +
	       Psi4*(shift_yL*(gxzL*P_rad0xl + gyzL*P_rad0yl + gzzL*P_rad0zl) + (shift_zL*(gxyL*P_rad0xl + gyyL*P_rad0yl + gyzL*P_rad0zl)) )+
	       SQR(Psi4)*(gxyL*gxzL*P_radxxl + gyyL*gyzL*P_radyyl + gyzL*gzzL*P_radzzl +
			  (gxyL*gyzL + gxzL*gyyL)*P_radxyl + (gxyL*gzzL + gyzL*gxzL)*P_radxzl + (gyyL*gzzL + gyzL*gyzL)*P_radyzl);
	     
	     
	     if (recom_rad){
	       tau_rad[index] = alpn1*alpn1*Psi6*(E_radl*u0L*u0L+2.0*F_rad0l*u0L+P_rad00l);
	       S_rad_x[index] = alpn1*Psi6*(E_radl*u0L*AUX.u_x + F_rad0l*AUX.u_x + F_rad_x * u0L + P_rad00l*shift_xL + Psi4*(P_rad0xl*gxxL + P_rad0yl*gxyL + P_rad0zl*gxzL));
	       S_rad_y[index] = alpn1*Psi6*(E_radl*u0L*AUX.u_y + F_rad0l*AUX.u_y + F_rad_y * u0L + P_rad00l*shift_yL + Psi4*(P_rad0xl*gxyL + P_rad0yl*gyyL + P_rad0zl*gyzL));
	       S_rad_z[index] = alpn1*Psi6*(E_radl*u0L*AUX.u_z + F_rad0l*AUX.u_z + F_rad_z * u0L + P_rad00l*shift_zL + Psi4*(P_rad0xl*gxzL + P_rad0yl*gyzL + P_rad0zl*gzzL));
	     }



	     //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++            
	     // Just a blank run of newton solver to check the bug                                                                                            
	     //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++                                                            
	     

             //      double P_radl, F_radxl, F_radyl, F_radzl, F_rad0l;                                                                                                      
	       //Apply atm floor of rad primitives:                                                                                                                                                                                                                                                                                                                                                                       
	     /*             if (tau_radl < 0.0)                                                                                                                                                                                 
	       {                                                                                                                                                          
		 recom_rad = true;                                                                                                                 
		 E_radl = Erad_atm_cut;                                              
		 F_radxl = 0.0;                                                                                                                                    
		 F_radyl = 0.0;                                                                                              
		 F_radzl = 0.0;                                                                                                                                                  
		 F_rad0l = 0.0;                                                                                                                                   
		 P_radl = E_radl/3.0;
	       }                                                                                                                                                        
             else{                                             
               double au0 = au0m1 + 1.0;                                                           
	       double F_rad0_denom = alpn1 * Psi6 *(2.0*au0*au0 + 1.0)/(4.0*au0*au0 - 1.0);                                                                                                                    
               double F_rad0_num = 4.0*au0*(1.0 - au0*au0)/(4.0*au0*au0 - 1.0)*tau_radl       
		 +( (shiftxL*u0L + uxL)*S_rad_xl                                                                                                            
                    +(shiftyL*u0L + uyL)*S_rad_yl                                                                                                                                          
                    +(shiftzL*u0L + uzL)*S_rad_zl);                                                                                                  
	       
               F_rad0l = F_rad0_num/F_rad0_denom;                                                                      
                                                                                                             
               E_radl = 3.0*(tau_radl/Psi6 - 2.0*alpn1*au0*F_rad0l)/(4.0*au0*au0 - 1.0);                                                                                                            
                                                                                                                                               
                                                                                                                                                                                            
                                                                                                                                               
               P_radl = E_radl/3.0;                                                                                                                                                                    
                                                                                                                                               
               F_radxl = ((gupxx_phys*S_rad_xl +                                                                                              
                           gupxy_phys*S_rad_yl +                                                                                                                                                    
			   gupxz_phys*S_rad_zl))/(alpn1*Psi6*u0L)-                                                                                   
                 (E_radl+P_radl)*u0L*(vx[index]+shiftxL) -                                                                                         
		 2.0*(F_rad0l * shiftxL) - F_rad0l*vx[index];                                                
                                                                                                                              
                                                                                                                                                                                   
                                                                     
               F_radyl = ((gupxy_phys*S_rad_xl +                                                                                         
                           gupyy_phys*S_rad_yl +                                                                                                           
			   gupyz_phys*S_rad_zl))/(alpn1*Psi6*u0L)-                                                                                                     
		 (E_radl+P_radl)*u0L*(vy[index]+shiftyL) -                                                                                                                
		 2.0*(F_rad0l * shiftyL) - F_rad0l*vy[index];                                                                                                                                 
                                                                                                                                                                                                  
                                                                                                                              
               F_radzl = ((gupxz_phys*S_rad_xl +                                                                                                               
                           gupyz_phys*S_rad_yl +                                                                                                                                             
			   gupzz_phys*S_rad_zl))/(alpn1*Psi6*u0L)-                                                                                                                                
		 (E_radl+P_radl)*u0L*(vz[index]+shiftzL) -                                                                               
                 2.0*(F_rad0l * shiftzL) - F_rad0l*vz[index];                                                                          
             }                                                                                                                                                                                  
                                                                                                 
                                                                                                        
             double shift_xl = Psi4 * (shiftxL*gxxL + shiftyL*gxyL +shiftzL*gxzL);                                                                                          
             double shift_yl = Psi4 * (shiftxL*gxyL + shiftyL*gyyL +shiftzL*gyzL);        
             double shift_zl = Psi4 * (shiftxL*gxzL + shiftyL*gyzL +shiftzL*gzzL);                                                             
                                                                                                                                               
                                                                                                                                                                                                     
                                                                                                                                               
             //Apply Psi6threshold fix on E_radl, and F_radl:                                                                                                                                        
                                                                                                                                               
             if(Psi6 > Psi6threshold) {                                                                                                                                         
               recom_rad = true;                                                                                                                                            
               double Erad_horiz_cap = 2.0*Erad_atm_cut;                     
               if(E_radl < Erad_horiz_cap){                                                                                                                
                 E_radl = Erad_horiz_cap;                                                                                                                                                              
                 P_radl = E_radl/3.0;                                                                                                                                    
                 F_radxl = 0.0;                                                                                                                                                           
                 F_radyl = 0.0;                                                                                                                                
                 F_radzl = 0.0;                                                                                                                                                                
                 F_rad0l = 0.0;                                                                                                                                                                       
               }                                                                                                                   
             }                                                                                                                                                                                         
	     E_rad[index] = E_radl;                                                                
	     F_radx[index] = F_radxl;                                                                                                                                         
	     F_rady[index] = F_radyl;                                                                                                                                                     
	     F_radz[index] = F_radzl;                                                                                                                                          
             //      F_rad0[index] = F_rad0l;                                                                                                                              
             //      F[index] = sqrt(Psi4*(gxxL*SQR(F_radxl)+gyyL*SQR(F_radyl)+gzzL*SQR(F_radzl) + 2.0*( gxyL*F_radxl*F_radyl + gxzL*F_radxl*F_radzl + gyzL*F_radyl*F_radzl))  )     
                                                                      
             double temp_rad = alpn1*u0L;                                                                                                                                                             
             double temp_rad1 = temp_rad*temp_rad*(E_radl+P_radl) - P_radl + 2.0*SQR(alpn1)*u0L*F_rad0l;                                                                                       
             //F_rad_\alpha = g_{0\alpha} F_rad0 + g_{z\alpha} F_radz + g_{y\alpha} F_rady + g_{z\alpha} F_radz                                                            
             double F_rad_xl = Psi4 * (gxxL * F_radxl + gxyL * F_radyl + gxzL * F_radzl) + shift_xl*F_rad0l;                                                                             
             double F_rad_yl = Psi4 * (gxyL * F_radxl + gyyL * F_radyl + gyzL * F_radzl) + shift_yl*F_rad0l;                                                                                  
             double F_rad_zl = Psi4 * (gxzL * F_radxl + gyzL * F_radyl + gzzL * F_radzl) + shift_zl*F_rad0l;                                                                                           
                                                                                                                                                                                                   
             double v_xl = Psi4 * (vx[index]*gxxL + vy[index]*gxyL +vz[index]*gxzL);                                                                                 
             double v_yl = Psi4 * (vx[index]*gxyL + vy[index]*gyyL +vz[index]*gyzL);   
	     double v_zl = Psi4 * (vx[index]*gxzL + vy[index]*gyzL +vz[index]*gzzL);             
	     
	     rho[index]   = rho[index] + temp_rad1;                                                                         
	     Sx[index] = Sx[index] + temp_rad*( ( (E_rad[index]+E_rad[index]/3.0)*u0L + F_rad0l) * (shift_xl + v_xl) + F_rad_xl); 
	     Sy[index] = Sy[index] + temp_rad*( ( (E_rad[index]+E_rad[index]/3.0)*u0L + F_rad0l) * (shift_yl + v_yl) + F_rad_yl);   
	     Sz[index] = Sz[index] + temp_rad*( ( (E_rad[index]+E_rad[index]/3.0)*u0L + F_rad0l) * (shift_zl + v_zl) + F_rad_zl); 

             Sxx[index] = Sxx[index] + (E_rad[index]+E_rad[index]/3.0)*SQR(u0L*(shift_xl+v_xl)) + 2.0*F_rad_xl*u0L*(shift_xl+v_xl) + Psi4 * E_rad[index]/3.0 * gxxL;
             Syy[index] = Syy[index] + (E_rad[index]+E_rad[index]/3.0)*SQR(u0L*(shift_yl+v_yl)) + 2.0*F_rad_yl*u0L*(shift_yl+v_yl) + Psi4 * E_rad[index]/3.0 * gyyL;
             Szz[index] = Szz[index] + (E_rad[index]+E_rad[index]/3.0)*SQR(u0L*(shift_zl+v_zl)) + 2.0*F_rad_yl*u0L*(shift_zl+v_zl) + Psi4 * E_rad[index]/3.0 * gzzL;

             Sxy[index] = Sxy[index] + (E_rad[index]+E_rad[index]/3.0)*SQR(u0L)*(shift_xl+v_xl)*(shift_yl+v_yl) + u0L*(F_rad_xl*(shift_yl+v_yl)+F_rad_yl*(shift_xl+v_xl)) + Psi4 * E_rad[index]/3.0 * gxyL;          
	     Sxz[index] = Sxz[index] + (E_rad[index]+E_rad[index]/3.0)*SQR(u0L)*(shift_xl+v_xl)*(shift_zl+v_zl) + u0L*(F_rad_xl*(shift_zl+v_zl)+F_rad_zl*(shift_xl+v_xl)) + Psi4 * E_rad[index]/3.0 * gxzL;
             Syz[index] = Syz[index] + (E_rad[index]+E_rad[index]/3.0)*SQR(u0L)*(shift_yl+v_yl)*(shift_zl+v_zl) + u0L*(F_rad_yl*(shift_zl+v_zl)+F_rad_zl*(shift_yl+v_yl)) + Psi4 * E_rad[index]/3.0 * gyzL;


             if (recom_rad){
               tau_rad[index] = SQR(alpn1)*Psi6*(E_radl*SQR(u0L) + 2.0*F_rad0l*u0L + P_radl*SQR(u0L)) - Psi6*P_radl;
               S_rad_x[index] = alpn1*Psi6*((E_radl+P_radl)*u0L*u_xL + F_rad0l*u_xL + F_rad_xl*u0L); 
               S_rad_x[index] = alpn1*Psi6*((E_radl+P_radl)*u0L*u_yL + F_rad0l*u_yL + F_rad_yl*u0L);
               S_rad_z[index] = alpn1*Psi6*((E_radl+P_radl)*u0L*u_zL + F_rad0l*u_zL + F_rad_zl*u0L);
	     }

	     //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
             // Just a blank run of newton solver to check the bug                                                                                                                                        
             //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++      	

	     */
      }
  
	     
  gettimeofday(&end, NULL);
  
  seconds  = end.tv_sec  - start.tv_sec;
  useconds = end.tv_usec - start.tv_usec;

  mtime = ((seconds) * 1000 + useconds/1000.0) + 0.5;  


  printf("Primitives radiation solver: %f solutions/second\n",(iend-istart)*(jend-jstart)*(kend-kstart) / ((double)mtime/1000.0));

}
 








#define NP 5
//#define MAXITS 50
#define TOLF 1.0e-10
#define TOLMIN 1.0e-14
#define TOLX 3.0e-16

//#define STPMX 10000.0
//#define STPMX 100.0
#define FREERETURN {return;}


void newt2_cpp_rad_pr(double x[5],struct auxarray_rad &aux,
	       void (*funcv1)(int &n,double *x,double *fvec,struct auxarray_rad &aux),
	       void (*fdjac1)
	       (int &n,double *x,struct auxarray_rad &aux, double *fvec,int &np,double fjac[][5]),
	       int &n,bool &check,
	       int *indx,double *g,double *p,double *xold,double *fvec,int &MAXITS,double &STPMX)
{
  void lubksb_rad_pr(double a[][5], int n, int *indx, double b[]);
  void ludcmp_rad_pr(double a[][5], int n, int *indx, double *d);
  int i,its,j;
  double d,den,f,fold,stpmax,sum,temp,test,fjac[5][5];


  //printf ("INSIDE newt2_cpp_rad, before fmin_rad, fvec[1]=%e, fvec[2]=%e,fvec[3]=%e, fvec[4]=%e, ux=x[1]=%e \n", fvec[1], fvec[2],fvec[3], fvec[4], x[1]);
 
  f=fmin_rad_pr(x,aux,fvec,funcv1,n);

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
  stpmax=STPMX*max_val_rad_pr(sqrt(sum),(double)n);
  for (its=1;its<=MAXITS;its++) {
    int np=NP;
    fdjac1(n,x,aux,fvec,np,fjac);

    for (i=1;i<=n;i++) {
      for (sum=0.0,j=1;j<=n;j++) sum += fjac[j][i]*fvec[j];
      g[i]=sum;
    }

    for (i=1;i<=n;i++) xold[i]=x[i];
    fold=f;
    for (i=1;i<=n;i++) p[i] = -fvec[i];

    // The two following lines solve the equation fjac * x = p, where p[i] = -fvec[i] = -\delta f(i).
    // After two functions are called, the vector p(i) is replaced by the new solution x(i), used as the Newton directon of new solution. 
    ludcmp_rad_pr(fjac,n,indx,&d);
    lubksb_rad_pr(fjac,n,indx,p);
    
    
    lnsrch_rad_pr(n,xold,fold,g,p,x,&f,stpmax,check,fmin_rad_pr,funcv1,aux,fvec);

    test=0.0;
    for (i=1;i<=n;i++)
      if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
    if (test < TOLF) {
      check=0;
      FREERETURN
	}
    if (check) {
      test=0.0;
      den=max_val_rad_pr(f,0.5*n);
      for (i=1;i<=n;i++) {
	temp=fabs(g[i])*max_val_rad_pr(fabs(x[i]),1.0)/den;
	if (temp > test) test=temp;
      }
      check=(test < TOLMIN ? 1 : 0);
      if(check!=0) printf("BAD CHECK %d: %e\t%e\t%e\t%e\t%e\t%e\n",check,test,TOLMIN,fabs(x[1]),fabs(x[2]),fabs(x[3]),fabs(x[4]) );
      FREERETURN
	}
    test=0.0;
    for (i=1;i<=n;i++) {
      temp=(fabs(x[i]-xold[i]))/max_val_rad_pr(fabs(x[i]),1.0);
      if (temp > test) test=temp;
    }
    if (test < TOLX) FREERETURN
		       }
  check=1;
  //if(check!=0) printf("BAD STUFF %d: %e\t%e.  ITS=%d\n",check,test,TOLMIN,its);
  //printf("MAXITS exceeded in newt\n");
  //exit(0);
}


double fmin_rad_pr(double *x,struct auxarray_rad &aux,double *fvec,
		   void (*nrfuncv)(int &n,double *x,double *fvec,struct auxarray_rad &aux),
		 int &n)
{
  int i;
  double sum;

  (*nrfuncv)(n,x,fvec,aux);
  for (sum=0.0,i=1;i<=n;i++) sum += SQR(fvec[i]);
  return 0.5*sum;
}







#define ALF 1.0e-4
#define TOLX2 1.0e-16

void lnsrch_rad_pr(int n, double xold[], double fold, double g[], double p[], double x[5],
		 double *f, double stpmax, bool &check, 
		 double (*fmin)
		 (double *x,struct auxarray_rad &aux,double *fvec,
		  void (*nrfuncv)(int &n,double *x,double *fvec,struct auxarray_rad &aux),
		  int &n),
		 void (*nrfuncv)(int &n,double *x,double *fvec,struct auxarray_rad &aux),
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
  if (slope >= 0.0) printf("OUCH!  Roundoff problem in lnsrch.  Note that this function has been updated to correct a bug in Numerical Recipes 2.06.  This is version 2.08, with the 2.06->2.08 diff obtained from http://www.numerical-recipes.com/upgrade/upgrade-208.html  x[1-4] = %e,%e,%e,%e \n", xold[1],xold[2],xold[3],xold[4]);
  test=0.0;
  for (i=1;i<=n;i++) {
    temp=fabs(p[i])/max_val_rad_pr(fabs(xold[i]),1.0);
    if (temp > test) test=temp;
  }
  alamin=TOLX2/test;
  alam=1.0;
  for (;;) {
    for (i=1;i<=n;i++) x[i]=xold[i]+alam*p[i];
    *f=(*fmin)(x,aux,fvec,nrfuncv,n);
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
    alam=max_val_rad_pr(tmplam,0.1*alam);
    //    alam=FMAX(tmplam,0.1*alam);
  }

}
#undef ALF
#undef TOLX2



void lubksb_rad_pr(double a[][5], int n, int *indx, double b[])
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
void ludcmp_rad_pr(double a[][5], int n, int *indx, double *d)
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



void function_rad_pr(int &n,double *x,double *fvec,struct auxarray_rad &aux){
  double F_rad_x    = x[1];
  double F_rad_y    = x[2];
  double F_rad_z    = x[3];
  double E_rad     = x[4];


  double gijuiuj = aux.gupxx_phys*SQR(aux.u_x) + 2.0*aux.gupxy_phys*aux.u_x*aux.u_y +
    2.0*aux.gupxz_phys*aux.u_x*aux.u_z + aux.gupyy_phys*SQR(aux.u_y) + 2.0*aux.gupyz_phys*aux.u_y*aux.u_z +
    aux.gupzz_phys*SQR(aux.u_z);
  double au0m1 = gijuiuj/( 1.0+sqrt(1.0+gijuiuj) );
  //if (aux.rho_s < 0.0) au0m1 = gijuiuj/( 1.0-sqrt(1.0+gijuiuj) );
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
  if (FaFa <= 0.0 || E_rad <= aux.zeta_cut){
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



void jacobian_rad_pr
(int &n,double *x,struct auxarray_rad &aux, double *fvec,int &np, double fjac[][5]){
  double F_rad_x    = x[1];
  double F_rad_y    = x[2];
  double F_rad_z    = x[3];
  double E_rad     = x[4];

  double gijuiuj = aux.gupxx_phys*SQR(aux.u_x) + 2.0*aux.gupxy_phys*aux.u_x*aux.u_y +
    2.0*aux.gupxz_phys*aux.u_x*aux.u_z + aux.gupyy_phys*SQR(aux.u_y) + 2.0*aux.gupyz_phys*aux.u_y*aux.u_z +
    aux.gupzz_phys*SQR(aux.u_z);
  double au0m1 = gijuiuj/( 1.0+sqrt(1.0+gijuiuj) );
  //  if (aux.rho_s < 0.0) au0m1 = gijuiuj/( 1.0-sqrt(1.0+gijuiuj) );
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
  
  
  double FaFa_inv;
  if (FaFa <= 0.0 || E_rad <= aux.zeta_cut){
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



void function_rad_pr_font_fix(int &n,double *x,double *fvec,struct auxarray_rad &aux){

    double F_rad_x = x[1]*(SQR(x[1]) + 1.0);
    double F_rad_y = x[2]*(SQR(x[2]) + 1.0);
    double F_rad_z = x[3]*(SQR(x[3]) + 1.0);
    //    double E_rad = x[4];


    double gijuiuj = aux.gupxx_phys*SQR(aux.u_x) + 2.0*aux.gupxy_phys*aux.u_x*aux.u_y +
      2.0*aux.gupxz_phys*aux.u_x*aux.u_z + aux.gupyy_phys*SQR(aux.u_y) + 2.0*aux.gupyz_phys*aux.u_y*aux.u_z +
      aux.gupzz_phys*SQR(aux.u_z);
    double au0m1 = gijuiuj/( 1.0+sqrt(1.0+gijuiuj) );
    //    if (aux.rho_s < 0.0) au0m1 = gijuiuj/( 1.0-sqrt(1.0+gijuiuj) );
    double u0 = (au0m1+1.0)/aux.alpn1;

    //    double rho_b =  fabs(aux.rho_s/(aux.alpn1*aux.Psi6*u0));
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
    if (FaFa <= 0.0 || E_rad <= aux.zeta_cut){
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
    
}



void jacobian_rad_pr_font_fix(int &n,double *x,struct auxarray_rad &aux, double *fvec,int &np, double fjac[][5]){
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
  //  if (aux.rho_s < 0.0) au0m1 = gijuiuj/( 1.0-sqrt(1.0+gijuiuj) );
  double u0 = (au0m1+1.0)/aux.alpn1;
  double u0_inv = aux.alpn1/(au0m1+1.0);

  //  double rho_b =  fabs(aux.rho_s/(aux.alpn1*aux.Psi6*u0));
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
  if (FaFa <= 0.0 || E_rad <= aux.zeta_cut){
    FaFa_inv == 0.0;
  }else{
    FaFa_inv = 1.0/FaFa;
  }
 
 
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



extern "C" void CCTK_FCALL CCTK_FNAME(primitive_vars_rad_pr)
  (int *ext,int *nghostzones, 
   double *tau_rad, double *S_rad_x, double *S_rad_y, double *S_rad_z,
   double *Sx, double *Sy, double *Sz, double *rho,
   double *Sxx, double *Sxy, double *Sxz, double *Syy, double *Syz, double *Szz,
   double *E_rad, double *F_radx, double *F_rady, double *F_radz, 
   double *phi, double *alpha, double *shiftx, double *shifty, double *shiftz, double *gxx, double *gxy, double *gxz, double *gyy, double *gyz, double *gzz,
   double *gupxx, double *gupxy, double *gupxz, double *gupyy, double *gupyz, double *gupzz,
   double *vx, double *vy, double *vz, double *u0,
   const cGH **cctkGH,int &ignore_ghostzones, int &repairs_rad_needed, double &Psi6threshold, double &Erad_atm_cut) {
  primitive_vars_rad_pr(ext,nghostzones,
			     tau_rad,S_rad_x,S_rad_y,S_rad_z,
			     Sx, Sy, Sz, rho,
			     Sxx, Sxy, Sxz, Syy, Syz, Szz,
			     E_rad,F_radx,F_rady,F_radz,
			     phi, alpha, shiftx, shifty, shiftz, gxx, gxy, gxz, gyy, gyz, gzz,
			     gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,
			     vx, vy, vz, u0,
			     *cctkGH, ignore_ghostzones, repairs_rad_needed, Psi6threshold, Erad_atm_cut);
}
