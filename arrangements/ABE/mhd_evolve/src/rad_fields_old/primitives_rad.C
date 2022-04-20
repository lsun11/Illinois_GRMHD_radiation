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
  double rho_s;
};


void newt2_cpp_rad(double x[],struct auxarray_rad &aux,
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
void lnsrch_rad(int n, double xold[], double fold, double g[], double p[], double x[],
		 double *f, double stpmax, bool &check, 
		 double (*fmin)
		 (double *x,struct auxarray_rad &aux,double *fvec,
		  void (*nrfuncv)(int &n,double *x,double *fvec, struct auxarray_rad &aux),
		  int &n),
		 void (*nrfuncv)(int &n,double *x,double *fvec,struct auxarray_rad &aux),
		 struct auxarray_rad &aux,double *fvec);
double fmin_rad(double *x,struct auxarray_rad &aux,double *fvec,
		 void (*nrfuncv)(int &n,double *x,double *fvec,struct auxarray_rad &aux),
		 int &n);

void compute_M1_prim(double &Pij,double &Fi,double &Fj,double &Fksq,double &E,double &gij,double &ui, double &uj, double &chi, double &psim4);
//double fasterpow_rad(double inputvar,double inputpow);

extern "C" void CCTK_FCALL CCTK_FNAME(primitive_vars_rad_cpp)
  (int *ext,int *nghostzones, double *X, double *Y, double *Z, 
   double *rho_star,
   double *tau_rad, double *S_rad_x, double *S_rad_y, double *S_rad_z,
   double *Sx, double *Sy, double *Sz, double *rho, 
   double *Sxx, double *Sxy, double *Sxz, double *Syy, double *Syz, double *Szz,
   double *E_rad,  double *F_radx, double *F_rady, double *F_radz, 
   double *P_radxx, double *P_radyy, double *P_radzz, double *P_radxy, double *P_radxz, double *P_radyz,
   double *phi, double *alpha, double *shiftx, double *shifty, double *shiftz, double *gxx, double *gxy, double *gxz, double *gyy, double *gyz, double *gzz,
   double *gupxx, double *gupxy, double *gupxz, double *gupyy, double *gupyz, double *gupzz,
   double *vx, double *vy, double *vz, double *u0,
   double *failure_tracker,const cGH **cctkGH,int &ignore_ghostzones, int &repairs_needed);


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
 double *E_rad, double *F_radx, double *F_rady, double *F_radz,
 double *P_radxx, double *P_radyy, double *P_radzz, double *P_radxy, double *P_radxz,double *P_radyz,
 double *phi, double *alpha, double *shiftx, double *shifty, double *shiftz, double *gxx, double *gxy, double *gxz, double *gyy, double *gyz, double *gzz,
 double *gupxx, double *gupxy, double *gupxz, double *gupyy, double *gupyz, double *gupzz,
 double *vx, double *vy, double *vz, double *u0,
 double *failure_tracker,const cGH *cctkGH,int &ignore_ghostzones, int &repairs_needed) {

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

  // Initialize the failure_tracker array
    #pragma omp parallel for
    for(int i=0;i<ext[0];i++) for(int j=0;j<ext[1];j++) for(int k=0;k<ext[2];k++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    failure_tracker[index] = 0.0;
    }
  repairs_needed = 0;

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

        if (compute_primitives==1) { 
	   int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

	   int MAXITS_beforefontfix=50;
	   double STPMX_beforefontfix=20.0; // <== faster than 100.0, when maxits set to 50

	   double au0m1;

	   //Structure that holds auxilliary data for Newton-Raphson solver
	   struct auxarray_rad AUX;
	   
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


           double tau_radl = tau_rad[index];
           double S_rad_xl = S_rad_x[index];
           double S_rad_yl = S_rad_y[index];
           double S_rad_zl = S_rad_z[index];


	   double u0L = u0[index];
	   double F_rad_x = (gxxL*F_radx[index] + 
	   		gxyL*F_rady[index] + 
	   		gxzL*F_radz[index])*Psi4;
	   double F_rad_y = (gxyL*F_radx[index] + 
	   		gyyL*F_rady[index] + 
	   		gyzL*F_radz[index])*Psi4;
	   double F_rad_z = (gxzL*F_radx[index] + 
	   		gyzL*F_rady[index] + 
	   		gzzL*F_radz[index])*Psi4;
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
	   
           double rho_s = rho_star[index];
          

	   double shift_xL = Psi4*(shiftxL*gxxL + shiftyL*gxyL + shiftzL*gxzL);
	   double shift_yL = Psi4*(shiftxL*gxyL + shiftyL*gyyL + shiftzL*gyzL); 
	   double shift_zL = Psi4*(shiftxL*gxzL + shiftyL*gyzL + shiftzL*gzzL); 

	   double beta2 = shiftxL*shift_xL + shiftyL*shift_yL + shiftzL*shift_zL;  
	   double udotbeta = u0L*(vx[index]*shift_xL + vy[index]*shift_yL + vz[index]*shift_zL);
	   double g_00L =beta2-alpn1*alpn1;
	   double u_0L = g_00L*u0L + udotbeta;
	  

	     int nn = 4;
	     double UU[5];

	     double E_radl;

	     UU[1] = F_rad_x;
	     UU[2] = F_rad_y;
	     UU[3] = F_rad_z;
	     UU[4] = E_rad[index];

	     AUX.rho_s = rho_s;
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


	     UUguess[1] = UU[1];
	     UUguess[2] = UU[2];
	     UUguess[3] = UU[3];
	     UUguess[4] = UU[4];

	    	     
	     bool check;
	     int newt_indx[5];
	     double newt_g[5],newt_p[5],newt_xold[5],newt_fvec[5];

	     newt2_cpp_rad(UU,AUX,function_rad,jacobian_rad,nn,check, 
		       newt_indx,newt_g,newt_p,newt_xold,newt_fvec,MAXITS_beforefontfix,STPMX_beforefontfix);
	   
	     int MAXITS_taustildefix=MAXITS_beforefontfix*10; //  <-- Try harder when doing a Font fix; we don't want to fail here!!!
	     double STPMX_taustildefix=STPMX_beforefontfix*10; //  <-- Try harder when doing a Font fix; we don't want to fail here!!!
	    
	    
	       while(check || UU[4] < 0.0) {
		 check = false;
		 newt2_cpp_rad(UU,AUX,function_rad,jacobian_rad,nn,check,
			   newt_indx,newt_g,newt_p,newt_xold,newt_fvec,MAXITS_taustildefix,STPMX_taustildefix);
		 MAXITS_taustildefix*=10;  // <-- Try even harder!  Don't give up!
		 //if(MAXITS_taustildefix>=MAXITS_beforefontfix*1e6) break;
		 if(MAXITS_taustildefix>=MAXITS_beforefontfix*1e4) break;
	       }
	     
	     bool recom = false;

	     //****************************************************************
	     //                          FONT FIX
	     // Impose Font fix when Newt-Raph fails (check = true) or when it
	     //  gives negative eps [UU[4] < 0]
	     //****************************************************************
	     if(check || UU[4] < 0.0) {
               
	       double UU_font_fix[4];

	       UU_font_fix[1] = 1.0;
	       UU_font_fix[2] = 1.0;
	       UU_font_fix[3] = 1.0;
	       // Following if statement will make the Fixed: line compatible with DAGH version in axisymmetry
	       //                 if(i > 1 && k > 1) then
	       //count = count + 1;
	       //                 }
	       nn = 3;
	         
	       check = false;
	       int MAXITS_fontfix=MAXITS_beforefontfix*10; //  <-- Try harder when doing a Font fix; we don't want to fail here!!!
	       double STPMX_fontfix=STPMX_beforefontfix*10; //  <-- Try harder when doing a Font fix; we don't want to fail here!!!
	       newt2_cpp_rad(UU_font_fix,AUX,function_rad_font_fix,jacobian_rad_font_fix,nn,check,
			 newt_indx,newt_g,newt_p,newt_xold,newt_fvec,MAXITS_fontfix,STPMX_fontfix);
	       while(check) {
		 check = false;
		 newt2_cpp_rad(UU_font_fix,AUX,function_rad_font_fix,jacobian_rad_font_fix,nn,check,
			   newt_indx,newt_g,newt_p,newt_xold,newt_fvec,MAXITS_fontfix,STPMX_fontfix);
		 MAXITS_fontfix*=10;  // <-- Try even harder!  Don't give up!
		 if(MAXITS_fontfix>=MAXITS_beforefontfix*1e6) break;
		 //if(MAXITS_fontfix>=MAXITS_beforefontfix*1e4) break;
	       }
	       if(check) {
	         printf("ERROR: FONT FIX (secondary solver) JUST FAILED\n");
	         printf("Problem at (x,y,z) = %e %e %e, %d %d %d\n", X[index],Y[index],Z[index],i,j,k);
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
                       
	         //FIXME: FAILURE_TRACKER ARRAY NOT SET:
	         failure_tracker[index] = 1.0;
	         //repairs_needed = true;
   	         repairs_needed = 1;
	         // Set everything to 0.0 before calculating 
	         // these quantities from averages.
	         F_rad_x         = 0.0;
	         F_rad_y         = 0.0;
	         F_rad_z         = 0.0;
	         recom        = true;
	       } else {
	         //if the Font fix worked, do the following:                 
	         F_rad_x = UU_font_fix[1]*(SQR(UU_font_fix[1]) + 1.0);
	         F_rad_y = UU_font_fix[2]*(SQR(UU_font_fix[2]) + 1.0);
	         F_rad_z = UU_font_fix[3]*(SQR(UU_font_fix[3]) + 1.0);
	         recom = true;    
	       }
	       //************************************************************************************************************** 
	     } else {
	       // Inversion worked without the Font fix!  Now set the primitives:
	       F_rad_x = UU[1];
	       F_rad_y = UU[2];
	       F_rad_z = UU[3];
	       E_radl  = UU[4];
	     }
	     

	   //Add radiation terms to the source.
	     double F_rad_0l = - (F_rad_x * AUX.ux + F_rad_y * AUX.uy + F_rad_z * AUX.uz)/u0L;
	     double F_rad0l = (-F_rad_0l+F_rad_x*AUX.shiftx+F_rad_y*AUX.shifty+F_rad_z*AUX.shiftz)/SQR(AUX.alpn1);
	     double F_radxl = F_rad_0l*AUX.shiftx/SQR(AUX.alpn1) + F_rad_x*AUX.gupxx_phys + F_rad_y*AUX.gupxy_phys +F_rad_z*AUX.gupxz_phys;
	     double F_radyl = F_rad_0l*AUX.shifty/SQR(AUX.alpn1) + F_rad_x*AUX.gupxy_phys + F_rad_y*AUX.gupyy_phys +F_rad_z*AUX.gupyz_phys;
	     double F_radzl = F_rad_0l*AUX.shiftz/SQR(AUX.alpn1) + F_rad_x*AUX.gupxz_phys + F_rad_y*AUX.gupyz_phys +F_rad_z*AUX.gupzz_phys;
	     double Fksq = F_rad_x*F_radxl +  F_rad_y*F_radyl +  F_rad_z*F_radzl;

	     double zeta = (F_rad_0l*F_rad0l + F_rad_x*F_radxl +  F_rad_y*F_radyl +  F_rad_z*F_radzl )/SQR(E_radl);
	     double chi = 1/3.0 + SQR(zeta)*(6-2*zeta+6*SQR(zeta))/15.0;

	     double P_radxxl,P_radyyl,P_radzzl,P_radxyl,P_radxzl,P_radyzl;
	     // Since we use gupij_phys = gupij * psim4, we set psim4 = 1.0 in this function.
	     double psim4_temp = 1.0;
	     compute_M1_prim(P_radxxl, F_radxl, F_radxl, Fksq, E_radl, AUX.gupxx_phys, AUX.ux, AUX.ux, chi, psim4_temp);
	     compute_M1_prim(P_radyyl, F_radyl, F_radyl, Fksq, E_radl, AUX.gupyy_phys, AUX.uy, AUX.uy, chi, psim4_temp);
	     compute_M1_prim(P_radzzl, F_radzl, F_radzl, Fksq, E_radl, AUX.gupzz_phys, AUX.uz, AUX.uz, chi, psim4_temp);
	     compute_M1_prim(P_radxyl, F_radxl, F_radyl, Fksq, E_radl, AUX.gupxy_phys, AUX.ux, AUX.uy, chi, psim4_temp);
	     compute_M1_prim(P_radxzl, F_radxl, F_radzl, Fksq, E_radl, AUX.gupxz_phys, AUX.ux, AUX.uz, chi, psim4_temp);	       
             compute_M1_prim(P_radyzl, F_radyl, F_radzl, Fksq, E_radl, AUX.gupyz_phys, AUX.uy, AUX.uz, chi, psim4_temp);
  

	     double P_rad0xl = - (P_radxxl * AUX.u_x + P_radxyl * AUX.u_y + P_radxzl * AUX.u_z)/AUX.u_0;
	     double P_rad0yl = - (P_radxyl * AUX.u_x + P_radyyl * AUX.u_y + P_radyzl * AUX.u_z)/AUX.u_0;
	     double P_rad0zl = - (P_radxzl * AUX.u_x + P_radyzl * AUX.u_y + P_radzzl * AUX.u_z)/AUX.u_0;

	     double P_rad00l = - (P_rad0xl * AUX.u_x + P_rad0yl * AUX.u_y + P_rad0zl * AUX.u_z)/AUX.u_0;


	   E_rad[index] = E_radl;
	   F_radx[index] = F_radxl;
	   F_rady[index] = F_radyl;
	   F_radz[index] = F_radzl;

	   P_radxx[index] = P_radxxl;
	   P_radyy[index] = P_radyyl;
	   P_radzz[index] = P_radzzl;
	   P_radxy[index] = P_radxyl;
	   P_radxz[index] = P_radxzl;
	   P_radyz[index] = P_radyzl;

	   rho[index] = rho[index] + SQR(alpn1)*(E_radl*SQR(u0L) + 2.0*F_rad0l*u0L + P_rad00l);
	   
	   Sx[index] = Sx[index] + alpn1*(u0L*E_radl*AUX.u_x + F_rad0l*AUX.u_x + u0L*F_rad_x + 
					  P_rad00l*shift_xL + Psi4*(P_rad0xl*gxxL + P_rad0yl*gxyL + P_rad0zl*gxzL));
					 
	   Sy[index] = Sy[index] + alpn1*(u0L*E_radl*AUX.u_y + F_rad0l*AUX.u_y + u0L*F_rad_y +
                                          P_rad00l*shift_yL + Psi4*(P_rad0xl*gxyL + P_rad0yl*gyyL + P_rad0zl*gyzL));

	   Sz[index] = Sz[index] + alpn1*(u0L*E_radl*AUX.u_z + F_rad0l*AUX.u_z + u0L*F_rad_z +
                                          P_rad00l*shift_zL + Psi4*(P_rad0xl*gxzL + P_rad0yl*gyzL + P_rad0zl*gzzL));


	  
	   Sxx[index] = Sxx[index] + E_radl*AUX.u_x*AUX.u_x + 2*F_rad_x*AUX.u_x + 
	                             SQR(shift_xL)*P_rad00l + shift_xL*2.0*(gxxL*P_rad0xl+gxyL*P_rad0yl+gxzL*P_rad0zl) + 
	                             + SQR(Psi4)*( SQR(gxxL)*P_radxxl + SQR(gxyL)*P_radyyl + SQR(gxzL)*P_radzzl +
				     2.0*(gxxL*gxyL*P_radxyl+gxxL*gxzL*P_radxzl+gxyL*gxzL*P_radyzl) );

	   Syy[index] = Syy[index] + E_radl*AUX.u_y*AUX.u_y + 2*F_rad_y*AUX.u_y +
	                             SQR(shift_yL)*P_rad00l + shift_yL*2.0*(gxyL*P_rad0xl+gyyL*P_rad0yl+gyzL*P_rad0zl) + 
	                             SQR(Psi4)*( SQR(gxyL)*P_radxxl + SQR(gyyL)*P_radyyl + SQR(gyzL)*P_radzzl +
				     2.0*(gxyL*gyyL*P_radxyl+gxyL*gyzL*P_radxzl+gyyL*gyzL*P_radyzl) );
	    
	   Szz[index] = Szz[index] + E_radl*AUX.u_z*AUX.u_z + 2*F_rad_z*AUX.u_z +
	                             SQR(shift_zL)*P_rad00l + shift_zL*2.0*(gxzL*P_rad0xl+gxzL*P_rad0yl+gzzL*P_rad0zl) + 
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
			             (gxxL*gxzL + gxyL*gyzL)*P_radxyl + (gxxL*gzzL + gxzL*gxzL)*P_radxzl + (gxyL*gzzL + gxzL*gyzL)*P_radyzl);

	   Syz[index] = Syz[index] + E_radl*AUX.u_y*AUX.u_z + F_rad_y*AUX.u_z + F_rad_z*AUX.u_y +
                                     shift_yL*shift_zL*P_rad00l +
	                             Psi4*(shift_yL*(gxzL*P_rad0xl + gyzL*P_rad0yl + gzzL*P_rad0zl) + (shift_zL*(gxyL*P_rad0xl + gyyL*P_rad0yl + gyzL*P_rad0zl)) )+
	                             SQR(Psi4)*(gxyL*gxzL*P_radxxl + gyyL*gyzL*P_radyyl + gyzL*gzzL*P_radzzl +
			             (gxyL*gyzL + gxyL*gyzL)*P_radxyl + (gxyL*gzzL + gyzL*gxzL)*P_radxzl + (gyyL*gzzL + gyzL*gyzL)*P_radyzl);
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
#define TOLF 1.0e-8
#define TOLMIN 1.0e-12
#define TOLX 3.0e-16
/*
#define TOLF 1.0e-13
#define TOLMIN 1.0e-15
#define TOLX 1.0e-13
*/
//#define STPMX 10000.0
//#define STPMX 100.0
#define FREERETURN {return;}
//#define FREERETURN {free_dmatrix_newt(fjac,1,NP,1,NP);return;}

void newt2_cpp_rad(double x[],struct auxarray_rad &aux,
	       void (*funcv1)(int &n,double *x,double *fvec,struct auxarray_rad &aux),
	       void (*fdjac1)
	       (int &n,double *x,struct auxarray_rad &aux, double *fvec,int &np,double fjac[][5]),
	       int &n,bool &check,
	       int *indx,double *g,double *p,double *xold,double *fvec,int &MAXITS,double &STPMX)
{
  void lubksb_rad(double a[][5], int n, int *indx, double b[]);
  void ludcmp_rad(double a[][5], int n, int *indx, double *d);
  int i,its,j;
  double d,den,f,fold,stpmax,sum,temp,test,fjac[5][5];


  //printf ("INSIDE newt2_cpp_rad, before fmin_rad, fvec[1]=%e, fvec[2]=%e,fvec[3]=%e, fvec[4]=%e, ux=x[1]=%e \n", fvec[1], fvec[2],fvec[3], fvec[4], x[1]);
 

  f=fmin_rad(x,aux,fvec,funcv1,n);


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
    ludcmp_rad(fjac,n,indx,&d);
    lubksb_rad(fjac,n,indx,p);
    
    
    lnsrch_rad(n,xold,fold,g,p,x,&f,stpmax,check,fmin_rad,funcv1,aux,fvec);

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
      if(check!=0) printf("BAD CHECK %d: %e\t%e\n",check,test,TOLMIN);
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
  //printf("MAXITS exceeded in newt\n");
  //exit(0);
}

double fmin_rad(double *x,struct auxarray_rad &aux,double *fvec,
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
//#define TOLX 1.0e-7
#define TOLX 3.0e-16
//#define TOLX 1.0e-15

void lnsrch_rad(int n, double xold[], double fold, double g[], double p[], double x[],
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

  //  printf("START lnsrch_rad!!!!!\n");

  check=0;
  for (sum=0.0,i=1;i<=n;i++) sum += p[i]*p[i];
  sum=sqrt(sum);
  if (sum > stpmax)
    for (i=1;i<=n;i++) p[i] *= stpmax/sum;
  for (slope=0.0,i=1;i<=n;i++)
    slope += g[i]*p[i];
  if (slope >= 0.0) printf("OUCH!  Roundoff problem in lnsrch.  Note that this function has been updated to correct a bug in Numerical Recipes 2.06.  This is version 2.08, with the 2.06->2.08 diff obtained from http://www.numerical-recipes.com/upgrade/upgrade-208.html\n");
  test=0.0;
  for (i=1;i<=n;i++) {
    temp=fabs(p[i])/max_val_rad(fabs(xold[i]),1.0);
    if (temp > test) test=temp;
  }
  alamin=TOLX/test;
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
    alam=max_val_rad(tmplam,0.1*alam);
    //    alam=FMAX(tmplam,0.1*alam);
  }

}
#undef ALF
#undef TOLX

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
	for(int www=0;www<10000;www++) {
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
  double F_radx = F_rad_0*aux.shiftx/SQR(aux.alpn1) + F_rad_x*aux.gupxx_phys + F_rad_y*aux.gupxy_phys +F_rad_z*aux.gupxz_phys;
  double F_rady = F_rad_0*aux.shifty/SQR(aux.alpn1) + F_rad_x*aux.gupxy_phys + F_rad_y*aux.gupyy_phys +F_rad_z*aux.gupyz_phys;
  double F_radz = F_rad_0*aux.shiftz/SQR(aux.alpn1) + F_rad_x*aux.gupxz_phys + F_rad_y*aux.gupyz_phys +F_rad_z*aux.gupzz_phys;

  double FaFa = F_radx*F_rad_x + F_rady*F_rad_y + F_radz*F_rad_z + F_rad0*F_rad_0; // F_rad^alpha*F_rad_alpha
  double FkFk = F_radx*F_rad_x + F_rady*F_rad_y + F_radz*F_rad_z;

  double FaFa_o_E2 = FaFa/SQR(E_rad);
  double A = FaFa_o_E2*(3.0 - sqrt(FaFa_o_E2) + 3.0*FaFa_o_E2)/5.0;

  double hxx = aux.gupxx_phys + SQR(aux.ux);
  double hyy = aux.gupyy_phys + SQR(aux.uy);
  double hzz = aux.gupzz_phys + SQR(aux.uz);
  double hxy = aux.gupxy_phys + aux.ux*aux.uy;
  double hxz = aux.gupxz_phys + aux.ux*aux.uz;
  double hyz = aux.gupyz_phys + aux.uy*aux.uz;

  
  double Bxx = A * (SQR(F_radx)/FkFk - hxx/3.0) + hxx/3.0;
  double Byy = A * (SQR(F_rady)/FkFk - hyy/3.0) + hyy/3.0;
  double Bzz = A * (SQR(F_radz)/FkFk - hzz/3.0) + hzz/3.0;
  double Bxy = A * (F_radx*F_rady/FkFk - hxy/3.0) + hxy/3.0;
  double Bxz = A * (F_radx*F_radz/FkFk - hxz/3.0) + hxz/3.0;
  double Byz = A * (F_rady*F_radz/FkFk - hyz/3.0) + hyz/3.0;

  double Pxx = E_rad * Bxx;
  double Pyy = E_rad * Byy;
  double Pzz = E_rad * Bzz;
  double Pxy = E_rad * Bxy;
  double Pxz = E_rad * Bxz;
  double Pyz = E_rad * Byz;

  double uiujPij = Pxx*SQR(aux.u_x) + Pyy*SQR(aux.u_y) + Pzz*SQR(aux.u_z) +
    2.0*( Pxy*aux.u_x*aux.u_y +  Pxz*aux.u_x*aux.u_z + Pyz*aux.u_y*aux.u_z) ;
  
  double hx_x = 1.0 + aux.ux * aux.u_x ;
  double hy_y = 1.0 + aux.uy * aux.u_y ;
  double hz_z = 1.0 + aux.uz * aux.u_z ;
  double hx_y = aux.ux*aux.u_y;
  double hx_z = aux.ux*aux.u_z;
  double hy_z = aux.uy*aux.u_z;
  double hy_x = aux.uy*aux.u_x;
  double hz_x = aux.uz*aux.u_x;
  double hz_y = aux.uz*aux.u_y;  

  double Bx_x = A * (F_radx*F_rad_x/FkFk - hx_x/3.0) + hx_x/3.0;
  double By_y = A * (F_rady*F_rad_y/FkFk - hy_y/3.0) + hy_y/3.0;
  double Bz_z = A * (F_radz*F_rad_z/FkFk - hz_z/3.0) + hz_z/3.0;
  double Bx_y = A * (F_radx*F_rad_y/FkFk - hx_y/3.0)- hx_y/3.0;
  double Bx_z = A * (F_radx*F_rad_z/FkFk - hx_z/3.0)- hx_z/3.0;
  double By_z = A * (F_rady*F_rad_z/FkFk - hy_z/3.0)- hy_z/3.0;
  double By_x = A * (F_rady*F_rad_x/FkFk - hy_x/3.0)- hy_x/3.0;
  double Bz_x = A * (F_radz*F_rad_x/FkFk - hz_x/3.0)- hz_x/3.0;
  double Bz_y = A * (F_radz*F_rad_y/FkFk - hz_y/3.0)- hz_y/3.0;


  double Px_x = E_rad * Bx_x;
  double Py_y = E_rad * By_y;
  double Pz_z = E_rad * Bz_z;
  double Px_y = E_rad * Bx_y;
  double Px_z = E_rad * Bx_z;
  double Py_z = E_rad * By_z;
  double Py_x = E_rad * By_x;
  double Pz_x = E_rad * Bz_x;
  double Pz_y = E_rad * Bz_y;


  double uiPi_x = aux.ux*Px_x + aux.uy*Py_x + aux.uz*Pz_x;
  double uiPi_y = aux.ux*Px_y + aux.uy*Py_y + aux.uz*Pz_y;
  double uiPi_z = aux.ux*Px_z + aux.uy*Py_z + aux.uz*Pz_z;

  //
  // fvec(1): Eq. for S_rad_x; fvec[2]: Eq. for S_rad_y; fvec[3]: Eq. for S_rad_z; 
  // fvec[4]: Eq. for tau_rad.
  //
  fvec[1] = aux.alpn1*aux.Psi6*(E_rad*u0*aux.u_x + F_rad0*aux.u_x + F_rad_x*u0 - uiPi_x/aux.u_0) - aux.S_rad_x;
  fvec[2] = aux.alpn1*aux.Psi6*(E_rad*u0*aux.u_y + F_rad0*aux.u_y + F_rad_y*u0 - uiPi_y/aux.u_0) - aux.S_rad_y;
  fvec[3] = aux.alpn1*aux.Psi6*(E_rad*u0*aux.u_z + F_rad0*aux.u_z + F_rad_z*u0 - uiPi_z/aux.u_0) - aux.S_rad_z;
  fvec[4] = SQR(aux.alpn1)*aux.Psi6*(E_rad*SQR(u0) + 2.0*F_rad0*u0 + uiujPij/SQR(aux.u_0)) - aux.tau_rad;

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
  double F_radx = F_rad_0*aux.shiftx/SQR(aux.alpn1) + F_rad_x*aux.gupxx_phys + F_rad_y*aux.gupxy_phys +F_rad_z*aux.gupxz_phys;
  double F_rady = F_rad_0*aux.shifty/SQR(aux.alpn1) + F_rad_x*aux.gupxy_phys + F_rad_y*aux.gupyy_phys +F_rad_z*aux.gupyz_phys;
  double F_radz = F_rad_0*aux.shiftz/SQR(aux.alpn1) + F_rad_x*aux.gupxz_phys + F_rad_y*aux.gupyz_phys +F_rad_z*aux.gupzz_phys;


  double FaFa = F_radx*F_rad_x + F_rady*F_rad_y + F_radz*F_rad_z + F_rad0*F_rad_0; // F_rad^alpha*F_rad_alpha                                  
  double FkFk = F_radx*F_rad_x + F_rady*F_rad_y + F_radz*F_rad_z;

  double FaFa_o_E2 = FaFa/SQR(E_rad);
  double A = FaFa_o_E2*(3.0 - sqrt(FaFa_o_E2) + 3.0*FaFa_o_E2)/5.0;
  // D is for Jacobian only!!!                                                                                                                  
  double D = (-3.0*FaFa/SQR(E_rad) + 2.0*pow(FaFa, 1.5)/pow(E_rad, 3.0) - 9.0*SQR(FaFa)/pow(E_rad,4.0))*0.2;

  double dF_0dF_x = -aux.ux*u0_inv;
  double dF_0dF_y = -aux.uy*u0_inv;
  double dF_0dF_z = -aux.uz*u0_inv;

  double dF0dF_x = (aux.ux*aux.shiftx/(SQR(aux.alpn1)*u0)-aux.gupxx_phys)*aux.u_x 
                 + (aux.ux*aux.shifty/(SQR(aux.alpn1)*u0)-aux.gupxy_phys)*aux.u_y
                 + (aux.ux*aux.shiftz/(SQR(aux.alpn1)*u0)-aux.gupxz_phys)*aux.u_z;
  double dF0dF_y = (aux.uy*aux.shiftx/(SQR(aux.alpn1)*u0)-aux.gupxy_phys)*aux.u_x
                 + (aux.uy*aux.shifty/(SQR(aux.alpn1)*u0)-aux.gupyy_phys)*aux.u_y
                 + (aux.uy*aux.shiftz/(SQR(aux.alpn1)*u0)-aux.gupyz_phys)*aux.u_z;
  double dF0dF_z = (aux.uz*aux.shiftx/(SQR(aux.alpn1)*u0)-aux.gupxz_phys)*aux.u_x
                 + (aux.uz*aux.shifty/(SQR(aux.alpn1)*u0)-aux.gupyz_phys)*aux.u_y
                 + (aux.uz*aux.shiftz/(SQR(aux.alpn1)*u0)-aux.gupzz_phys)*aux.u_z;

  double Cx = 0.2*(3.0/SQR(E_rad) - 1.5*sqrt(FaFa)/pow(E_rad, 3.0) + 6.0*FaFa/pow(E_rad,4.0))*(F_rad0*dF_0dF_x + F_rad_0*dF0dF_x + F_radx) ;
  double Cy = 0.2*(3.0/SQR(E_rad) - 1.5*sqrt(FaFa)/pow(E_rad, 3.0) + 6.0*FaFa/pow(E_rad,4.0))*(F_rad0*dF_0dF_y + F_rad_0*dF0dF_y + F_rady) ;
  double Cz = 0.2*(3.0/SQR(E_rad) - 1.5*sqrt(FaFa)/pow(E_rad, 3.0) + 6.0*FaFa/pow(E_rad,4.0))*(F_rad0*dF_0dF_z + F_rad_0*dF0dF_z + F_radz) ;

  double hxx = aux.gupxx_phys + SQR(aux.ux);
  double hyy = aux.gupyy_phys + SQR(aux.uy);
  double hzz = aux.gupzz_phys + SQR(aux.uz);
  double hxy = aux.gupxy_phys + aux.ux*aux.uy;
  double hxz = aux.gupxz_phys + aux.ux*aux.uz;
  double hyz = aux.gupyz_phys + aux.uy*aux.uz;


  double Bxx = A * (SQR(F_radx)/FkFk - hxx/3.0) + hxx/3.0;
  double Byy = A * (SQR(F_rady)/FkFk - hyy/3.0) + hyy/3.0;
  double Bzz = A * (SQR(F_radz)/FkFk - hzz/3.0) + hzz/3.0;
  double Bxy = A * (F_radx*F_rady/FkFk - hxy/3.0) + hxy/3.0;
  double Bxz = A * (F_radx*F_radz/FkFk - hxz/3.0) + hxz/3.0;
  double Byz = A * (F_rady*F_radz/FkFk - hyz/3.0) + hyz/3.0;

  double Pxx = E_rad * Bxx;
  double Pyy = E_rad * Byy;
  double Pzz = E_rad * Bzz;
  double Pxy = E_rad * Bxy;
  double Pxz = E_rad * Bxz;
  double Pyz = E_rad * Byz;


  double uiujhij = hxx*SQR(aux.u_x) + hyy*SQR(aux.u_y) + hzz*SQR(aux.u_z) +
    2.0*( hxy*aux.u_x*aux.u_y +  hxz*aux.u_x*aux.u_z + hyz*aux.u_y*aux.u_z) ;
  double uiujFiFj = SQR(aux.u_x*F_radx) + SQR(aux.u_y*F_rady) + SQR(aux.u_z*F_radz) +
    2.0*(F_radx*F_rady*aux.u_x*aux.u_y +  F_radx*F_radz*aux.u_x*aux.u_z + F_rady*F_radz*aux.u_y*aux.u_z) ;

 
  double hx_x = 1.0 + aux.ux * aux.u_x ;
  double hy_y = 1.0 + aux.uy * aux.u_y ;
  double hz_z = 1.0 + aux.uz * aux.u_z ;
  double hx_y = aux.ux*aux.u_y;
  double hx_z = aux.ux*aux.u_z;
  double hy_z = aux.uy*aux.u_z;
  double hy_x = aux.uy*aux.u_x;
  double hz_x = aux.uz*aux.u_x;
  double hz_y = aux.uz*aux.u_y;

  double Bx_x = E_rad * (F_radx*F_rad_x/FkFk - (1.0 + hx_x)/3.0);
  double By_y = E_rad * (F_rady*F_rad_y/FkFk - (1.0 + hy_y)/3.0);
  double Bz_z = E_rad * (F_radz*F_rad_z/FkFk - (1.0 + hz_z)/3.0);
  double Bx_y = E_rad * (F_radx*F_rad_y/FkFk - hx_y/3.0);
  double Bx_z = E_rad * (F_radx*F_rad_z/FkFk - hx_z/3.0);
  double By_z = E_rad * (F_rady*F_rad_z/FkFk - hy_z/3.0);
  double By_x = E_rad * (F_rady*F_rad_x/FkFk - hy_x/3.0);
  double Bz_x = E_rad * (F_radz*F_rad_x/FkFk - hz_x/3.0);
  double Bz_y = E_rad * (F_radz*F_rad_y/FkFk - hz_y/3.0);


  double Px_x = E_rad * Bx_x;
  double Py_y = E_rad * By_y;
  double Pz_z = E_rad * Bz_z;
  double Px_y = E_rad * Bx_y;
  double Px_z = E_rad * Bx_z;
  double Py_z = E_rad * By_z;
  double Py_x = E_rad * By_x;
  double Pz_x = E_rad * Bz_x;
  double Pz_y = E_rad * Bz_y;


  double uihi_x = aux.u_x*hx_x + aux.u_y*hy_x + aux.u_z*hz_x;
  double uihi_y = aux.u_x*hx_y + aux.u_y*hy_y + aux.u_z*hz_y;
  double uihi_z = aux.u_x*hx_z + aux.u_y*hy_z + aux.u_z*hz_z;

  double uiFiF_x = aux.u_x*F_radx*F_rad_x +aux.u_y*F_rady*F_rad_x +aux.u_z*F_radz*F_rad_x; 
  double uiFiF_y = aux.u_x*F_radx*F_rad_y +aux.u_y*F_rady*F_rad_y +aux.u_z*F_radz*F_rad_y;
  double uiFiF_z = aux.u_x*F_radx*F_rad_z +aux.u_y*F_rady*F_rad_z +aux.u_z*F_radz*F_rad_z;


  double duiujFiFjdFx = 2.0*(SQR(aux.u_x)*F_radx + aux.u_x*(aux.u_y*F_rady + aux.u_z*F_radz) ); 
  double duiujFiFjdFy = 2.0*(SQR(aux.u_y)*F_rady + aux.u_y*(aux.u_x*F_radx + aux.u_z*F_radz) );
  double duiujFiFjdFz = 2.0*(SQR(aux.u_z)*F_radz + aux.u_z*(aux.u_x*F_radx + aux.u_y*F_rady) );

  double vx = aux.ux/u0;
  double vy = aux.uy/u0;
  double vz = aux.uz/u0;

  // f(1) = S_rad_x; f(2) = S_rad_y; f(3) = S_rad_z; f(4) = tau_rad;
  // x[1] = F_rad_x; x[2] = F_rad_y; x[3] = F_rad_z; x[4] = E_rad;
  // fjac(i,j) = partial f(i) / partial x(j) 
  //

  fjac[4][4] = SQR(aux.alpn1)*aux.Psi6*(SQR(u0) + (uiujhij*(1.0-D)/3.0 + uiujFiFj*D/FkFk)/SQR(aux.u_0) );


  fjac[1][4] = aux.alpn1*aux.Psi6*(u0*aux.u_x - (uihi_x*(1.0-D)/3.0 + uiFiF_x*D/FkFk)/aux.u_0);

  fjac[2][4] = aux.alpn1*aux.Psi6*(u0*aux.u_y - (uihi_y*(1.0-D)/3.0 + uiFiF_y*D/FkFk)/aux.u_0);

  fjac[3][4] = aux.alpn1*aux.Psi6*(u0*aux.u_z - (uihi_z*(1.0-D)/3.0 + uiFiF_z*D/FkFk)/aux.u_0);


  fjac[4][1] = SQR(aux.alpn1)*aux.Psi6*(2.0*u0*dF0dF_x + E_rad/SQR(aux.u_0) * (A*(duiujFiFjdFx/FkFk - uiujFiFj/SQR(FkFk)*F_radx) +  (uiujFiFj/FkFk-uiujhij/3.0) * Cx));

  fjac[4][2] = SQR(aux.alpn1)*aux.Psi6*(2.0*u0*dF0dF_y + E_rad/SQR(aux.u_0) * (A*(duiujFiFjdFy/FkFk - uiujFiFj/SQR(FkFk)*F_rady) +  (uiujFiFj/FkFk-uiujhij/3.0) * Cy));

  fjac[4][3] = SQR(aux.alpn1)*aux.Psi6*(2.0*u0*dF0dF_z + E_rad/SQR(aux.u_0) * (A*(duiujFiFjdFz/FkFk - uiujFiFj/SQR(FkFk)*F_radz) +  (uiujFiFj/FkFk-uiujhij/3.0) * Cz));



  fjac[1][1] = aux.alpn1*aux.Psi6*(u0 + aux.u_x*dF0dF_x + E_rad*( Cx*(F_rad0*F_rad_x/FkFk - u0*aux.u_x/3.0) + A*(F_rad0*F_rad_x/SQR(FkFk)*F_radx +
													       (F_rad_x*(vx + aux.shiftx)/SQR(aux.alpn1) + F_rad0))));					   
  fjac[2][2] = aux.alpn1*aux.Psi6*(u0 + aux.u_y*dF0dF_y + E_rad*( Cy*(F_rad0*F_rad_y/FkFk - u0*aux.u_y/3.0) + A*(F_rad0*F_rad_y/SQR(FkFk)*F_rady +
														(F_rad_y*(vy + aux.shifty)/SQR(aux.alpn1) + F_rad0))));
  fjac[3][3] = aux.alpn1*aux.Psi6*(u0 + aux.u_z*dF0dF_z + E_rad*( Cz*(F_rad0*F_rad_z/FkFk - u0*aux.u_z/3.0) + A*(F_rad0*F_rad_z/SQR(FkFk)*F_radz +
														(F_rad_z*(vz + aux.shiftz)/SQR(aux.alpn1) + F_rad0))));
							

  fjac[1][2] = aux.alpn1*aux.Psi6*(aux.u_x*dF0dF_y + E_rad*( Cy*(F_rad0*F_rad_x/FkFk - u0*aux.u_x/3.0) + A*(F_rad0*F_rad_x/SQR(FkFk)*F_rady +
													   F_rad_x*(vy + aux.shifty)/SQR(aux.alpn1) )));

  fjac[1][3] = aux.alpn1*aux.Psi6*(aux.u_x*dF0dF_z + E_rad*( Cz*(F_rad0*F_rad_x/FkFk - u0*aux.u_x/3.0) + A*(F_rad0*F_rad_x/SQR(FkFk)*F_radz +
                                                                                                           F_rad_x*(vz + aux.shiftz)/SQR(aux.alpn1) )));


  fjac[2][1] = aux.alpn1*aux.Psi6*(aux.u_y*dF0dF_x + E_rad*( Cx*(F_rad0*F_rad_y/FkFk - u0*aux.u_y/3.0) + A*(F_rad0*F_rad_y/SQR(FkFk)*F_radx +
                                                                                                           F_rad_y*(vx + aux.shiftx)/SQR(aux.alpn1) )));


  fjac[3][1] = aux.alpn1*aux.Psi6*(aux.u_z*dF0dF_x + E_rad*( Cx*(F_rad0*F_rad_z/FkFk - u0*aux.u_z/3.0) + A*(F_rad0*F_rad_z/SQR(FkFk)*F_radx +
													    F_rad_z*(vx + aux.shiftx)/SQR(aux.alpn1) )));


  fjac[2][3] = aux.alpn1*aux.Psi6*(aux.u_y*dF0dF_z + E_rad*( Cz*(F_rad0*F_rad_y/FkFk - u0*aux.u_y/3.0) + A*(F_rad0*F_rad_y/SQR(FkFk)*F_radz +
													    F_rad_y*(vz + aux.shiftz)/SQR(aux.alpn1) )));
	
	
  fjac[3][2] = aux.alpn1*aux.Psi6*(aux.u_z*dF0dF_y + E_rad*( Cy*(F_rad0*F_rad_z/FkFk - u0*aux.u_z/3.0) + A*(F_rad0*F_rad_z/SQR(FkFk)*F_rady +
													    F_rad_z*(vy + aux.shifty)/SQR(aux.alpn1) )));



}



void function_rad_font_fix(int &n,double *x,double *fvec,struct auxarray_rad &aux){

    double F_rad_x = x[1]*(SQR(x[1]) + 1.0);
    double F_rad_y = x[2]*(SQR(x[2]) + 1.0);
    double F_rad_z = x[3]*(SQR(x[3]) + 1.0);
    double E_rad = x[4];

    double gijuiuj = aux.gupxx_phys*SQR(aux.u_x) + 2.0*aux.gupxy_phys*aux.u_x*aux.u_y +
      2.0*aux.gupxz_phys*aux.u_x*aux.u_z + aux.gupyy_phys*SQR(aux.u_y) + 2.0*aux.gupyz_phys*aux.u_y*aux.u_z +
      aux.gupzz_phys*SQR(aux.u_z);
    double au0m1 = gijuiuj/( 1.0+sqrt(1.0+gijuiuj) );
    if (aux.rho_s < 0.0) au0m1 = gijuiuj/( 1.0-sqrt(1.0+gijuiuj) );
    double u0 = (au0m1+1.0)/aux.alpn1;


    double F_rad_0 = - (F_rad_x * aux.ux + F_rad_y * aux.uy + F_rad_z * aux.uz)/u0;
    double F_rad0 = (-F_rad_0+F_rad_x*aux.shiftx+F_rad_y*aux.shifty+F_rad_z*aux.shiftz)/SQR(aux.alpn1);
    double F_radx = F_rad_0*aux.shiftx/SQR(aux.alpn1) + F_rad_x*aux.gupxx_phys + F_rad_y*aux.gupxy_phys +F_rad_z*aux.gupxz_phys;
    double F_rady = F_rad_0*aux.shifty/SQR(aux.alpn1) + F_rad_x*aux.gupxy_phys + F_rad_y*aux.gupyy_phys +F_rad_z*aux.gupyz_phys;
    double F_radz = F_rad_0*aux.shiftz/SQR(aux.alpn1) + F_rad_x*aux.gupxz_phys + F_rad_y*aux.gupyz_phys +F_rad_z*aux.gupzz_phys;


    double FaFa = F_radx*F_rad_x + F_rady*F_rad_y + F_radz*F_rad_z + F_rad0*F_rad_0; // F_rad^alpha*F_rad_alpha                                                                           
    double FkFk = F_radx*F_rad_x + F_rady*F_rad_y + F_radz*F_rad_z;

    double FaFa_o_E2 = FaFa/SQR(E_rad);
    double A = FaFa_o_E2*(3.0 - sqrt(FaFa_o_E2) + 3.0*FaFa_o_E2)/5.0;

    double hxx = aux.gupxx_phys + SQR(aux.ux);
    double hyy = aux.gupyy_phys + SQR(aux.uy);
    double hzz = aux.gupzz_phys + SQR(aux.uz);
    double hxy = aux.gupxy_phys + aux.ux*aux.uy;
    double hxz = aux.gupxz_phys + aux.ux*aux.uz;
    double hyz = aux.gupyz_phys + aux.uy*aux.uz;

    double Bxx = A * (SQR(F_radx)/FkFk - hxx/3.0) + hxx/3.0;
    double Byy = A * (SQR(F_rady)/FkFk - hyy/3.0) + hyy/3.0;
    double Bzz = A * (SQR(F_radz)/FkFk - hzz/3.0) + hzz/3.0;
    double Bxy = A * (F_radx*F_rady/FkFk - hxy/3.0) + hxy/3.0;
    double Bxz = A * (F_radx*F_radz/FkFk - hxz/3.0) + hxz/3.0;
    double Byz = A * (F_rady*F_radz/FkFk - hyz/3.0) + hyz/3.0;

    double Pxx = E_rad * Bxx;
    double Pyy = E_rad * Byy;
    double Pzz = E_rad * Bzz;
    double Pxy = E_rad * Bxy;
    double Pxz = E_rad * Bxz;
    double Pyz = E_rad * Byz;

    double uiujPij = Pxx*SQR(aux.u_x) + Pyy*SQR(aux.u_y) + Pzz*SQR(aux.u_z) +
      2.0*( Pxy*aux.u_x*aux.u_y +  Pxz*aux.u_x*aux.u_z + Pyz*aux.u_y*aux.u_z) ;

    double hx_x = 1.0 + aux.ux * aux.u_x ;
    double hy_y = 1.0 + aux.uy * aux.u_y ;
    double hz_z = 1.0 + aux.uz * aux.u_z ;
    double hx_y = aux.ux*aux.u_y;
    double hx_z = aux.ux*aux.u_z;
    double hy_z = aux.uy*aux.u_z;
    double hy_x = aux.uy*aux.u_x;
    double hz_x = aux.uz*aux.u_x;
    double hz_y = aux.uz*aux.u_y;

    double Bx_x = A * (F_radx*F_rad_x/FkFk - hx_x/3.0) + hx_x/3.0;
    double By_y = A * (F_rady*F_rad_y/FkFk - hy_y/3.0) + hy_y/3.0;
    double Bz_z = A * (F_radz*F_rad_z/FkFk - hz_z/3.0) + hz_z/3.0;
    double Bx_y = A * (F_radx*F_rad_y/FkFk - hx_y/3.0)- hx_y/3.0;
    double Bx_z = A * (F_radx*F_rad_z/FkFk - hx_z/3.0)- hx_z/3.0;
    double By_z = A * (F_rady*F_rad_z/FkFk - hy_z/3.0)- hy_z/3.0;
    double By_x = A * (F_rady*F_rad_x/FkFk - hy_x/3.0)- hy_x/3.0;
    double Bz_x = A * (F_radz*F_rad_x/FkFk - hz_x/3.0)- hz_x/3.0;
    double Bz_y = A * (F_radz*F_rad_y/FkFk - hz_y/3.0)- hz_y/3.0;


    double Px_x = E_rad * Bx_x;
    double Py_y = E_rad * By_y;
    double Pz_z = E_rad * Bz_z;
    double Px_y = E_rad * Bx_y;
    double Px_z = E_rad * Bx_z;
    double Py_z = E_rad * By_z;
    double Py_x = E_rad * By_x;
    double Pz_x = E_rad * Bz_x;
    double Pz_y = E_rad * Bz_y;

    double uiPi_x = aux.ux*Px_x + aux.uy*Py_x + aux.uz*Pz_x;
    double uiPi_y = aux.ux*Px_y + aux.uy*Py_y + aux.uz*Pz_y;
    double uiPi_z = aux.ux*Px_z + aux.uy*Py_z + aux.uz*Pz_z;


    //                                                             
    // fvec(1): Eq. for S_rad_x; fvec[2]: Eq. for S_rad_y; fvec[3]: Eq. for S_rad_z;
    // fvec[4]: Eq. for tau_rad.                    
    //         
    fvec[1] = aux.alpn1*aux.Psi6*(E_rad*u0*aux.u_x + F_rad0*aux.u_x + F_rad_x*u0 - uiPi_x/aux.u_0) - aux.S_rad_x;
    fvec[2] = aux.alpn1*aux.Psi6*(E_rad*u0*aux.u_y + F_rad0*aux.u_y + F_rad_y*u0 - uiPi_y/aux.u_0) - aux.S_rad_y;
    fvec[3] = aux.alpn1*aux.Psi6*(E_rad*u0*aux.u_z + F_rad0*aux.u_z + F_rad_z*u0 - uiPi_z/aux.u_0) - aux.S_rad_z;
    fvec[4] = SQR(aux.alpn1)*aux.Psi6*(E_rad*SQR(u0) + 2.0*F_rad0*u0 + uiujPij/SQR(aux.u_0)) - aux.tau_rad;

}



void jacobian_rad_font_fix(int &n,double *x,struct auxarray_rad &aux, double *fvec,int &np, double fjac[][5]){
  //
  double F_rad_x  = x[1]*(SQR(x[1]) + 1.0);
  double F_rad_y  = x[2]*(SQR(x[2]) + 1.0);
  double F_rad_z  = x[3]*(SQR(x[3]) + 1.0);
  double E_rad = x[4];
  double fac[4];
  fac[1] = 3.0*SQR(x[1]) + 1.0;
  fac[2] = 3.0*SQR(x[2]) + 1.0;
  fac[3] = 3.0*SQR(x[3]) + 1.0;

  double gijuiuj = aux.gupxx_phys*SQR(aux.u_x) + 2.0*aux.gupxy_phys*aux.u_x*aux.u_y +
    2.0*aux.gupxz_phys*aux.u_x*aux.u_z + aux.gupyy_phys*SQR(aux.u_y) + 2.0*aux.gupyz_phys*aux.u_y*aux.u_z +
    aux.gupzz_phys*SQR(aux.u_z);
  double au0m1 = gijuiuj/( 1.0+sqrt(1.0+gijuiuj) );
  if (aux.rho_s < 0.0) au0m1 = gijuiuj/( 1.0-sqrt(1.0+gijuiuj) );
  double u0 = (au0m1+1.0)/aux.alpn1;
  double u0_inv = aux.alpn1/(au0m1+1.0);

  double F_rad_0 = - (F_rad_x * aux.ux + F_rad_y * aux.uy + F_rad_z * aux.uz)/u0;
  double F_rad0 = (-F_rad_0+F_rad_x*aux.shiftx+F_rad_y*aux.shifty+F_rad_z*aux.shiftz)/SQR(aux.alpn1);
  double F_radx = F_rad_0*aux.shiftx/SQR(aux.alpn1) + F_rad_x*aux.gupxx_phys + F_rad_y*aux.gupxy_phys +F_rad_z*aux.gupxz_phys;
  double F_rady = F_rad_0*aux.shifty/SQR(aux.alpn1) + F_rad_x*aux.gupxy_phys + F_rad_y*aux.gupyy_phys +F_rad_z*aux.gupyz_phys;
  double F_radz = F_rad_0*aux.shiftz/SQR(aux.alpn1) + F_rad_x*aux.gupxz_phys + F_rad_y*aux.gupyz_phys +F_rad_z*aux.gupzz_phys;

  double FaFa = F_radx*F_rad_x + F_rady*F_rad_y + F_radz*F_rad_z + F_rad0*F_rad_0; // F_rad^alpha*F_rad_alpha                                                                           
  double FkFk = F_radx*F_rad_x + F_rady*F_rad_y + F_radz*F_rad_z;

  double FaFa_o_E2 = FaFa/SQR(E_rad);
  double A = FaFa_o_E2*(3.0 - sqrt(FaFa_o_E2) + 3.0*FaFa_o_E2)/5.0;
  // D is for Jacobian only!!!                                                                                        
  double D = (-3.0*FaFa/SQR(E_rad) + 2.0*pow(FaFa, 1.5)/pow(E_rad, 3.0) - 9.0*SQR(FaFa)/pow(E_rad,4.0))*0.2;

  double dF_0dF_x = -aux.ux*u0_inv;
  double dF_0dF_y = -aux.uy*u0_inv;
  double dF_0dF_z = -aux.uz*u0_inv;

  double dF0dF_x = (aux.ux*aux.shiftx/(SQR(aux.alpn1)*u0)-aux.gupxx_phys)*aux.u_x
    + (aux.ux*aux.shifty/(SQR(aux.alpn1)*u0)-aux.gupxy_phys)*aux.u_y
    + (aux.ux*aux.shiftz/(SQR(aux.alpn1)*u0)-aux.gupxz_phys)*aux.u_z;
  double dF0dF_y = (aux.uy*aux.shiftx/(SQR(aux.alpn1)*u0)-aux.gupxy_phys)*aux.u_x
    + (aux.uy*aux.shifty/(SQR(aux.alpn1)*u0)-aux.gupyy_phys)*aux.u_y
    + (aux.uy*aux.shiftz/(SQR(aux.alpn1)*u0)-aux.gupyz_phys)*aux.u_z;
  double dF0dF_z = (aux.uz*aux.shiftx/(SQR(aux.alpn1)*u0)-aux.gupxz_phys)*aux.u_x
    + (aux.uz*aux.shifty/(SQR(aux.alpn1)*u0)-aux.gupyz_phys)*aux.u_y
    + (aux.uz*aux.shiftz/(SQR(aux.alpn1)*u0)-aux.gupzz_phys)*aux.u_z;

  double Cx = 0.2*(2.0/SQR(E_rad) - 1.5*sqrt(FaFa)/pow(E_rad, 3.0) + 6.0*FaFa/pow(E_rad,4.0))*(F_rad0*dF_0dF_x + F_rad_0*dF0dF_x + F_radx) ;
  double Cy = 0.2*(2.0/SQR(E_rad) - 1.5*sqrt(FaFa)/pow(E_rad, 3.0) + 6.0*FaFa/pow(E_rad,4.0))*(F_rad0*dF_0dF_y + F_rad_0*dF0dF_y + F_rady) ;
  double Cz = 0.2*(2.0/SQR(E_rad) - 1.5*sqrt(FaFa)/pow(E_rad, 3.0) + 6.0*FaFa/pow(E_rad,4.0))*(F_rad0*dF_0dF_z + F_rad_0*dF0dF_z + F_radz) ;

  double hxx = aux.gupxx_phys + SQR(aux.ux);
  double hyy = aux.gupyy_phys + SQR(aux.uy);
  double hzz = aux.gupzz_phys + SQR(aux.uz);
  double hxy = aux.gupxy_phys + aux.ux*aux.uy;
  double hxz = aux.gupxz_phys + aux.ux*aux.uz;
  double hyz = aux.gupyz_phys + aux.uy*aux.uz;

  double Bxx = E_rad * (SQR(F_radx)/FkFk - hxx/3.0);
  double Byy = E_rad * (SQR(F_rady)/FkFk - hyy/3.0);
  double Bzz = E_rad * (SQR(F_radz)/FkFk - hzz/3.0);
  double Bxy = E_rad * (F_radx*F_rady/FkFk - hxy/3.0);
  double Bxz = E_rad * (F_radx*F_radz/FkFk - hxz/3.0);
  double Byz = E_rad * (F_rady*F_radz/FkFk - hyz/3.0);

  double uiujhij = hxx*SQR(aux.u_x) + hyy*SQR(aux.u_y) + hzz*SQR(aux.u_z) +
    2.0*( hxy*aux.u_x*aux.u_y +  hxz*aux.u_x*aux.u_z + hyz*aux.u_y*aux.u_z) ;
  double uiujFiFj = SQR(aux.u_x*F_radx) + SQR(aux.u_y*F_rady) + SQR(aux.u_z*F_radz) +
    2.0*(F_radx*F_rady*aux.u_x*aux.u_y +  F_radx*F_radz*aux.u_x*aux.u_z + F_rady*F_radz*aux.u_y*aux.u_z) ;


  double hx_x = 1.0 + aux.ux * aux.u_x ;
  double hy_y = 1.0 + aux.uy * aux.u_y ;
  double hz_z = 1.0 + aux.uz * aux.u_z ;
  double hx_y = aux.ux*aux.u_y;
  double hx_z = aux.ux*aux.u_z;
  double hy_z = aux.uy*aux.u_z;
  double hy_x = aux.uy*aux.u_x;
  double hz_x = aux.uz*aux.u_x;
  double hz_y = aux.uz*aux.u_y;

  double Bx_x = E_rad * (F_radx*F_rad_x/FkFk - (1.0 + hx_x)/3.0);
  double By_y = E_rad * (F_rady*F_rad_y/FkFk - (1.0 + hy_y)/3.0);
  double Bz_z = E_rad * (F_radz*F_rad_z/FkFk - (1.0 + hz_z)/3.0);
  double Bx_y = E_rad * (F_radx*F_rad_y/FkFk - hx_y/3.0);
  double Bx_z = E_rad * (F_radx*F_rad_z/FkFk - hx_z/3.0);
  double By_z = E_rad * (F_rady*F_rad_z/FkFk - hy_z/3.0);
  double By_x = E_rad * (F_rady*F_rad_x/FkFk - hy_x/3.0);
  double Bz_x = E_rad * (F_radz*F_rad_x/FkFk - hz_x/3.0);
  double Bz_y = E_rad * (F_radz*F_rad_y/FkFk - hz_y/3.0);


  double Px_x = E_rad * Bx_x;
  double Py_y = E_rad * By_y;
  double Pz_z = E_rad * Bz_z;
  double Px_y = E_rad * Bx_y;
  double Px_z = E_rad * Bx_z;
  double Py_z = E_rad * By_z;
  double Py_x = E_rad * By_x;
  double Pz_x = E_rad * Bz_x;
  double Pz_y = E_rad * Bz_y;


  double uihi_x = aux.u_x*hx_x + aux.u_y*hy_x + aux.u_z*hz_x;
  double uihi_y = aux.u_x*hx_y + aux.u_y*hy_y + aux.u_z*hz_y;
  double uihi_z = aux.u_x*hx_z + aux.u_y*hy_z + aux.u_z*hz_z;

  double uiFiF_x = aux.u_x*F_radx*F_rad_x +aux.u_y*F_rady*F_rad_x +aux.u_z*F_radz*F_rad_x;
  double uiFiF_y = aux.u_x*F_radx*F_rad_y +aux.u_y*F_rady*F_rad_y +aux.u_z*F_radz*F_rad_y;
  double uiFiF_z = aux.u_x*F_radx*F_rad_z +aux.u_y*F_rady*F_rad_z +aux.u_z*F_radz*F_rad_z;


  double duiujFiFjdFx = 2.0*(SQR(aux.u_x)*F_radx + aux.u_x*(aux.u_y*F_rady + aux.u_z*F_radz) );
  double duiujFiFjdFy = 2.0*(SQR(aux.u_y)*F_rady + aux.u_y*(aux.u_x*F_radx + aux.u_z*F_radz) );
  double duiujFiFjdFz = 2.0*(SQR(aux.u_z)*F_radz + aux.u_z*(aux.u_x*F_radx + aux.u_y*F_rady) );

  double vx = aux.ux/u0;
  double vy = aux.uy/u0;
  double vz = aux.uz/u0;

  fjac[1][1] = aux.alpn1*aux.Psi6*(u0 + aux.u_x*dF0dF_x + E_rad*( Cx*(F_rad0*F_rad_x/FkFk - u0*aux.u_x/3.0) + A*(F_rad0*F_rad_x/SQR(FkFk)*F_radx +
														 (F_rad_x*(vx + aux.shiftx)/SQR(aux.alpn1) + F_rad0))));
  fjac[2][2] = aux.alpn1*aux.Psi6*(u0 + aux.u_y*dF0dF_y + E_rad*( Cy*(F_rad0*F_rad_y/FkFk - u0*aux.u_y/3.0) + A*(F_rad0*F_rad_y/SQR(FkFk)*F_rady +
														 (F_rad_y*(vy + aux.shifty)/SQR(aux.alpn1) + F_rad0))));
  fjac[3][3] = aux.alpn1*aux.Psi6*(u0 + aux.u_z*dF0dF_z + E_rad*( Cz*(F_rad0*F_rad_z/FkFk - u0*aux.u_z/3.0) + A*(F_rad0*F_rad_z/SQR(FkFk)*F_radz +
														 (F_rad_z*(vz + aux.shiftz)/SQR(aux.alpn1) + F_rad0))));


  fjac[1][2] = aux.alpn1*aux.Psi6*(aux.u_x*dF0dF_y + E_rad*( Cy*(F_rad0*F_rad_x/FkFk - u0*aux.u_x/3.0) + A*(F_rad0*F_rad_x/SQR(FkFk)*F_rady +
													    F_rad_x*(vy + aux.shifty)/SQR(aux.alpn1) )));

  fjac[1][3] = aux.alpn1*aux.Psi6*(aux.u_x*dF0dF_z + E_rad*( Cz*(F_rad0*F_rad_x/FkFk - u0*aux.u_x/3.0) + A*(F_rad0*F_rad_x/SQR(FkFk)*F_radz +
													    F_rad_x*(vz + aux.shiftz)/SQR(aux.alpn1) )));


  fjac[2][1] = aux.alpn1*aux.Psi6*(aux.u_y*dF0dF_x + E_rad*( Cx*(F_rad0*F_rad_y/FkFk - u0*aux.u_y/3.0) + A*(F_rad0*F_rad_y/SQR(FkFk)*F_radx +
													    F_rad_y*(vx + aux.shiftx)/SQR(aux.alpn1) )));


  fjac[3][1] = aux.alpn1*aux.Psi6*(aux.u_z*dF0dF_x + E_rad*( Cx*(F_rad0*F_rad_z/FkFk - u0*aux.u_z/3.0) + A*(F_rad0*F_rad_z/SQR(FkFk)*F_radx +
                                                                                                            F_rad_z*(vx + aux.shiftx)/SQR(aux.alpn1) )));


  fjac[2][3] = aux.alpn1*aux.Psi6*(aux.u_y*dF0dF_z + E_rad*( Cz*(F_rad0*F_rad_y/FkFk - u0*aux.u_y/3.0) + A*(F_rad0*F_rad_y/SQR(FkFk)*F_radz +
                                                                                                            F_rad_y*(vz + aux.shiftz)/SQR(aux.alpn1) )));


  fjac[3][2] = aux.alpn1*aux.Psi6*(aux.u_z*dF0dF_y + E_rad*( Cy*(F_rad0*F_rad_z/FkFk - u0*aux.u_z/3.0) + A*(F_rad0*F_rad_z/SQR(FkFk)*F_rady +
                                                                                                            F_rad_z*(vy + aux.shifty)/SQR(aux.alpn1) )));
}





extern "C" void CCTK_FCALL CCTK_FNAME(primitive_vars_rad_cpp)
  (int *ext,int *nghostzones, double *X, double *Y, double *Z, 
   double *rho_star,
   double *tau_rad, double *S_rad_x, double *S_rad_y, double *S_rad_z,
   double *Sx, double *Sy, double *Sz, double *rho,
   double *Sxx, double *Sxy, double *Sxz, double *Syy, double *Syz, double *Szz,
   double *E_rad, double *F_radx, double *F_rady, double *F_radz,
   double *P_radxx, double *P_radyy, double *P_radzz, double *P_radxy, double *P_radxz,double *P_radyz,
   double *phi, double *alpha, double *shiftx, double *shifty, double *shiftz, double *gxx, double *gxy, double *gxz, double *gyy, double *gyz, double *gzz,
   double *gupxx, double *gupxy, double *gupxz, double *gupyy, double *gupyz, double *gupzz,
   double *vx, double *vy, double *vz, double *u0,
   double *failure_tracker,const cGH **cctkGH,int &ignore_ghostzones, int &repairs_needed) {
  primitive_vars_rad_cpp(ext,nghostzones, X, Y, Z,
			     rho_star,
			     tau_rad,S_rad_x,S_rad_y,S_rad_z,
			     Sx, Sy, Sz, rho,
			     Sxx, Sxy, Sxz, Syy, Syz, Szz,
			     E_rad,F_radx,F_rady,F_radz,
			     P_radxx, P_radyy, P_radzz, P_radxy, P_radxz, P_radyz,
			     phi, alpha, shiftx, shifty, shiftz, gxx, gxy, gxz, gyy, gyz, gzz,
			     gupxx, gupxy, gupxz, gupyy, gupyz, gupzz,
			     vx, vy, vz, u0, 
			     failure_tracker, *cctkGH,ignore_ghostzones, repairs_needed);
}

void compute_M1_prim(double &Pij,double &Fi,double &Fj,double &Fksq,double &E,double &gupij,double &ui, double &uj, double &chi, double &psim4) {
  double Pij_thin = Fi*Fj/Fksq*E;
  double Pij_thick = (psim4*gupij+ui*uj)*E/3.0;
  Pij = Pij_thin*(3.0*chi -1.0)/2.0 + Pij_thick*1.5*(1.0-chi);
}
