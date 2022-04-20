//-----------------------------------------------------------------------------
// routines for evolving the conducting fluid
// in the presence of a magnetic field
//-----------------------------------------------------------------------------
// Compute the flux for advecting rho_star, tau (Font's energy variable), 
//  and S_i .
//-----------------------------------------------------------------------------
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "cctk.h"
#include "primitives_solver_header.h"

#define FAC 0.99
#define SQR(x) ((x) * (x))

inline void find_cp_cm(double &cplus,double &cminus,double &v02,double &u0,
		       double &vi,double &lapse,double &shifti,double &psim4,double &gupii);

extern "C" void CCTK_FCALL mhdfluxnew_hybrid_
  (int *flux_direction, const cGH **cctkGH,int *cctk_lsh,int *nghostzones,
   int *enable_HARM_energyvariable,
   double *X,double *dX,double *dY,double *dZ,
   double *f_rho_i, double *f_rhoYe_i, double *f_tau_i, 
   double *f_Sx_i, double *f_Sy_i, double *f_Sz_i, 
   double *Pr, double *Pl, 
   double *rho_br, double *rho_bl, double *Y_er, double *Y_el,
   double *T_fluidl, double *T_fluidr,
   double *vxr, double *vxl, double *vyr, double *vyl, double *vzr, double *vzl, 
   double *Bxr, double *Bxl, double *Byr, double *Byl, double *Bzr, double *Bzl, 
   double *v02r, double *v02l, 
   double *gupxx_f, double *gupyy_f, double *gupzz_f, 
   double *cmax, double *cmin, 
   double *alpha_f, 
   double *betax_f, double *betay_f, double *betaz_f, 
   double *gxx_f, double *gxy_f, double *gxz_f, double *gyy_f, double *gyz_f, double *gzz_f, 
   double *phi_f, int *Symmetry, 
   int *enable_OS_collapse, double *rho_b_atm,
   int *neos, int *ergo_star, double *ergo_sigma, double *rho_tab, double *P_tab, double *eps_tab, double *k_tab, double *gamma_tab,double *gamma_th, 
   int *use_central_scheme_instead_of_hll, int *compute_microphysics);

double fasterpow_mhdflux(double inputvar,double inputpow);
void compute_eps_th(double &T_fluid, double &rho_b, double &eps_thermal);

extern "C" void mhdfluxnew_hybrid(int flux_direction, const cGH *cctkGH,int *cctk_lsh,int *nghostzones,
				  int enable_HARM_energyvariable,
				  double *X,double dX,double dY,double dZ,
				  double *f_rho_i,  double *f_rhoYe_i, double *f_tau_i, 
				  double *f_Sx_i, double *f_Sy_i, double *f_Sz_i, 
				  double *Pr, double *Pl, 
				  double *rho_br, double *rho_bl, double *Y_er, double *Y_el, 
				  double *T_fluidl, double *T_fluidr,
				  double *vxr, double *vxl, double *vyr, double *vyl, double *vzr, double *vzl, 
				  double *Bxr, double *Bxl, double *Byr, double *Byl, double *Bzr, double *Bzl, 
				  double *v02r, double *v02l, 
				  double *gupxx_f, double *gupyy_f, double *gupzz_f, 
				  double *cmax, double *cmin, 
				  double *alpha_f, 
				  double *betax_f, double *betay_f, double *betaz_f, 
				  double *gxx_f, double *gxy_f, double *gxz_f, double *gyy_f, double *gyz_f, double *gzz_f, 
				  double *phi_f, int Symmetry, 
				  int enable_OS_collapse,double rho_b_atm,
				  int neos, int ergo_star, double ergo_sigma, double *rho_tab, double *P_tab, double *eps_tab, double *k_tab, double *gamma_tab,double gamma_th,
				  int use_central_scheme_instead_of_hll, int compute_microphysics) {
  
  int AXISYM = 4;

  double f1os4pi = 1.0/sqrt(4.0*M_PI);

  int jstart = 0;
  int jend = cctk_lsh[1];
  
  if(Symmetry==4) {
    jstart++;
    jend--;
  }
  //  printf ("Inside mhdfluxnew_hybrid.C \n");
  //  printf("NGHOSTZONES = %d %d %d\n",nghostzones[0],nghostzones[1],nghostzones[2]);

  //printf ("neos = %d", neos);
  
#pragma omp parallel for
  for(int k=0;k<cctk_lsh[2];k++) for(int j=jstart;j<jend;j++) for(int i=0;i<cctk_lsh[0];i++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

    double alpha_fL = alpha_f[index];
    double al = 1.0 + alpha_fL;

    double phi_fL = phi_f[index];
    double psi2 = exp(2.0 * phi_fL);
    double psi4 = psi2*psi2;
    double psim4 = 1.0/(psi2*psi2);

    double gxx_fL = gxx_f[index];
    double gxy_fL = gxy_f[index];
    double gxz_fL = gxz_f[index];
    double gyy_fL = gyy_f[index];
    double gyz_fL = gyz_f[index];
    double gzz_fL = gzz_f[index];

    double gupxx_fL = gupxx_f[index];
    double gupyy_fL = gupyy_f[index];
    double gupzz_fL = gupzz_f[index];

    double betax_fL = betax_f[index];
    double betay_fL = betay_f[index];
    double betaz_fL = betaz_f[index];

    double vxrL = vxr[index];
    double vyrL = vyr[index];
    double vzrL = vzr[index];

    double vxlL = vxl[index];
    double vylL = vyl[index];
    double vzlL = vzl[index];

    // Compute al*u0-1
    double er = psi4*(gxx_fL*SQR(vxrL + betax_fL) +
		      2.0*gxy_fL*(vxrL + betax_fL)*(vyrL + betay_fL) +
		      2.0*gxz_fL*(vxrL + betax_fL)*(vzrL + betaz_fL) +
		      gyy_fL*SQR(vyrL + betay_fL) +
		      2.0*gyz_fL*(vyrL + betay_fL)*(vzrL + betaz_fL) +
		      gzz_fL*SQR(vzrL + betaz_fL) )/al/al;

    // *** Check for superluminal velocity ***
    if (er > 1.0) {
      double sqrtFACoer_r = sqrt(FAC/er);
      vxrL = (vxrL + betax_fL)*sqrtFACoer_r-betax_fL;
      vyrL = (vyrL + betay_fL)*sqrtFACoer_r-betay_fL;
      vzrL = (vzrL + betaz_fL)*sqrtFACoer_r-betaz_fL;
      er = FAC;
    }
    // ***************************************
    double el = sqrt(1.0-er);
    double au0r1 = er/el/(1.0+el);
    double u0r = (au0r1+1.0)/al;

    er = psi4*(gxx_fL*SQR(vxlL + betax_fL) +
	       2.0*gxy_fL*(vxlL + betax_fL)*(vylL + betay_fL) +
	       2.0*gxz_fL*(vxlL + betax_fL)*(vzlL + betaz_fL) +
	       gyy_fL*SQR(vylL + betay_fL) +
	       2.0*gyz_fL*(vylL + betay_fL)*(vzlL + betaz_fL) +
	       gzz_fL*SQR(vzlL + betaz_fL) )/al/al;

    // *** Check for superluminal velocity ***
    if (er > 1.0) {
      double sqrtFACoer_l = sqrt(FAC/er);
      vxlL = (vxlL + betax_fL)*sqrtFACoer_l-betax_fL;
      vylL = (vylL + betay_fL)*sqrtFACoer_l-betay_fL;
      vzlL = (vzlL + betaz_fL)*sqrtFACoer_l-betaz_fL;
      er = FAC;
    }



    // ***************************************
    el = sqrt(1.0-er);
    double au0l1 = er/el/(1.0+el);
    double u0l = (au0l1+1.0)/al;

    double rho_brL = rho_br[index];
    double T_fluidrL = T_fluidr[index];
    if (rho_brL < rho_b_atm) {rho_brL = rho_b_atm;}
    double rho_brLinv = 1.0/rho_brL;
    double PrL = Pr[index];
    double v02lL = v02l[index];
    double cplusl,dPcold_drho,hr1,eps_th;
    double P_over_rho = 1.0e-4;

  
    //Test:
    double eps_thermal_bhns = 1.5;

    
    if (rho_brL <= rho_tab[0]) {
     
      if (rho_brL == 0.0) {
	dPcold_drho = 0.0;
	eps_th = 0.0;	
	hr1 = 0.0;
      } else {
	v02lL = k_tab[0]*fasterpow_mhdflux(rho_brL,gamma_tab[0]);           // v02lL -> P_cold
	dPcold_drho = gamma_tab[0]*v02lL*rho_brLinv;
	cplusl = v02lL*rho_brLinv/(gamma_tab[0]-1.0);       // cplusl -> eps_cold
	if (compute_microphysics==0){
	  eps_th = (PrL - v02lL)/(gamma_th-1.0)*rho_brLinv;
	}
	else{
	  //	  eps_th = (eps_thermal_bhns -1.0) * cplusl;
	  compute_eps_th(T_fluidrL, rho_brL, eps_th);
	}
	}      
	hr1 = cplusl + eps_th + PrL*rho_brLinv;
      }
      

    for(int nn=1;nn<neos;nn++) {
      if (rho_brL <= rho_tab[nn]  &&  rho_brL > rho_tab[nn-1]) {
	v02lL = k_tab[nn]*fasterpow_mhdflux(rho_brL,gamma_tab[nn]);         // v02lL -> P_cold
	cplusl = eps_tab[nn-1] + (v02lL*rho_brLinv - P_tab[nn-1]/rho_tab[nn-1])/(gamma_tab[nn]-1.0);      // cplusl -> eps_cold
	dPcold_drho = gamma_tab[nn]*v02lL*rho_brLinv;
	if (compute_microphysics==0){
	  eps_th = (PrL - v02lL)/(gamma_th-1.0)*rho_brLinv;
	}
	else{
	  //	  eps_th = (eps_thermal_bhns -1.0) * cplusl;
	  compute_eps_th(T_fluidrL, rho_brL, eps_th);
	}
	hr1 = cplusl + eps_th + PrL*rho_brLinv;
      }
    }

    if (rho_brL > rho_tab[neos-1]) {
      if (ergo_star == 0){
	v02lL = k_tab[neos]*fasterpow_mhdflux(rho_brL,gamma_tab[neos]);    // v02lL -> P_cold
	dPcold_drho = gamma_tab[neos]*v02lL*rho_brLinv;
	cplusl = eps_tab[neos-1] + (v02lL*rho_brLinv - P_tab[neos-1]/rho_tab[neos-1])/(gamma_tab[neos]-1.0); // cplusl -> eps_cold
	if (compute_microphysics==0){
	  eps_th = (PrL - v02lL)/(gamma_th-1.0)*rho_brLinv;
	}
	else{
	  //	  eps_th = (eps_thermal_bhns -1.0) * cplusl;
	  compute_eps_th(T_fluidrL, rho_brL, eps_th);
	} 
      }
     else
       {
	 //printf("Ergo Star!");
	 v02lL = (ergo_sigma* (1+eps_tab[neos-1]+P_tab[neos-1]/rho_tab[neos-1])/fasterpow_mhdflux(rho_tab[neos-1],ergo_sigma)* fasterpow_mhdflux(rho_brL, ergo_sigma+1) + P_tab[neos-1] - ergo_sigma*((1+eps_tab[neos-1])*rho_tab[neos-1]) )/(ergo_sigma+1);
	cplusl = ((1+eps_tab[neos-1]+P_tab[neos-1]/rho_tab[neos-1])/fasterpow_mhdflux(rho_tab[neos-1],ergo_sigma) * fasterpow_mhdflux(rho_brL, ergo_sigma+1) - P_tab[neos-1] + ergo_sigma*((1+eps_tab[neos-1])*rho_tab[neos-1]) )/((ergo_sigma+1)*rho_brL)-1;
	dPcold_drho = ergo_sigma*(1+eps_tab[neos-1]+P_tab[neos-1]/rho_tab[neos-1])/fasterpow_mhdflux(rho_tab[neos-1],ergo_sigma)*rho_brL;
	if (compute_microphysics==0){
	  eps_th = (PrL - v02lL)/(gamma_th-1.0)*rho_brLinv;
	}else{
	  //	  eps_th = (eps_thermal_bhns -1.0) * cplusl;
	  compute_eps_th(T_fluidrL, rho_brL, eps_th);
	}
       }
      hr1 = cplusl + eps_th + PrL*rho_brLinv;
    }
    

    double BxrL = Bxr[index];
    double ByrL = Byr[index];
    double BzrL = Bzr[index];
 
    // Here sba = b^a defined in Gammie's paper
    double sbtr = f1os4pi*psi4*u0r/al*( gxx_fL*BxrL*(vxrL+betax_fL) + 
					gxy_fL*( BxrL*(vyrL+betay_fL) + ByrL*(vxrL+betax_fL) ) + 
					gxz_fL*( BxrL*(vzrL+betaz_fL) + BzrL*(vxrL+betax_fL) ) + 
					gyy_fL*ByrL*(vyrL+betay_fL) + 
					gyz_fL*( ByrL*(vzrL+betaz_fL) + BzrL*(vyrL+betay_fL) ) + 
					gzz_fL*BzrL*(vzrL+betaz_fL) );
    double sbxr = f1os4pi*BxrL/(al*u0r) + sbtr*vxrL;
    double sbyr = f1os4pi*ByrL/(al*u0r) + sbtr*vyrL;
    double sbzr = f1os4pi*BzrL/(al*u0r) + sbtr*vzrL;
    // Compute b^2
    double sb2r = -SQR(al*sbtr) + psi4*( gxx_fL*SQR(sbxr+betax_fL*sbtr) + 
					 2.0*gxy_fL*(sbxr+betax_fL*sbtr)*(sbyr+betay_fL*sbtr) + 
					 2.0*gxz_fL*(sbxr+betax_fL*sbtr)*(sbzr+betaz_fL*sbtr) + 
					 gyy_fL*SQR(sbyr+betay_fL*sbtr) + 
					 2.0*gyz_fL*(sbyr+betay_fL*sbtr)*(sbzr+betaz_fL*sbtr) + 
					 gzz_fL*SQR(sbzr+betaz_fL*sbtr) );

    // Compute v02r
    double cminusl;
    double v02rL = v02r[index];

    double cplusRl;
    if (rho_brL > 0.0) {
      cplusRl = (dPcold_drho + gamma_th*(gamma_th-1.0)*eps_th)/(1.0+hr1);
      cminusl = sb2r/(sb2r + rho_brL*(1.0+hr1));
      v02rL = cminusl + cplusRl*(1.0-cminusl);
    } else {
      v02rL = sb2r/(sb2r+1.e-300);
    }

    //    if(cplusRl==0 || v02rL==0) { printf(" Right: cplusRl=%e, cplusl=%e,cminusl=%e,v02rL=%e,eps_th=%e,hr1=%e, gamma_th=%e, rho_brLinv=%e, PrL=%e \n", cplusRl, cplusl,cminusl,v02rL,eps_th,hr1,gamma_th,rho_brLinv,PrL); }



    double hl1;
    double rho_blL = rho_bl[index];
    double T_fluidlL = T_fluidl[index];
    if (rho_blL <= rho_b_atm) {rho_blL = rho_b_atm;}
    double rho_blLinv = 1.0/rho_blL;
    double PlL = Pl[index];


    if (rho_blL <= rho_tab[0]) {
      if (rho_blL == 0.0) {
        dPcold_drho = 0.0;
        eps_th = 0.0;
        hl1 = 0.0;
      } else {
        v02lL = k_tab[0]*fasterpow_mhdflux(rho_blL,gamma_tab[0]);           // v02lL -> P_cold                                       
	dPcold_drho = gamma_tab[0]*v02lL*rho_blLinv;
	cplusl = v02lL*rho_blLinv/(gamma_tab[0]-1.0);       // cplusl -> eps_cold                                                     
	if (compute_microphysics==0){
	  eps_th = (PlL - v02lL)/(gamma_th-1.0)*rho_blLinv;        
	}else{
	  //	  eps_th = (eps_thermal_bhns -1.0) * cplusl;
	  compute_eps_th(T_fluidlL, rho_blL, eps_th);
       	}
	hl1 = cplusl + eps_th + PlL*rho_blLinv;
      }
    }


    for(int nn=1;nn<neos;nn++) {
      if (rho_blL <= rho_tab[nn]  &&  rho_blL > rho_tab[nn-1]) {
	v02lL = k_tab[nn]*fasterpow_mhdflux(rho_blL,gamma_tab[nn]);         // v02lL -> P_cold
	cplusl = eps_tab[nn-1] + (v02lL*rho_blLinv - P_tab[nn-1]/rho_tab[nn-1])/(gamma_tab[nn]-1.0);      // cplusl -> eps_cold
	dPcold_drho = gamma_tab[nn]*v02lL*rho_blLinv;
	if (compute_microphysics==0){
	  eps_th = (PlL - v02lL)/(gamma_th-1.0)*rho_blLinv;
	}else{
	  //	  eps_th = (eps_thermal_bhns -1.0) * cplusl;
	  compute_eps_th(T_fluidlL, rho_blL, eps_th);
	}
     	hl1 = cplusl + eps_th + PlL*rho_blLinv;
      }
    }
      if (rho_blL > rho_tab[neos-1]) {
	if (ergo_star == 0){
	  v02lL = k_tab[neos]*fasterpow_mhdflux(rho_blL,gamma_tab[neos]);    // v02lL -> P_cold  
	  dPcold_drho = gamma_tab[neos]*v02lL*rho_blLinv;
	  cplusl = eps_tab[neos-1] + (v02lL*rho_blLinv - P_tab[neos-1]/rho_tab[neos-1])/(gamma_tab[neos]-1.0); // cplusl -> eps_cold
	  if (compute_microphysics==0){
	    eps_th = (PlL - v02lL)/(gamma_th-1.0)*rho_blLinv;
	  }else{
	    //	    eps_th = (eps_thermal_bhns -1.0) * cplusl;
	    compute_eps_th(T_fluidlL, rho_blL, eps_th);
	  }
	}
     else
       {
	 //printf("Ergo Star! ---2 ");
	v02lL = (ergo_sigma* (1+eps_tab[neos-1]+P_tab[neos-1]/rho_tab[neos-1])/fasterpow_mhdflux(rho_tab[neos-1],ergo_sigma)* fasterpow_mhdflux(rho_blL, ergo_sigma+1) + P_tab[neos-1] - ergo_sigma*((1+eps_tab[neos-1])*rho_tab[neos-1]) )/(ergo_sigma+1);
	cplusl = ((1+eps_tab[neos-1]+P_tab[neos-1]/rho_tab[neos-1])/fasterpow_mhdflux(rho_tab[neos-1],ergo_sigma) * fasterpow_mhdflux(rho_blL, ergo_sigma+1) - P_tab[neos-1] + ergo_sigma*((1+eps_tab[neos-1])*rho_tab[neos-1]) )/((ergo_sigma+1)*rho_blL)-1;
        dPcold_drho = ergo_sigma*(1+eps_tab[neos-1]+P_tab[neos-1]/rho_tab[neos-1])/fasterpow_mhdflux(rho_tab[neos-1],ergo_sigma)*rho_blL;
	if (compute_microphysics==0){
	  eps_th = (PlL - v02lL)/(gamma_th-1.0)*rho_blLinv;
	}else{
	  //	  eps_th = (eps_thermal_bhns -1.0) * cplusl;
	  compute_eps_th(T_fluidlL, rho_blL, eps_th);
	}
       }
	hl1 = cplusl + eps_th + PlL*rho_blLinv;
    }
    



    double BxlL = Bxl[index];
    double BylL = Byl[index];

    double BzlL = Bzl[index];

    // Here sba = b^a defined in Gammie's paper
    double sbtl = f1os4pi*psi4*u0l/al*( gxx_fL*BxlL*(vxlL+betax_fL) + 
					gxy_fL*( BxlL*(vylL+betay_fL) + BylL*(vxlL+betax_fL) ) + 
					gxz_fL*( BxlL*(vzlL+betaz_fL) + BzlL*(vxlL+betax_fL) ) + 
					gyy_fL*BylL*(vylL+betay_fL) + 
					gyz_fL*( BylL*(vzlL+betaz_fL) + BzlL*(vylL+betay_fL) ) + 
					gzz_fL*BzlL*(vzlL+betaz_fL) );
    double sbxl = f1os4pi*BxlL/(al*u0l) + sbtl*vxlL;
    double sbyl = f1os4pi*BylL/(al*u0l) + sbtl*vylL;
    double sbzl = f1os4pi*BzlL/(al*u0l) + sbtl*vzlL;
    // Compute b^2
    double sb2l = -SQR(al*sbtl) + psi4*( gxx_fL*SQR(sbxl+betax_fL*sbtl) + 
					 2.0*gxy_fL*(sbxl+betax_fL*sbtl)*(sbyl+betay_fL*sbtl) + 
					 2.0*gxz_fL*(sbxl+betax_fL*sbtl)*(sbzl+betaz_fL*sbtl) + 
					 gyy_fL*SQR(sbyl+betay_fL*sbtl) + 
					 2.0*gyz_fL*(sbyl+betay_fL*sbtl)*(sbzl+betaz_fL*sbtl) + 
					 gzz_fL*SQR(sbzl+betaz_fL*sbtl) );



    // Compute v02lL
    double cplusLl;
    if (rho_blL > 0.0) {
      cplusLl = (dPcold_drho + gamma_th*(gamma_th-1.0)*eps_th)/(1.0+hl1);
      cminusl = sb2l/(sb2l + rho_blL*(1.0+hl1));
      v02lL = cminusl + cplusLl*(1.0-cminusl);
    } else {
      v02lL = sb2l/(sb2l+1.e-300);
    }
    
    
    if(cplusLl==0 || v02lL==0) { printf("Left: cplusLl=%e, cplusl=%e,cminusl=%e,v02lL=%e,eps_th=%e,hl1=%e, gamma_th=%e, rho_blL=%e, PlL=%e, gamma_tab[0]=%e, \n",cplusLl, cplusl,cminusl,v02lL,eps_th,hl1,gamma_th,rho_blL,PlL, gamma_tab[0]); }

    double cplusr,cminusr;
    if (flux_direction==1) {
      find_cp_cm(cplusl,cminusl,v02lL,u0l,
		 vxlL,alpha_fL,betax_fL,psim4,gupxx_fL);
      //if(cplusl==0 || cminusl==0) { printf("cplusl=%e,cminusl=%e,v02lL=%e,u0l=%e,vxlL=%e,alpha_fL=%e, gupxx_fL=%e \n",cplusl,cminusl,v02lL,u0l,vxlL,alpha_fL,gupxx_fL); }
      find_cp_cm(cplusr,cminusr,v02rL,u0r,
		 vxrL,alpha_fL,betax_fL,psim4,gupxx_fL);
    } else if (flux_direction==2) {
      find_cp_cm(cplusl,cminusl,v02lL,u0l,
		 vylL,alpha_fL,betay_fL,psim4,gupyy_fL);
      find_cp_cm(cplusr,cminusr,v02rL,u0r,
		 vyrL,alpha_fL,betay_fL,psim4,gupyy_fL);
    } else if (flux_direction==3){
      find_cp_cm(cplusl,cminusl,v02lL,u0l,
		 vzlL,alpha_fL,betaz_fL,psim4,gupzz_fL);
      find_cp_cm(cplusr,cminusr,v02rL,u0r,
		 vzrL,alpha_fL,betaz_fL,psim4,gupzz_fL);
    }


    // cminL = max(0.0,cplusl,cplusr)
    double cmaxL = 0.0;
    if(cmaxL < cplusl) cmaxL = cplusl;
    if(cmaxL < cplusr) cmaxL = cplusr;

    // cminL = -min(0.0,cminusl,cminusr)
    double cminL = 0.0;
    if(cminL > cminusl) cminL = cminusl;
    if(cminL > cminusr) cminL = cminusr;
    cminL *= -1;

    if (use_central_scheme_instead_of_hll==1) {
       double clax = cmaxL;
       if (cminL > clax) clax = cminL;
       cmaxL = clax; cminL = clax;
    }


    double psi6 = psi2*psi4;

    //*********************************************************************
    // density flux
    //*********************************************************************
    double Y_erL = Y_er[index];
    double Y_elL = Y_el[index];

    double rhor = al*psi6*rho_brL*u0r;                   
    double rhol = al*psi6*rho_blL*u0l; 

    
    double Fr,Fl;

    if (flux_direction==1) {
      Fr = al*psi6*rho_brL*vxrL*u0r;
      Fl = al*psi6*rho_blL*vxlL*u0l;
    } else if (flux_direction==2) {
      Fr = al*psi6*rho_brL*vyrL*u0r;
      Fl = al*psi6*rho_blL*vylL*u0l;
    } else if (flux_direction==3) {
      Fr = al*psi6*rho_brL*vzrL*u0r;
      Fl = al*psi6*rho_blL*vzlL*u0l;
    }
    

    // HLL step for rho:
    f_rho_i[index] = (cminL*Fr + cmaxL*Fl - cminL*cmaxL*(rhor-rhol) )/(cmaxL + cminL);

    if(isnan(f_rho_i[index])){
	printf("f_rho_i is nan, cminL=%e, cmaxL=%e, cplusl=%e, cplusr=%e, cminusl=%e, cminusr=%e, v02rL=%e, v02lL=%e \n", cminL, cmaxL, cplusl, cplusr, cminusl, cminusr, v02rL, v02lL);
      }

    double rhoYer = al*psi6*rho_brL*Y_erL*u0r;                                                                                                         
    double rhoYel = al*psi6*rho_blL*Y_elL*u0l;


    //    double Fr,Fl;
    if (flux_direction==1) {
      Fr = al*psi6*rho_brL*Y_erL*vxrL*u0r;
      Fl = al*psi6*rho_blL*Y_elL*vxlL*u0l;
    } else if (flux_direction==2) {
      Fr = al*psi6*rho_brL*Y_erL*vyrL*u0r;
      Fl = al*psi6*rho_blL*Y_elL*vylL*u0l;
    } else if (flux_direction==3) {
      Fr = al*psi6*rho_brL*Y_erL*vzrL*u0r;
      Fl = al*psi6*rho_blL*Y_elL*vzlL*u0l;
    }

    // HLL step for rho:                                                                                                                                                         
    f_rhoYe_i[index] = (cminL*Fr + cmaxL*Fl - cminL*cmaxL*(rhoYer-rhoYel) )/(cmaxL + cminL);

    if(isnan(f_rho_i[index])){
      printf("f_rho_i is nan, cminL=%e, cmaxL=%e, cplusl=%e, cplusr=%e, cminusl=%e, cminusr=%e, v02rL=%e, v02lL=%e \n", cminL, cmaxL, cplusl, cplusr, cminusl, cminusr, v02rL, \
	     v02lL);
    }


    //*********************************************************************
    // energy flux
    //*********************************************************************

    double u_xl = u0l*psi4*( gxx_fL*(vxlL+betax_fL) + gxy_fL*(vylL+betay_fL) + 
			     gxz_fL*(vzlL+betaz_fL) );
    double u_yl = u0l*psi4*( gxy_fL*(vxlL+betax_fL) + gyy_fL*(vylL+betay_fL) + 
			     gyz_fL*(vzlL+betaz_fL) );
    double u_zl = u0l*psi4*( gxz_fL*(vxlL+betax_fL) + gyz_fL*(vylL+betay_fL) + 
			     gzz_fL*(vzlL+betaz_fL) );

    double u_xr = u0r*psi4*( gxx_fL*(vxrL+betax_fL) + gxy_fL*(vyrL+betay_fL) + 
			     gxz_fL*(vzrL+betaz_fL) );
    double u_yr = u0r*psi4*( gxy_fL*(vxrL+betax_fL) + gyy_fL*(vyrL+betay_fL) + 
			     gyz_fL*(vzrL+betaz_fL) );
    double u_zr = u0r*psi4*( gxz_fL*(vxrL+betax_fL) + gyz_fL*(vyrL+betay_fL) + 
			     gzz_fL*(vzrL+betaz_fL) );

    double sb_xl = psi4*( gxx_fL*(sbxl+sbtl*betax_fL) + gxy_fL*(sbyl+sbtl*betay_fL) + 
			  gxz_fL*(sbzl+sbtl*betaz_fL) );
    double sb_yl = psi4*( gxy_fL*(sbxl+sbtl*betax_fL) + gyy_fL*(sbyl+sbtl*betay_fL) + 
			  gyz_fL*(sbzl+sbtl*betaz_fL) );
    double sb_zl = psi4*( gxz_fL*(sbxl+sbtl*betax_fL) + gyz_fL*(sbyl+sbtl*betay_fL) + 
			  gzz_fL*(sbzl+sbtl*betaz_fL) );

    double sb_xr = psi4*( gxx_fL*(sbxr+sbtr*betax_fL) + gxy_fL*(sbyr+sbtr*betay_fL) + 
			  gxz_fL*(sbzr+sbtr*betaz_fL) );
    double sb_yr = psi4*( gxy_fL*(sbxr+sbtr*betax_fL) + gyy_fL*(sbyr+sbtr*betay_fL) + 
			  gyz_fL*(sbzr+sbtr*betaz_fL) );
    double sb_zr = psi4*( gxz_fL*(sbxr+sbtr*betax_fL) + gyz_fL*(sbyr+sbtr*betay_fL) + 
			  gzz_fL*(sbzr+sbtr*betaz_fL) );

    if (enable_HARM_energyvariable==0) {
     
      el = (au0l1+hl1+au0l1*hl1)*rhol-psi6*PlL;
      if (flux_direction==1) {
	Fl = el*vxlL + PlL*psi6*(vxlL+betax_fL) + psi6*( SQR(al*u0l) *sb2l*vxlL 
							 + 0.5*betax_fL*sb2l - SQR(al)*sbtl*sbxl );
      } else if (flux_direction==2) {
	Fl = el*vylL + PlL*psi6*(vylL+betay_fL) + psi6*( SQR(al*u0l) *sb2l*vylL 
							 + 0.5*betay_fL*sb2l - SQR(al)*sbtl*sbyl );
      } else if (flux_direction==3){
	Fl = el*vzlL + PlL*psi6*(vzlL+betaz_fL) + psi6*( SQR(al*u0l) *sb2l*vzlL 
							 + 0.5*betaz_fL*sb2l - SQR(al)*sbtl*sbzl );
      }

      el = el + psi6*(sb2l*SQR(al*u0l) - 0.5*sb2l - SQR(al*sbtl));
      er = (au0r1+hr1+au0r1*hr1)*rhor-psi6*PrL;

      if (flux_direction==1) {
	Fr = er*vxrL + PrL*psi6*(vxrL+betax_fL) + psi6*( SQR(al*u0r) *sb2r*vxrL 
							 + 0.5*betax_fL*sb2r - SQR(al)*sbtr*sbxr );
      } else if (flux_direction==2) {
	Fr = er*vyrL + PrL*psi6*(vyrL+betay_fL) + psi6*( SQR(al*u0r) *sb2r*vyrL 
							 + 0.5*betay_fL*sb2r - SQR(al)*sbtr*sbyr );
      } else if (flux_direction==3){
	Fr = er*vzrL + PrL*psi6*(vzrL+betaz_fL) + psi6*( SQR(al*u0r) *sb2r*vzrL 
							 + 0.5*betaz_fL*sb2r - SQR(al)*sbtr*sbzr );
      }
      er = er + psi6*(sb2r*SQR(al*u0r) - 0.5*sb2r - SQR(al*sbtr));

    } else {
     
    // Temporarily store 1+u_0r and 1+u_0l in er and el
      el = -alpha_fL - al*au0l1 + u_xl*betax_fL + u_yl*betay_fL + u_zl*betaz_fL;
      er = -alpha_fL - al*au0r1 + u_xr*betax_fL + u_yr*betay_fL + u_zr*betaz_fL;
      // Temporarily store b_tr and b_tl to cplusr and cplusl
      cplusl = -al*al*sbtl + sb_xl*betax_fL + sb_yl*betay_fL + sb_zl*betaz_fL;
      cplusr = -al*al*sbtr + sb_xr*betax_fL + sb_yr*betay_fL + sb_zr*betaz_fL;
      // Now compute energy flux
      if (flux_direction==1) {
	Fl = -rhol*vxlL*(el + (el-1.0)*hl1) - al*psi6*( sb2l*u0l*vxlL*(el-1.0) - sbxl*cplusl );
	Fr = -rhor*vxrL*(er + (er-1.0)*hr1) - al*psi6*( sb2r*u0r*vxrL*(er-1.0) - sbxr*cplusr );
      } else if (flux_direction==2) {
	Fl = -rhol*vylL*(el + (el-1.0)*hl1) - al*psi6*( sb2l*u0l*vylL*(el-1.0) - sbyl*cplusl );
	Fr = -rhor*vyrL*(er + (er-1.0)*hr1) - al*psi6*( sb2r*u0r*vyrL*(er-1.0) - sbyr*cplusr );
      } else if (flux_direction==3){
	Fl = -rhol*vzlL*(el + (el-1.0)*hl1) - al*psi6*( sb2l*u0l*vzlL*(el-1.0) - sbzl*cplusl );
	Fr = -rhor*vzrL*(er + (er-1.0)*hr1) - al*psi6*( sb2r*u0r*vzrL*(er-1.0) - sbzr*cplusr );
      }
      // Finally, compute taur and taul and store them in el and er
      el = -rhol*(el + (el-1.0)*hl1)-al*psi6*(PlL + (u0l*(el-1.0)+0.5)*sb2l  
					      -sbtl*cplusl );
      er = -rhor*(er + (er-1.0)*hr1)-al*psi6*(PrL + (u0r*(er-1.0)+0.5)*sb2r 
					      -sbtr*cplusr );
    } // if(enable_HARM_energyvariable == 0)

    // HLL step for tau:
    f_tau_i[index] = (cminL*Fr + cmaxL*Fl - cminL*cmaxL*(er-el) )/(cmaxL + cminL);

    /*
    if (isnan(f_tau_i[index])){
      printf("inside mhdfluxnew-hybrid.C, f_tau_i[index] is nan \n");
      printf("cmax=%e, cmin=%e, Fr=%e, Fl=%e, er=%e, el=%e \n", cmax, cmin, Fr, Fl, er, el);
    }
    */

    if (Symmetry==AXISYM  &&  flux_direction==1) {
      double X_fL = X[index] - 0.5*dX;
      f_rho_i[index] = f_rho_i[index] * X_fL;
      f_rhoYe_i[index] =  f_rhoYe_i[index] * X_fL , 
      f_tau_i[index] = f_tau_i[index] * X_fL;
    }

    //*********************************************************************
    // Flux for S_j
    //*********************************************************************

    // Flux for S_x
    if (flux_direction==1) {
      Fl = al*psi6*( (rho_blL*(1.0+hl1)+sb2l)*u0l*vxlL*u_xl  
		     + PlL+0.5*sb2l - sbxl*sb_xl );
      Fr = al*psi6*( (rho_brL*(1.0+hr1)+sb2r)*u0r*vxrL*u_xr  
		     + PrL+0.5*sb2r - sbxr*sb_xr );
    } else if (flux_direction==2) {
      Fl = al*psi6*( (rho_blL*(1.0+hl1)+sb2l)*u0l*vylL*u_xl  
		     - sbyl*sb_xl );
      Fr = al*psi6*( (rho_brL*(1.0+hr1)+sb2r)*u0r*vyrL*u_xr  
		     - sbyr*sb_xr );
    } else if (flux_direction==3){
      Fl = al*psi6*( (rho_blL*(1.0+hl1)+sb2l)*u0l*vzlL*u_xl  
		     - sbzl*sb_xl );
      Fr = al*psi6*( (rho_brL*(1.0+hr1)+sb2r)*u0r*vzrL*u_xr  
		     - sbzr*sb_xr );
    }

    el = al*psi6*( (rho_blL*(1.0+hl1)+sb2l)*u0l*u_xl - sbtl*sb_xl );
    er = al*psi6*( (rho_brL*(1.0+hr1)+sb2r)*u0r*u_xr - sbtr*sb_xr );
  
    // HLL step for Sx:
    f_Sx_i[index] = (cminL*Fr + cmaxL*Fl - cminL*cmaxL*(er-el) )/(cmaxL + cminL);

    // Flux for S_y
    if (flux_direction==1) {
      Fl = al*psi6*( (rho_blL*(1.0+hl1)+sb2l)*u0l*vxlL*u_yl  
		     - sbxl*sb_yl );
      Fr = al*psi6*( (rho_brL*(1.0+hr1)+sb2r)*u0r*vxrL*u_yr  
		     - sbxr*sb_yr );
    } else if (flux_direction==2) {
      Fl = al*psi6*( (rho_blL*(1.0+hl1)+sb2l)*u0l*vylL*u_yl  
		     + PlL+0.5*sb2l - sbyl*sb_yl );
      Fr = al*psi6*( (rho_brL*(1.0+hr1)+sb2r)*u0r*vyrL*u_yr  
		     + PrL+0.5*sb2r - sbyr*sb_yr );
    } else if (flux_direction==3){
      Fl = al*psi6*( (rho_blL*(1.0+hl1)+sb2l)*u0l*vzlL*u_yl  
		     - sbzl*sb_yl );
      Fr = al*psi6*( (rho_brL*(1.0+hr1)+sb2r)*u0r*vzrL*u_yr  
		     - sbzr*sb_yr );
    }
    el = al*psi6*( (rho_blL*(1.0+hl1)+sb2l)*u0l*u_yl - sbtl*sb_yl );
    er = al*psi6*( (rho_brL*(1.0+hr1)+sb2r)*u0r*u_yr - sbtr*sb_yr );

    // HLL step for Sy:
    f_Sy_i[index] = (cminL*Fr + cmaxL*Fl - cminL*cmaxL*(er-el) )/(cmaxL + cminL);

    // Flux for S_z
    if (flux_direction==1) {
      Fl = al*psi6*( (rho_blL*(1.0+hl1)+sb2l)*u0l*vxlL*u_zl  
		     - sbxl*sb_zl );
      Fr = al*psi6*( (rho_brL*(1.0+hr1)+sb2r)*u0r*vxrL*u_zr  
		     - sbxr*sb_zr );
    } else if (flux_direction==2) {
      Fl = al*psi6*( (rho_blL*(1.0+hl1)+sb2l)*u0l*vylL*u_zl  
		     - sbyl*sb_zl );
      Fr = al*psi6*( (rho_brL*(1.0+hr1)+sb2r)*u0r*vyrL*u_zr  
		     - sbyr*sb_zr );
    } else if (flux_direction==3){
      Fl = al*psi6*( (rho_blL*(1.0+hl1)+sb2l)*u0l*vzlL*u_zl  
		     + PlL+0.5*sb2l - sbzl*sb_zl );
      Fr = al*psi6*( (rho_brL*(1.0+hr1)+sb2r)*u0r*vzrL*u_zr  
		     + PrL+0.5*sb2r - sbzr*sb_zr );
    }
    el = al*psi6*( (rho_blL*(1.0+hl1)+sb2l)*u0l*u_zl - sbtl*sb_zl );
    //if(i==4 && j==1 && k==1) printf("el: %.15e\t%.15e\t%.15e\n%.15e\t%.15e\n",rho_blL,hl1,sb2l,u0l,u_zl);
    er = al*psi6*( (rho_brL*(1.0+hr1)+sb2r)*u0r*u_zr - sbtr*sb_zr );
                                                                               
    // HLL step for Sz:
    f_Sz_i[index] = (cminL*Fr + cmaxL*Fl - cminL*cmaxL*(er-el) )/(cmaxL + cminL);
    //if(i==4 && j==1 && k==1) printf("FLUXES: %.15e\t%.15e\t%.15e\n%.15e\t%.15e\n",cminL,Fr,cmaxL,er,el);

    // Gotta be careful about x-direction in axisymmetry:
    if (Symmetry==AXISYM  &&  flux_direction==1) {
      double X_fL = X[index] - 0.5*dX;
      f_Sx_i[index] = f_Sx_i[index] * X_fL;
      f_Sy_i[index] = f_Sy_i[index] * SQR(X_fL);
      f_Sz_i[index] = f_Sz_i[index] * X_fL;
    }
    /*
      if(i==9 && j==1 && k==9) printf("FLUXES: %.15e %.15e %.15e %.15e %.15e\n",f_rho_i[index],f_tau_i[index],f_Sx_i[index],f_Sy_i[index],f_Sz_i[index]);

      vxr[index] = vxrL;
      vyr[index] = vyrL;
      vzr[index] = vzrL;

      vxl[index] = vxlL;
      vyl[index] = vylL;
      vzl[index] = vzlL;
    */

    cmax[index] = cmaxL;
    cmin[index] = cminL;
      }
}

inline void find_cp_cm(double &cplus,double &cminus,double &v02,double &u0,
		       double &vi,double &lapse,double &shifti,double &psim4,double &gupii) {


  //Find cplus, cminus:
  double a = SQR(u0) * (1.0-v02) + v02/SQR(1.0+lapse);
  double b = 2.0* ( shifti/SQR(1.0+lapse) * v02 - SQR(u0) * vi * (1.0-v02) );
  double c = SQR(u0*vi) * (1.0-v02) - v02 * ( psim4*gupii -
					      SQR(shifti/(1.0+lapse)) );
  double detm = b*b - 4.0*a*c;
  if(detm < 0.0) detm = 0.0;
  detm = sqrt(detm);

  //cplus = (detm-b)/2.0/a;
  //cminus = -(detm+b)/2.0/a;
  cplus = 0.5*(detm-b)/a;
  cminus = -0.5*(detm+b)/a;
  if (cplus < cminus) {
    double cp = cminus;
    cminus = cplus;
    cplus = cp;
  }
  /* 
  if (isnan(cplus)||isnan(cminus)){
    printf ("inside find_cp_cm, NaN found!!! v02=%e, u0=%e, detm = %e, a = %e, b=%e, c=%e \n", v02, u0, detm, a, b, c);
    }*/
}

double fasterpow_mhdflux(double inputvar,double inputpow) {
  if(inputpow==2.0) return SQR(inputvar);
  return pow(inputvar,inputpow);
}



extern "C" void CCTK_FCALL mhdfluxnew_hybrid_
  (int *flux_direction, const cGH **cctkGH,int *cctk_lsh,int *nghostzones,
   int *enable_HARM_energyvariable,
   double *X,double *dX,double *dY,double *dZ,
   double *f_rho_i,  double *f_rhoYe_i, double *f_tau_i, 
   double *f_Sx_i, double *f_Sy_i, double *f_Sz_i, 
   double *Pr, double *Pl, 
   double *rho_br, double *rho_bl, double *Y_er, double *Y_el,
   double *T_fluidl, double *T_fluidr,
   double *vxr, double *vxl, double *vyr, double *vyl, double *vzr, double *vzl, 
   double *Bxr, double *Bxl, double *Byr, double *Byl, double *Bzr, double *Bzl, 
   double *v02r, double *v02l, 
   double *gupxx_f, double *gupyy_f, double *gupzz_f, 
   double *cmax, double *cmin, 
   double *alpha_f, 
   double *betax_f, double *betay_f, double *betaz_f, 
   double *gxx_f, double *gxy_f, double *gxz_f, double *gyy_f, double *gyz_f, double *gzz_f, 
   double *phi_f, int *Symmetry, 
   int *enable_OS_collapse,double *rho_b_atm,
   int *neos, int *ergo_star, double *ergo_sigma, double *rho_tab, double *P_tab, double *eps_tab, double *k_tab, double *gamma_tab,double *gamma_th, 
   int *use_central_scheme_instead_of_hll,int *compute_microphysics) 
{
  mhdfluxnew_hybrid(*flux_direction, *cctkGH,cctk_lsh,nghostzones,
		    *enable_HARM_energyvariable,
		    X,*dX,*dY,*dZ,
		    f_rho_i, f_rhoYe_i, f_tau_i, 
		    f_Sx_i, f_Sy_i, f_Sz_i, 
		    Pr, Pl, 
		    rho_br, rho_bl, Y_er, Y_el,
		    T_fluidl, T_fluidr,
		    vxr, vxl, vyr, vyl, vzr, vzl, 
		    Bxr, Bxl, Byr, Byl, Bzr, Bzl, 
		    v02r, v02l, 
		    gupxx_f, gupyy_f, gupzz_f, 
		    cmax, cmin, 
		    alpha_f, 
		    betax_f, betay_f, betaz_f, 
		    gxx_f, gxy_f, gxz_f, gyy_f, gyz_f, gzz_f, 
		    phi_f, *Symmetry, 
		    *enable_OS_collapse, *rho_b_atm,
		    *neos, *ergo_star, *ergo_sigma,rho_tab, P_tab, eps_tab, k_tab, gamma_tab,*gamma_th, 
		    *use_central_scheme_instead_of_hll,*compute_microphysics);
}
