 //---------------------------
// compute radiation flux
//---------------------------
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "cctk.h"

#define SQR(x) ((x) * (x))
#define FAC 0.99
#define AXISYM 4

extern "C" void CCTK_FCALL CCTK_FNAME(flux_rad_cpp)
  (int *flux_direction, const cGH **cctkGH,int *ext,double *X,
   double *tau_rad_flux, 
   double *S_radx_flux, double *S_rady_flux, double *S_radz_flux,
   double *E_radr, double *E_radl, 
   double *F_rad0r, double *F_rad0l, 
   double *F_radxr, double *F_radxl,
   double *F_radyr, double *F_radyl,
   double *F_radzr, double *F_radzl,
   double *P_radr, double *P_radl, 
   double *vxr, double *vxl, double *vyr, double *vyl, double *vzr, double *vzl,
   double *gxx_f, double *gxy_f, double *gxz_f, double *gyy_f, double *gyz_f, double *gzz_f,
   double *cmax,double *cmin,
   double *betax_f, double *betay_f, double *betaz_f,
   double *alpha_f, double *phi_f, 
   int &pow_axi,int &Symmetry);

void flux_hll_cpp_frad(int *ext,double &Ur,double &Ul,double &Fr,double &Fl,double &F,double &cmax,double &cmin);

extern "C" void flux_rad_cpp(int flux_direction, const cGH *cctkGH,int *ext,double *X,
			     double *tau_rad_flux,
			     double *S_radx_flux, double *S_rady_flux, double *S_radz_flux,
			     double *E_radr, double *E_radl,
			     double *F_rad0r, double *F_rad0l,
			     double *F_radxr, double *F_radxl,
			     double *F_radyr, double *F_radyl,
			     double *F_radzr, double *F_radzl,
			     double *P_radr, double *P_radl, 
			     double *vxr, double *vxl, double *vyr, double *vyl, double *vzr, double *vzl,
			     double *gxx_f, double *gxy_f, double *gxz_f, double *gyy_f, double *gyz_f, double *gzz_f,
			     double *cmax,double *cmin,
			     double *betax_f, double *betay_f, double *betaz_f,
			     double *alpha_f, double *phi_f,
			     int &pow_axi,int &Symmetry) {

// Auxiliary arrays
// real*8, dimension(ext(1),ext(2),ext(3))        :: Fr,Fl, X_f
//

  int istart = 0;
  int jstart = 0;
  int kstart = 0;
  int iend = ext[0];
  int jend = ext[1];
  int kend = ext[2];


  // Find flux for tau_rad (alpha^2*psi6*R^{0i})

#pragma omp parallel for
    for(int k=kstart;k<kend;k++)
      for(int j=jstart;j<jend;j++)
	for(int i=istart;i<iend;i++) {
	  int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	  //printf("hi %d %d %d %e %e %e %e\n",i,j,k,vir[index],Bjr[index],vjr[index],Bir[index]);

	  double alpha_fL = alpha_f[index];
	  double al = 1.0 + alpha_fL;
	  double phi_fL = phi_f[index];
	  double psi2 = exp(2.0 * phi_fL);
	  double psi4 = psi2*psi2; 
	  double psi6 = psi4*psi2;
	  
	  double gxx_fL = gxx_f[index];
	  double gxy_fL = gxy_f[index];
	  double gxz_fL = gxz_f[index];
	  double gyy_fL = gyy_f[index];
	  double gyz_fL = gyz_f[index];
	  double gzz_fL = gzz_f[index];


	  double betax_fL = betax_f[index];
	  double betay_fL = betay_f[index];
	  double betaz_fL = betaz_f[index];

	  double vxrL = vxr[index];
	  double vyrL = vyr[index];
	  double vzrL = vzr[index];

	  double vxlL = vxl[index];
	  double vylL = vyl[index];  
	  double vzlL = vzl[index];  


	  double E_radrL = E_radr[index];
	  double E_radlL = E_radl[index];
	  double F_rad0rL = F_rad0r[index];
          double F_rad0lL = F_rad0l[index];
	  double F_radxrL = F_radxr[index];
          double F_radxlL = F_radxl[index];
	  double F_radyrL = F_radyr[index];
          double F_radylL = F_radyl[index];
          double F_radzrL = F_radzr[index];
          double F_radzlL = F_radzl[index];
	  double P_radrL = P_radr[index];
          double P_radlL = P_radl[index];


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

	  double el = sqrt(1.0-er);
	  double au0r1 = er/el/(1.0+el);
	  double u0r = (au0r1+1.0)/al;

	  double uxr = u0r*vxrL; 
	  double uyr = u0r*vyrL;
	  double uzr = u0r*vzrL;

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

	  el = sqrt(1.0-er);
	  double au0l1 = er/el/(1.0+el);
	  double u0l = (au0l1+1.0)/al;

	  double uxl = u0l*vxlL;
          double uyl = u0l*vylL;
          double uzl = u0l*vzlL;


	  double Fr, Fl;
	  if (flux_direction==1) {
	    Fr = al*al*psi6*((E_radrL + P_radrL)*u0r*uxr + F_rad0rL*uxr + u0r*F_radxrL);
	    Fl = al*al*psi6*((E_radlL + P_radlL)*u0l*uxl + F_rad0lL*uxl + u0l*F_radxlL);
	  } else  if (flux_direction==2) {
	    Fr = al*al*psi6*((E_radrL + P_radrL)*u0r*uyr + F_rad0rL*uyr + u0r*F_radyrL);
	    Fl = al*al*psi6*((E_radlL + P_radlL)*u0l*uyl + F_rad0lL*uyl + u0l*F_radylL);
	  } else if  (flux_direction==3) {
	    Fr = al*al*psi6*((E_radrL + P_radrL)*u0r*uzr + F_rad0rL*uzr + u0r*F_radzrL);
	    Fl = al*al*psi6*((E_radlL + P_radlL)*u0l*uzl + F_rad0lL*uzl + u0l*F_radzlL);
	  }

	  double tau_radr =  al*al*psi6*((E_radrL + P_radrL)*u0r*u0r + 2*F_rad0rL*u0r);
	  double tau_radl =  al*al*psi6*((E_radlL + P_radlL)*u0l*u0l + 2*F_rad0lL*u0l);

	  /*
	  if(isnan(Fr) || isnan(Fl)) {
	    printf("BadFr: %d %d %d %e %e %e %e\n",i,j,k,vir[index],Bjr[index],vjr[index],Bir[index]);
	    exit(1);
	  }
	  */

	  flux_hll_cpp_frad(ext,tau_radr,tau_radl,Fr,Fl,tau_rad_flux[index],cmax[index],cmin[index]);


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

	  double F_rad_xl = psi4*( gxx_fL*(F_radxlL+F_rad0lL*betax_fL) + gxy_fL*(F_radylL+F_rad0lL*betay_fL) +
				  gxz_fL*(F_radzlL+F_rad0lL*betaz_fL) );
	  double F_rad_yl = psi4*( gxy_fL*(F_radxlL+F_rad0lL*betax_fL) + gyy_fL*(F_radylL+F_rad0lL*betay_fL) +
				  gyz_fL*(F_radzlL+F_rad0lL*betaz_fL) );
	  double F_rad_zl = psi4*( gxz_fL*(F_radxlL+F_rad0lL*betax_fL) + gyz_fL*(F_radylL+F_rad0lL*betay_fL) +
				  gzz_fL*(F_radzlL+F_rad0lL*betaz_fL) );	
	 
	  double F_rad_xr = psi4*( gxx_fL*(F_radxrL+F_rad0rL*betax_fL) + gxy_fL*(F_radyrL+F_rad0rL*betay_fL) +
				  gxz_fL*(F_radzrL+F_rad0rL*betaz_fL) );
          double F_rad_yr = psi4*( gxy_fL*(F_radxrL+F_rad0rL*betax_fL) + gyy_fL*(F_radyrL+F_rad0rL*betay_fL) +
				  gyz_fL*(F_radzrL+F_rad0rL*betaz_fL) );
	  double F_rad_zr = psi4*( gxz_fL*(F_radxrL+F_rad0rL*betax_fL) + gyz_fL*(F_radyrL+F_rad0rL*betay_fL) +
                                  gzz_fL*(F_radzrL+F_rad0rL*betaz_fL) );
	
	  // Find flux for S_rad_i (alpha*psi6*R^j_i                                                                                                                                                                 //= alpha*psi6*(4/3*E_rad*u^j*u_i + F_rad^j*u_i + F_rad_i*u^j)  ) 
	  double Fxr, Fxl;
	  double Fyr, Fyl;
	  double Fzr, Fzl;

	  if (flux_direction==1) {
            Fxr = al*psi6*((E_radrL + P_radrL)*u_xr*uxr + F_radxrL*u_xr + uxr*F_rad_xr);
            Fxl = al*psi6*((E_radlL + P_radlL)*u_xl*uxl + F_radxlL*u_xl + uxl*F_rad_xl);
	    Fyr = al*psi6*((E_radrL + P_radrL)*u_yr*uxr + F_radxrL*u_yr + uxr*F_rad_yr);
            Fyl = al*psi6*((E_radlL + P_radlL)*u_yl*uxl + F_radxlL*u_yl + uxl*F_rad_yl);
	    Fzr = al*psi6*((E_radrL + P_radrL)*u_zr*uxr + F_radxrL*u_zr + uxr*F_rad_zr);
            Fzl = al*psi6*((E_radlL + P_radlL)*u_zl*uxl + F_radxlL*u_zl + uxl*F_rad_zl);          
	  } else  if (flux_direction==2) {
            Fxr = al*psi6*((E_radrL + P_radrL)*u_xr*uyr + F_radyrL*u_xr + uyr*F_rad_xr);
            Fxl = al*psi6*((E_radlL + P_radlL)*u_xl*uyl + F_radylL*u_xl + uyl*F_rad_xl);
            Fyr = al*psi6*((E_radrL + P_radrL)*u_yr*uyr + F_radyrL*u_yr + uyr*F_rad_yr);
            Fyl = al*psi6*((E_radlL + P_radlL)*u_yl*uyl + F_radylL*u_yl + uyl*F_rad_yl);
            Fzr = al*psi6*((E_radrL + P_radrL)*u_zr*uyr + F_radyrL*u_zr + uyr*F_rad_zr);
            Fzl = al*psi6*((E_radlL + P_radlL)*u_zl*uyl + F_radylL*u_zl + uyl*F_rad_zl);
	  } else if  (flux_direction==3) {
	    Fxr = al*psi6*((E_radrL + P_radrL)*u_xr*uzr + F_radzrL*u_xr + uzr*F_rad_xr);
            Fxl = al*psi6*((E_radlL + P_radlL)*u_xl*uzl + F_radzlL*u_xl + uzl*F_rad_xl);
            Fyr = al*psi6*((E_radrL + P_radrL)*u_yr*uzr + F_radzrL*u_yr + uzr*F_rad_yr);
            Fyl = al*psi6*((E_radlL + P_radlL)*u_yl*uzl + F_radzlL*u_yl + uzl*F_rad_yl);
            Fzr = al*psi6*((E_radrL + P_radrL)*u_zr*uzr + F_radzrL*u_zr + uzr*F_rad_zr);
            Fzl = al*psi6*((E_radlL + P_radlL)*u_zl*uzl + F_radzlL*u_zl + uzl*F_rad_zl);
	  } 

          double S_radxr =  al*psi6*((E_radrL + P_radrL)*u0r*u_xr + F_rad0rL*u_xr + u0r*F_rad_xr);
          double S_radxl =  al*psi6*((E_radlL + P_radlL)*u0l*u_xl + F_rad0lL*u_xl + u0l*F_rad_xl);
          flux_hll_cpp_frad(ext,S_radxr,S_radxl,Fxr,Fxl,S_radx_flux[index],cmax[index],cmin[index]);

	  double S_radyr =  al*psi6*((E_radrL + P_radrL)*u0r*u_yr + F_rad0rL*u_yr + u0r*F_rad_yr);
          double S_radyl =  al*psi6*((E_radlL + P_radlL)*u0l*u_yl + F_rad0lL*u_yl + u0l*F_rad_yl);
          flux_hll_cpp_frad(ext,S_radyr,S_radyl,Fyr,Fyl,S_rady_flux[index],cmax[index],cmin[index]);

	  double S_radzr =  al*psi6*((E_radrL + P_radrL)*u0r*u_zr + F_rad0rL*u_zr + u0r*F_rad_zr);
          double S_radzl =  al*psi6*((E_radlL + P_radlL)*u0l*u_zl + F_rad0lL*u_zl + u0l*F_rad_zl);
          flux_hll_cpp_frad(ext,S_radzr,S_radzl,Fzr,Fzl,S_radz_flux[index],cmax[index],cmin[index]);

	}
 
    if (Symmetry==AXISYM && flux_direction==1) {
      double dX2 = 0.5*(X[CCTK_GFINDEX3D(cctkGH,1,0,0)]-X[CCTK_GFINDEX3D(cctkGH,0,0,0)]);

#pragma omp parallel for
      for(int k=kstart;k<kend;k++)
	for(int j=jstart;j<jend;j++)
	  for(int i=istart;i<iend;i++){
	    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	    double X_fL = X[index] - 0.5*dX2;
      
	    tau_rad_flux[index] *= X_fL;
	    S_radx_flux[index] *= X_fL;
	    S_rady_flux[index] *= SQR(X_fL);
	    S_radz_flux[index] *= X_fL;
	  }
    }

}




extern "C" void CCTK_FCALL CCTK_FNAME(flux_rad_cpp)
  (int *flux_direction, const cGH **cctkGH,int *ext,double *X,
   double *tau_rad_flux,
   double *S_radx_flux, double *S_rady_flux, double *S_radz_flux,
   double *E_radr, double *E_radl,
   double *F_rad0r, double *F_rad0l,
   double *F_radxr, double *F_radxl,
   double *F_radyr, double *F_radyl,
   double *F_radzr, double *F_radzl,
   double *P_radr, double *P_radl,
   double *vxr, double *vxl, double *vyr, double *vyl, double *vzr, double *vzl,
   double *gxx_f, double *gxy_f, double *gxz_f, double *gyy_f, double *gyz_f, double *gzz_f,
   double *cmax,double *cmin,
   double *betax_f, double *betay_f, double *betaz_f,
   double *alpha_f, double *phi_f,
   int &pow_axi,int &Symmetry)

{
  flux_rad_cpp
    (*flux_direction, *cctkGH,ext,X,
     tau_rad_flux,
     S_radx_flux, S_rady_flux, S_radz_flux,
     E_radr, E_radl, 
     F_rad0r, F_rad0l, 
     F_radxr,F_radxl, 
     F_radyr, F_radyl,
     F_radzr, F_radzl,
     P_radr,P_radl,
     vxr,vxl,vyr,vyl,vzr,vzl,
     gxx_f, gxy_f, gxz_f, gyy_f, gyz_f, gzz_f,
     cmax,cmin,
     betax_f, betay_f, betaz_f, 
     alpha_f, phi_f, 
     pow_axi,Symmetry);
}

void flux_hll_cpp_frad(int *ext,double &Ur,double &Ul,double &Fr,double &Fl,double &F,double &cmax,double &cmin) {
  F = (cmin*Fr + cmax*Fl - cmin*cmax*(Ur-Ul) )/(cmax + cmin);
  //We find that cmax==cmin==0 when we are in the ghostzone sometimes when B fields are being evolved.
  //   In the above line of code, this means that F=NaN, which leads to instability.
  //if(cmax+cmin==0) F=0;
}
