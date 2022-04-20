//The procedure of computing cmax/cmin is the same for radiation variables and for others
//(we assume the dispersion relation is independent of radiation), so no need to change.
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "cctk.h"

#define SQR(x) ((x) * (x))
#define FACTOR 0.99

extern "C" void CCTK_FCALL CCTK_FNAME(compute_cmax_cmin_hybrid_rad_cpp)
  (const cGH **cctkGH,int *ext,double *v02r,double *v02l,double *cmax,double *cmin,double *rho_br,double *rho_bl,
   double *Pr,double *Pl,double *vxr,double *vxl,double *vyr,double *vyl,double *vzr,double *vzl,
   double *Bxr,double *Bxl,double *Byr,double *Byl,double *Bzr,double *Bzl,
   double *alpha,double *betax,double *betay,double *betaz,double *phi,
   double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
   double *gupxx,double *gupyy,double *gupzz,
   int &m,int &neos, int &ergo_star, double &ergo_sigma, double *rho_tab,double *P_tab, double *eps_tab, double *k_tab,
   double *gamma_tab,double &gamma_th,
   double *cplusr,double *cminusr,double *u0,double *cplusl,double *cminusl);

double fasterpow_cmaxcmin(double inputvar,double inputpow);

void find_cplus_cminus_cpp(const cGH *cctkGH,int *ext,double *cplus,double *cminus,double *v02,double *u0, 
			   double *vi,double *lapse,double *shifti,double *phi,double *gupii);

void find_v02_hybrid_cpp(const cGH *cctkGH,int *ext,double *v02,double *rho_b,double *P,
			 double *u0,double *vx,double *vy,double *vz, 
			 int &neos, int &ergo_star, double &ergo_sigma, 
			 double *rho_tab,double *P_tab,double *eps_tab,double *k_tab,double *gamma_tab,double &gamma_th, 
			 double *sqBx,double *sqBy,double *sqBz,
			 double *lapse,double *shiftx,double *shifty,double *shiftz,double * phi,
			 double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz);

//--------------------------------------------------!
// cmin and cmax				    !
//--------------------------------------------------!
void compute_cmax_cmin_hybrid_rad_cpp(const cGH *cctkGH,int *ext,double *v02r,double *v02l,double *cmax,double *cmin,double *rho_br,double *rho_bl,
				  double *Pr,double *Pl,double *vxr,double *vxl,double *vyr,double *vyl,double *vzr,double *vzl,
				  double *Bxr,double *Bxl,double *Byr,double *Byl,double *Bzr,double *Bzl,
				  double *alpha,double *betax,double *betay,double *betaz,double *phi,
				  double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
				  double *gupxx,double *gupyy,double *gupzz,
				  int &m,int &neos, int &ergo_star, double &ergo_sigma, 
				  double *rho_tab,double *P_tab, double *eps_tab, double *k_tab,
				  double *gamma_tab,double &gamma_th,
				  double *cplusr,double *cminusr,double *u0,double *cplusl,double *cminusl) {

  printf("hi inside compute_cmax_cmin_hybrid_RAD_cpp m=%d\n",m);

  find_v02_hybrid_cpp(cctkGH,ext,v02r,rho_br,Pr,u0,vxr,vyr,vzr,
		      neos, ergo_star, ergo_sigma, rho_tab, P_tab, eps_tab, k_tab, gamma_tab, gamma_th,
		      Bxr,Byr,Bzr, alpha,betax,betay,betaz, phi,gxx,gxy,gxz,gyy,gyz,gzz);
    
  if (m==1) {
    find_cplus_cminus_cpp(cctkGH,ext,cplusr,cminusr,v02r,u0,
			  vxr,alpha,betax,phi,gupxx);
  } else if (m==2) {
    find_cplus_cminus_cpp(cctkGH,ext,cplusr,cminusr,v02r,u0,
			  vyr,alpha,betay,phi,gupyy);
  } else {
    find_cplus_cminus_cpp(cctkGH,ext,cplusr,cminusr,v02r,u0,
			  vzr,alpha,betaz,phi,gupzz);
  }
  
  find_v02_hybrid_cpp(cctkGH,ext,v02l,rho_bl,Pl,u0,vxl,vyl,vzl,
		      neos, ergo_star, ergo_sigma, rho_tab, P_tab, eps_tab, k_tab, gamma_tab, gamma_th,
		      Bxl,Byl,Bzl, alpha,betax,betay,betaz, phi,gxx,gxy,gxz,gyy,gyz,gzz);
  
  if (m==1) {
    find_cplus_cminus_cpp(cctkGH,ext,cplusl,cminusl,v02l,u0,
			  vxl,alpha,betax,phi,gupxx);
  } else if (m==2) {
    find_cplus_cminus_cpp(cctkGH,ext,cplusl,cminusl,v02l,u0,
			  vyl,alpha,betay,phi,gupyy);
  } else {
    find_cplus_cminus_cpp(cctkGH,ext,cplusl,cminusl,v02l,u0,
			  vzl,alpha,betaz,phi,gupzz);
  }
 
  /* Set up variables used in the grid loop for the physical grid points */
  int istart = 0;
  int jstart = 0;
  int kstart = 0;
  int iend = ext[0];
  int jend = ext[1];
  int kend = ext[2];

#pragma omp parallel for
  for(int k=kstart;k<kend;k++)
    for(int j=jstart;j<jend;j++)
      for(int i=istart;i<iend;i++) {
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	cmax[index] = 0.0;
	if(cplusl[index]>0.0)         cmax[index] = cplusl[index];
	if(cplusr[index]>cmax[index]) cmax[index] = cplusr[index];

	cmin[index] = 0.0;
	if(cminusl[index]<0.0)         cmin[index] = cminusl[index];
	if(cminusr[index]<cmin[index]) cmin[index] = cminusr[index];
	cmin[index]*=-1.0;

       }
}

//----------------------------------------------------//
// Find maximum and minimum velocities along x^i-axis //
//----------------------------------------------------//
//
void find_cplus_cminus_cpp(const cGH *cctkGH,int *ext,double *cplus,double *cminus,double *v02,double *u0, 
			   double *vi,double *lapse,double *shifti,double *phi,double *gupii) {

  /* Set up variables used in the grid loop for the physical grid points */
  int istart = 0;
  int jstart = 0;
  int kstart = 0;
  int iend = ext[0];
  int jend = ext[1];
  int kend = ext[2];

#pragma omp parallel for
  for(int k=kstart;k<kend;k++)
    for(int j=jstart;j<jend;j++)
      for(int i=istart;i<iend;i++) {
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
	
	double u0L = u0[index];
	double v02L = v02[index];
	double lapseL = lapse[index];
	double viL = vi[index];
	double phiL = phi[index];
	double shiftiL = shifti[index];
	double gupiiL = gupii[index];

	double a = SQR(u0L) * (1.0-v02L) + v02L/SQR(1.0+lapseL);
	double b = 2.0* ( shiftiL/SQR(1.0+lapseL) * v02L - SQR(u0L) * viL * (1.0-v02L) );
	double c = SQR(u0L*viL) * (1.0-v02L) - v02L * ( exp(-4.0*phiL)*gupiiL - 
							SQR(shiftiL/(1.0+lapseL)) );

	double detm=0;
	double determinant = b*b - 4.0*a*c;
	if(determinant > 0) detm = sqrt(determinant);
	//double detm = sqrt( max(b*b - 4.0*a*c,0.0) );
	cplus[index] =   0.5*(detm-b)/a;
	cminus[index] = -0.5*(detm+b)/a;
	
	if (cplus[index] < cminus[index]) {
	  double cp = cminus[index];
	  cminus[index] = cplus[index];
	  cplus[index] = cp;
	}
      }
}

//-----------------------------------------------------------------//
// Calculate an estimate of maximum velocity in fluid's rest frame //
// for hybrid EOS						  //
//-----------------------------------------------------------------//
//
void find_v02_hybrid_cpp(const cGH *cctkGH,int *ext,double *v02,double *rho_b,double *P,
			 double *u0,double *vx,double *vy,double *vz, 
			 int &neos, int &ergo_star, double &ergo_sigma, 
			 double *rho_tab,double *P_tab,double *eps_tab,double *k_tab,double *gamma_tab,double &gamma_th, 
			 double *sqBx,double *sqBy,double *sqBz,
			 double *lapse,double *shiftx,double *shifty,double *shiftz,double * phi,
			 double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz) {
  /* Set up variables used in the grid loop for the physical grid points */
  int istart = 0;
  int jstart = 0;
  int kstart = 0;
  int iend = ext[0];
  int jend = ext[1];
  int kend = ext[2];

#pragma omp parallel for
  for(int k=kstart;k<kend;k++)
    for(int j=jstart;j<jend;j++)
      for(int i=istart;i<iend;i++) {
	int index = CCTK_GFINDEX3D(cctkGH,i,j,k);

	double rho_bL = rho_b[index];
	double PL = P[index];

	double vxL = vx[index];
	double vyL = vy[index];
	double vzL = vz[index];
	
	double sqBxL = sqBx[index];
	double sqByL = sqBy[index];
	double sqBzL = sqBz[index];

	double lapseL = lapse[index];

	double shiftxL = shiftx[index];
	double shiftyL = shifty[index];
	double shiftzL = shiftz[index];

	double Psi4 = exp(4.0*phi[index]);
	double Psim12 = 1.0/(Psi4*Psi4*Psi4);

	double gxxL = gxx[index];
	double gxyL = gxy[index];
	double gxzL = gxz[index];
	double gyyL = gyy[index];
	double gyzL = gyz[index];
	double gzzL = gzz[index];

	double h,dPcold_drho,eps_th,va2,cs2,b2;

	if(rho_bL <= rho_tab[0]) {
	  if (rho_bL <= 0.0) {
	    v02[index] = 0.0;                                   // v02[index] -> P_cold
	    dPcold_drho = 0.0;
	    eps_th = 0.0;
	    h = 1.0;
	  } else {
	    v02[index] = k_tab[0]*fasterpow_cmaxcmin(rho_bL,gamma_tab[0]);           // v02[index] -> P_cold
	    va2 = v02[index]/rho_bL/(gamma_tab[0]-1.0);          // va2 -> eps_cold
	    dPcold_drho = gamma_tab[0]*v02[index]/rho_bL;
	    eps_th = (PL-v02[index])/(gamma_th-1.0)/rho_bL;
	    h = 1.0 + eps_th + va2 + PL/rho_bL;
	  }
	}
 
      for(int ii=1;ii<neos;ii++) {
	if (rho_bL <= rho_tab[ii] && rho_bL > rho_tab[ii-1]) {
	  v02[index] = k_tab[ii]*fasterpow_cmaxcmin(rho_bL,gamma_tab[ii]);	       // v02[index] -> P_cold
	  va2 = eps_tab[ii-1] + (v02[index]/rho_bL - P_tab[ii-1]/rho_tab[ii-1])/ 
	     (gamma_tab[ii]-1.0);                   // va2 -> eps_cold
	  dPcold_drho = gamma_tab[ii]*v02[index]/rho_bL;
	  eps_th = (PL-v02[index])/(gamma_th-1.0)/rho_bL;
	  h = 1.0 + eps_th + va2 + PL/rho_bL;
          }
        }

 
	if(rho_bL > rho_tab[neos-1]) {
	  if (ergo_star == 0){
	  v02[index] = k_tab[neos]*fasterpow_cmaxcmin(rho_bL,gamma_tab[neos]);   // v02[index] -> P_cold
	  va2 = eps_tab[neos-1] + (v02[index]/rho_bL - P_tab[neos-1]/rho_tab[neos-1])/ 
	    (gamma_tab[neos]-1.0);               // va2 -> eps_cold
	  dPcold_drho = gamma_tab[neos]*v02[index]/rho_bL;
	  eps_th = (PL-v02[index])/(gamma_th-1.0)/rho_bL;
	  }
	  else {
	  v02[index] = (ergo_sigma* (1+eps_tab[neos-1]+P_tab[neos-1]/rho_tab[neos-1])/fasterpow_cmaxcmin(rho_tab[neos-1],ergo_sigma)* fasterpow_cmaxcmin(rho_bL, ergo_sigma+1) + P_tab[neos-1] - ergo_sigma*((1+eps_tab[neos-1])*rho_tab[neos-1]) )/(ergo_sigma+1);
          va2 = ((1+eps_tab[neos-1]+P_tab[neos-1]/rho_tab[neos-1])/fasterpow_cmaxcmin(rho_tab[neos-1],ergo_sigma) * fasterpow_cmaxcmin(rho_bL, ergo_sigma+1) - P_tab[neos-1] + ergo_sigma*((1+eps_tab[neos-1])*rho_tab[neos-1]) )/((ergo_sigma+1)*rho_bL)-1;
	  dPcold_drho = gamma_tab[neos]*v02[index]/rho_bL;
	  eps_th = (PL-v02[index])/(gamma_th-1.0)/rho_bL;
	  }
	 h = 1.0 + eps_th + va2 + PL/rho_bL;
	}
      

	u0[index] = Psi4*(gxxL*SQR(vxL + shiftxL) +
			  2.0*gxyL*(vxL + shiftxL)*(vyL + shiftyL) +        
			  2.0*gxzL*(vxL + shiftxL)*(vzL + shiftzL) +        
			  gyyL*SQR(vyL + shiftyL) +                          
			  2.0*gyzL*(vyL + shiftyL)*(vzL + shiftzL) +        
			  gzzL*SQR(vzL + shiftzL) )/SQR(1.0+lapseL);
	double u0L = u0[index];
	// *** Check for superluminal velocity ***
	if(u0L > 1.0) {
	  vx[index] = (vxL + shiftxL)*sqrt(FACTOR/u0L)-shiftxL;
	  vy[index] = (vyL + shiftyL)*sqrt(FACTOR/u0L)-shiftyL;
	  vz[index] = (vzL + shiftzL)*sqrt(FACTOR/u0L)-shiftzL;
	  u0[index] = SQR(FACTOR);
	  u0L = u0[index];
	  vxL = vx[index];
	  vyL = vy[index];
	  vzL = vz[index];
	}
	u0[index] = 1.0/(1.0+lapseL)/sqrt(1.0-u0L);
	u0L = u0[index];
 
	double sbt = u0L*Psi4*( gxxL*(vxL+shiftxL)*sqBxL +
				gxyL*( (vxL+shiftxL)*sqByL + (vyL+shiftyL)*sqBxL ) +
				gxzL*( (vxL+shiftxL)*sqBzL + (vzL+shiftzL)*sqBxL ) +
				gyyL*(vyL+shiftyL)*sqByL +
				gyzL*( (vyL+shiftyL)*sqBzL + (vzL+shiftzL)*sqByL ) +
				gzzL*(vzL+shiftzL)*sqBzL );
	double sbx = sqBxL/u0L + sbt*vxL;
	double sby = sqByL/u0L + sbt*vyL;
	double sbz = sqBzL/u0L + sbt*vzL;

	b2 = - SQR((1.0+lapseL)*sbt) + Psi4*(gxxL*SQR(sbx + shiftxL*sbt) +
					     2.0*gxyL*(sbx + shiftxL*sbt)*(sby + shiftyL*sbt) +
					     2.0*gxzL*(sbx + shiftxL*sbt)*(sbz + shiftzL*sbt) +
					     gyyL*SQR(sby + shiftyL*sbt) +
					     2.0*gyzL*(sby + shiftyL*sbt)*(sbz + shiftzL*sbt) +
					     gzzL*SQR(sbz + shiftzL*sbt) );
	
	b2 = b2/(4.0*M_PI)*Psim12/SQR(SQR(1.0+lapseL));

	if(rho_bL > 0.0) {
	  cs2 = (dPcold_drho + gamma_th*(gamma_th-1.0)*eps_th)/h;
	  va2 = b2/(b2 + rho_bL*h);
	  v02[index] = va2+cs2*(1.0-va2);
	} else {
	  v02[index] = b2/(b2+1e-300);
	}
      }
}

double fasterpow_cmaxcmin(double inputvar,double inputpow) {
  if(inputpow==2.0) return SQR(inputvar);
  return pow(inputvar,inputpow);
}

extern "C" void CCTK_FCALL CCTK_FNAME(compute_cmax_cmin_hybrid_rad_cpp)
  (const cGH **cctkGH,int *ext,double *v02r,double *v02l,double *cmax,double *cmin,double *rho_br,double *rho_bl,
   double *Pr,double *Pl,double *vxr,double *vxl,double *vyr,double *vyl,double *vzr,double *vzl,
   double *Bxr,double *Bxl,double *Byr,double *Byl,double *Bzr,double *Bzl,
   double *alpha,double *betax,double *betay,double *betaz,double *phi,
   double *gxx,double *gxy,double *gxz,double *gyy,double *gyz,double *gzz,
   double *gupxx,double *gupyy,double *gupzz,
   int &m,int &neos, int &ergo_star, double &ergo_sigma, double *rho_tab,double *P_tab, double *eps_tab, double *k_tab,
   double *gamma_tab,double &gamma_th,
   double *cplusr,double *cminusr,double *u0,double *cplusl,double *cminusl)
{
  compute_cmax_cmin_hybrid_rad_cpp(*cctkGH,ext,v02r,v02l,cmax,cmin,rho_br,rho_bl,
			       Pr,Pl,vxr,vxl,vyr,vyl,vzr,vzl,
			       Bxr,Bxl,Byr,Byl,Bzr,Bzl,
			       alpha,betax,betay,betaz,phi,
			       gxx,gxy,gxz,gyy,gyz,gzz,
			       gupxx,gupyy,gupzz,
			       m,neos, ergo_star, ergo_sigma, rho_tab,P_tab, eps_tab, k_tab,
			       gamma_tab,gamma_th,
			       cplusr,cminusr,u0,cplusl,cminusl);
}

