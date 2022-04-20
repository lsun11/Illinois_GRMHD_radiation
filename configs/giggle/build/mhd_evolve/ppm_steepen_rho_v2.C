//-----------------------------------------------------------------//
// steepen rho_b
// Note from YT: This subroutine is a major modification to 
//  the DAGH's version, and I believe there are bugs in 
//  the DAGH's version.
//-----------------------------------------------------------------//
//
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "cctk.h"

#define AXISYM 4
#define ONESIXTH 0.16666666666666666666666666666666666666666666
#define SQR(x) ((x) * (x))

// standard Colella-Woodward parameters:
//    K0 = 0.1d0, eta1 = 20.0, eta2 = 0.05, epsilon = 0.01d0
// our parameters:
//parameter(K0 = 0.1d0, eta1 = 20.0, eta2 = 0.05, epsilon = 0.01d0)
#define K0      0.1
#define ETA1   20.0
#define ETA2    0.05
#define EPSILON 0.01

// The following values are Font-like:
// parameter(K0 = 1.0, eta1 = 5.0, eta2 = 0.05, epsilon = 0.1d0)
//

extern "C" void CCTK_FCALL ppm_steepen_rho_v2_cpp_
  (const cGH **cctkGH,int *ext,double *rho,double *d1rho,double *d2rho,double *rhor_ppm,double *rhol_ppm,
   double *rhor_mc,double *rhol_mc,double *P,
   double *rho_tab,double *P_tab,double *eps_tab,double *k_tab,double *gamma_tab,double &gamma_th,  
   int &neos,int &ergo_star, double &ergo_sigma,int &m,int &Symmetry, int &enable_OS_collapse);

double max_ppmsteepenrho(double val1, double val2);
double min_ppmsteepenrho(double val1, double val2);
double fasterpow_ppmsteepenrho(double inputvar,double inputpow);

void compute_pcold_epscold_steepenrho_cpp(double &rhob, double &P_cold, double &eps_cold, 
					  int &neos,int &ergo_star, double &ergo_sigma,double *rho_tab,double *P_tab,double *eps_tab,double *k_tab,double *gamma_tab, int &enable_OS_collapse);

extern "C" void ppm_steepen_rho_v2_cpp(const cGH *cctkGH,int *ext,double *rho,double *d1rho,double *d2rho,double *rhor_ppm,double *rhol_ppm,
				       double *rhor_mc,double *rhol_mc,double *P,
				       double *rho_tab,double *P_tab,double *eps_tab,double *k_tab,double *gamma_tab,double &gamma_th,  
				       int &neos,int &ergo_star, double &ergo_sigma,int &m,int &Symmetry, int &enable_OS_collapse) {
  // 
  int imin=0,jmin=0,kmin=0;
  int imax=ext[0];
  int jmax=ext[1];
  int kmax=ext[2];
  
  if (Symmetry==AXISYM) { 
    jmin = 1;
    jmax = 1;
  }
  
  //Fill all points except for upper boundaries in m-direction.
  //It should not be necessary to fill lower and upper boundaries, since the values will be overwritten at the next sync call.
  // Besides, the values at the boundaries aren't correct.
  if(m==1) { imin++; imax--; }
  if(m==2) { jmin++; jmax--; }
  if(m==3) { kmin++; kmax--; }
#pragma omp parallel for
  for(int k=kmin;k<kmax;k++) for(int j=jmin;j<jmax;j++) for(int i=imin;i<imax;i++) {
    int index = CCTK_GFINDEX3D(cctkGH,i,j,k);
    int indexm1,indexp1;
    if(m==1) {
      indexp1 = CCTK_GFINDEX3D(cctkGH,i+1,j,k);
      indexm1 = CCTK_GFINDEX3D(cctkGH,i-1,j,k);
    } else if(m==2) {
      indexp1 = CCTK_GFINDEX3D(cctkGH,i,j+1,k);
      indexm1 = CCTK_GFINDEX3D(cctkGH,i,j-1,k);
    } else {
      indexp1 = CCTK_GFINDEX3D(cctkGH,i,j,k+1);
      indexm1 = CCTK_GFINDEX3D(cctkGH,i,j,k-1);
    }

    double gamma_cold = gamma_tab[neos];
    int exit_do = 0;
    for(int n=0;n<neos;n++) {
      if (rho[index] < rho_tab[n]) { 
	gamma_cold = gamma_tab[n];
	exit_do = 1;
      }
      if (exit_do == 1) n=neos+100;
    }

    double P_cold,eps_cold;
    compute_pcold_epscold_steepenrho_cpp(rho[index], P_cold, eps_cold, 
					 neos,ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,k_tab,gamma_tab,enable_OS_collapse);
    // Compute effective Gamma = (partial P / partial rho0)_s /(P/rho0)
    double Gamma = gamma_th + (gamma_cold-gamma_th)*P_cold/P[index];
    double contact_disc_check = Gamma*K0*fabs(rho[indexp1]-rho[indexm1])*  
      min_ppmsteepenrho(P[indexp1],P[indexm1])  
      -fabs(P[indexp1]-P[indexm1])*min_ppmsteepenrho(rho[indexp1],rho[indexm1]);
    double second_dir_check = -d2rho[indexp1]*d2rho[indexm1];
    double relative_change_check = fabs(2.0*d1rho[index]) 
      - EPSILON*min_ppmsteepenrho(rho[indexp1],rho[indexm1]);
    if(contact_disc_check >= 0.0 && second_dir_check >= 0.0 
       && relative_change_check >= 0.0) {
      double eta_tilde=0.0;
      if (fabs(d1rho[index]) > 0.0) { 
	eta_tilde = -ONESIXTH*(d2rho[indexp1]-d2rho[indexm1])/(2.0*d1rho[index]);
      }
      double eta = max_ppmsteepenrho(0.0,min_ppmsteepenrho(ETA1*(eta_tilde - ETA2),1.0));
      // rhol_mc is defined as "Ul" here:
      //      Ul[index] = U[indexm1] + 0.5*nabla_U[indexm1]; <- nabla_U is the MC slope-limited \delta U
      rhol_ppm[index] = rhol_ppm[index]*(1.0-eta) 
	+ rhol_mc[index]*eta; 
      // rhor_mc is defined as "Ur" here:
      //      Ur[index] = U[index]   - 0.5*nabla_U[index]; <- nabla_U is the MC slope-limited \delta U
      // Based on this line of code, Ur[index] = a_j - \delta_m a_j / 2. (cf. Eq. 65 in Marti & Muller's "PPM Method for 1D Relativistic Hydro." paper)
      //    So: Ur[indexp1] = a_{j+1} - \delta_m a_{j+1} / 2. This is why we have rhor_mc[indexp1]
      rhor_ppm[index] = rhor_ppm[index]*(1.0-eta) 
	+ rhor_mc[indexp1]*eta;
    }
  }
}

//The following is just copied from primitives_hybrid.C, but with the function name changed.  
//   We do this so that the compiler will inline it into the above loops.  I.e., to make
//   the code go FASTER!
void compute_pcold_epscold_steepenrho_cpp(double &rhob, double &P_cold, double &eps_cold, 
					  int &neos, int &ergo_star, double &ergo_sigma, double *rho_tab,double *P_tab,double *eps_tab,double *k_tab,double *gamma_tab,int &enable_OS_collapse) {
  
  bool exit_do;
  int i = 0;

  if (enable_OS_collapse == 1){
    
    P_cold = k_tab[0]*fasterpow_ppmsteepenrho(rhob,1.0);
    eps_cold = P_cold/rhob/(gamma_tab[1]-1.0);
    /*
    P_cold = 0.0;
    eps_cold = 0.0;*/
  }

  else{
  exit_do = 0;
  while(exit_do==0) {
    if (rhob <= rho_tab[i]) {
      exit_do = 1;
      P_cold = k_tab[i]*fasterpow_ppmsteepenrho(rhob,gamma_tab[i]);
      if (i==0) {
        if (rhob != 0.0) {
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
    if (ergo_star==0){
      P_cold = k_tab[neos]*fasterpow_ppmsteepenrho(rhob,gamma_tab[neos]);
      eps_cold = eps_tab[neos-1] + (P_cold/rhob - P_tab[neos-1]/rho_tab[neos-1])/(gamma_tab[neos]-1.0);
    }
    else {
      P_cold = ((ergo_sigma* (1+eps_tab[neos-1]+P_tab[neos-1]/rho_tab[neos-1])/fasterpow_ppmsteepenrho(rho_tab[neos-1],ergo_sigma)) * fasterpow_ppmsteepenrho(rhob, ergo_sigma+1) + P_tab[neos-1] - ergo_sigma*\
		((1+eps_tab[neos-1])*rho_tab[neos-1]))/(ergo_sigma+1);
      eps_cold = ((1+eps_tab[neos-1]+P_tab[neos-1]/rho_tab[neos-1])/fasterpow_ppmsteepenrho(rho_tab[neos-1],ergo_sigma) * fasterpow_ppmsteepenrho(rhob, ergo_sigma+1) - P_tab[neos-1] + ergo_sigma*((1+eps_tab[\
																							     neos-1])*rho_tab[neos-1]) )/((ergo_sigma+1.0)*rhob)-1.0;
    }
  }
  }
}

double max_ppmsteepenrho(double val1, double val2) {
  if(val1>val2) return val1;
  return val2;
}
double min_ppmsteepenrho(double val1, double val2) {
  if(val1<val2) return val1;
  return val2;
}

double fasterpow_ppmsteepenrho(double inputvar,double inputpow) {
  if(inputpow==2.0) return SQR(inputvar);
  return pow(inputvar,inputpow);
}

extern "C" void CCTK_FCALL ppm_steepen_rho_v2_cpp_
  (const cGH **cctkGH,int *ext,double *rho,double *d1rho,double *d2rho,double *rhor_ppm,double *rhol_ppm,
   double *rhor_mc,double *rhol_mc,double *P,
   double *rho_tab,double *P_tab,double *eps_tab,double *k_tab,double *gamma_tab,double &gamma_th,  
   int &neos, int &ergo_star, double &ergo_sigma, int &m,int &Symmetry, int &enable_OS_collapse)
{
  ppm_steepen_rho_v2_cpp(*cctkGH,ext,rho,d1rho,d2rho,rhor_ppm,rhol_ppm,
			 rhor_mc,rhol_mc,P,
			 rho_tab,P_tab,eps_tab,k_tab,gamma_tab,gamma_th,  
			 neos,ergo_star,ergo_sigma,m,Symmetry,enable_OS_collapse);
}
