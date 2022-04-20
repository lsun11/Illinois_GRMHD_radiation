//---------------------------
// Conversion between P and P_thermal
// Used when the parameter reconstruct_Pthermal_instead_of_P == 1
//---------------------------
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "cctk.h"

extern void compute_pcold_epscold_cpp(double &rhob, double &P_cold, double &eps_cold,
				      int &neos, int &ergo_star, double &ergo_sigma, double *rho_tab, double *P_tab,
                               double *eps_tab, double *k_tab,
				      double *gamma_tab, int &enable_OS_collapse);

extern "C" void CCTK_FCALL convert_p_to_pthermal_cpp_
  (const cGH **cctkGH,int *ext, double *rho_b, double *P, 
   double *P_th, int &neos, int &ergo_star, double &ergo_sigma, double *rho_tab, double *P_tab,
   double *eps_tab,double *k_tab,double *gamma_tab, int &enable_OS_collapse);

extern "C" void convert_P_to_Pthermal_cpp(const cGH *cctkGH, int *ext,
				double *rho_b, double *P, double *P_th, 
				int &neos, int &ergo_star, double &ergo_sigma,
				double *rho_tab, double *P_tab,
				double *eps_tab, double *k_tab, 
					  double *gamma_tab, int &enable_OS_collapse)
{

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

          double rho_i = rho_b[index];
	  double P_cold, eps_cold;

	  compute_pcold_epscold_cpp(rho_i, P_cold, eps_cold,neos,ergo_star,ergo_sigma,rho_tab,P_tab,eps_tab,k_tab,gamma_tab, enable_OS_collapse);

	  P_th[index] = P[index] - P_cold;

	}	  
}

extern "C" void CCTK_FCALL convert_p_to_pthermal_cpp_
  (const cGH **cctkGH,int *ext, double *rho_b, double *P,
   double *P_th, int &neos, int &ergo_star, double &ergo_sigma, double *rho_tab, double *P_tab,
   double *eps_tab,double *k_tab,double *gamma_tab, int &enable_OS_collapse)
{
  convert_P_to_Pthermal_cpp(*cctkGH,ext,rho_b,P,P_th,neos, ergo_star, ergo_sigma, rho_tab,P_tab,eps_tab,
			    k_tab,gamma_tab, enable_OS_collapse);
}

// Note: Input gridfunction Pr and Pl as P_thermal, returns P as output.
extern "C" void CCTK_FCALL convert_pthermal_to_p_cpp_
  (const cGH **cctkGH,int *ext, double *rho_br, double *rho_bl, double *Pr,
   double *Pl, int &neos, int &ergo_star, double &ergo_sigma, double *rho_tab, double *P_tab,
   double *eps_tab,double *k_tab,double *gamma_tab, int &enable_OS_collapse);

extern "C" void convert_Pthermal_to_P_cpp(const cGH *cctkGH, int *ext,
                                double *rho_br, double *rho_bl, 
				double *Pr, double *Pl, int &neos,int &ergo_star, double &ergo_sigma,
                                double *rho_tab, double *P_tab,
                                double *eps_tab, double *k_tab,
					  double *gamma_tab, int &enable_OS_collapse)
{

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

          double rhor_i = rho_br[index];
	  double rhol_i = rho_bl[index];
          double P_cold, eps_cold;
          compute_pcold_epscold_cpp(rhor_i, P_cold, eps_cold,neos,
				    ergo_star, ergo_sigma,rho_tab,P_tab,eps_tab,k_tab,gamma_tab, enable_OS_collapse);
          Pr[index] += P_cold;
	  compute_pcold_epscold_cpp(rhol_i, P_cold, eps_cold,neos,
				    ergo_star,ergo_sigma,rho_tab,P_tab,eps_tab,k_tab,gamma_tab, enable_OS_collapse);
          Pl[index] += P_cold;

        }
}

extern "C" void CCTK_FCALL convert_pthermal_to_p_cpp_
  (const cGH **cctkGH,int *ext,double *rho_br, double *rho_bl, 
   double *Pr, double *Pl, int &neos, int &ergo_star, double &ergo_sigma, double *rho_tab, double *P_tab,
   double *eps_tab,double *k_tab,double *gamma_tab,int &enable_OS_collapse)
{
  convert_Pthermal_to_P_cpp(*cctkGH,ext,rho_br,rho_bl,Pr,Pl,neos,ergo_star,ergo_sigma,
			    rho_tab,P_tab,eps_tab,k_tab,gamma_tab, enable_OS_collapse);
}
