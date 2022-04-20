#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <sys/time.h>
#include "cctk.h"

#include "primitives_solver_header.h"
#define SQR(x) ((x) * (x))
double max(double val1,double val2);
double max(double val1, double val2) {
  if(val1>val2) return val1;
  return val2;
}
double min(double val1,double val2);
double min(double val1, double val2) {
  if(val1<val2) return val1;
  return val2;
}

void compute_P_T_microphys(double &P, double &T_fluid, double &P_cold,  double &eps, double &eps_cold, double &rho_b);
void compute_dP_nucldrho_microphys(double &dP_nucldrho, double &T_fluid, double &rho_b);
double find_T(double &T, double &eng_thermal_cgs, double &rad_const, double& n_nucleon, double& N_nu, double& T_scal);
void solve_T(double &Ta, double &Tb, double &T, double &eng_thermal_cgs, double &rad_const, double& n_nucleon, double& N_nu, double& T_scal);

double find_T_w(double &T, double &w_lhs_cgs, double &rad_const, double& n_nucleon, double& N_nu, double& T_scal);
void solve_T_w(double &Ta, double &Tb, double &T, double &w_lhs_cgs, double &rad_const, double& n_nucleon, double& N_nu, double& T_scal);
void compute_P_T_microphys_w(double &P, double &w, double &T_fluid, double &P_cold, double &eps_cold, double &rho_b);
void compute_opacity_emissivity(double &rho_star, double &al, double &Psi6, double &u0, double &kappa_a, double &kappa_s, double &eta_gf, double &kappa_a_nue, double &kappa_s_nue, double &eta_gf_nue, double &kappa_a_nux, double &kappa_s_nux, double &eta_gf_nux, double &T_fluid, double &chi_rad, double &chi_rad_nue, double &Y_e, double &optd, double &eta_nue, double &T_fluid_cgs_atm, int &microphysics_scheme, int &rad_fix);

void compute_eps_T_microphys (double &eps, double &T_fluid, double &eps_cold,  double &P, double &P_cold, double &rho_b);
double find_T_P(double &T, double &P_thermal_cgs, double &rad_const, double& n_nucleon, double& N_nu, double& T_scal);
void solve_T_P(double &Ta, double &Tb,double &T, double &P_thermal_cgs, double &rad_const, double& n_nucleon, double& N_nu, double& T_scal);

void compute_eps_th(double &T_fluid, double &rho_b, double &eps_thermal);



void compute_opacity_emissivity(double &rho_star, double &al, double &Psi6, double &u0, double &kappa_a, double &kappa_s, double &eta_gf, double &kappa_a_nue, double &kappa_s_nue, double &eta_gf_nue, double &kappa_a_nux, double &kappa_s_nux, double &eta_gf_nux, double &T_fluid, double &chi_rad, double &chi_rad_nue, double &Y_e, double &optd, double &eta_nue, double &T_fluid_cgs_atm, int &microphysics_scheme, int &rad_fix){

  //set nux_fac = 4.0 if rad_fix is 1
  double nux_fac;
  if (rad_fix==1){
    nux_fac = 4.0;
  }else{
    nux_fac = 1.0;
  }


  double T_fluid_cgs = max(1.0, T_fluid/(kb_cgs * G_cgs) * pow(c_cgs,4.0)*pow(k2c, 1.0));
  
  if (T_fluid_cgs > T_fluid_cgs_atm){

    double rho_b_cgs = rho_star/(al*Psi6*u0) * SQR(c_cgs)/G_cgs * pow(k2c, -2.0);
    double n_nucleon = rho_b_cgs/m_n_cgs;
    
    double T_o_T11_inv = pow(T_fluid_cgs/1.0e11, -1.0);
    double T_o_T11 = T_fluid_cgs/1.0e11;
    double rho_o_rho14 = rho_b_cgs/1.0e14;
    double rho_o_rho14_2o3 = pow(rho_b_cgs/1.0e14, 2.0/3.0);
  

       
    //chemical potential
    double eta_e = 3.505*rho_o_rho14_2o3*T_o_T11_inv;    
    double m_eta_e = eta_e;// chemical potential (e+) = chemical potential (e-) - 2*restmass of electron ~ chemical potential (e-)   
    double eta_n = 3.592*rho_o_rho14_2o3*T_o_T11_inv;
    double eta_p = 0.1128*SQR(rho_o_rho14_2o3)*T_o_T11_inv;
    double Q = 0.15*T_o_T11_inv; //(1.2935MeV/kbT)
    

    //    double eta_nue = 0.0;
    // eta_nue is a ouput variable
    double eta_nue_temp = eta_e + eta_p - eta_n - Q; 
    if ( eta_nue_temp < -20.0) eta_nue_temp = -20.0;
    if ( eta_nue_temp > 20.0) eta_nue_temp = 20.0;
    eta_nue = eta_nue_temp;

    double eta_anue = -1.0*eta_nue;
    double eta_nux = 0.0;
    
    
    //Fermi Functions
    double F5oF3 = Fermi5(eta_nue)/Fermi3(eta_nue);
    if (isnan(F5oF3))
      {
	printf("in compute_T_microphysics, F5oF3 is nan!!!! \n");
	printf("Fermi5(eta_nue), Fermi3(eta_nue) = %e,%e \n", Fermi5(eta_nue),Fermi3(eta_nue));
	printf("eta_nue, eta_e, eta_p, eta_n, Q = %e,%e,%e,%e,%e \n", eta_nue, eta_e, eta_p, eta_n, Q);
      }
    double F5oF3_a = Fermi5(eta_anue)/Fermi3(eta_anue);
    double F5oF3_x = Fermi5(eta_nux)/Fermi3(eta_nux);

    double F5oF4 = Fermi5(eta_nue)/Fermi4(eta_nue);
    double F5oF4_a = Fermi5(eta_anue)/Fermi4(eta_anue);
    
    //Ruffert et al (A15,A16) 
    //    double e_frac = 1.0/(1.0 + exp(3.505 * rho_o_rho14_2o3 * T_o_T11_inv - F5oF4));
    double ep_frac = 1.0/(1.0 + exp(- 3.505 * rho_o_rho14_2o3 * T_o_T11_inv - F5oF4_a));
    
    // set e_frac to 1 since it becomes too small in some regions. 
    double e_frac = 1.0/(1.0 + exp(3.505 * rho_o_rho14_2o3 * T_o_T11_inv - F5oF4)); 
    //    double e_frac = 1.0;


    double ka_gf_cgs, ks_gf_cgs, eta_gf_cgs, ka_gf_nue_cgs, ks_gf_nue_cgs, eta_gf_nue_cgs, ka_gf_nux_cgs, ks_gf_nux_cgs, eta_gf_nux_cgs;
    if (microphysics_scheme==0){
      // Simplified scheme, set nue nux terms to zero!
      kappa_a_nue = 0.0;
      kappa_s_nue = 0.0;
      eta_gf_nue = 0.0;
      kappa_a_nux = 0.0;
      kappa_s_nux = 0.0;
      eta_gf_nux = 0.0;
      //      Y_e = 0.0;

      double Yn = 1.0;
      double Yp = 0.0;

      double Ynp = 1.0;
      double Ypn = 0.0;

      //Ruffert et al (A11)
      double ka_nue_n_cgs =  min(0.0001, rho_b_cgs * k0 * sigma_0 * SQR(kb_cgs*T_fluid_cgs/(m_e_cgs * c_cgs * c_cgs))* F5oF3/ m_n_cgs * Ynp * e_frac);      
      if(isnan(ka_nue_n_cgs)){
	printf("ka_nue_n_cgs is nan, rho_b_cgs, k0, sigma_0, kb_cgs, T_fluid_cgs, m_e_cgs, F5oF3, m_n_cgs, Ynp, e_frac= %e,%e,%e,%e, %e,%e,%e,%e, %e,%e \n", rho_b_cgs, k0, sigma_0, kb_cgs, T_fluid_cgs, m_e_cgs, F5oF3, m_n_cgs, Ynp, e_frac);
      }
      //Ruffert et al (A12) 
      double ka_anue_p_cgs = min(0.0001, rho_b_cgs * k0 * sigma_0 * SQR(kb_cgs*T_fluid_cgs/(m_e_cgs * c_cgs * c_cgs)) * F5oF3_a/ m_n_cgs * Ypn * ep_frac);
      
      ka_gf_cgs = ka_nue_n_cgs + ka_anue_p_cgs;
      
      double Ynn = Yn/(1.0 + 2.0/3.0 * 3.592 * rho_o_rho14_2o3 * T_o_T11_inv);
      double Ypp = Yp/(1.0 + 2.0/3.0 * 1.1282 * 0.1 * SQR(rho_o_rho14_2o3) * T_o_T11_inv);
      //Ruffert et al (A6)
      //    double ks_n_cgs =  rho_b_cgs * k1 * sigma_0 * SQR(kb_cgs*T_fluid_cgs/(m_e_cgs * c_cgs * c_cgs)) * F5oF3/ m_n_cgs * Ynn;
      double ks_n_cgs;
      if (ka_nue_n_cgs == 0){
	ks_n_cgs = 0.0; 
      }else{
	ks_n_cgs = ka_nue_n_cgs / (k0 * Ynp * e_frac) * Ynn;
      }
      //    double ks_p_cgs =  rho_b_cgs * k2 * sigma_0 * SQR(kb_cgs*T_fluid_cgs/(m_e_cgs * c_cgs * c_cgs)) * F5oF3/ m_n_cgs * Ypp;
      double ks_p_cgs;
      if (ka_anue_p_cgs == 0){
	ks_p_cgs = 0.0;
      }else{
	ks_p_cgs = ka_anue_p_cgs / (k0 * Ypn * ep_frac) * Ypp * F5oF3/F5oF3_a;
      }
      double ks_gf_cgs = ks_n_cgs + ks_p_cgs;
      
      if (isnan(ks_gf_cgs)){
	printf ("ks_gf_cgs is nan, ks_n_cgs, ks_p_cgs, k0, e_frac, ep_frac = %e,%e,%e,%e,%e \n", ks_n_cgs, ks_p_cgs, k0, e_frac, ep_frac);
	printf ("ka_nue_n_cgs, ka_anue_p_cgs, F5oF3, Ynn, Ypp = %e,%e,%e,%e,%e \n", ka_nue_n_cgs, ka_anue_p_cgs, F5oF3, Ynn, Ypp);
	printf (" F5oF3_a, F5oF3_x, eta_nue = %e,%e,%e \n", F5oF3_a, F5oF3_x, eta_nue);
	printf ("T_fluid_cgs, rho_b_cgs = %e, %e \n", T_fluid_cgs, rho_b_cgs);
      }

      // Integrated Plank's law over frequnecy. ( photon: x1  neutrino: x7N/8 for N species) unit: erg cm^-3                   
      double B_T, B_T_nux;
      if (rad_fix==1){
	B_T = rad_const_cgs*pow(T_fluid_cgs,4.0) * 7.0/16.0;	
      }
      else
	{
	  B_T = rad_const_cgs*pow(T_fluid_cgs,4.0) * 7.0/8.0;
	}
      B_T_nux = nux_fac * B_T;

      //beta decay emissivity (LTE, using Kirchoff's law) unit: erg cm^-3 s^-1                            
      double eta_beta_LTE = ka_gf_cgs * B_T * c_cgs;
      double eta_gf_cgs = eta_beta_LTE;
    }
    else{
      // ===========================================================================================================================================//
      // =========================================================!!!!!!FULL SCHEME!!!!!!===========================================================//
      // ===========================================================================================================================================// 
      double Yn = 1.0 - Y_e;
      double Yp = Y_e;
      double Ynp = (2.0*Y_e - 1.0)/(exp( (1.1282 * 0.1 * rho_o_rho14_2o3 -3.592) * rho_o_rho14_2o3 * T_o_T11_inv ) - 1.0);
      double Ypn = Ynp * exp((1.1282 * 0.1 * rho_o_rho14_2o3 -3.592) * rho_o_rho14_2o3 * T_o_T11_inv);

      //Absorption:
      //Ruffert et al (A11) (electron neutrino) absrobed by neutron
      double ka_nue_n_cgs =  min(0.0001, rho_b_cgs * k0 * sigma_0 * SQR(kb_cgs*T_fluid_cgs/(m_e_cgs * c_cgs * c_cgs))* F5oF3/ m_n_cgs * Ynp * e_frac);

      if (isnan(ka_nue_n_cgs)){
        printf ("ka_nue_n_cgs is nan, Ynp, k0, e_frac, ep_frac = %e,%e,%e,%e \n", Ynp, k0, e_frac, ep_frac);
	printf ("Y_e, rho_o_rho14_2o3, T_o_T11_inv, F5oF3 = %e, %e, %e, %e \n", Y_e, rho_o_rho14_2o3, T_o_T11_inv, F5oF3);
        printf ("T_fluid_cgs, rho_b_cgs = %e, %e \n", T_fluid_cgs, rho_b_cgs);
      }

      //Ruffert et al (A12) (electron antineutrino) absrobed by proton
      double ka_anue_p_cgs = min(0.0001, rho_b_cgs * k0 * sigma_0 * SQR(kb_cgs*T_fluid_cgs/(m_e_cgs * c_cgs * c_cgs)) * F5oF3_a/ m_n_cgs * Ypn * ep_frac);

      ka_gf_cgs = ka_anue_p_cgs;
      ka_gf_nue_cgs = ka_nue_n_cgs;
      ka_gf_nux_cgs = 0.0; // Set heavy lepton term to zero here.
      //===========================================================================================================================================//
      //Scattering: Ruffert et al (A6) for electron neutrino and electron anitneutrino, zero for oterh species
      double Ynn = Yn/(1.0 + 2.0/3.0 * 3.592 * rho_o_rho14_2o3 * T_o_T11_inv);
      double Ypp = Yp/(1.0 + 2.0/3.0 * 1.1282 * 0.1 * SQR(rho_o_rho14_2o3) * T_o_T11_inv);
      double ks_nue_n_cgs, ks_anue_n_cgs, ks_nux_n_cgs;
      double ks_nue_p_cgs, ks_anue_p_cgs, ks_nux_p_cgs; 


      if (k0 * Ynp * e_frac == 0){
        ks_anue_n_cgs = 0.0;
	ks_anue_p_cgs = 0.0;
      }else{
        ks_anue_n_cgs = ka_nue_n_cgs / (k0 * Ynp * e_frac) * Ynn *k1;
	ks_anue_p_cgs = ka_nue_n_cgs / (k0 * Ynp * e_frac) * Ypp *k2;
      }


      if (k0 * Ypn * ep_frac == 0){
        ks_nue_n_cgs = 0.0;
        ks_nue_p_cgs = 0.0;
      }else{
        ks_nue_n_cgs = ka_anue_p_cgs / (k0 * Ypn * ep_frac) * Ynn * k1;
	ks_nue_p_cgs = ka_anue_p_cgs / (k0 * Ypn * ep_frac) * Ypp * k2;
      }

      ks_nux_n_cgs = ks_anue_n_cgs * F5oF3_x/F5oF3_a;
      ks_nux_p_cgs = ks_anue_p_cgs * F5oF3_x/F5oF3_a;

      ks_gf_cgs = min(0.0001, ks_anue_n_cgs + ks_anue_p_cgs);
      ks_gf_nue_cgs = min(0.0001, ks_nue_n_cgs + ks_nue_p_cgs);
      ks_gf_nux_cgs = min(0.0001, ks_nux_n_cgs + ks_nux_p_cgs);


      if(ks_gf_cgs > 1.0e-2 || ks_gf_nue_cgs > 1.0e-2 || ks_gf_nux_cgs > 1.0e-2)
	{
	  printf("Inside compute_T_microphysics.C, ks_gf_cgs,ks_gf_nue_cgs, ks_gf_nux_cgs=%e,%e,%e \n", ks_gf_cgs,ks_gf_nue_cgs,ks_gf_nux_cgs);
	  printf("ks_nue_n_cgs, ks_anue_n_cgs, ks_nux_n_cgs, ks_nue_p_cgs, ks_anue_p_cgs, ks_nux_p_cgs =%e,%e,%e,%e,%e,%e \n", ks_nue_n_cgs, ks_anue_n_cgs, ks_nux_n_cgs, ks_nue_p_cgs, ks_anue_p_cgs, ks_nux_p_cgs);
	  printf("Ynp,Ypn,Ynn,Ypp,e_frac,ep_frac,F5oF3_a,F5oF3_x=%e,%e,%e,%e,%e,%e,%e,%e \n",Ynp,Ypn,Ynn,Ypp,e_frac,ep_frac,F5oF3_a,F5oF3_x);
	}
      //===========================================================================================================================================// 
      //Emission:
      //Ruffert et al (B5)                                                                                      
      double eps_pos = eightpi_o_hc3 * pow(kb_cgs*T_fluid_cgs, 4.0) * Fermi3 (m_eta_e);
      double eps_mi = eightpi_o_hc3 * pow(kb_cgs*T_fluid_cgs, 4.0) * Fermi3 (eta_e);
      //Ruffert et al (B6)                                                                                                                                  
      double eps_tilde_pos = eightpi_o_hc3 * pow(kb_cgs*T_fluid_cgs, 5.0) * Fermi4 (m_eta_e);
      double eps_tilde_mi = eightpi_o_hc3 * pow(kb_cgs*T_fluid_cgs, 5.0) * Fermi4 (eta_e);
      //Ruffert et al (B7)                                                                                                                                    
      double eps_st_pos = eightpi_o_hc3 * pow(kb_cgs*T_fluid_cgs, 6.0) * Fermi5 (m_eta_e);
      double eps_st_mi = eightpi_o_hc3 * pow(kb_cgs*T_fluid_cgs, 6.0) * Fermi5 (eta_e);
    
    
      //beta decay emissivity (free) unit: erg s^-1 cm^-3 Ruffert et al (B15)     
      //Ruffert et al (B3,B4)                                                                                                                             
      double nue_frac = 1.0/(1.0 + exp(eta_nue - Fermi5(eta_e)/Fermi4(eta_e)));
      double anue_frac = 1.0/(1.0 + exp(eta_anue - Fermi5(m_eta_e)/Fermi4(m_eta_e)));
    
      //Ruffert et al (B1) (electron neutrino) emission of electron neutrino p + e^- --> n + nue               
      double eta_nue_free = rho_b_cgs * k0/2.0 * sigma_0 * c_cgs/SQR((m_e_cgs * c_cgs * c_cgs))/ m_n_cgs  * Ypn * eps_st_pos * nue_frac;
      //Ruffert et al (B2) (electron antineutrino) emission of electron antineutrino n + e^+ --> p + anue                                                                       
      double eta_anue_free = rho_b_cgs * k0/2.0 * sigma_0 * c_cgs/SQR((m_e_cgs * c_cgs * c_cgs))/ m_n_cgs  * Ynp * eps_st_mi * anue_frac;
    
      // Integrated Plank's law over frequnecy. ( photon: x1  neutrino: x7N/8 for N species) unit: erg cm^-3
      double B_T, B_T_nux;
      if (rad_fix==1){
	B_T = rad_const_cgs*pow(T_fluid_cgs,4.0) * 7.0/16.0;
      }
      else
        {
          B_T = rad_const_cgs*pow(T_fluid_cgs,4.0) * 7.0/8.0;
        }
      B_T_nux = nux_fac * B_T;

      //beta decay emissivity (LTE, using Kirchoff's law) unit: erg cm^-3 s^-1
      // (electron antineutrino) (A12) * B_T
      double eta_beta_LTE = ka_anue_p_cgs * B_T * c_cgs;
      // (electron neutrino) (A11) * B_T 
      double eta_beta_LTE_nue = ka_nue_n_cgs * B_T * c_cgs;

      //beta decay emissivity (total): intepolation between LTE and free                                                                     
      // NEED TO FIX CHI_RAD = 0 ISSUE!!!!
      double chi_radL = chi_rad;
      if (chi_radL < 1.0/3.0) {
	chi_radL = 1.0/3.0;
      }else if (chi_radL > 1.0) {
	chi_radL = 1.0;
      }
      else{ // test: set chi to 1 for some other cases such as NAN.
	chi_radL = 1.0;
      }

      double chi_rad_nueL = chi_rad_nue;
      if (chi_rad_nueL < 1.0/3.0) {
        chi_rad_nueL = 1.0/3.0;
      }else if (chi_rad_nueL > 1.0) {
        chi_rad_nueL = 1.0;
      }
      else{
	chi_rad_nueL = 1.0;
      }
      // [beta-process: electron antineutrino]
      double eta_beta = 1.5*(1.0-chi_radL) * eta_beta_LTE + (3.0*chi_radL -1.0)/2.0 * eta_anue_free;
      // [beta-process: electron neutrino]
      double eta_beta_nue = 1.5*(1.0-chi_rad_nueL) * eta_beta_LTE_nue + (3.0*chi_rad_nueL -1.0)/2.0 * eta_nue_free;
      if (isnan(eta_beta_nue)) printf("Inside compute_T_microphysics.C, eta_beta_nue is nan, chi_rad_nueL, eta_beta_LTE_nue, eta_nue_free=%e,%e,%e \n", chi_rad_nueL, eta_beta_LTE_nue, eta_nue_free);
      // [beta-process: heavy lepton neutrino]
      double eta_beta_nux = 0.0;	
      
      
      // Pair annihilation emissivity unit: erg s^-1 cm^-3 Ruffert et al (B16)                                                                    
      double nue_frac_ee = 1.0/(1.0 + exp(eta_nue - 0.5 * ( Fermi4(eta_e)/Fermi3(eta_e) + Fermi4(m_eta_e)/Fermi3(m_eta_e) ) ));
      double anue_frac_ee = 1.0/(1.0 + exp(eta_anue - 0.5 * ( Fermi4(eta_e)/Fermi3(eta_e) + Fermi4(m_eta_e)/Fermi3(m_eta_e) )));
      double nux_frac_ee = 1.0/(1.0 + exp(eta_nux - 0.5 * ( Fermi4(eta_e)/Fermi3(eta_e) + Fermi4(m_eta_e)/Fermi3(m_eta_e) ) ));

      //(B8) [pair: electron antineutrino]
      double eta_ee = 2.0*2.3432/36.0 * sigma_0 * c_cgs/SQR((m_e_cgs * c_cgs * c_cgs)) * 0.5*(eps_tilde_mi*eps_pos + eps_tilde_pos*eps_mi) * anue_frac_ee * nue_frac_ee;
      //(B8) [pair: electron neutrino]  
      double eta_ee_nue = 2.0*2.3432/36.0 * sigma_0 * c_cgs/SQR((m_e_cgs * c_cgs * c_cgs)) * 0.5*(eps_tilde_mi*eps_pos + eps_tilde_pos*eps_mi) * anue_frac_ee * nue_frac_ee;
      //(B10) [pair: heavy lepton neutrino]  
      double eta_ee_nux = nux_fac * 0.05591 * sigma_0 * c_cgs/SQR((m_e_cgs * c_cgs * c_cgs)) * 0.5*(eps_tilde_mi*eps_pos + eps_tilde_pos*eps_mi) * SQR(nux_frac_ee);    
      
           
      // Plasmon annihilation emissivity unit: erg s^-1 cm^-3 Ruffert et al (B17)                                                                                                  
      double gamma_plas = 5.565e-2 * sqrt( (SQR(M_PI) + 3.0 * SQR(eta_e))/3.0 );
      double nue_frac_plas = 1.0/(1.0 + exp(eta_nue -  (1.0  + 0.5 *  SQR(gamma_plas)/(1.0+gamma_plas)) ));
      double anue_frac_plas = 1.0/(1.0 + exp(eta_anue -  (1.0  + 0.5 *  SQR(gamma_plas)/(1.0+gamma_plas)) ));
      double Coe_plas = 1305.2854674928328; //pi^3/(3*alpha^\star)*C_v^2 in (B11)
      double Coe_plas2 = pow(gamma_plas, 6.0) *exp(-gamma_plas)*(1.0+gamma_plas);//gamma^6*exp(-gamma)*(1+gamma) in (B11/B12)
      double Coe_plas3 = 9.064482413144674; //4pi^3/(3*alpha^\star)*(C_v-1)^2 in (B12)   
      double nux_frac_plas = 1.0/(1.0 + exp(eta_nux -  (1.0  + 0.5 *  SQR(gamma_plas)/(1.0+gamma_plas)) ));

      //(B11) [plasmon: electron antineutrino]
      double eta_plasmon = Coe_plas * sigma_0 * c_cgs/SQR((m_e_cgs * c_cgs * c_cgs)) * pow (kb_cgs*T_fluid_cgs, 8.0) * pow( hbar_cgs*c_cgs*2.0*M_PI,-6.0) * Coe_plas2 * nue_frac_plas * anue_frac_plas;
      //(B11) [plasmon: electron neutrino]
      double eta_plasmon_nue = Coe_plas * sigma_0 * c_cgs/SQR((m_e_cgs * c_cgs * c_cgs)) * pow (kb_cgs*T_fluid_cgs, 8.0) * pow( hbar_cgs*c_cgs*2.0*M_PI,-6.0) * Coe_plas2 * nue_frac_plas * anue_frac_plas;
      //(B12) [plasmon: heavy lepton neutrino]
      double eta_plasmon_nux = nux_fac * Coe_plas3 * sigma_0 * c_cgs/SQR((m_e_cgs * c_cgs * c_cgs)) * pow (kb_cgs*T_fluid_cgs, 8.0) * pow( hbar_cgs*c_cgs*2.0*M_PI,-6.0) * Coe_plas2 * SQR(nux_frac_plas);

          
      // Neucleon bremsstrahlung emissivity unit: erg s^-1 cm^-3 in arXiv:astro-ph/9905132 Eqn.(49) ( set zeta = 0.5 )                                           
      double XNsq = SQR(Yn) + SQR(Yp) + 28.0/3.0 * Yn * Yp;
      if (isnan(XNsq)) {printf ("XNsq is nan!!!! Yn, Yp, Y_e = %e, %e, %e \n", Yn, Yp, Y_e);}
      // Same for all three species //
      // [Bremsstrahlung: electron antineutrino]
      double eta_brem = 1.04e30 * 0.5 * XNsq * SQR(rho_o_rho14) * pow(T_o_T11, 5.5);
      // [Bremsstrahlung: electron neutrino]
      double eta_brem_nue = 1.04e30 * 0.5 * XNsq * SQR(rho_o_rho14) * pow(T_o_T11, 5.5);
      // [Bremsstrahlung: heavy lepton neutrino]
      double eta_brem_nux = nux_fac * 1.04e30 * 0.5 * XNsq * SQR(rho_o_rho14) * pow(T_o_T11, 5.5);
      //===========================================================================================================================================//  
      
      // Compute the final sum of ka and eta
      // Here we sum ee + plasmon + brem in order to get the corresponding inverse opacitiy
      double eta_sum = eta_ee + eta_plasmon + eta_brem;
      double eta_sum_nue = eta_ee_nue + eta_plasmon_nue + eta_brem_nue;
      double eta_sum_nux = eta_ee_nux + eta_plasmon_nux + eta_brem_nux;
      
      // Total ka = ka_beta + ka_sum, where ka_sum comes from Kirchoff's law using eta_sum
      double ka_sum = eta_sum/B_T/c_cgs;
      double ka_sum_nue = eta_sum_nue/B_T/c_cgs;
      double ka_sum_nux = eta_sum_nux/B_T_nux/c_cgs;

      double ka_gf_cgs_temp = ka_gf_cgs + ka_sum;
      double ka_gf_nue_cgs_temp = ka_gf_nue_cgs + ka_sum_nue;
      double ka_gf_nux_cgs_temp = ka_gf_nux_cgs + ka_sum_nux;

      // Set a cap on ka
      ka_gf_cgs = min(0.0001, ka_gf_cgs_temp);
      ka_gf_nue_cgs= min(0.0001, ka_gf_nue_cgs_temp);
      ka_gf_nux_cgs= min(0.0001, ka_gf_nux_cgs_temp);

      // Now set the final emissivity using eta = eta_beta + eta_sum
      eta_gf_cgs = eta_beta + eta_sum;
      eta_gf_nue_cgs = eta_beta_nue + eta_sum_nue;
      eta_gf_nux_cgs = eta_beta_nux + eta_sum_nux;

      if (isnan(ka_gf_nue_cgs) || isnan(eta_gf_nue_cgs))
	{
	  printf("Inside compute_T_microphysics.C, ka_gf_nue_cgs, eta_gf_nue_cgs= %e, %e \n", ka_gf_nue_cgs, eta_gf_nue_cgs);
	  printf("ka_gf_nue_cgs, ka_sum_nue, eta_beta_nue, eta_sum_nue = %e,%e,%e,%e \n", ka_gf_nue_cgs, ka_sum_nue, eta_beta_nue, eta_sum_nue);
	  printf("ka_nue_n_cgs = %e \n", ka_nue_n_cgs);
	  printf("eta_beta_LTE_nue, eta_nue_free, eta_ee_nue, eta_plasmon_nue, eta_brem_nue = %e,%e,%e,%e, %e \n", eta_beta_LTE_nue, eta_nue_free, eta_ee_nue, eta_plasmon_nue, eta_brem_nue);
	}
      
      // TEST!!!!!!!!!!!!!!!!!!!!!!!!!
      //    double eta_gf_cgs = eta_beta_LTE    
      /*
	if(eta_gf_cgs > 1.0e41) {
	printf("eta_gf_cgs > 1.0e41, eta_beta, eta_pair, eta_plasmon, eta_brem = %e,%e,%e,%e \n", eta_beta,eta_pair,eta_plasmon,eta_brem);
	printf("eta_beta_LTE, eta_beta_free, eta_ee_nue, eta_ee_nux, eta_plasmon_nue, eta_plasmon_nux = %e, %e, %e, %e, %e, %e \n", eta_beta_LTE, eta_beta_free, eta_ee_nue, eta_ee_nux, eta_plasmon_nue, eta_plasmon_nux);
	printf ("T_fluid, T_fluid_cgs, T_o_T11, T_o_T11_inv, rho_o_rho14, rho_o_rho14_2o3, F5oF3, F5oF3_a, F5oF4, F5oF4_a = %e, %e, %e, %e, %e, %e, %e, %e, %e, %e \n", T_fluid, T_fluid_cgs, T_o_T11, T_o_T11_inv, rho_o_rho14, rho_o_rho14_2o3, F5oF3, F5oF3_a, F5oF4, F5oF4_a);
	printf ("chi_rad, XNsq, gamma_plas, nue_frac_plas, anue_frac_plas, Coe_plas2, nux_frac_plas = %e,%e,%e,%e, %e,%e,%e\n", chi_rad, XNsq, gamma_plas, nue_frac_plas, anue_frac_plas, Coe_plas2, nux_frac_plas);
	printf ("B_T, ka_gf_cgs, ks_gf_cgs = %e, %e, %e \n", B_T, ka_gf_cgs, ks_gf_cgs);      
	printf ("ka_nue_n_cgs, ka_anue_p_cgs, e_frac, ep_frac = %e,%e,%e,%e \n", ka_nue_n_cgs, ka_anue_p_cgs, e_frac, ep_frac);
	
	printf ("rho_b_cgs, sigma_0 * SQR(kb_cgs*T_fluid_cgs/(m_e_cgs * c_cgs * c_cgs)) = %e, %e \n", rho_b_cgs, sigma_0 * SQR(kb_cgs*T_fluid_cgs/(m_e_cgs * c_cgs * c_cgs)) );
	printf("----------------------------------------------------------------------------------------------------------------------------------------------------------\n");
	
	}
      */ 
	
    }
    
    // convert kappas from cm^-1 to km^-1                                                                          
    kappa_a = max(0.0, ka_gf_cgs * pow(c2k, -1.0));
    kappa_s = max(0.0, ks_gf_cgs * pow(c2k, -1.0));
    kappa_a_nue = max(0.0, ka_gf_nue_cgs * pow(c2k, -1.0));
    kappa_s_nue = max(0.0, ks_gf_nue_cgs * pow(c2k, -1.0));
    kappa_a_nux = max(0.0, ka_gf_nux_cgs * pow(c2k, -1.0));
    kappa_s_nux = max(0.0, ks_gf_nux_cgs * pow(c2k, -1.0));    
    // convert eta from erg s^-1 cm^-3 to km^-3                                                                                             
    eta_gf = max(0.0, eta_gf_cgs * G_cgs * pow (c_cgs, -5.0) * pow(c2k, -3.0));
    eta_gf_nue = max(0.0, eta_gf_nue_cgs * G_cgs * pow (c_cgs, -5.0) * pow(c2k, -3.0));
    eta_gf_nux = max(0.0, eta_gf_nux_cgs * G_cgs * pow (c_cgs, -5.0) * pow(c2k, -3.0));

    if ( fabs(kappa_a) > 10.0 || fabs(kappa_a_nue) > 10.0 || fabs(kappa_a_nux) > 10.0 ||  fabs(eta_gf) > 10.0 || fabs(eta_gf_nue) > 10.0 || fabs(eta_gf_nux) > 10.0)
      {
	printf("Inside compute_T_microphysics.C, kappa_a, kappa_a_nue, kappa_a_nux = %e,%e,%e \n", kappa_a, kappa_a_nue, kappa_a_nux);
	printf("ka_gf_cgs, ks_gf_cgs, ka_gf_nue_cgs, ks_gf_nue_cgs, ka_gf_nux_cgs, ks_gf_nux_cgs = %e,%e,%e,%e,%e,%e \n", ka_gf_cgs, ks_gf_cgs, ka_gf_nue_cgs, ks_gf_nue_cgs, ka_gf_nux_cgs, ks_gf_nux_cgs);
	printf("eta_gf, eta_gf_nue, eta_gf_nux = %e,%e,%e \n", eta_gf, eta_gf_nue, eta_gf_nux);
	printf("eta_gf_cgs, eta_gf_nue_cgs, eta_gf_nux_cgs = %e,%e,%e \n", eta_gf_cgs, eta_gf_nue_cgs, eta_gf_nux_cgs);
      }

  }
  else{
    kappa_a = 0.0;
    kappa_s = 0.0;
    eta_gf = 0.0;
    kappa_a_nue = 0.0;
    kappa_s_nue = 0.0;
    eta_gf_nue = 0.0;
    kappa_a_nux = 0.0;
    kappa_s_nux = 0.0;
    eta_gf_nux = 0.0;
  }
  
}



void compute_eps_th(double &T_fluid, double &rho_b, double &eps_thermal){
  double rho_b_cgs = fasterpow_prim(c_cgs,2.0) / G_cgs * rho_b * fasterpow_prim(k2c, -2.0);
  double n_nucleon = rho_b_cgs/m_n_cgs;
  double T_fluid_cgs = T_fluid/(kb_cgs * G_cgs) * fasterpow_prim(c_cgs,4.0) * fasterpow_prim(k2c, 1.0);

  double eng_thermal_cgs = 1.5*n_nucleon*kb_cgs*T_fluid_cgs + (11.0/4.0)*rad_const_cgs * fasterpow_prim(T_fluid_cgs, 4.0);
  double eng_thermal = eng_thermal_cgs * G_cgs * fasterpow_prim(c_cgs,-4.0) * fasterpow_prim(c2k, -2.0);

  eps_thermal = eng_thermal/rho_b;

}




extern "C" void CCTK_FCALL CCTK_FNAME(compute_P_T_microphys)
(double &P, double &T_fluid, double &P_cold,  double &eps, double &eps_cold, double &rho_b);


void compute_P_T_microphys(double &P, double &T_fluid, double &P_cold,  double &eps, double &eps_cold, double &rho_b){

  //|||| Now include thermal and nuclear parts   ||||//                                                                                                   
  // Note that the code unit quantities are powers of km, we need to convert them to powers of cm first!!!///                                               
  bool recom = false;
  double eps_thermal = eps - eps_cold;
  
  /*  if (eps_thermal < 1.0e-8){
    eps_thermal = fasterpow_prim(rho_b,0.357)*13.0;
    //Need to recompute eps_tot if eps_thermal is revalued.
    eps = eps_thermal + eps_cold;
  }
  */

  double eng_thermal = rho_b*eps_thermal;
  double eng_thermal_cgs = fasterpow_prim(c_cgs,4.0) / G_cgs * eng_thermal * fasterpow_prim(k2c, -2.0);
  double rho_b_cgs = fasterpow_prim(c_cgs,2.0) / G_cgs * rho_b * fasterpow_prim(k2c, -2.0);
  double n_nucleon = rho_b_cgs/m_n_cgs;
  double N_nu = 3.0;
  double T_fluid_cgs;
  double T_scal= 1.0e0;
  //double T_scal = 1.0;                                                                                             
  //  printf("Begin compute_P_T_microphys, T_fluid!!! T_fluid=%e \n", T_fluid);                                                              

  double T_fluid_cgs_floor = 1.0e1;

  //check if the thermal part if positive:                                                                                                            
  if (eps_thermal <= 0.0){
    //printf("eps_thermal = %e \n", eps_thermal);                                                                                        
    T_fluid_cgs = T_fluid_cgs_floor;
  }
  else{
    // If we use this part rad_const_cgs is set to be 7.5646e-15 erg cm^-3 K^-4                                                                   
    // Solve for T_fluid                                                                                                                       
    // As 1.5nkBT >> 11/4*aT^4 in most cases, make the initial guess of Tcgs = eng_thermal_cgs/(1.5*n*kB) and then solve for Tcgs.                       
    // If T_fluid is postive (calculated in the last iteration), use it as the guess, otherwise  use Tcgs = eng_thermal_cgs/(1.5*n*kB).
    // T guess needs to be adjusted for high temperature, low density region (rare).                                                                     
    
    //double threshold = fasterpow_prim(eng_thermal_cgs/(kb_cgs*0.21512), 0.75);
    double threshold = fasterpow_prim(eng_thermal_cgs,0.75)*0.8585*fasterpow_prim(rad_const_cgs,0.25)/kb_cgs;
    double T_fluid_cgs_guess;
    if (n_nucleon < threshold){
      //T_fluid_cgs_guess= eng_thermal_cgs/(1.5 * kb_cgs * n_nucleon);   
      T_fluid_cgs = fasterpow_prim(4.0*eng_thermal_cgs/rad_const_cgs/11.0, 0.25);
      T_fluid_cgs_guess = fasterpow_prim(4.0*eng_thermal_cgs/rad_const_cgs/11.0, 0.25);
    }
    else{
      //T_fluid_cgs = fasterpow_prim(4.0*eng_thermal_cgs/rad_const_cgs/11.0, 0.25);     
      T_fluid_cgs = eng_thermal_cgs/(1.5 * kb_cgs * n_nucleon);
      T_fluid_cgs_guess= eng_thermal_cgs/(1.5 * kb_cgs * n_nucleon);
    } 
 
    //The final T should be smaller than T_guess since T_guess only consists of one of the two parts of eng_thermal.
    // Therefore, no matter we set F1(T_guess) or F2(T_guess) to eng_thermal, F1(T_guess) + F2(T_guess) > eng_thermal and thus real T is smaller. 
  double Ta = T_fluid_cgs*0.8;
  double Tb = T_fluid_cgs*1.3;
  
  double rad_constl = rad_const_cgs;
  //  if (T_fluid_cgs == 0.0) printf("inside compute_P_T_microphys, T_fluid_cgs is 0!!! T_fluid=%e \n", T_fluid);                                
  solve_T(Ta, Tb, T_fluid_cgs, eng_thermal_cgs, rad_constl, n_nucleon, N_nu, T_scal);

  // WE DO A CHECK AFTER GETTING THE TEMPERATURE!!!!
  double E1 = (11.0/4.0 * rad_const_cgs)* fasterpow_prim(T_fluid_cgs, 4.0);
  double E2 = 1.5*n_nucleon*T_fluid_cgs*kb_cgs;
  double E_ratio = (E1 + E2)/eng_thermal_cgs;
  if (T_fluid_cgs > T_fluid_cgs_floor && E_ratio > 1.05){
    printf("after solve_T, E_ratio too large, is %.16e, E1=%.16e, E2=%.16e, eng_thermal_cgs=%.16e, T_fluid_cgs=%.16e, T_fluid_cgs_guess=%.16e, n_nucleon=%.16e \n", E_ratio, E1, E2, eng_thermal_cgs, T_fluid_cgs, T_fluid_cgs_guess, n_nucleon);
  }
  
  //set T_fluid_cgs floor (10 K)
  if (T_fluid_cgs < T_fluid_cgs_floor){
    T_fluid_cgs = T_fluid_cgs_floor;    
  }

  //set T_fluid_cgs cap (10^12 K ~ 100 MeV)   
  if (T_fluid_cgs > 1.0e11* T_fluid_cgs_floor){
    T_fluid_cgs = 1.0e11*T_fluid_cgs_floor;
  }
  }

  //  Here fasterpow_prim(c2k, 1.0) means converting from cm^1 to km^1 (Temperature scales as [L]);                             
  T_fluid = T_fluid_cgs * kb_cgs * G_cgs / fasterpow_prim(c_cgs,4.0) * fasterpow_prim(c2k, 1.0);
 
  //  double Fac2 = min_val(1.0,  2.0*fasterpow_prim((3.0*M_PI*M_PI), 1.0/3.0)*m_n_cgs/(18.0*hbar_cgs*hbar_cgs)*fasterpow_prim(n_nucleon, -2.0/3.0)*kb_cgs*T_fluid_cgs  );       
  double Fac2 = 1.0;
  double P_nucl_cgs = (5.0/3.0-1.0)*(1.5*n_nucleon*kb_cgs*T_fluid_cgs*Fac2);
  double P_nucl = fasterpow_prim(c_cgs, -4.0) * fasterpow_prim(G_cgs, 1.0) * P_nucl_cgs * fasterpow_prim(c2k, -2.0);

  //double P_rad_m_cgs = (3.0 + 7.0 * N_nu/8.0) * rad_const_cgs * fasterpow_prim(T_fluid_cgs, 4.0)/3.0;                                                
  double P_rad_ph_cgs = (11.0/12.0)*rad_const_cgs * fasterpow_prim(T_fluid_cgs, 4.0);
  double P_rad_ph = fasterpow_prim(c_cgs, -4.0) * fasterpow_prim(G_cgs, 1.0) * P_rad_ph_cgs * fasterpow_prim(c2k, -2.0);

  //  double P_thermal =  E_rad/3.0 + P_nucl;                                                                                                               
  double P_thermal = P_rad_ph + P_nucl;
  P = P_cold + P_thermal;
  
  if (T_fluid > 1.0e-30){
    printf("Inside compute_P_T_microphys, T_fluid is too large, is %.16e, eps=%.16e, eps_cold=%.16e, eps_thermal=%.16e, P=%.16e, P_cold=%.16e, P_thermal=%.16e, T_fluid_cgs=%.16e, n_nucleon=%.16e \n", T_fluid, eps, eps_cold, eps_thermal, P, P_cold, P_thermal,T_fluid_cgs, n_nucleon);
    printf("eng_thermal_cgs = %.16e, rho_b_cgs = %.16e, rho_b = %.16e \n", eng_thermal_cgs, rho_b_cgs, rho_b);
  }

}


extern "C" void CCTK_FCALL CCTK_FNAME(compute_P_T_microphys)
  (double &P, double &T_fluid, double &P_cold,  double &eps, double &eps_cold, double &rho_b) {

  compute_P_T_microphys(P,T_fluid,P_cold,eps,eps_cold,rho_b);

}




//|||| Done with getting thermal and nuclear parts   ||||//  

void compute_dP_nucldrho_microphys(double &dP_nucldrho, double &T_fluid, double &rho_b)                                                         
{                                                                           
  double dP_nucldrho_cgs;       
  //  double T_fluid_cgs = T_fluid * fasterpow_prim(c_cgs, 4.0) * M_sun_geo / (G_cgs * kb_cgs);                   
  double T_fluid_cgs = T_fluid * fasterpow_prim(c_cgs,4.0)/(G_cgs*kb_cgs) * fasterpow_prim(k2c, 1.0);                            
  double rho_b_cgs = fasterpow_prim(c_cgs,2.0) / G_cgs * rho_b * fasterpow_prim(k2c, -2.0);                          
  double m2_o_rho5 = fasterpow_prim(m_n_cgs*m_n_cgs/fasterpow_prim(rho_b_cgs,5.0) , 1.0/3.0);               
  double kbT = T_fluid_cgs*kb_cgs;                                                                                
  double T0_limits = fasterpow_prim((3.0*M_PI*M_PI), 1.0/3.0)*fasterpow_prim(m_n_cgs, 2.0/3.0)/(18.0*hbar_cgs*hbar_cgs)*fasterpow_prim(rho_b_cgs, -2.0/3.0)*SQR(kbT);  
                                                                                                                                                   
  //dP_nucldrho_cgs = (5.0/3.0-1.0)*min_val(1.5*kbT/m_n_cgs,  T0_limits);                                   
  dP_nucldrho_cgs = (5.0/3.0-1.0)*1.5*kbT/m_n_cgs;                                                                                                     
  dP_nucldrho = dP_nucldrho_cgs /( c_cgs * c_cgs);            
}               



double find_T(double &T, double &eng_thermal_cgs, double &rad_const, double& n_nucleon, double& N_nu, double& T_scal){
  //  double Fac = min_val(1.0,  2.0*fasterpow_prim((3.0*M_PI*M_PI), 1.0/3.0)*m_n_cgs/(18.0*hbar_cgs*hbar_cgs)*fasterpow_prim(n_nucleon, -2.0/3.0)*kb_cgs*T );                   
  double Fac = 1.0;
  double diff = (11.0/4.0 * rad_const)* fasterpow_prim(T*T_scal, 4.0) + 1.5*n_nucleon*T*T_scal*kb_cgs*Fac - eng_thermal_cgs;
  //  if (diff == 0) {printf("Inside find_T, diff is 0!!!, T=%.16e, (11.0/4.0 * rad_const)* fasterpow_prim(T*T_scal, 4.0)=%.16e, 1.5*n_nucleon*T*T_scal*kb_cgs*Fac=%.16e,  eng_thermal_cgs=%.16e \n", T, (11.0/4.0 * rad_const)* fasterpow_prim(T*T_scal, 4.0), 1.5*n_nucleon*T*T_scal*kb_cgs, eng_thermal_cgs);}
  return diff/eng_thermal_cgs;
  //return diff;                                                                                                                                                                  
}


void solve_T(double &Ta, double &Tb,
             double &T, double &eng_thermal_cgs, double &rad_const, double& n_nucleon, double& N_nu, double& T_scal)
{
  double fa, fb, fc;
  int bisec = 1;
  double T_init =T;
  double tiny = 1.0;
  int iter=0;
  int iter1=0;

  //If T_cgs_guess is smaller than 1K, set T_cgs to 1K, don't do the bisection!
  if (T_init < tiny){
    T = tiny;
    bisec = 0;
  }
  else{ 
    Ta /=T_scal;
    Tb /=T_scal;

    // Obtain the result fa = F1(T_guess_a) + F2(T_guess_a) - eng_thermal_cgs & fb = F1(T_guess_b) + F2(T_guess_b) - eng_thermal_cgs
    // In principle, fb >  F1(T_guess) + F2(T_guess) - eng_thermal_cgs > 0, and we want fa < 0
    fa = find_T(Ta, eng_thermal_cgs, rad_const, n_nucleon, N_nu, T_scal);
    fb = find_T(Tb, eng_thermal_cgs, rad_const, n_nucleon, N_nu, T_scal);
 
    while (fa *fb >= 0)
      {
	if (iter > 20) {//We fail to converge to find a proper Ta -- Tb range!!!!
	  T = T_init;
	  bisec = 0;
	  printf("break loop: Ta = %e, Tb =%e, Ta_init= %e, Tb_init=%e, eng_thermal_cgs = %e, fa=%e, fb=%e \n", Ta, Tb, T_init*0.8, T_init*1.05, eng_thermal_cgs, fa, fb);
	  break;
	}
	Ta *= 0.5;
	Tb *= 1.5;
	fa = find_T(Ta, eng_thermal_cgs, rad_const, n_nucleon, N_nu, T_scal);
	fb = find_T(Tb, eng_thermal_cgs, rad_const, n_nucleon, N_nu, T_scal);
	iter += 1;
      } 
    
    /*    if (T_init > 1.0e10){
    printf ("Done finding boundary!!! Ta=%e,Tb=%e,fa=%e,fb=%e,iter=%d,eng_thermal_cgs=%e \n", Ta,Tb,fa,fb,iter,eng_thermal_cgs);
    }*/
    //Now we find the proper fa < eng_thermal_cgs and fb > eng_thermal_cgs, start bisection to find T.
 
    if (bisec){
      T = Ta;
      while ( fabs(Tb-Ta)/(Tb+Ta) >= 1.0e-12)
	{
	  // Find middle point                                                                                              
	  T = (Ta+Tb)/2;
	  fc = find_T(T, eng_thermal_cgs, rad_const, n_nucleon, N_nu, T_scal);
	  // Check if middle point is root                                                                                      
	  if (fc == 0.0)
	    break;
	  // Decide the side to repeat the steps                                                       
	  else if (fc*fa < 0)
	    Tb = T;
	  else
	    Ta = T;
	  
	  iter1 +=1;
	  
	  if (iter1 > 50){
	    // Failed to converge after 50 iteration, do 50 more!!!!!!!
	    //	 printf("inside solve_T, over iteration, fabs(Tb-Ta)/(Tb+Ta)=%e \n", fabs(Tb-Ta)/(Tb+Ta));
	    int iter2 = 0;
	    while (fabs(Tb-Ta)/(Tb+Ta) >= 1.0e-12)
	      {
		// Find middle point                                                                                                                         
		T = (Ta+Tb)/2;
		fc = find_T(T, eng_thermal_cgs, rad_const, n_nucleon, N_nu, T_scal);
		// Check if middle point is root                                                                               
		if (fc == 0.0)
		  break;
		// Decide the side to repeat the steps                                                                                                           
		else if (fc*fa < 0)
		  Tb = T;
		else
		  Ta = T;
		
		iter2 +=1;
	      }
	    if (iter2 > 50){
	      printf("inside solve_T, over iteration AGAIN, fabs(Tb-Ta)/(Tb+Ta)=%e \n", fabs(Tb-Ta)/(Tb+Ta));      
	      break;
	    }	 
	  }

	}
 
      T *=T_scal;
      Ta *=T_scal;
      Tb *=T_scal;  

    } // End of bisection       
  } // End of T > T_tiny
  /*  if (T_init > 1.0e10){
    printf ("At the end of solve_T, T=%e, T_init=%e, bisec=%d, iter1=%d,  Ta=%e, Tb=%e, abs(Tb-Ta)/(Tb+Ta)=%e \n", T, T_init,bisec, iter1, Ta,Tb,abs(Tb-Ta)/(Tb+Ta));
    printf ("Tb-Ta = %e, fabs(Tb-Ta)=%e, Tb+Ta=%e, fabs(Tb-Ta)/(Tb+Ta)=%e  \n", Tb-Ta, fabs(Tb-Ta), (Tb+Ta), fabs(Tb-Ta)/(Tb+Ta) );
    }*/
} 


//!!!!!!!!!!!!!!!!!!!!! The following three functions are used only when primitives_solver = 11.  i.e. the harm P.S. !!!!!!!!!!!!!!!!!!!!!!                                       
//w_lhs = w - rho*eps_cold - P_cold - rho_b = rho_b*eps_th + P_th                                                                                                                 
void compute_P_T_microphys_w(double &P, double &w, double &T_fluid, double &P_cold, double &eps_cold, double &rho_b){
  double rho_b_cgs = fasterpow_prim(c_cgs,2.0) / G_cgs * rho_b * fasterpow_prim(k2c, -2.0);
  double n_nucleon = rho_b_cgs/m_n_cgs;
  double N_nu = 3.0;
  double T_fluid_cgs;

  double w_lhs = w - rho_b*eps_cold - P_cold -rho_b;

  /*
  if (w_lhs < 1.0e-19){
    w_lhs = fasterpow_prim(rho_b, 1.357)*19.5;
    //Need to recompute w if w_lhs is revalued.                                   
    w = w_lhs + rho_b*eps_cold + P_cold +rho_b;
  }
  */

  double w_lhs_cgs = fasterpow_prim(c_cgs,4.0) / G_cgs *  w_lhs * fasterpow_prim(k2c, -2.0);
  double T_scal = 1.0e0;

  double T_fluid_cgs_floor = 1.0e1;
  double T_fluid_cgs_guess;
  if(w_lhs <=0.0){
    T_fluid_cgs = T_fluid_cgs_floor;
    //    T_fluid = T_fluid_cgs_floor * kb_cgs * G_cgs / fasterpow_prim(c_cgs,4.0) * fasterpow_prim(c2k, 1.0);
    //    P = P_cold;
  }
  else{
    //In most cases, 11/3 aT^4 << 2.5 nkBT , so use w_lhs/2.5kB as Tcgs guess.   
    double threshold = fasterpow_prim(w_lhs_cgs,0.75)*0.553513*fasterpow_prim(rad_const_cgs,0.25)/kb_cgs;
    if (n_nucleon < threshold){
      T_fluid_cgs = fasterpow_prim(3.0*w_lhs_cgs/rad_const_cgs/11.0, 0.25);
      T_fluid_cgs_guess = fasterpow_prim(3.0*w_lhs_cgs/rad_const_cgs/11.0, 0.25);
    }
    else{
      T_fluid_cgs = w_lhs_cgs/ (2.5 * kb_cgs * n_nucleon);
      T_fluid_cgs_guess = w_lhs_cgs/ (2.5 * kb_cgs * n_nucleon);
    }

    double Ta = T_fluid_cgs*0.9;
    double Tb = T_fluid_cgs*1.1;

    double rad_constl = rad_const_cgs;
    solve_T_w(Ta, Tb, T_fluid_cgs, w_lhs_cgs, rad_constl, n_nucleon, N_nu, T_scal);
    // WE DO A CHECK AFTER GETTING THE TEMPERATURE!!!!                                                                                                                        
  
    double E1 = (11.0/3.0 * rad_const_cgs)* fasterpow_prim(T_fluid_cgs, 4.0);
    double E2 = 2.5*n_nucleon*T_fluid_cgs*kb_cgs;
    double E_ratio = (E1 + E2)/w_lhs_cgs;
    if (T_fluid_cgs > T_fluid_cgs_floor && E_ratio > 1.05){
      printf("after solve_T_w, E_ratio too large, is %.16e, E1=%.16e, E2=%.16e, w_lhs_cgs=%.16e, T_fluid_cgs=%.16e, T_fluid_cgs_guess=%.16e, n_nucleon=%.16e \n", E_ratio, E1, E2, w_lhs_cgs, T_fluid_cgs, T_fluid_cgs_guess, n_nucleon);
    }
    
    if (T_fluid_cgs < T_fluid_cgs_floor){
      T_fluid_cgs = T_fluid_cgs_floor;
    }
    
    //set T_fluid_cgs cap (10^12 K ~ 100 MeV)
    if (T_fluid_cgs > 1.0e11* T_fluid_cgs_floor){
      T_fluid_cgs = 1.0e11*T_fluid_cgs_floor;
    }
  }


    //  T_fluid = T_fluid_cgs * kb_cgs/(SQR(c_cgs) * m_n_cgs);                                                                             
    T_fluid = T_fluid_cgs * kb_cgs * G_cgs / fasterpow_prim(c_cgs,4.0) * fasterpow_prim(c2k, 1.0);
    double Fac2 = 1.0;
    double P_nucl_cgs = (5.0/3.0-1.0)*(1.5*n_nucleon*kb_cgs*T_fluid_cgs*Fac2);
    double P_nucl = fasterpow_prim(c_cgs, -4.0) * fasterpow_prim(G_cgs, 1.0) * P_nucl_cgs * fasterpow_prim(c2k, -2.0);
    double P_rad_ph_cgs = (11.0/12.0)*rad_const_cgs * fasterpow_prim(T_fluid_cgs, 4.0);
    double P_rad_ph = fasterpow_prim(c_cgs, -4.0) * fasterpow_prim(G_cgs, 1.0) * P_rad_ph_cgs * fasterpow_prim(c2k, -2.0);
    
    double P_thermal =  P_rad_ph + P_nucl;
    P = P_cold + P_thermal;

    if (T_fluid > 1.0e-30){
      printf("Inside compute_P_T_microphys_w, T_fluid is too large, is %.16e, eps_cold=%.16e, P=%.16e, P_cold=%.16e, P_thermal=%.16e,T_fluid_cgs=%.16e, n_nucleon=%.16e \n", T_fluid, eps_cold, P, P_cold, P_thermal,T_fluid_cgs, n_nucleon);
      printf("w_lhs_cgs = %.16e, rho_b_cgs = %.16e, rho_b = %.16e \n", w_lhs_cgs, rho_b_cgs, rho_b);
    }
   
}



double find_T_w(double &T, double &w_lhs_cgs, double &rad_const, double& n_nucleon, double& N_nu, double &T_scal){
  /*  double gamma_fermi = 5.0/3.0;                                                                                                     
  double Fac = min_val(1.0,  2.0*fasterpow_prim((3.0*M_PI*M_PI), 1.0/3.0)*m_n_cgs/(18.0*hbar_cgs*hbar_cgs)*fasterpow_prim(n_nucleon, -2.0/3.0)*kb_cgs*T );     
  */
  // w_lhs = w - rho_b*eps_cold - P_cold = eng_th + P_th = (11/4 aT^4 + 3/2*nkT) + (nkT + 11/12aT^4) = 11/3 aT^4 + 5/2*nkT                                                                              
  double Fac=1.0;
  double diff =  11.0/3.0 * rad_const * fasterpow_prim(T*T_scal, 4.0) + 2.5*n_nucleon*T*T_scal*kb_cgs*Fac - w_lhs_cgs;
  return  diff/w_lhs_cgs;
}

void solve_T_w(double &Ta, double &Tb,
	       double &T, double &w_lhs_cgs, double &rad_const, double& n_nucleon, double& N_nu, double &T_scal){
  double fa, fb, fc;
  double gamma_fermi = 5.0/3.0;
  int bisec = 1;
  double T_init = T;
  double tiny = 1.0;

  //If T_cgs_guess is smaller than 1K, set T_cgs to 1K, don't do the bisection!    
  if (T_init < tiny){
    T = tiny;
    bisec = 0;
  }
  else{
    Ta /=T_scal;
    Tb /=T_scal;

    fa = find_T_w(Ta, w_lhs_cgs, rad_const, n_nucleon, N_nu, T_scal);
    fb = find_T_w(Tb, w_lhs_cgs, rad_const, n_nucleon, N_nu, T_scal);
  
  int iter=0;

  while (fa * fb >= 0)
    {
      if (iter > 20) { //We fail to converge to find a proper Ta -- Tb range!!!! 
        T = T_init;
        bisec = 0;
	printf("break loop Solve_T_w: Ta = %e, Tb =%e, Ta_init= %e, Tb_init=%e, w_lhs_cgs = %e, fa=%e, fb=%e \n", Ta, Tb, T_init*0.8, T_init*1.05, w_lhs_cgs, fa, fb);
        break;
      }
      /*                                                                                                                                                            
      double Fac_a_w = min_val(1.0,  2.0*fasterpow_prim((3.0*M_PI*M_PI), 1.0/3.0)*m_n_cgs/(18.0*hbar_cgs*hbar_cgs)*fasterpow_prim(n_nucleon, -2.0/3.0)*kb_cgs*Ta );       
      double RHS_a_w = E_rad/3.0 + 3.0 * rad_const * fasterpow_prim(Ta, 4.0) + gamma_fermi*1.5*n_nucleon*T*kb_cgs*Fac_a_w;                                                  
      double Fac_b_w = min_val(1.0,  2.0*fasterpow_prim((3.0*M_PI*M_PI), 1.0/3.0)*m_n_cgs/(18.0*hbar_cgs*hbar_cgs)*fasterpow_prim(n_nucleon, -2.0/3.0)*kb_cgs*Tb );              
      double RHS_b_w = E_rad/3.0 + 3.0 * rad_const * fasterpow_prim(Tb, 4.0) + gamma_fermi*1.5*n_nucleon*T*kb_cgs*Fac_b_w;                                                       
      printf("You have not assumed right a and b, Ta=%e, Tb=%e, fa=%e, fb=%e, iter=%d \n", Ta, Tb, fa, fb, iter);                                                                
      printf("w_lhs_cgs=%e, RHS_a_w=%e, RHS_b_w=%e \n", w_lhs_cgs, RHS_a_w, RHS_b_w);                                                                                            
      printf("Fac_a_w=%e, rad_const=%e, fasterpow_prim(Ta, 4.0)=%e, n_nucleon=%e\n ", Fac_a_w, rad_const, fasterpow_prim(Ta, 4.0), n_nucleon);                     
      */
     // Expamd the gussing range                                                           
                                                                                                 
      Ta *= 0.5;
      Tb *= 1.5;
      fa = find_T_w(Ta, w_lhs_cgs, rad_const, n_nucleon, N_nu,  T_scal);
      fb = find_T_w(Tb, w_lhs_cgs, rad_const, n_nucleon, N_nu,  T_scal);
      iter +=1;
    }

  if (bisec){
    T = Ta;
    int iter = 0;
    while (fabs(Tb-Ta)/(Tb+Ta) >= 1.0e-12)
      {
        // Find middle point                                                                                                                        
        T = (Ta+Tb)/2;
        fc = find_T_w(T, w_lhs_cgs, rad_const, n_nucleon, N_nu, T_scal);
        // Check if middle point is root                                                                                                                        
        if (fc == 0.0)
          break;
        // Decide the side to repeat the steps                                                                                                  
        else if (fc*fa < 0)
          Tb = T;
        else
          Ta = T;
	
	iter +=1;
	
	if (iter > 50){
	  // Failed to converge after 50 iteration, do 50 more!!!!!!!                                                                                      
	  printf("inside solve_T_w, over iteration, fabs(Tb-Ta)/(Tb+Ta)=%e \n", fabs(Tb-Ta)/(Tb+Ta));                                                                                                            
	  int iter2 = 0;
	  while (fabs(Tb-Ta)/(Tb+Ta) >= 1.0e-12)
	    {
	      // Find middle point                                                                                                  
	      T = (Ta+Tb)/2;
	      fc = find_T_w(T, w_lhs_cgs, rad_const, n_nucleon, N_nu, T_scal);
	      // Check if middle point is root                                                                                                    
	      if (fc == 0.0)
		break;
	      // Decide the side to repeat the steps                                                                                     
	      else if (fc*fa < 0)
		Tb = T;
	      else
		Ta = T;

	      iter2 +=1;
	    }
	  if (iter2 > 50){
	    printf("inside solve_T_w, over iteration AGAIN, fabs(Tb-Ta)/(Tb+Ta)=%e \n", fabs(Tb-Ta)/(Tb+Ta));                                   
	    break;
	  }	 	  
	}
      } 
   
    T *= T_scal;
    Ta *= T_scal;
    Tb *=T_scal;

  } // End of bisection
  } // End of T > T_tiny 

}



extern "C" void CCTK_FCALL CCTK_FNAME(compute_eps_T_microphys)
  (double &eps, double &T_fluid, double &eps_cold,  double &P, double &P_cold, double &rho_b);


void compute_eps_T_microphys (double &eps, double &T_fluid, double &eps_cold,  double &P, double &P_cold, double &rho_b){

  double P_thermal = P - P_cold;

  /*
  if (P_thermal < 1.0e-19){
    P_thermal = fasterpow_prim(rho_b, 1.357)*6.5;
    //Need to recompute P if P_thermal is revalued.                                                                           
    P = P_thermal + P_cold;
  }
  */

  double P_thermal_cgs = P_thermal * fasterpow_prim(c_cgs, 4.0) * fasterpow_prim(G_cgs, -1.0) * fasterpow_prim(k2c, -2.0);
  double rho_b_cgs = fasterpow_prim(c_cgs,2.0) / G_cgs * rho_b * fasterpow_prim(k2c, -2.0);
  double n_nucleon = rho_b_cgs/m_n_cgs;
  double N_nu = 3.0;
  double T_fluid_cgs;
  double T_scal= 1.0e0;

  double T_fluid_cgs_floor = 1.0e1;

  //check if the thermal part if positive:                                                         
  
  if (P_thermal <= 0.0){
    T_fluid_cgs = T_fluid_cgs_floor;
    //    T_fluid = T_fluid_cgs_floor * kb_cgs * G_cgs / fasterpow_prim(c_cgs,4.0) * fasterpow_prim(c2k, 1.0);
    //    eps = eps_cold;
  }
  else{
    double threshold = fasterpow_prim(P_thermal_cgs,0.75)*0.978482042633557*fasterpow_prim(rad_const_cgs,0.25)/kb_cgs;
    double T_fluid_cgs_guess;

    if (n_nucleon < threshold){
      T_fluid_cgs = fasterpow_prim(12.0*P_thermal_cgs/rad_const_cgs/11.0, 0.25);
      T_fluid_cgs_guess = fasterpow_prim(12.0*P_thermal_cgs/rad_const_cgs/11.0, 0.25);
    }
    else{
      T_fluid_cgs = P_thermal_cgs/(1.0 * kb_cgs * n_nucleon);
      T_fluid_cgs_guess= P_thermal_cgs/(1.0 * kb_cgs * n_nucleon);
    }

    double Ta = T_fluid_cgs*0.9;
    double Tb = T_fluid_cgs*1.1;

    double rad_constl = rad_const_cgs;

    solve_T_P(Ta, Tb, T_fluid_cgs, P_thermal_cgs, rad_constl, n_nucleon, N_nu, T_scal);

    // WE DO A CHECK AFTER GETTING THE TEMPERATURE!!!!                                                                                               
    double E1 = (11.0/12.0 * rad_const_cgs)* fasterpow_prim(T_fluid_cgs, 4.0);
    double E2 = 1.0*n_nucleon*T_fluid_cgs*kb_cgs;
    double E_ratio = (E1 + E2)/P_thermal_cgs;
    if (T_fluid_cgs > T_fluid_cgs_floor && E_ratio > 1.05){
      printf("after solve_T_P, E_ratio too large, is %.16e, E1=%.16e, E2=%.16e, P_thermal_cgs=%.16e, T_fluid_cgs=%.16e, T_fluid_cgs_guess=%.16e, n_nucleon=%.16e \n", E_ratio, E1, E2, P_thermal_cgs, T_fluid_cgs, T_fluid_cgs_guess, n_nucleon);
    }

    //set T_fluid_cgs floor 
    if (T_fluid_cgs < T_fluid_cgs_floor){
      T_fluid_cgs = T_fluid_cgs_floor;
    }

    //set T_fluid_cgs cap (10^12 K ~ 100 MeV) 
    if (T_fluid_cgs > 1.0e11* T_fluid_cgs_floor){
      T_fluid_cgs = 1.0e11*T_fluid_cgs_floor;
    }
  }

    T_fluid = T_fluid_cgs * kb_cgs * G_cgs / fasterpow_prim(c_cgs,4.0) * fasterpow_prim(c2k, 1.0);
    double eng_nucl_cgs = 1.5*n_nucleon*kb_cgs*T_fluid_cgs;
    double eng_rad_ph_cgs = (11.0/4.0)*rad_const_cgs * fasterpow_prim(T_fluid_cgs, 4.0);
    double eng_thermal_cgs = eng_nucl_cgs + eng_rad_ph_cgs;

    double eng_thermal = eng_thermal_cgs * G_cgs * fasterpow_prim(c_cgs,-4.0) * fasterpow_prim(c2k, -2.0);
    double eps_thermal = eng_thermal/rho_b;
    eps = eps_thermal + eps_cold;

    if (T_fluid > 1.0e-30){
      printf("Inside compute_eps_T_microphys, T_fluid is too large, is %.16e, epsd=%.16e, eps_cold=%.16e, eps_thermal=%.16e, P=%.16e, P_cold=%.16e, P_thermal=%.16e, P_thermal_cgs=%.16e, T_fluid_cgs=%.16e, n_nucleon=%.16e \n", T_fluid, eps, eps_cold, eps_thermal, P, P_cold, P_thermal, P_thermal_cgs, T_fluid_cgs, n_nucleon);
      printf("eng_thermal_cgs = %.16e, rho_b_cgs = %.16e, rho_b = %.16e \n", eng_thermal_cgs, rho_b_cgs, rho_b);
    }
}

extern "C" void CCTK_FCALL CCTK_FNAME(compute_eps_T_microphys)
  (double &eps, double &T_fluid, double &eps_cold, double &P, double &P_cold, double &rho_b) {
  compute_eps_T_microphys(eps,T_fluid,eps_cold,P,P_cold,rho_b);
}



double find_T_P(double &T, double &P_thermal_cgs, double &rad_const, double& n_nucleon, double& N_nu, double& T_scal){
  double Fac = 1.0;
  double diff = (11.0/12.0 * rad_const)* fasterpow_prim(T*T_scal, 4.0) + 1.0*n_nucleon*T*T_scal*kb_cgs*Fac - P_thermal_cgs;
  return diff/P_thermal_cgs;
}




void solve_T_P(double &Ta, double &Tb,
             double &T, double &P_thermal_cgs, double &rad_const, double& n_nucleon, double& N_nu, double& T_scal)
{
  double fa, fb, fc;
  int bisec = 1;
  double T_init =T;
  double tiny = 1.0;

  //If T_cgs_guess is smaller than 1K, set T_cgs to 1K, don't do the bisection!                                                                                             
  if (T_init < tiny){
    T = tiny;
    bisec = 0;
  }
  else{
    Ta /=T_scal;
    Tb /=T_scal;

    fa = find_T_P(Ta, P_thermal_cgs, rad_const, n_nucleon, N_nu, T_scal);
    fb = find_T_P(Tb, P_thermal_cgs, rad_const, n_nucleon, N_nu, T_scal);

    int iter=0;

    while (fa *fb >= 0)
      {
        if (iter > 20) {//We fail to converge to find a proper Ta -- Tb range!!!!                                                                      
          T = T_init;
          bisec = 0;
          printf("break loop solve_T_P: Ta = %e, Tb =%e, Ta_init= %e, Tb_init=%e, P_thermal_cgs = %e, fa=%e, fb=%e \n", Ta, Tb, T_init*0.8, T_init*1.05, P_thermal_cgs, fa, fb);
          break;
        }
	Ta *= 0.5;
        Tb *= 1.5;
        fa = find_T_P(Ta, P_thermal_cgs, rad_const, n_nucleon, N_nu, T_scal);
        fb = find_T_P(Tb, P_thermal_cgs, rad_const, n_nucleon, N_nu, T_scal);
        iter += 1;
        //    printf ("iter = %d \n", iter);                                                                          
      }

    if (bisec){
      T = Ta;
      int iter = 0;
      while (fabs(Tb-Ta)/(Tb+Ta) >= 1.0e-12)
        {
	  T = (Ta+Tb)/2;
          fc = find_T_P(T, P_thermal_cgs, rad_const, n_nucleon, N_nu, T_scal);
	  if (fc == 0.0)
            break;
	  else if (fc*fa < 0)
            Tb = T;
          else
            Ta = T;

          iter +=1;

	  if (iter > 50){
	    int iter2 = 0;
            while (fabs(Tb-Ta)/(Tb+Ta) >= 1.0e-12)
              {
		T = (Ta+Tb)/2;
                fc = find_T_P(T, P_thermal_cgs, rad_const, n_nucleon, N_nu, T_scal);
		if (fc == 0.0)
                  break;
		else if (fc*fa < 0)
                  Tb = T;
                else
                  Ta = T;

                iter2 +=1;
              }
            if (iter2 > 50){
              printf("inside solve_T, over iteration AGAIN, fabs(Tb-Ta)/(Tb+Ta)=%e \n", fabs(Tb-Ta)/(Tb+Ta));
              break;
            }
          }
        }

      T *=T_scal;
      Ta *=T_scal;
      Tb *=T_scal;
    }
  }
}
