#define c_cgs 29979245800.0 // cm s^-1                                                                                                    
#define G_cgs 6.6743015e-8 // erg cm g^-2                                                                                                                                
#define hbar_cgs 1.05457266e-27 // erg s                                                                                                                                   
#define kb_cgs 1.380649e-16 //erg K^-1                                                                                                                    
#define m_n_cgs 1.674927485e-24 // g                                                          

#define m_e_cgs 9.10938e-28 //g
                            
#define rad_const_cgs 7.5646e-15 // erg cm^-3 K^-4     
#define k0 1.365    // (1+3*alpha^2)/4 * C_v                                                              
#define k1 0.3671875 // (1+5*alpha^2)/24 
#define k2 0.3257875 // (4(C_v-1)^2 + 5*alpha^2)/24  
                                                                                
#define M_sun_geo 1.0
#define sigma_0 1.76e-44 //cm^2
#define k2c 100000.0
#define c2k 0.00001
#define eightpi_o_hc3 3.2017829478627814e-51 // erg^-3 cm^-3
#define sigma_cgs 5.670373e-5 // erg s^-1 cm^-2 K^-4

void compute_pcold_epscold_cpp(double &rhob, double &P_cold, double &eps_cold,
				       int &neos, int &ergo_star, double &ergo_sigma, double *rho_tab,double *P_tab,double *eps_tab,
			       double *k_tab,double *gamma_tab, int &enable_OS_collapse);

double fasterpow_prim(double inputvar,double inputpow);

void compute_T_fluid_cpp(double &X, double &Y, double &Z, int&num_CO, 
			 double &xNS1, double &yNS1, double &xNS2, double &yNS2, double &rhob, double &rho_star, double &P,
			 double &M_B, double &rad_T_fac, double &rad_T_cutoff, double &rad_T_pow, double &rad_T_floor, 
			 double &rho_b_atm, double &T_fluid);

void compute_T_fluid_shock_cpp2(double &P, double &rho_b, double &T_fluid);
void compute_T_fluid_OS(double &E_rad, double &rad_const, double &T_fluid);

void flux_hll_cpp_frad(int *ext,double &Ur,double &Ul,double &Fr,double &Fl,double &F,double &cmax,double &cmin);



// opacities, emissivities computation
void compute_M1(double &Pij,double &Fi,double &Fj,double &Fasq,double &E,double &gupij, double &shifti, double &shiftj, double &lapse, double &ui, double &uj, double &chi, double &psim4, double &Erad_atm_cut);


//T_fluid finders
double find_T(double &T, double &eng_thermal_cgs, double &rad_const, double& n_nucleon, double& N_nu, double& T_scal);
void solve_T(double &Ta, double &Tb, double &T, double &eng_thermal_cgs, double &rad_const, double& n_nucleon, double& N_nu, double& T_scal);
double min_val(double val1,double val2);
void compute_P_T_microphys(double &P, double &T_fluid, double &P_cold,  double &eps, double &eps_cold, double &rho_b);
void compute_dP_nucldrho_microphys(double &dP_nucldrho, double &T_fluid, double &rho_b);

double find_T_w(double &T, double &w_lhs_cgs, double &rad_const, double& n_nucleon, double& N_nu, double &T_scal);
void solve_T_w(double &Ta, double &Tb, double &T, double &w_lhs_cgs, double &rad_const, double& n_nucleon, double& N_nu, double &T_scal);
void compute_P_T_microphys_w(double &P, double &w, double &T_fluid, double &P_cold, double &eps_cold, double &rho_b);
void compute_opacity_emissivity(double &rho_star, double &al, double &Psi6, double &u0, double &kappa_a, double &kappa_s, double &eta_gf, double &kappa_a_nue, double &kappa_s_nue, double &eta_gf_nue, double &kappa_a_nux, double &kappa_s_nux, double &eta_gf_nux, double &T_fluid, double &chi_rad, double &chi_rad_nue, double &Y_e, double &optd, double &eta_nue, double &T_fluid_cgs_atm, int &microphysics_scheme, int &rad_fix);

void compute_eps_T_microphys (double &eps, double &T_fluid, double &eps_cold,  double &P, double &P_cold, double &rho_b);
double find_T_P(double &T, double &P_thermal_cgs, double &rad_const, double& n_nucleon, double& N_nu, double& T_scal);
void solve_T_P(double &Ta, double &Tb, double &T, double &P_thermal_cgs, double &rad_const, double& n_nucleon, double& N_nu, double& T_scal);

void compute_eps_th(double &T_fluid, double &rho_b, double &eps_thermal);

//Fermi_functions
double Fermi3p(double &eta);
double Fermi3m(double &eta);
double Fermi3(double &eta);

double Fermi4p(double &eta);
double Fermi4m(double &eta);
double Fermi4(double &eta);

double Fermi5p(double &eta);
double Fermi5m(double &eta);
double Fermi5(double &eta);
