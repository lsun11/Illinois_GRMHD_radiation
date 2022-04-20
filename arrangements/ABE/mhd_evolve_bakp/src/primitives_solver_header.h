void compute_pcold_epscold_cpp(double &rhob, double &P_cold, double &eps_cold,
				       int &neos, int &ergo_star, double &ergo_sigma, double *rho_tab,double *P_tab,double *eps_tab,
				       double *k_tab,double *gamma_tab);

double fasterpow_prim(double inputvar,double inputpow);


void compute_T_fluid_cpp(double &X, double &Y, double &Z, int&num_CO, 
			 double &xNS1, double &yNS1, double &xNS2, double &yNS2, double &rhob, double &rho_star, double &P,
			 double &M_B, double &rad_T_fac, double &rad_T_cutoff, double &rad_T_pow, double &rad_T_floor, 
			 double &rho_b_atm, double &T_fluid);


void compute_T_fluid_shock_cpp(double &rho_star, double &P,double &M_B, double &rad_T_floor,double &T_fluid);
