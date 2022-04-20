#define DECLARE_CCTK_PARAMETERS \
CCTK_REAL  excision_radius&&\
CCTK_REAL  run_time&&\
CCTK_INT Symmetry&&\
CCTK_INT bssn_enable&&\
CCTK_INT cowling_enable&&\
CCTK_INT excision_enable&&\
CCTK_INT fisheye_enable&&\
CCTK_INT iter_count&&\
CCTK_INT number_of_mol_ministeps&&\
CCTK_INT rot_metric&&\
CCTK_INT trA_detg_enforce&&\
COMMON /cctk_params_global/excision_radius,run_time,Symmetry,bssn_enable,cowling_enable,excision_enable,fisheye_enable,iter_count,number_of_mol_ministeps,rot_metric,trA_detg_enforce&&\
CCTK_REAL  CoMExcision_Radius&&\
CCTK_REAL  M_ADM&&\
CCTK_REAL  initial_monopole_value&&\
CCTK_REAL  inner_lum_rad_ratio&&\
CCTK_REAL  lum_outer_rad&&\
CCTK_REAL  radius_esc1&&\
CCTK_REAL  radius_esc2&&\
CCTK_REAL  radius_esc3&&\
CCTK_REAL  radius_esc4&&\
CCTK_REAL  rhob_cutoff&&\
CCTK_REAL  rhosurf&&\
CCTK_REAL  rhosurf_rmax&&\
CCTK_REAL  rhosurf_rmin&&\
CCTK_REAL  u0sch&&\
CCTK_INT N_rad_ray&&\
CCTK_INT Nphi_points&&\
CCTK_INT Ntheta_points&&\
CCTK_INT const_rad_surf_diagnostics&&\
CCTK_INT drho_dtau_calc_enable&&\
CCTK_INT escape_mass_diag&&\
CCTK_INT luminosity_diagnostics&&\
COMMON /diagnostics_mhdrest/CoMExcision_Radius,M_ADM,initial_monopole_value,inner_lum_rad_ratio,lum_outer_rad,radius_esc1,radius_esc2,radius_esc3,radius_esc4,rhob_cutoff,rhosurf,rhosurf_rmax,rhosurf_rmin,u0sch,N_rad_ray,Nphi_points,Ntheta_points,const_rad_surf_diagnostics,drho_dtau_calc_enable,escape_mass_diag,luminosity_diagnostics&&\
CCTK_REAL  CCTKH5&&\
CCTK_REAL  CCTKH6&&\
CCTK_REAL  CCTKH7&&\
CCTK_REAL  CCTKH8&&\
CCTK_REAL  binary_orb_freq&&\
CCTK_REAL  CCTKH11&&\
CCTK_REAL  CCTKH12&&\
CCTK_REAL  CCTKH13&&\
CCTK_REAL  CCTKH14&&\
CCTK_REAL  CCTKH15&&\
CCTK_REAL  CCTKH16&&\
CCTK_REAL  rah1_initial&&\
CCTK_REAL  rah2_initial&&\
CCTK_REAL  CCTKH17&&\
CCTK_REAL  xbh1_initial&&\
CCTK_REAL  xbh2_initial&&\
CCTK_REAL  CCTKH19&&\
CCTK_INT CCTKH0&&\
CCTK_INT CCTKH1&&\
CCTK_INT CCTKH2&&\
CCTK_INT CCTKH3&&\
CCTK_INT CCTKH4&&\
CCTK_INT CCTKH9&&\
CCTK_INT CCTKH10&&\
CCTK_INT CCTKH18&&\
CCTK_INT CCTKH20&&\
COMMON /BSSNrest/CCTKH5,CCTKH6,CCTKH7,CCTKH8,binary_orb_freq,CCTKH11,CCTKH12,CCTKH13,CCTKH14,CCTKH15,CCTKH16,rah1_initial,rah2_initial,CCTKH17,xbh1_initial,xbh2_initial,CCTKH19,CCTKH0,CCTKH1,CCTKH2,CCTKH3,CCTKH4,CCTKH9,CCTKH10,CCTKH18,CCTKH20&&\
CCTK_REAL  CCTKH21&&\
CCTK_REAL  dcostheta&&\
CCTK_REAL  ddrbddr&&\
CCTK_REAL  dphi&&\
CCTK_REAL  drbdr&&\
CCTK_REAL  CCTKH26&&\
CCTK_REAL  rbr&&\
CCTK_REAL  CCTKH30&&\
CCTK_REAL  surf_radius&&\
CCTK_INT Compute_VolIntegrands_Every&&\
CCTK_INT N_phi&&\
CCTK_INT N_theta&&\
CCTK_INT WhichIntegral&&\
CCTK_INT CCTKH22&&\
CCTK_INT CCTKH23&&\
CCTK_INT CCTKH24&&\
CCTK_INT CCTKH25&&\
CCTK_INT CCTKH27&&\
CCTK_INT ntot&&\
CCTK_INT num_BHs&&\
CCTK_INT CCTKH28&&\
CCTK_INT CCTKH29&&\
CCTK_INT CCTKH31&&\
CCTK_INT sym_factor&&\
COMMON /DIAGNOSTICS_VACUUMrest/CCTKH21,dcostheta,ddrbddr,dphi,drbdr,CCTKH26,rbr,CCTKH30,surf_radius,Compute_VolIntegrands_Every,N_phi,N_theta,WhichIntegral,CCTKH22,CCTKH23,CCTKH24,CCTKH25,CCTKH27,ntot,num_BHs,CCTKH28,CCTKH29,CCTKH31,sym_factor&&\
CCTK_REAL  CCTKH32&&\
CCTK_REAL  CCTKH33&&\
CCTK_REAL  CCTKH34&&\
CCTK_REAL  CCTKH35&&\
CCTK_REAL  CCTKH45&&\
CCTK_REAL  CCTKH46&&\
CCTK_REAL  CCTKH47&&\
CCTK_REAL  CCTKH48&&\
CCTK_REAL  CCTKH49&&\
CCTK_REAL  CCTKH50&&\
CCTK_REAL  CCTKH51&&\
CCTK_REAL  CCTKH52&&\
CCTK_STRING  bitant_plane&&\
CCTK_STRING  domain&&\
CCTK_STRING  CCTKH36&&\
CCTK_STRING  CCTKH37&&\
CCTK_STRING  CCTKH44&&\
CCTK_INT CCTKH38&&\
CCTK_INT CCTKH39&&\
CCTK_INT CCTKH40&&\
CCTK_INT CCTKH41&&\
CCTK_INT CCTKH42&&\
CCTK_INT CCTKH43&&\
COMMON /GRIDrest/CCTKH32,CCTKH33,CCTKH34,CCTKH35,CCTKH45,CCTKH46,CCTKH47,CCTKH48,CCTKH49,CCTKH50,CCTKH51,CCTKH52,bitant_plane,domain,CCTKH36,CCTKH37,CCTKH44,CCTKH38,CCTKH39,CCTKH40,CCTKH41,CCTKH42,CCTKH43&&\
CCTK_REAL  CCTKH55&&\
CCTK_REAL  CCTKH56&&\
CCTK_REAL  CCTKH57&&\
CCTK_STRING  slicing_type&&\
CCTK_INT CCTKH53&&\
CCTK_INT CCTKH54&&\
CCTK_INT CCTKH58&&\
CCTK_INT CCTKH59&&\
CCTK_INT CCTKH60&&\
CCTK_INT CCTKH61&&\
CCTK_INT CCTKH62&&\
COMMON /LAPSErest/CCTKH55,CCTKH56,CCTKH57,slicing_type,CCTKH53,CCTKH54,CCTKH58,CCTKH59,CCTKH60,CCTKH61,CCTKH62&&\
CCTK_REAL  CCTKH65&&\
CCTK_REAL  CCTKH66&&\
CCTK_REAL  K_poly&&\
CCTK_REAL  CCTKH69&&\
CCTK_REAL  CCTKH71&&\
CCTK_REAL  P_fact&&\
CCTK_REAL  CCTKH72&&\
CCTK_REAL  CCTKH73&&\
CCTK_REAL  Sym_Bz&&\
CCTK_REAL  CCTKH75&&\
CCTK_REAL  CCTKH79&&\
CCTK_REAL  CCTKH80&&\
CCTK_REAL  CCTKH83&&\
CCTK_REAL  CCTKH90&&\
CCTK_REAL  ergo_sigma&&\
CCTK_REAL  CCTKH92&&\
CCTK_REAL  gamma_th&&\
CCTK_REAL  CCTKH94&&\
CCTK_REAL  CCTKH96&&\
CCTK_REAL  CCTKH97&&\
CCTK_REAL  CCTKH98&&\
CCTK_REAL  CCTKH99&&\
CCTK_REAL  CCTKH101&&\
CCTK_REAL  pfloor&&\
CCTK_REAL  CCTKH105&&\
CCTK_REAL  CCTKH106&&\
CCTK_REAL  CCTKH107&&\
CCTK_REAL  CCTKH108&&\
CCTK_REAL  CCTKH109&&\
CCTK_REAL  CCTKH111&&\
CCTK_REAL  CCTKH112&&\
CCTK_REAL  CCTKH115&&\
CCTK_REAL  CCTKH116&&\
CCTK_REAL  rho_b_atm&&\
CCTK_REAL  rho_b_max&&\
CCTK_REAL  rho_fact&&\
CCTK_REAL  CCTKH119&&\
CCTK_REAL  CCTKH120&&\
CCTK_REAL  CCTKH121&&\
CCTK_REAL  tau_atm&&\
CCTK_REAL  tau_fact&&\
CCTK_INT CCTKH63&&\
CCTK_INT CCTKH64&&\
CCTK_INT CCTKH67&&\
CCTK_INT CCTKH68&&\
CCTK_INT CCTKH70&&\
CCTK_INT CCTKH74&&\
CCTK_INT CCTKH76&&\
CCTK_INT CCTKH77&&\
CCTK_INT CCTKH78&&\
CCTK_INT CCTKH81&&\
CCTK_INT constrained_transport_scheme&&\
CCTK_INT CCTKH82&&\
CCTK_INT em_evolve_enable&&\
CCTK_INT CCTKH84&&\
CCTK_INT CCTKH85&&\
CCTK_INT enable_OS_collapse&&\
CCTK_INT CCTKH86&&\
CCTK_INT CCTKH87&&\
CCTK_INT CCTKH88&&\
CCTK_INT CCTKH89&&\
CCTK_INT eps_flag&&\
CCTK_INT ergo_star&&\
CCTK_INT CCTKH91&&\
CCTK_INT CCTKH93&&\
CCTK_INT hyperbolic_divergence_cleaning_centered_differencing&&\
CCTK_INT hyperbolic_divergence_cleaning_enable&&\
CCTK_INT CCTKH95&&\
CCTK_INT CCTKH100&&\
CCTK_INT neos&&\
CCTK_INT CCTKH102&&\
CCTK_INT CCTKH103&&\
CCTK_INT CCTKH104&&\
CCTK_INT CCTKH110&&\
CCTK_INT rad_evolve_enable&&\
CCTK_INT CCTKH113&&\
CCTK_INT CCTKH114&&\
CCTK_INT CCTKH117&&\
CCTK_INT CCTKH118&&\
CCTK_INT CCTKH122&&\
CCTK_INT CCTKH123&&\
CCTK_INT CCTKH124&&\
CCTK_INT CCTKH125&&\
COMMON /MHD_EVOLVErest/CCTKH65,CCTKH66,K_poly,CCTKH69,CCTKH71,P_fact,CCTKH72,CCTKH73,Sym_Bz,CCTKH75,CCTKH79,CCTKH80,CCTKH83,CCTKH90,ergo_sigma,CCTKH92,gamma_th,CCTKH94,CCTKH96,CCTKH97,CCTKH98,CCTKH99,CCTKH101,pfloor,CCTKH105,CCTKH106,CCTKH107,CCTKH108,CCTKH109,CCTKH111,CCTKH112,CCTKH115,CCTKH116,rho_b_atm,rho_b_max,rho_fact,CCTKH119,CCTKH120,CCTKH121,tau_atm,tau_fact,CCTKH63,CCTKH64,CCTKH67,CCTKH68,CCTKH70,CCTKH74,CCTKH76,CCTKH77,CCTKH78,CCTKH81,constrained_transport_scheme,CCTKH82,em_evolve_enable,CCTKH84,CCTKH85,enable_OS_collapse,CCTKH86,CCTKH87,CCTKH88,CCTKH89,eps_flag,ergo_star,CCTKH91,CCTKH93,hyperbolic_divergence_cleaning_centered_differencing,hyperbolic_divergence_cleaning_enable,CCTKH95,CCTKH100,neos,CCTKH102,CCTKH103,CCTKH104,CCTKH110,rad_evolve_enable,CCTKH113,CCTKH114,CCTKH117,CCTKH118,CCTKH122,CCTKH123,CCTKH124,CCTKH125&&\
CCTK_INT Spatial_Gauge&&\
CCTK_INT CCTKH126&&\
CCTK_INT CCTKH127&&\
COMMON /SHIFTrest/Spatial_Gauge,CCTKH126,CCTKH127&&\

