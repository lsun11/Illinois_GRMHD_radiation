#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL Erad_atm_cut;
  CCTK_REAL Erad_cut;
  CCTK_REAL K_poly;
  CCTK_REAL M_B;
  CCTK_REAL OS_T_Rmax;
  CCTK_REAL P_fact;
  CCTK_REAL P_max;
  CCTK_REAL Psi6threshold;
  CCTK_REAL Sym_Bz;
  CCTK_REAL T_fluid_cgs_atm;
  CCTK_REAL c_h_default;
  CCTK_REAL c_p_default;
  CCTK_REAL dxvec;
  CCTK_REAL eps_thermal_bhns;
  CCTK_REAL ergo_sigma;
  CCTK_REAL gamma_OS;
  CCTK_REAL gamma_th;
  CCTK_REAL lambda_line;
  CCTK_REAL mhd_bigp;
  CCTK_REAL mhd_mbh;
  CCTK_REAL mhd_mdot;
  CCTK_REAL mhd_r_crit;
  CCTK_REAL min_BH_radius;
  CCTK_REAL pfloor;
  CCTK_REAL r_ex_line;
  CCTK_REAL rad_T_cutoff;
  CCTK_REAL rad_T_fac;
  CCTK_REAL rad_T_floor;
  CCTK_REAL rad_T_pow;
  CCTK_REAL rad_emissivity_abs;
  CCTK_REAL rad_emissivity_sct;
  CCTK_REAL rad_opacity_abs;
  CCTK_REAL rad_opacity_sct;
  CCTK_REAL rho_b_atm;
  CCTK_REAL rho_b_max;
  CCTK_REAL rho_fact;
  CCTK_REAL rhos_max;
  CCTK_REAL sdots_o_rhot;
  CCTK_REAL t_cool;
  CCTK_REAL tau_atm;
  CCTK_REAL tau_fact;
  CCTK_INT A_BC_rad_mag_bondi_enable;
  CCTK_INT EM_BC;
  CCTK_INT MHD_MaxNumConstrainedVars;
  CCTK_INT MHD_MaxNumEvolvedVars;
  CCTK_INT Matter_BC;
  CCTK_INT Reconstruction;
  CCTK_INT adm_ppm_b;
  CCTK_INT allow_negative_eps_th;
  CCTK_INT artificial_cooling_enable;
  CCTK_INT compute_microphysics;
  CCTK_INT constrained_transport_scheme;
  CCTK_INT cooling_in_St_eq;
  CCTK_INT em_evolve_enable;
  CCTK_INT em_gauge;
  CCTK_INT enable_HARM_energyvariable;
  CCTK_INT enable_OS_collapse;
  CCTK_INT enable_disk_em_flux_induction;
  CCTK_INT enable_primitives_disk;
  CCTK_INT enable_shocktest_primitive_mode;
  CCTK_INT enable_trace_field_line;
  CCTK_INT eps_flag;
  CCTK_INT ergo_star;
  CCTK_INT force_font_fix_fail;
  CCTK_INT horizon_enforce_rho_profile;
  CCTK_INT hyperbolic_divergence_cleaning_centered_differencing;
  CCTK_INT hyperbolic_divergence_cleaning_enable;
  CCTK_INT m;
  CCTK_INT microphysics_scheme;
  CCTK_INT neos;
  CCTK_INT nrhovec;
  CCTK_INT primitives_debug;
  CCTK_INT primitives_solver;
  CCTK_INT rad_closure_scheme;
  CCTK_INT rad_evolve_enable;
  CCTK_INT rad_fix;
  CCTK_INT rad_fourforce_enable;
  CCTK_INT reconstruct_Bitildes_instead_of_Bis;
  CCTK_INT reconstruct_Pthermal_instead_of_P;
  CCTK_INT tau_stildefix_enable;
  CCTK_INT use_HARM_primitives;
  CCTK_INT use_central_scheme_instead_of_hll;
  CCTK_INT use_new_code;
} RESTRICTED_MHD_EVOLVE_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_RESTRICTED_MHD_EVOLVE_STRUCT_PARAMS \
  CCTK_REAL const Erad_atm_cut = RESTRICTED_MHD_EVOLVE_STRUCT.Erad_atm_cut; \
  CCTK_REAL const Erad_cut = RESTRICTED_MHD_EVOLVE_STRUCT.Erad_cut; \
  CCTK_REAL const K_poly = RESTRICTED_MHD_EVOLVE_STRUCT.K_poly; \
  CCTK_REAL const M_B = RESTRICTED_MHD_EVOLVE_STRUCT.M_B; \
  CCTK_REAL const OS_T_Rmax = RESTRICTED_MHD_EVOLVE_STRUCT.OS_T_Rmax; \
  CCTK_REAL const P_fact = RESTRICTED_MHD_EVOLVE_STRUCT.P_fact; \
  CCTK_REAL const P_max = RESTRICTED_MHD_EVOLVE_STRUCT.P_max; \
  CCTK_REAL const Psi6threshold = RESTRICTED_MHD_EVOLVE_STRUCT.Psi6threshold; \
  CCTK_REAL const Sym_Bz = RESTRICTED_MHD_EVOLVE_STRUCT.Sym_Bz; \
  CCTK_REAL const T_fluid_cgs_atm = RESTRICTED_MHD_EVOLVE_STRUCT.T_fluid_cgs_atm; \
  CCTK_REAL const c_h_default = RESTRICTED_MHD_EVOLVE_STRUCT.c_h_default; \
  CCTK_REAL const c_p_default = RESTRICTED_MHD_EVOLVE_STRUCT.c_p_default; \
  CCTK_REAL const dxvec = RESTRICTED_MHD_EVOLVE_STRUCT.dxvec; \
  CCTK_REAL const eps_thermal_bhns = RESTRICTED_MHD_EVOLVE_STRUCT.eps_thermal_bhns; \
  CCTK_REAL const ergo_sigma = RESTRICTED_MHD_EVOLVE_STRUCT.ergo_sigma; \
  CCTK_REAL const gamma_OS = RESTRICTED_MHD_EVOLVE_STRUCT.gamma_OS; \
  CCTK_REAL const gamma_th = RESTRICTED_MHD_EVOLVE_STRUCT.gamma_th; \
  CCTK_REAL const lambda_line = RESTRICTED_MHD_EVOLVE_STRUCT.lambda_line; \
  CCTK_REAL const mhd_bigp = RESTRICTED_MHD_EVOLVE_STRUCT.mhd_bigp; \
  CCTK_REAL const mhd_mbh = RESTRICTED_MHD_EVOLVE_STRUCT.mhd_mbh; \
  CCTK_REAL const mhd_mdot = RESTRICTED_MHD_EVOLVE_STRUCT.mhd_mdot; \
  CCTK_REAL const mhd_r_crit = RESTRICTED_MHD_EVOLVE_STRUCT.mhd_r_crit; \
  CCTK_REAL const min_BH_radius = RESTRICTED_MHD_EVOLVE_STRUCT.min_BH_radius; \
  CCTK_REAL const pfloor = RESTRICTED_MHD_EVOLVE_STRUCT.pfloor; \
  CCTK_REAL const r_ex_line = RESTRICTED_MHD_EVOLVE_STRUCT.r_ex_line; \
  CCTK_REAL const rad_T_cutoff = RESTRICTED_MHD_EVOLVE_STRUCT.rad_T_cutoff; \
  CCTK_REAL const rad_T_fac = RESTRICTED_MHD_EVOLVE_STRUCT.rad_T_fac; \
  CCTK_REAL const rad_T_floor = RESTRICTED_MHD_EVOLVE_STRUCT.rad_T_floor; \
  CCTK_REAL const rad_T_pow = RESTRICTED_MHD_EVOLVE_STRUCT.rad_T_pow; \
  CCTK_REAL const rad_emissivity_abs = RESTRICTED_MHD_EVOLVE_STRUCT.rad_emissivity_abs; \
  CCTK_REAL const rad_emissivity_sct = RESTRICTED_MHD_EVOLVE_STRUCT.rad_emissivity_sct; \
  CCTK_REAL const rad_opacity_abs = RESTRICTED_MHD_EVOLVE_STRUCT.rad_opacity_abs; \
  CCTK_REAL const rad_opacity_sct = RESTRICTED_MHD_EVOLVE_STRUCT.rad_opacity_sct; \
  CCTK_REAL const rho_b_atm = RESTRICTED_MHD_EVOLVE_STRUCT.rho_b_atm; \
  CCTK_REAL const rho_b_max = RESTRICTED_MHD_EVOLVE_STRUCT.rho_b_max; \
  CCTK_REAL const rho_fact = RESTRICTED_MHD_EVOLVE_STRUCT.rho_fact; \
  CCTK_REAL const rhos_max = RESTRICTED_MHD_EVOLVE_STRUCT.rhos_max; \
  CCTK_REAL const sdots_o_rhot = RESTRICTED_MHD_EVOLVE_STRUCT.sdots_o_rhot; \
  CCTK_REAL const t_cool = RESTRICTED_MHD_EVOLVE_STRUCT.t_cool; \
  CCTK_REAL const tau_atm = RESTRICTED_MHD_EVOLVE_STRUCT.tau_atm; \
  CCTK_REAL const tau_fact = RESTRICTED_MHD_EVOLVE_STRUCT.tau_fact; \
  CCTK_INT const A_BC_rad_mag_bondi_enable = RESTRICTED_MHD_EVOLVE_STRUCT.A_BC_rad_mag_bondi_enable; \
  CCTK_INT const EM_BC = RESTRICTED_MHD_EVOLVE_STRUCT.EM_BC; \
  CCTK_INT const MHD_MaxNumConstrainedVars = RESTRICTED_MHD_EVOLVE_STRUCT.MHD_MaxNumConstrainedVars; \
  CCTK_INT const MHD_MaxNumEvolvedVars = RESTRICTED_MHD_EVOLVE_STRUCT.MHD_MaxNumEvolvedVars; \
  CCTK_INT const Matter_BC = RESTRICTED_MHD_EVOLVE_STRUCT.Matter_BC; \
  CCTK_INT const Reconstruction = RESTRICTED_MHD_EVOLVE_STRUCT.Reconstruction; \
  CCTK_INT const adm_ppm_b = RESTRICTED_MHD_EVOLVE_STRUCT.adm_ppm_b; \
  CCTK_INT const allow_negative_eps_th = RESTRICTED_MHD_EVOLVE_STRUCT.allow_negative_eps_th; \
  CCTK_INT const artificial_cooling_enable = RESTRICTED_MHD_EVOLVE_STRUCT.artificial_cooling_enable; \
  CCTK_INT const compute_microphysics = RESTRICTED_MHD_EVOLVE_STRUCT.compute_microphysics; \
  CCTK_INT const constrained_transport_scheme = RESTRICTED_MHD_EVOLVE_STRUCT.constrained_transport_scheme; \
  CCTK_INT const cooling_in_St_eq = RESTRICTED_MHD_EVOLVE_STRUCT.cooling_in_St_eq; \
  CCTK_INT const em_evolve_enable = RESTRICTED_MHD_EVOLVE_STRUCT.em_evolve_enable; \
  CCTK_INT const em_gauge = RESTRICTED_MHD_EVOLVE_STRUCT.em_gauge; \
  CCTK_INT const enable_HARM_energyvariable = RESTRICTED_MHD_EVOLVE_STRUCT.enable_HARM_energyvariable; \
  CCTK_INT const enable_OS_collapse = RESTRICTED_MHD_EVOLVE_STRUCT.enable_OS_collapse; \
  CCTK_INT const enable_disk_em_flux_induction = RESTRICTED_MHD_EVOLVE_STRUCT.enable_disk_em_flux_induction; \
  CCTK_INT const enable_primitives_disk = RESTRICTED_MHD_EVOLVE_STRUCT.enable_primitives_disk; \
  CCTK_INT const enable_shocktest_primitive_mode = RESTRICTED_MHD_EVOLVE_STRUCT.enable_shocktest_primitive_mode; \
  CCTK_INT const enable_trace_field_line = RESTRICTED_MHD_EVOLVE_STRUCT.enable_trace_field_line; \
  CCTK_INT const eps_flag = RESTRICTED_MHD_EVOLVE_STRUCT.eps_flag; \
  CCTK_INT const ergo_star = RESTRICTED_MHD_EVOLVE_STRUCT.ergo_star; \
  CCTK_INT const force_font_fix_fail = RESTRICTED_MHD_EVOLVE_STRUCT.force_font_fix_fail; \
  CCTK_INT const horizon_enforce_rho_profile = RESTRICTED_MHD_EVOLVE_STRUCT.horizon_enforce_rho_profile; \
  CCTK_INT const hyperbolic_divergence_cleaning_centered_differencing = RESTRICTED_MHD_EVOLVE_STRUCT.hyperbolic_divergence_cleaning_centered_differencing; \
  CCTK_INT const hyperbolic_divergence_cleaning_enable = RESTRICTED_MHD_EVOLVE_STRUCT.hyperbolic_divergence_cleaning_enable; \
  CCTK_INT const m = RESTRICTED_MHD_EVOLVE_STRUCT.m; \
  CCTK_INT const microphysics_scheme = RESTRICTED_MHD_EVOLVE_STRUCT.microphysics_scheme; \
  CCTK_INT const neos = RESTRICTED_MHD_EVOLVE_STRUCT.neos; \
  CCTK_INT const nrhovec = RESTRICTED_MHD_EVOLVE_STRUCT.nrhovec; \
  CCTK_INT const primitives_debug = RESTRICTED_MHD_EVOLVE_STRUCT.primitives_debug; \
  CCTK_INT const primitives_solver = RESTRICTED_MHD_EVOLVE_STRUCT.primitives_solver; \
  CCTK_INT const rad_closure_scheme = RESTRICTED_MHD_EVOLVE_STRUCT.rad_closure_scheme; \
  CCTK_INT const rad_evolve_enable = RESTRICTED_MHD_EVOLVE_STRUCT.rad_evolve_enable; \
  CCTK_INT const rad_fix = RESTRICTED_MHD_EVOLVE_STRUCT.rad_fix; \
  CCTK_INT const rad_fourforce_enable = RESTRICTED_MHD_EVOLVE_STRUCT.rad_fourforce_enable; \
  CCTK_INT const reconstruct_Bitildes_instead_of_Bis = RESTRICTED_MHD_EVOLVE_STRUCT.reconstruct_Bitildes_instead_of_Bis; \
  CCTK_INT const reconstruct_Pthermal_instead_of_P = RESTRICTED_MHD_EVOLVE_STRUCT.reconstruct_Pthermal_instead_of_P; \
  CCTK_INT const tau_stildefix_enable = RESTRICTED_MHD_EVOLVE_STRUCT.tau_stildefix_enable; \
  CCTK_INT const use_HARM_primitives = RESTRICTED_MHD_EVOLVE_STRUCT.use_HARM_primitives; \
  CCTK_INT const use_central_scheme_instead_of_hll = RESTRICTED_MHD_EVOLVE_STRUCT.use_central_scheme_instead_of_hll; \
  CCTK_INT const use_new_code = RESTRICTED_MHD_EVOLVE_STRUCT.use_new_code; \
  enum { \
      dummy_RESTRICTED_MHD_EVOLVE_STRUCT_Erad_atm_cut = sizeof( Erad_atm_cut ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_Erad_cut = sizeof( Erad_cut ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_K_poly = sizeof( K_poly ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_M_B = sizeof( M_B ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_OS_T_Rmax = sizeof( OS_T_Rmax ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_P_fact = sizeof( P_fact ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_P_max = sizeof( P_max ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_Psi6threshold = sizeof( Psi6threshold ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_Sym_Bz = sizeof( Sym_Bz ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_T_fluid_cgs_atm = sizeof( T_fluid_cgs_atm ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_c_h_default = sizeof( c_h_default ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_c_p_default = sizeof( c_p_default ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_dxvec = sizeof( dxvec ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_eps_thermal_bhns = sizeof( eps_thermal_bhns ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_ergo_sigma = sizeof( ergo_sigma ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_gamma_OS = sizeof( gamma_OS ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_gamma_th = sizeof( gamma_th ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_lambda_line = sizeof( lambda_line ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_mhd_bigp = sizeof( mhd_bigp ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_mhd_mbh = sizeof( mhd_mbh ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_mhd_mdot = sizeof( mhd_mdot ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_mhd_r_crit = sizeof( mhd_r_crit ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_min_BH_radius = sizeof( min_BH_radius ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_pfloor = sizeof( pfloor ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_r_ex_line = sizeof( r_ex_line ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_rad_T_cutoff = sizeof( rad_T_cutoff ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_rad_T_fac = sizeof( rad_T_fac ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_rad_T_floor = sizeof( rad_T_floor ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_rad_T_pow = sizeof( rad_T_pow ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_rad_emissivity_abs = sizeof( rad_emissivity_abs ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_rad_emissivity_sct = sizeof( rad_emissivity_sct ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_rad_opacity_abs = sizeof( rad_opacity_abs ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_rad_opacity_sct = sizeof( rad_opacity_sct ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_rho_b_atm = sizeof( rho_b_atm ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_rho_b_max = sizeof( rho_b_max ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_rho_fact = sizeof( rho_fact ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_rhos_max = sizeof( rhos_max ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_sdots_o_rhot = sizeof( sdots_o_rhot ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_t_cool = sizeof( t_cool ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_tau_atm = sizeof( tau_atm ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_tau_fact = sizeof( tau_fact ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_A_BC_rad_mag_bondi_enable = sizeof( A_BC_rad_mag_bondi_enable ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_EM_BC = sizeof( EM_BC ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_MHD_MaxNumConstrainedVars = sizeof( MHD_MaxNumConstrainedVars ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_MHD_MaxNumEvolvedVars = sizeof( MHD_MaxNumEvolvedVars ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_Matter_BC = sizeof( Matter_BC ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_Reconstruction = sizeof( Reconstruction ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_adm_ppm_b = sizeof( adm_ppm_b ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_allow_negative_eps_th = sizeof( allow_negative_eps_th ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_artificial_cooling_enable = sizeof( artificial_cooling_enable ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_compute_microphysics = sizeof( compute_microphysics ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_constrained_transport_scheme = sizeof( constrained_transport_scheme ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_cooling_in_St_eq = sizeof( cooling_in_St_eq ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_em_evolve_enable = sizeof( em_evolve_enable ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_em_gauge = sizeof( em_gauge ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_enable_HARM_energyvariable = sizeof( enable_HARM_energyvariable ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_enable_OS_collapse = sizeof( enable_OS_collapse ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_enable_disk_em_flux_induction = sizeof( enable_disk_em_flux_induction ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_enable_primitives_disk = sizeof( enable_primitives_disk ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_enable_shocktest_primitive_mode = sizeof( enable_shocktest_primitive_mode ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_enable_trace_field_line = sizeof( enable_trace_field_line ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_eps_flag = sizeof( eps_flag ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_ergo_star = sizeof( ergo_star ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_force_font_fix_fail = sizeof( force_font_fix_fail ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_horizon_enforce_rho_profile = sizeof( horizon_enforce_rho_profile ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_hyperbolic_divergence_cleaning_centered_differencing = sizeof( hyperbolic_divergence_cleaning_centered_differencing ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_hyperbolic_divergence_cleaning_enable = sizeof( hyperbolic_divergence_cleaning_enable ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_m = sizeof( m ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_microphysics_scheme = sizeof( microphysics_scheme ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_neos = sizeof( neos ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_nrhovec = sizeof( nrhovec ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_primitives_debug = sizeof( primitives_debug ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_primitives_solver = sizeof( primitives_solver ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_rad_closure_scheme = sizeof( rad_closure_scheme ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_rad_evolve_enable = sizeof( rad_evolve_enable ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_rad_fix = sizeof( rad_fix ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_rad_fourforce_enable = sizeof( rad_fourforce_enable ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_reconstruct_Bitildes_instead_of_Bis = sizeof( reconstruct_Bitildes_instead_of_Bis ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_reconstruct_Pthermal_instead_of_P = sizeof( reconstruct_Pthermal_instead_of_P ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_tau_stildefix_enable = sizeof( tau_stildefix_enable ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_use_HARM_primitives = sizeof( use_HARM_primitives ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_use_central_scheme_instead_of_hll = sizeof( use_central_scheme_instead_of_hll ) \
    , dummy_RESTRICTED_MHD_EVOLVE_STRUCT_use_new_code = sizeof( use_new_code ) \
  }; \

