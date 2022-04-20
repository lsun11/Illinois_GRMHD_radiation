#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL Newton_tol;
  CCTK_REAL TP_Extend_Radius;
  CCTK_REAL TP_Tiny;
  CCTK_REAL TP_epsilon;
  CCTK_REAL adm_tol;
  CCTK_REAL center_offset[3];
  CCTK_REAL initial_lapse_psi_exponent;
  CCTK_REAL moncrief_radius_GW[100];
  CCTK_REAL par_P_minus[3];
  CCTK_REAL par_P_plus[3];
  CCTK_REAL par_S_minus[3];
  CCTK_REAL par_S_plus[3];
  CCTK_REAL par_b;
  CCTK_REAL par_m_minus;
  CCTK_REAL par_m_plus;
  CCTK_REAL target_M_minus;
  CCTK_REAL target_M_plus;
  const char * grid_setup_method;
  const char * initial_lapse;
  CCTK_INT Newton_maxit;
  CCTK_INT do_initial_debug_output;
  CCTK_INT do_residuum_debug_output;
  CCTK_INT give_bare_mass;
  CCTK_INT keep_u_around;
  CCTK_INT moncrief_gw_num_radii;
  CCTK_INT multiply_old_lapse;
  CCTK_INT npoints_A;
  CCTK_INT npoints_B;
  CCTK_INT npoints_phi;
  CCTK_INT rescale_sources;
  CCTK_INT schedule_in_ADMBase_InitialData;
  CCTK_INT solve_momentum_constraint;
  CCTK_INT swap_xz;
  CCTK_INT use_external_initial_guess;
  CCTK_INT use_sources;
  CCTK_INT verbose;
} PRIVATE_TWOPUNCTURESAEI_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_TWOPUNCTURESAEI_STRUCT_PARAMS \
  CCTK_REAL const Newton_tol = PRIVATE_TWOPUNCTURESAEI_STRUCT.Newton_tol; \
  CCTK_REAL const TP_Extend_Radius = PRIVATE_TWOPUNCTURESAEI_STRUCT.TP_Extend_Radius; \
  CCTK_REAL const TP_Tiny = PRIVATE_TWOPUNCTURESAEI_STRUCT.TP_Tiny; \
  CCTK_REAL const TP_epsilon = PRIVATE_TWOPUNCTURESAEI_STRUCT.TP_epsilon; \
  CCTK_REAL const adm_tol = PRIVATE_TWOPUNCTURESAEI_STRUCT.adm_tol; \
  CCTK_REAL const * const center_offset = PRIVATE_TWOPUNCTURESAEI_STRUCT.center_offset; \
  CCTK_REAL const initial_lapse_psi_exponent = PRIVATE_TWOPUNCTURESAEI_STRUCT.initial_lapse_psi_exponent; \
  CCTK_REAL const * const moncrief_radius_GW = PRIVATE_TWOPUNCTURESAEI_STRUCT.moncrief_radius_GW; \
  CCTK_REAL const * const par_P_minus = PRIVATE_TWOPUNCTURESAEI_STRUCT.par_P_minus; \
  CCTK_REAL const * const par_P_plus = PRIVATE_TWOPUNCTURESAEI_STRUCT.par_P_plus; \
  CCTK_REAL const * const par_S_minus = PRIVATE_TWOPUNCTURESAEI_STRUCT.par_S_minus; \
  CCTK_REAL const * const par_S_plus = PRIVATE_TWOPUNCTURESAEI_STRUCT.par_S_plus; \
  CCTK_REAL const par_b = PRIVATE_TWOPUNCTURESAEI_STRUCT.par_b; \
  CCTK_REAL const par_m_minus = PRIVATE_TWOPUNCTURESAEI_STRUCT.par_m_minus; \
  CCTK_REAL const par_m_plus = PRIVATE_TWOPUNCTURESAEI_STRUCT.par_m_plus; \
  CCTK_REAL const target_M_minus = PRIVATE_TWOPUNCTURESAEI_STRUCT.target_M_minus; \
  CCTK_REAL const target_M_plus = PRIVATE_TWOPUNCTURESAEI_STRUCT.target_M_plus; \
  const char * const grid_setup_method = PRIVATE_TWOPUNCTURESAEI_STRUCT.grid_setup_method; \
  const char * const initial_lapse = PRIVATE_TWOPUNCTURESAEI_STRUCT.initial_lapse; \
  CCTK_INT const Newton_maxit = PRIVATE_TWOPUNCTURESAEI_STRUCT.Newton_maxit; \
  CCTK_INT const do_initial_debug_output = PRIVATE_TWOPUNCTURESAEI_STRUCT.do_initial_debug_output; \
  CCTK_INT const do_residuum_debug_output = PRIVATE_TWOPUNCTURESAEI_STRUCT.do_residuum_debug_output; \
  CCTK_INT const give_bare_mass = PRIVATE_TWOPUNCTURESAEI_STRUCT.give_bare_mass; \
  CCTK_INT const keep_u_around = PRIVATE_TWOPUNCTURESAEI_STRUCT.keep_u_around; \
  CCTK_INT const moncrief_gw_num_radii = PRIVATE_TWOPUNCTURESAEI_STRUCT.moncrief_gw_num_radii; \
  CCTK_INT const multiply_old_lapse = PRIVATE_TWOPUNCTURESAEI_STRUCT.multiply_old_lapse; \
  CCTK_INT const npoints_A = PRIVATE_TWOPUNCTURESAEI_STRUCT.npoints_A; \
  CCTK_INT const npoints_B = PRIVATE_TWOPUNCTURESAEI_STRUCT.npoints_B; \
  CCTK_INT const npoints_phi = PRIVATE_TWOPUNCTURESAEI_STRUCT.npoints_phi; \
  CCTK_INT const rescale_sources = PRIVATE_TWOPUNCTURESAEI_STRUCT.rescale_sources; \
  CCTK_INT const schedule_in_ADMBase_InitialData = PRIVATE_TWOPUNCTURESAEI_STRUCT.schedule_in_ADMBase_InitialData; \
  CCTK_INT const solve_momentum_constraint = PRIVATE_TWOPUNCTURESAEI_STRUCT.solve_momentum_constraint; \
  CCTK_INT const swap_xz = PRIVATE_TWOPUNCTURESAEI_STRUCT.swap_xz; \
  CCTK_INT const use_external_initial_guess = PRIVATE_TWOPUNCTURESAEI_STRUCT.use_external_initial_guess; \
  CCTK_INT const use_sources = PRIVATE_TWOPUNCTURESAEI_STRUCT.use_sources; \
  CCTK_INT const verbose = PRIVATE_TWOPUNCTURESAEI_STRUCT.verbose; \
  enum { \
      dummy_PRIVATE_TWOPUNCTURESAEI_STRUCT_Newton_tol = sizeof( Newton_tol ) \
    , dummy_PRIVATE_TWOPUNCTURESAEI_STRUCT_TP_Extend_Radius = sizeof( TP_Extend_Radius ) \
    , dummy_PRIVATE_TWOPUNCTURESAEI_STRUCT_TP_Tiny = sizeof( TP_Tiny ) \
    , dummy_PRIVATE_TWOPUNCTURESAEI_STRUCT_TP_epsilon = sizeof( TP_epsilon ) \
    , dummy_PRIVATE_TWOPUNCTURESAEI_STRUCT_adm_tol = sizeof( adm_tol ) \
    , dummy_PRIVATE_TWOPUNCTURESAEI_STRUCT_center_offset = sizeof( center_offset ) \
    , dummy_PRIVATE_TWOPUNCTURESAEI_STRUCT_initial_lapse_psi_exponent = sizeof( initial_lapse_psi_exponent ) \
    , dummy_PRIVATE_TWOPUNCTURESAEI_STRUCT_moncrief_radius_GW = sizeof( moncrief_radius_GW ) \
    , dummy_PRIVATE_TWOPUNCTURESAEI_STRUCT_par_P_minus = sizeof( par_P_minus ) \
    , dummy_PRIVATE_TWOPUNCTURESAEI_STRUCT_par_P_plus = sizeof( par_P_plus ) \
    , dummy_PRIVATE_TWOPUNCTURESAEI_STRUCT_par_S_minus = sizeof( par_S_minus ) \
    , dummy_PRIVATE_TWOPUNCTURESAEI_STRUCT_par_S_plus = sizeof( par_S_plus ) \
    , dummy_PRIVATE_TWOPUNCTURESAEI_STRUCT_par_b = sizeof( par_b ) \
    , dummy_PRIVATE_TWOPUNCTURESAEI_STRUCT_par_m_minus = sizeof( par_m_minus ) \
    , dummy_PRIVATE_TWOPUNCTURESAEI_STRUCT_par_m_plus = sizeof( par_m_plus ) \
    , dummy_PRIVATE_TWOPUNCTURESAEI_STRUCT_target_M_minus = sizeof( target_M_minus ) \
    , dummy_PRIVATE_TWOPUNCTURESAEI_STRUCT_target_M_plus = sizeof( target_M_plus ) \
    , dummy_PRIVATE_TWOPUNCTURESAEI_STRUCT_grid_setup_method = sizeof( grid_setup_method ) \
    , dummy_PRIVATE_TWOPUNCTURESAEI_STRUCT_initial_lapse = sizeof( initial_lapse ) \
    , dummy_PRIVATE_TWOPUNCTURESAEI_STRUCT_Newton_maxit = sizeof( Newton_maxit ) \
    , dummy_PRIVATE_TWOPUNCTURESAEI_STRUCT_do_initial_debug_output = sizeof( do_initial_debug_output ) \
    , dummy_PRIVATE_TWOPUNCTURESAEI_STRUCT_do_residuum_debug_output = sizeof( do_residuum_debug_output ) \
    , dummy_PRIVATE_TWOPUNCTURESAEI_STRUCT_give_bare_mass = sizeof( give_bare_mass ) \
    , dummy_PRIVATE_TWOPUNCTURESAEI_STRUCT_keep_u_around = sizeof( keep_u_around ) \
    , dummy_PRIVATE_TWOPUNCTURESAEI_STRUCT_moncrief_gw_num_radii = sizeof( moncrief_gw_num_radii ) \
    , dummy_PRIVATE_TWOPUNCTURESAEI_STRUCT_multiply_old_lapse = sizeof( multiply_old_lapse ) \
    , dummy_PRIVATE_TWOPUNCTURESAEI_STRUCT_npoints_A = sizeof( npoints_A ) \
    , dummy_PRIVATE_TWOPUNCTURESAEI_STRUCT_npoints_B = sizeof( npoints_B ) \
    , dummy_PRIVATE_TWOPUNCTURESAEI_STRUCT_npoints_phi = sizeof( npoints_phi ) \
    , dummy_PRIVATE_TWOPUNCTURESAEI_STRUCT_rescale_sources = sizeof( rescale_sources ) \
    , dummy_PRIVATE_TWOPUNCTURESAEI_STRUCT_schedule_in_ADMBase_InitialData = sizeof( schedule_in_ADMBase_InitialData ) \
    , dummy_PRIVATE_TWOPUNCTURESAEI_STRUCT_solve_momentum_constraint = sizeof( solve_momentum_constraint ) \
    , dummy_PRIVATE_TWOPUNCTURESAEI_STRUCT_swap_xz = sizeof( swap_xz ) \
    , dummy_PRIVATE_TWOPUNCTURESAEI_STRUCT_use_external_initial_guess = sizeof( use_external_initial_guess ) \
    , dummy_PRIVATE_TWOPUNCTURESAEI_STRUCT_use_sources = sizeof( use_sources ) \
    , dummy_PRIVATE_TWOPUNCTURESAEI_STRUCT_verbose = sizeof( verbose ) \
  }; \

