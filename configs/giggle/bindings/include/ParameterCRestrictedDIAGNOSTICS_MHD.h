#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL CoMExcision_Radius;
  CCTK_REAL M_ADM;
  CCTK_REAL initial_monopole_value;
  CCTK_REAL inner_lum_rad_ratio;
  CCTK_REAL lum_outer_rad;
  CCTK_REAL radius_esc1;
  CCTK_REAL radius_esc2;
  CCTK_REAL radius_esc3;
  CCTK_REAL radius_esc4;
  CCTK_REAL rhob_cutoff;
  CCTK_REAL rhosurf;
  CCTK_REAL rhosurf_rmax;
  CCTK_REAL rhosurf_rmin;
  CCTK_REAL u0sch;
  CCTK_INT N_rad_ray;
  CCTK_INT Nphi_points;
  CCTK_INT Ntheta_points;
  CCTK_INT const_rad_surf_diagnostics;
  CCTK_INT drho_dtau_calc_enable;
  CCTK_INT escape_mass_diag;
  CCTK_INT luminosity_diagnostics;
} RESTRICTED_DIAGNOSTICS_MHD_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_RESTRICTED_DIAGNOSTICS_MHD_STRUCT_PARAMS \
  CCTK_REAL const CoMExcision_Radius = RESTRICTED_DIAGNOSTICS_MHD_STRUCT.CoMExcision_Radius; \
  CCTK_REAL const M_ADM = RESTRICTED_DIAGNOSTICS_MHD_STRUCT.M_ADM; \
  CCTK_REAL const initial_monopole_value = RESTRICTED_DIAGNOSTICS_MHD_STRUCT.initial_monopole_value; \
  CCTK_REAL const inner_lum_rad_ratio = RESTRICTED_DIAGNOSTICS_MHD_STRUCT.inner_lum_rad_ratio; \
  CCTK_REAL const lum_outer_rad = RESTRICTED_DIAGNOSTICS_MHD_STRUCT.lum_outer_rad; \
  CCTK_REAL const radius_esc1 = RESTRICTED_DIAGNOSTICS_MHD_STRUCT.radius_esc1; \
  CCTK_REAL const radius_esc2 = RESTRICTED_DIAGNOSTICS_MHD_STRUCT.radius_esc2; \
  CCTK_REAL const radius_esc3 = RESTRICTED_DIAGNOSTICS_MHD_STRUCT.radius_esc3; \
  CCTK_REAL const radius_esc4 = RESTRICTED_DIAGNOSTICS_MHD_STRUCT.radius_esc4; \
  CCTK_REAL const rhob_cutoff = RESTRICTED_DIAGNOSTICS_MHD_STRUCT.rhob_cutoff; \
  CCTK_REAL const rhosurf = RESTRICTED_DIAGNOSTICS_MHD_STRUCT.rhosurf; \
  CCTK_REAL const rhosurf_rmax = RESTRICTED_DIAGNOSTICS_MHD_STRUCT.rhosurf_rmax; \
  CCTK_REAL const rhosurf_rmin = RESTRICTED_DIAGNOSTICS_MHD_STRUCT.rhosurf_rmin; \
  CCTK_REAL const u0sch = RESTRICTED_DIAGNOSTICS_MHD_STRUCT.u0sch; \
  CCTK_INT const N_rad_ray = RESTRICTED_DIAGNOSTICS_MHD_STRUCT.N_rad_ray; \
  CCTK_INT const Nphi_points = RESTRICTED_DIAGNOSTICS_MHD_STRUCT.Nphi_points; \
  CCTK_INT const Ntheta_points = RESTRICTED_DIAGNOSTICS_MHD_STRUCT.Ntheta_points; \
  CCTK_INT const const_rad_surf_diagnostics = RESTRICTED_DIAGNOSTICS_MHD_STRUCT.const_rad_surf_diagnostics; \
  CCTK_INT const drho_dtau_calc_enable = RESTRICTED_DIAGNOSTICS_MHD_STRUCT.drho_dtau_calc_enable; \
  CCTK_INT const escape_mass_diag = RESTRICTED_DIAGNOSTICS_MHD_STRUCT.escape_mass_diag; \
  CCTK_INT const luminosity_diagnostics = RESTRICTED_DIAGNOSTICS_MHD_STRUCT.luminosity_diagnostics; \
  enum { \
      dummy_RESTRICTED_DIAGNOSTICS_MHD_STRUCT_CoMExcision_Radius = sizeof( CoMExcision_Radius ) \
    , dummy_RESTRICTED_DIAGNOSTICS_MHD_STRUCT_M_ADM = sizeof( M_ADM ) \
    , dummy_RESTRICTED_DIAGNOSTICS_MHD_STRUCT_initial_monopole_value = sizeof( initial_monopole_value ) \
    , dummy_RESTRICTED_DIAGNOSTICS_MHD_STRUCT_inner_lum_rad_ratio = sizeof( inner_lum_rad_ratio ) \
    , dummy_RESTRICTED_DIAGNOSTICS_MHD_STRUCT_lum_outer_rad = sizeof( lum_outer_rad ) \
    , dummy_RESTRICTED_DIAGNOSTICS_MHD_STRUCT_radius_esc1 = sizeof( radius_esc1 ) \
    , dummy_RESTRICTED_DIAGNOSTICS_MHD_STRUCT_radius_esc2 = sizeof( radius_esc2 ) \
    , dummy_RESTRICTED_DIAGNOSTICS_MHD_STRUCT_radius_esc3 = sizeof( radius_esc3 ) \
    , dummy_RESTRICTED_DIAGNOSTICS_MHD_STRUCT_radius_esc4 = sizeof( radius_esc4 ) \
    , dummy_RESTRICTED_DIAGNOSTICS_MHD_STRUCT_rhob_cutoff = sizeof( rhob_cutoff ) \
    , dummy_RESTRICTED_DIAGNOSTICS_MHD_STRUCT_rhosurf = sizeof( rhosurf ) \
    , dummy_RESTRICTED_DIAGNOSTICS_MHD_STRUCT_rhosurf_rmax = sizeof( rhosurf_rmax ) \
    , dummy_RESTRICTED_DIAGNOSTICS_MHD_STRUCT_rhosurf_rmin = sizeof( rhosurf_rmin ) \
    , dummy_RESTRICTED_DIAGNOSTICS_MHD_STRUCT_u0sch = sizeof( u0sch ) \
    , dummy_RESTRICTED_DIAGNOSTICS_MHD_STRUCT_N_rad_ray = sizeof( N_rad_ray ) \
    , dummy_RESTRICTED_DIAGNOSTICS_MHD_STRUCT_Nphi_points = sizeof( Nphi_points ) \
    , dummy_RESTRICTED_DIAGNOSTICS_MHD_STRUCT_Ntheta_points = sizeof( Ntheta_points ) \
    , dummy_RESTRICTED_DIAGNOSTICS_MHD_STRUCT_const_rad_surf_diagnostics = sizeof( const_rad_surf_diagnostics ) \
    , dummy_RESTRICTED_DIAGNOSTICS_MHD_STRUCT_drho_dtau_calc_enable = sizeof( drho_dtau_calc_enable ) \
    , dummy_RESTRICTED_DIAGNOSTICS_MHD_STRUCT_escape_mass_diag = sizeof( escape_mass_diag ) \
    , dummy_RESTRICTED_DIAGNOSTICS_MHD_STRUCT_luminosity_diagnostics = sizeof( luminosity_diagnostics ) \
  }; \

