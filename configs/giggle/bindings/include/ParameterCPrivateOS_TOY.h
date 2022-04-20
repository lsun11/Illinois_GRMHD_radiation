#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL E_over_rho;
  CCTK_REAL Fx_over_rho;
  CCTK_REAL Fy_over_rho;
  CCTK_REAL Fz_over_rho;
  CCTK_REAL M0_initial;
  CCTK_REAL M_OS;
  CCTK_REAL P_over_rho;
  CCTK_REAL Po4PiB;
  CCTK_REAL R_OS;
  CCTK_REAL opt_depth_a;
  CCTK_REAL opt_depth_s;
  CCTK_REAL particle_rad_cut;
  CCTK_REAL rounding;
  CCTK_REAL xc_OS;
  CCTK_REAL yc_OS;
  CCTK_REAL zc_OS;
  CCTK_INT narr;
} PRIVATE_OS_TOY_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_OS_TOY_STRUCT_PARAMS \
  CCTK_REAL const E_over_rho = PRIVATE_OS_TOY_STRUCT.E_over_rho; \
  CCTK_REAL const Fx_over_rho = PRIVATE_OS_TOY_STRUCT.Fx_over_rho; \
  CCTK_REAL const Fy_over_rho = PRIVATE_OS_TOY_STRUCT.Fy_over_rho; \
  CCTK_REAL const Fz_over_rho = PRIVATE_OS_TOY_STRUCT.Fz_over_rho; \
  CCTK_REAL const M0_initial = PRIVATE_OS_TOY_STRUCT.M0_initial; \
  CCTK_REAL const M_OS = PRIVATE_OS_TOY_STRUCT.M_OS; \
  CCTK_REAL const P_over_rho = PRIVATE_OS_TOY_STRUCT.P_over_rho; \
  CCTK_REAL const Po4PiB = PRIVATE_OS_TOY_STRUCT.Po4PiB; \
  CCTK_REAL const R_OS = PRIVATE_OS_TOY_STRUCT.R_OS; \
  CCTK_REAL const opt_depth_a = PRIVATE_OS_TOY_STRUCT.opt_depth_a; \
  CCTK_REAL const opt_depth_s = PRIVATE_OS_TOY_STRUCT.opt_depth_s; \
  CCTK_REAL const particle_rad_cut = PRIVATE_OS_TOY_STRUCT.particle_rad_cut; \
  CCTK_REAL const rounding = PRIVATE_OS_TOY_STRUCT.rounding; \
  CCTK_REAL const xc_OS = PRIVATE_OS_TOY_STRUCT.xc_OS; \
  CCTK_REAL const yc_OS = PRIVATE_OS_TOY_STRUCT.yc_OS; \
  CCTK_REAL const zc_OS = PRIVATE_OS_TOY_STRUCT.zc_OS; \
  CCTK_INT const narr = PRIVATE_OS_TOY_STRUCT.narr; \
  enum { \
      dummy_PRIVATE_OS_TOY_STRUCT_E_over_rho = sizeof( E_over_rho ) \
    , dummy_PRIVATE_OS_TOY_STRUCT_Fx_over_rho = sizeof( Fx_over_rho ) \
    , dummy_PRIVATE_OS_TOY_STRUCT_Fy_over_rho = sizeof( Fy_over_rho ) \
    , dummy_PRIVATE_OS_TOY_STRUCT_Fz_over_rho = sizeof( Fz_over_rho ) \
    , dummy_PRIVATE_OS_TOY_STRUCT_M0_initial = sizeof( M0_initial ) \
    , dummy_PRIVATE_OS_TOY_STRUCT_M_OS = sizeof( M_OS ) \
    , dummy_PRIVATE_OS_TOY_STRUCT_P_over_rho = sizeof( P_over_rho ) \
    , dummy_PRIVATE_OS_TOY_STRUCT_Po4PiB = sizeof( Po4PiB ) \
    , dummy_PRIVATE_OS_TOY_STRUCT_R_OS = sizeof( R_OS ) \
    , dummy_PRIVATE_OS_TOY_STRUCT_opt_depth_a = sizeof( opt_depth_a ) \
    , dummy_PRIVATE_OS_TOY_STRUCT_opt_depth_s = sizeof( opt_depth_s ) \
    , dummy_PRIVATE_OS_TOY_STRUCT_particle_rad_cut = sizeof( particle_rad_cut ) \
    , dummy_PRIVATE_OS_TOY_STRUCT_rounding = sizeof( rounding ) \
    , dummy_PRIVATE_OS_TOY_STRUCT_xc_OS = sizeof( xc_OS ) \
    , dummy_PRIVATE_OS_TOY_STRUCT_yc_OS = sizeof( yc_OS ) \
    , dummy_PRIVATE_OS_TOY_STRUCT_zc_OS = sizeof( zc_OS ) \
    , dummy_PRIVATE_OS_TOY_STRUCT_narr = sizeof( narr ) \
  }; \

