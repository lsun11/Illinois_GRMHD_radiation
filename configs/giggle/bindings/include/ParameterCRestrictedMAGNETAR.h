#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL Omega_Frame;
  CCTK_REAL P_deplete;
  CCTK_REAL betam1;
  CCTK_REAL lapse_regrid;
  CCTK_REAL magnetar_P_max;
  CCTK_REAL magnetar_rho_star_max;
  CCTK_REAL magnetar_tau_max;
  CCTK_REAL p_c;
  const char * mag_interp_operator;
  CCTK_INT binfile_checkpoint_iteration;
  CCTK_INT binfile_restart;
  CCTK_INT em_field_type;
} RESTRICTED_MAGNETAR_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_RESTRICTED_MAGNETAR_STRUCT_PARAMS \
  CCTK_REAL const Omega_Frame = RESTRICTED_MAGNETAR_STRUCT.Omega_Frame; \
  CCTK_REAL const P_deplete = RESTRICTED_MAGNETAR_STRUCT.P_deplete; \
  CCTK_REAL const betam1 = RESTRICTED_MAGNETAR_STRUCT.betam1; \
  CCTK_REAL const lapse_regrid = RESTRICTED_MAGNETAR_STRUCT.lapse_regrid; \
  CCTK_REAL const magnetar_P_max = RESTRICTED_MAGNETAR_STRUCT.magnetar_P_max; \
  CCTK_REAL const magnetar_rho_star_max = RESTRICTED_MAGNETAR_STRUCT.magnetar_rho_star_max; \
  CCTK_REAL const magnetar_tau_max = RESTRICTED_MAGNETAR_STRUCT.magnetar_tau_max; \
  CCTK_REAL const p_c = RESTRICTED_MAGNETAR_STRUCT.p_c; \
  const char * const mag_interp_operator = RESTRICTED_MAGNETAR_STRUCT.mag_interp_operator; \
  CCTK_INT const binfile_checkpoint_iteration = RESTRICTED_MAGNETAR_STRUCT.binfile_checkpoint_iteration; \
  CCTK_INT const binfile_restart = RESTRICTED_MAGNETAR_STRUCT.binfile_restart; \
  CCTK_INT const em_field_type = RESTRICTED_MAGNETAR_STRUCT.em_field_type; \
  enum { \
      dummy_RESTRICTED_MAGNETAR_STRUCT_Omega_Frame = sizeof( Omega_Frame ) \
    , dummy_RESTRICTED_MAGNETAR_STRUCT_P_deplete = sizeof( P_deplete ) \
    , dummy_RESTRICTED_MAGNETAR_STRUCT_betam1 = sizeof( betam1 ) \
    , dummy_RESTRICTED_MAGNETAR_STRUCT_lapse_regrid = sizeof( lapse_regrid ) \
    , dummy_RESTRICTED_MAGNETAR_STRUCT_magnetar_P_max = sizeof( magnetar_P_max ) \
    , dummy_RESTRICTED_MAGNETAR_STRUCT_magnetar_rho_star_max = sizeof( magnetar_rho_star_max ) \
    , dummy_RESTRICTED_MAGNETAR_STRUCT_magnetar_tau_max = sizeof( magnetar_tau_max ) \
    , dummy_RESTRICTED_MAGNETAR_STRUCT_p_c = sizeof( p_c ) \
    , dummy_RESTRICTED_MAGNETAR_STRUCT_mag_interp_operator = sizeof( mag_interp_operator ) \
    , dummy_RESTRICTED_MAGNETAR_STRUCT_binfile_checkpoint_iteration = sizeof( binfile_checkpoint_iteration ) \
    , dummy_RESTRICTED_MAGNETAR_STRUCT_binfile_restart = sizeof( binfile_restart ) \
    , dummy_RESTRICTED_MAGNETAR_STRUCT_em_field_type = sizeof( em_field_type ) \
  }; \

