#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL moncrief_radius_GW[11];
  CCTK_REAL p_atm;
  CCTK_REAL p_atm_index;
  CCTK_REAL rho_atm;
  CCTK_REAL rho_atm_index;
  CCTK_INT enable_alt_atmosphere;
  CCTK_INT moncrief_gw_num_radii;
} PRIVATE_MAGNETAR_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_MAGNETAR_STRUCT_PARAMS \
  CCTK_REAL const * const moncrief_radius_GW = PRIVATE_MAGNETAR_STRUCT.moncrief_radius_GW; \
  CCTK_REAL const p_atm = PRIVATE_MAGNETAR_STRUCT.p_atm; \
  CCTK_REAL const p_atm_index = PRIVATE_MAGNETAR_STRUCT.p_atm_index; \
  CCTK_REAL const rho_atm = PRIVATE_MAGNETAR_STRUCT.rho_atm; \
  CCTK_REAL const rho_atm_index = PRIVATE_MAGNETAR_STRUCT.rho_atm_index; \
  CCTK_INT const enable_alt_atmosphere = PRIVATE_MAGNETAR_STRUCT.enable_alt_atmosphere; \
  CCTK_INT const moncrief_gw_num_radii = PRIVATE_MAGNETAR_STRUCT.moncrief_gw_num_radii; \
  enum { \
      dummy_PRIVATE_MAGNETAR_STRUCT_moncrief_radius_GW = sizeof( moncrief_radius_GW ) \
    , dummy_PRIVATE_MAGNETAR_STRUCT_p_atm = sizeof( p_atm ) \
    , dummy_PRIVATE_MAGNETAR_STRUCT_p_atm_index = sizeof( p_atm_index ) \
    , dummy_PRIVATE_MAGNETAR_STRUCT_rho_atm = sizeof( rho_atm ) \
    , dummy_PRIVATE_MAGNETAR_STRUCT_rho_atm_index = sizeof( rho_atm_index ) \
    , dummy_PRIVATE_MAGNETAR_STRUCT_enable_alt_atmosphere = sizeof( enable_alt_atmosphere ) \
    , dummy_PRIVATE_MAGNETAR_STRUCT_moncrief_gw_num_radii = sizeof( moncrief_gw_num_radii ) \
  }; \

