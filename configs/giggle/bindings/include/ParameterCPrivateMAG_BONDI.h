#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL Bstrength;
  CCTK_REAL npoly;
  CCTK_REAL p_atm;
  CCTK_REAL r0;
  CCTK_REAL r_sonic;
  CCTK_REAL rho_atm;
  CCTK_INT N_phi;
  CCTK_INT N_theta;
  CCTK_INT enable_alt_atmosphere;
  CCTK_INT enable_movingbox;
  CCTK_INT puncture_id;
} PRIVATE_MAG_BONDI_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_MAG_BONDI_STRUCT_PARAMS \
  CCTK_REAL const Bstrength = PRIVATE_MAG_BONDI_STRUCT.Bstrength; \
  CCTK_REAL const npoly = PRIVATE_MAG_BONDI_STRUCT.npoly; \
  CCTK_REAL const p_atm = PRIVATE_MAG_BONDI_STRUCT.p_atm; \
  CCTK_REAL const r0 = PRIVATE_MAG_BONDI_STRUCT.r0; \
  CCTK_REAL const r_sonic = PRIVATE_MAG_BONDI_STRUCT.r_sonic; \
  CCTK_REAL const rho_atm = PRIVATE_MAG_BONDI_STRUCT.rho_atm; \
  CCTK_INT const N_phi = PRIVATE_MAG_BONDI_STRUCT.N_phi; \
  CCTK_INT const N_theta = PRIVATE_MAG_BONDI_STRUCT.N_theta; \
  CCTK_INT const enable_alt_atmosphere = PRIVATE_MAG_BONDI_STRUCT.enable_alt_atmosphere; \
  CCTK_INT const enable_movingbox = PRIVATE_MAG_BONDI_STRUCT.enable_movingbox; \
  CCTK_INT const puncture_id = PRIVATE_MAG_BONDI_STRUCT.puncture_id; \
  enum { \
      dummy_PRIVATE_MAG_BONDI_STRUCT_Bstrength = sizeof( Bstrength ) \
    , dummy_PRIVATE_MAG_BONDI_STRUCT_npoly = sizeof( npoly ) \
    , dummy_PRIVATE_MAG_BONDI_STRUCT_p_atm = sizeof( p_atm ) \
    , dummy_PRIVATE_MAG_BONDI_STRUCT_r0 = sizeof( r0 ) \
    , dummy_PRIVATE_MAG_BONDI_STRUCT_r_sonic = sizeof( r_sonic ) \
    , dummy_PRIVATE_MAG_BONDI_STRUCT_rho_atm = sizeof( rho_atm ) \
    , dummy_PRIVATE_MAG_BONDI_STRUCT_N_phi = sizeof( N_phi ) \
    , dummy_PRIVATE_MAG_BONDI_STRUCT_N_theta = sizeof( N_theta ) \
    , dummy_PRIVATE_MAG_BONDI_STRUCT_enable_alt_atmosphere = sizeof( enable_alt_atmosphere ) \
    , dummy_PRIVATE_MAG_BONDI_STRUCT_enable_movingbox = sizeof( enable_movingbox ) \
    , dummy_PRIVATE_MAG_BONDI_STRUCT_puncture_id = sizeof( puncture_id ) \
  }; \

