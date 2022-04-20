#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL Bxos4pi;
  CCTK_REAL P_in;
  CCTK_REAL P_out;
  CCTK_REAL disk_inner_radius;
  CCTK_REAL disk_outer_radius;
  CCTK_REAL npolyshock;
  CCTK_REAL rho_in;
  CCTK_REAL rho_out;
  CCTK_REAL rot_omega;
} PRIVATE_CYLINDRICAL2D_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_CYLINDRICAL2D_STRUCT_PARAMS \
  CCTK_REAL const Bxos4pi = PRIVATE_CYLINDRICAL2D_STRUCT.Bxos4pi; \
  CCTK_REAL const P_in = PRIVATE_CYLINDRICAL2D_STRUCT.P_in; \
  CCTK_REAL const P_out = PRIVATE_CYLINDRICAL2D_STRUCT.P_out; \
  CCTK_REAL const disk_inner_radius = PRIVATE_CYLINDRICAL2D_STRUCT.disk_inner_radius; \
  CCTK_REAL const disk_outer_radius = PRIVATE_CYLINDRICAL2D_STRUCT.disk_outer_radius; \
  CCTK_REAL const npolyshock = PRIVATE_CYLINDRICAL2D_STRUCT.npolyshock; \
  CCTK_REAL const rho_in = PRIVATE_CYLINDRICAL2D_STRUCT.rho_in; \
  CCTK_REAL const rho_out = PRIVATE_CYLINDRICAL2D_STRUCT.rho_out; \
  CCTK_REAL const rot_omega = PRIVATE_CYLINDRICAL2D_STRUCT.rot_omega; \
  enum { \
      dummy_PRIVATE_CYLINDRICAL2D_STRUCT_Bxos4pi = sizeof( Bxos4pi ) \
    , dummy_PRIVATE_CYLINDRICAL2D_STRUCT_P_in = sizeof( P_in ) \
    , dummy_PRIVATE_CYLINDRICAL2D_STRUCT_P_out = sizeof( P_out ) \
    , dummy_PRIVATE_CYLINDRICAL2D_STRUCT_disk_inner_radius = sizeof( disk_inner_radius ) \
    , dummy_PRIVATE_CYLINDRICAL2D_STRUCT_disk_outer_radius = sizeof( disk_outer_radius ) \
    , dummy_PRIVATE_CYLINDRICAL2D_STRUCT_npolyshock = sizeof( npolyshock ) \
    , dummy_PRIVATE_CYLINDRICAL2D_STRUCT_rho_in = sizeof( rho_in ) \
    , dummy_PRIVATE_CYLINDRICAL2D_STRUCT_rho_out = sizeof( rho_out ) \
    , dummy_PRIVATE_CYLINDRICAL2D_STRUCT_rot_omega = sizeof( rot_omega ) \
  }; \

