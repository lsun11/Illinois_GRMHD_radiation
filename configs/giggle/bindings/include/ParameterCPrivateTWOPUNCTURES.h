#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL bh_mass_minus;
  CCTK_REAL bh_mass_plus;
  CCTK_REAL bh_px_minus;
  CCTK_REAL bh_px_plus;
  CCTK_REAL bh_py_minus;
  CCTK_REAL bh_py_plus;
  CCTK_REAL bh_spin_minus;
  CCTK_REAL bh_spin_plus;
  CCTK_REAL excis_radius;
  CCTK_REAL half_binary_separation;
  CCTK_REAL moncrief_radius_GW[100];
  CCTK_REAL x_offset;
  CCTK_INT fill_excision_enable;
  CCTK_INT genID_cmdline_output_enable;
  CCTK_INT moncrief_gw_num_radii;
} PRIVATE_TWOPUNCTURES_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_TWOPUNCTURES_STRUCT_PARAMS \
  CCTK_REAL const bh_mass_minus = PRIVATE_TWOPUNCTURES_STRUCT.bh_mass_minus; \
  CCTK_REAL const bh_mass_plus = PRIVATE_TWOPUNCTURES_STRUCT.bh_mass_plus; \
  CCTK_REAL const bh_px_minus = PRIVATE_TWOPUNCTURES_STRUCT.bh_px_minus; \
  CCTK_REAL const bh_px_plus = PRIVATE_TWOPUNCTURES_STRUCT.bh_px_plus; \
  CCTK_REAL const bh_py_minus = PRIVATE_TWOPUNCTURES_STRUCT.bh_py_minus; \
  CCTK_REAL const bh_py_plus = PRIVATE_TWOPUNCTURES_STRUCT.bh_py_plus; \
  CCTK_REAL const bh_spin_minus = PRIVATE_TWOPUNCTURES_STRUCT.bh_spin_minus; \
  CCTK_REAL const bh_spin_plus = PRIVATE_TWOPUNCTURES_STRUCT.bh_spin_plus; \
  CCTK_REAL const excis_radius = PRIVATE_TWOPUNCTURES_STRUCT.excis_radius; \
  CCTK_REAL const half_binary_separation = PRIVATE_TWOPUNCTURES_STRUCT.half_binary_separation; \
  CCTK_REAL const * const moncrief_radius_GW = PRIVATE_TWOPUNCTURES_STRUCT.moncrief_radius_GW; \
  CCTK_REAL const x_offset = PRIVATE_TWOPUNCTURES_STRUCT.x_offset; \
  CCTK_INT const fill_excision_enable = PRIVATE_TWOPUNCTURES_STRUCT.fill_excision_enable; \
  CCTK_INT const genID_cmdline_output_enable = PRIVATE_TWOPUNCTURES_STRUCT.genID_cmdline_output_enable; \
  CCTK_INT const moncrief_gw_num_radii = PRIVATE_TWOPUNCTURES_STRUCT.moncrief_gw_num_radii; \
  enum { \
      dummy_PRIVATE_TWOPUNCTURES_STRUCT_bh_mass_minus = sizeof( bh_mass_minus ) \
    , dummy_PRIVATE_TWOPUNCTURES_STRUCT_bh_mass_plus = sizeof( bh_mass_plus ) \
    , dummy_PRIVATE_TWOPUNCTURES_STRUCT_bh_px_minus = sizeof( bh_px_minus ) \
    , dummy_PRIVATE_TWOPUNCTURES_STRUCT_bh_px_plus = sizeof( bh_px_plus ) \
    , dummy_PRIVATE_TWOPUNCTURES_STRUCT_bh_py_minus = sizeof( bh_py_minus ) \
    , dummy_PRIVATE_TWOPUNCTURES_STRUCT_bh_py_plus = sizeof( bh_py_plus ) \
    , dummy_PRIVATE_TWOPUNCTURES_STRUCT_bh_spin_minus = sizeof( bh_spin_minus ) \
    , dummy_PRIVATE_TWOPUNCTURES_STRUCT_bh_spin_plus = sizeof( bh_spin_plus ) \
    , dummy_PRIVATE_TWOPUNCTURES_STRUCT_excis_radius = sizeof( excis_radius ) \
    , dummy_PRIVATE_TWOPUNCTURES_STRUCT_half_binary_separation = sizeof( half_binary_separation ) \
    , dummy_PRIVATE_TWOPUNCTURES_STRUCT_moncrief_radius_GW = sizeof( moncrief_radius_GW ) \
    , dummy_PRIVATE_TWOPUNCTURES_STRUCT_x_offset = sizeof( x_offset ) \
    , dummy_PRIVATE_TWOPUNCTURES_STRUCT_fill_excision_enable = sizeof( fill_excision_enable ) \
    , dummy_PRIVATE_TWOPUNCTURES_STRUCT_genID_cmdline_output_enable = sizeof( genID_cmdline_output_enable ) \
    , dummy_PRIVATE_TWOPUNCTURES_STRUCT_moncrief_gw_num_radii = sizeof( moncrief_gw_num_radii ) \
  }; \

