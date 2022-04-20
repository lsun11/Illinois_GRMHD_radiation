#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL each_black_hole_mass;
  CCTK_REAL each_black_hole_momentum;
  CCTK_REAL each_black_hole_spin;
  CCTK_REAL excis_radius;
  CCTK_REAL moncrief_radius_GW[100];
  CCTK_REAL total_binary_separation;
  CCTK_INT fill_excision_enable;
  CCTK_INT genID_cmdline_output_enable;
  CCTK_INT moncrief_gw_num_radii;
} PRIVATE_BBHLORENE_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_BBHLORENE_STRUCT_PARAMS \
  CCTK_REAL const each_black_hole_mass = PRIVATE_BBHLORENE_STRUCT.each_black_hole_mass; \
  CCTK_REAL const each_black_hole_momentum = PRIVATE_BBHLORENE_STRUCT.each_black_hole_momentum; \
  CCTK_REAL const each_black_hole_spin = PRIVATE_BBHLORENE_STRUCT.each_black_hole_spin; \
  CCTK_REAL const excis_radius = PRIVATE_BBHLORENE_STRUCT.excis_radius; \
  CCTK_REAL const * const moncrief_radius_GW = PRIVATE_BBHLORENE_STRUCT.moncrief_radius_GW; \
  CCTK_REAL const total_binary_separation = PRIVATE_BBHLORENE_STRUCT.total_binary_separation; \
  CCTK_INT const fill_excision_enable = PRIVATE_BBHLORENE_STRUCT.fill_excision_enable; \
  CCTK_INT const genID_cmdline_output_enable = PRIVATE_BBHLORENE_STRUCT.genID_cmdline_output_enable; \
  CCTK_INT const moncrief_gw_num_radii = PRIVATE_BBHLORENE_STRUCT.moncrief_gw_num_radii; \
  enum { \
      dummy_PRIVATE_BBHLORENE_STRUCT_each_black_hole_mass = sizeof( each_black_hole_mass ) \
    , dummy_PRIVATE_BBHLORENE_STRUCT_each_black_hole_momentum = sizeof( each_black_hole_momentum ) \
    , dummy_PRIVATE_BBHLORENE_STRUCT_each_black_hole_spin = sizeof( each_black_hole_spin ) \
    , dummy_PRIVATE_BBHLORENE_STRUCT_excis_radius = sizeof( excis_radius ) \
    , dummy_PRIVATE_BBHLORENE_STRUCT_moncrief_radius_GW = sizeof( moncrief_radius_GW ) \
    , dummy_PRIVATE_BBHLORENE_STRUCT_total_binary_separation = sizeof( total_binary_separation ) \
    , dummy_PRIVATE_BBHLORENE_STRUCT_fill_excision_enable = sizeof( fill_excision_enable ) \
    , dummy_PRIVATE_BBHLORENE_STRUCT_genID_cmdline_output_enable = sizeof( genID_cmdline_output_enable ) \
    , dummy_PRIVATE_BBHLORENE_STRUCT_moncrief_gw_num_radii = sizeof( moncrief_gw_num_radii ) \
  }; \

