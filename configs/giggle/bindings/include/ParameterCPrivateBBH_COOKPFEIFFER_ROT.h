#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL BigM;
  CCTK_REAL excis_radius;
  CCTK_REAL moncrief_radius_GW[100];
  CCTK_REAL x_offset;
  CCTK_INT fill_excision_enable;
  CCTK_INT moncrief_gw_num_radii;
  CCTK_INT rotate_refinement_boxes;
  CCTK_INT use_cookpfeiffer_lapseshift;
} PRIVATE_BBH_COOKPFEIFFER_ROT_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_BBH_COOKPFEIFFER_ROT_STRUCT_PARAMS \
  CCTK_REAL const BigM = PRIVATE_BBH_COOKPFEIFFER_ROT_STRUCT.BigM; \
  CCTK_REAL const excis_radius = PRIVATE_BBH_COOKPFEIFFER_ROT_STRUCT.excis_radius; \
  CCTK_REAL const * const moncrief_radius_GW = PRIVATE_BBH_COOKPFEIFFER_ROT_STRUCT.moncrief_radius_GW; \
  CCTK_REAL const x_offset = PRIVATE_BBH_COOKPFEIFFER_ROT_STRUCT.x_offset; \
  CCTK_INT const fill_excision_enable = PRIVATE_BBH_COOKPFEIFFER_ROT_STRUCT.fill_excision_enable; \
  CCTK_INT const moncrief_gw_num_radii = PRIVATE_BBH_COOKPFEIFFER_ROT_STRUCT.moncrief_gw_num_radii; \
  CCTK_INT const rotate_refinement_boxes = PRIVATE_BBH_COOKPFEIFFER_ROT_STRUCT.rotate_refinement_boxes; \
  CCTK_INT const use_cookpfeiffer_lapseshift = PRIVATE_BBH_COOKPFEIFFER_ROT_STRUCT.use_cookpfeiffer_lapseshift; \
  enum { \
      dummy_PRIVATE_BBH_COOKPFEIFFER_ROT_STRUCT_BigM = sizeof( BigM ) \
    , dummy_PRIVATE_BBH_COOKPFEIFFER_ROT_STRUCT_excis_radius = sizeof( excis_radius ) \
    , dummy_PRIVATE_BBH_COOKPFEIFFER_ROT_STRUCT_moncrief_radius_GW = sizeof( moncrief_radius_GW ) \
    , dummy_PRIVATE_BBH_COOKPFEIFFER_ROT_STRUCT_x_offset = sizeof( x_offset ) \
    , dummy_PRIVATE_BBH_COOKPFEIFFER_ROT_STRUCT_fill_excision_enable = sizeof( fill_excision_enable ) \
    , dummy_PRIVATE_BBH_COOKPFEIFFER_ROT_STRUCT_moncrief_gw_num_radii = sizeof( moncrief_gw_num_radii ) \
    , dummy_PRIVATE_BBH_COOKPFEIFFER_ROT_STRUCT_rotate_refinement_boxes = sizeof( rotate_refinement_boxes ) \
    , dummy_PRIVATE_BBH_COOKPFEIFFER_ROT_STRUCT_use_cookpfeiffer_lapseshift = sizeof( use_cookpfeiffer_lapseshift ) \
  }; \

