#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL R_edge;
  CCTK_REAL betam1;
  CCTK_REAL moncrief_radius_GW[11];
  const char * OS_filename;
  CCTK_INT em_field_type;
  CCTK_INT moncrief_gw_num_radii;
  CCTK_INT narr;
  CCTK_INT set_lapse_one;
} PRIVATE_OS_COLLAPSE_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_OS_COLLAPSE_STRUCT_PARAMS \
  CCTK_REAL const R_edge = PRIVATE_OS_COLLAPSE_STRUCT.R_edge; \
  CCTK_REAL const betam1 = PRIVATE_OS_COLLAPSE_STRUCT.betam1; \
  CCTK_REAL const * const moncrief_radius_GW = PRIVATE_OS_COLLAPSE_STRUCT.moncrief_radius_GW; \
  const char * const OS_filename = PRIVATE_OS_COLLAPSE_STRUCT.OS_filename; \
  CCTK_INT const em_field_type = PRIVATE_OS_COLLAPSE_STRUCT.em_field_type; \
  CCTK_INT const moncrief_gw_num_radii = PRIVATE_OS_COLLAPSE_STRUCT.moncrief_gw_num_radii; \
  CCTK_INT const narr = PRIVATE_OS_COLLAPSE_STRUCT.narr; \
  CCTK_INT const set_lapse_one = PRIVATE_OS_COLLAPSE_STRUCT.set_lapse_one; \
  enum { \
      dummy_PRIVATE_OS_COLLAPSE_STRUCT_R_edge = sizeof( R_edge ) \
    , dummy_PRIVATE_OS_COLLAPSE_STRUCT_betam1 = sizeof( betam1 ) \
    , dummy_PRIVATE_OS_COLLAPSE_STRUCT_moncrief_radius_GW = sizeof( moncrief_radius_GW ) \
    , dummy_PRIVATE_OS_COLLAPSE_STRUCT_OS_filename = sizeof( OS_filename ) \
    , dummy_PRIVATE_OS_COLLAPSE_STRUCT_em_field_type = sizeof( em_field_type ) \
    , dummy_PRIVATE_OS_COLLAPSE_STRUCT_moncrief_gw_num_radii = sizeof( moncrief_gw_num_radii ) \
    , dummy_PRIVATE_OS_COLLAPSE_STRUCT_narr = sizeof( narr ) \
    , dummy_PRIVATE_OS_COLLAPSE_STRUCT_set_lapse_one = sizeof( set_lapse_one ) \
  }; \

