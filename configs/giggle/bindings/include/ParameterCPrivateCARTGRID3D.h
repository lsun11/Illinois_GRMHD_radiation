#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  const char * set_coordinate_ranges_on;
  CCTK_INT avoid_origin;
  CCTK_INT avoid_originx;
  CCTK_INT avoid_originy;
  CCTK_INT avoid_originz;
  CCTK_INT no_origin;
  CCTK_INT no_originx;
  CCTK_INT no_originy;
  CCTK_INT no_originz;
  CCTK_INT register_default_coordinate_systems;
} PRIVATE_CARTGRID3D_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_CARTGRID3D_STRUCT_PARAMS \
  const char * const set_coordinate_ranges_on = PRIVATE_CARTGRID3D_STRUCT.set_coordinate_ranges_on; \
  CCTK_INT const avoid_origin = PRIVATE_CARTGRID3D_STRUCT.avoid_origin; \
  CCTK_INT const avoid_originx = PRIVATE_CARTGRID3D_STRUCT.avoid_originx; \
  CCTK_INT const avoid_originy = PRIVATE_CARTGRID3D_STRUCT.avoid_originy; \
  CCTK_INT const avoid_originz = PRIVATE_CARTGRID3D_STRUCT.avoid_originz; \
  CCTK_INT const no_origin = PRIVATE_CARTGRID3D_STRUCT.no_origin; \
  CCTK_INT const no_originx = PRIVATE_CARTGRID3D_STRUCT.no_originx; \
  CCTK_INT const no_originy = PRIVATE_CARTGRID3D_STRUCT.no_originy; \
  CCTK_INT const no_originz = PRIVATE_CARTGRID3D_STRUCT.no_originz; \
  CCTK_INT const register_default_coordinate_systems = PRIVATE_CARTGRID3D_STRUCT.register_default_coordinate_systems; \
  enum { \
      dummy_PRIVATE_CARTGRID3D_STRUCT_set_coordinate_ranges_on = sizeof( set_coordinate_ranges_on ) \
    , dummy_PRIVATE_CARTGRID3D_STRUCT_avoid_origin = sizeof( avoid_origin ) \
    , dummy_PRIVATE_CARTGRID3D_STRUCT_avoid_originx = sizeof( avoid_originx ) \
    , dummy_PRIVATE_CARTGRID3D_STRUCT_avoid_originy = sizeof( avoid_originy ) \
    , dummy_PRIVATE_CARTGRID3D_STRUCT_avoid_originz = sizeof( avoid_originz ) \
    , dummy_PRIVATE_CARTGRID3D_STRUCT_no_origin = sizeof( no_origin ) \
    , dummy_PRIVATE_CARTGRID3D_STRUCT_no_originx = sizeof( no_originx ) \
    , dummy_PRIVATE_CARTGRID3D_STRUCT_no_originy = sizeof( no_originy ) \
    , dummy_PRIVATE_CARTGRID3D_STRUCT_no_originz = sizeof( no_originz ) \
    , dummy_PRIVATE_CARTGRID3D_STRUCT_register_default_coordinate_systems = sizeof( register_default_coordinate_systems ) \
  }; \

