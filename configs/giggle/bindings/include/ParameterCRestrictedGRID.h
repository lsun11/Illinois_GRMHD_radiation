#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL dx;
  CCTK_REAL dxyz;
  CCTK_REAL dy;
  CCTK_REAL dz;
  CCTK_REAL xmax;
  CCTK_REAL xmin;
  CCTK_REAL xyzmax;
  CCTK_REAL xyzmin;
  CCTK_REAL ymax;
  CCTK_REAL ymin;
  CCTK_REAL zmax;
  CCTK_REAL zmin;
  const char * bitant_plane;
  const char * domain;
  const char * quadrant_direction;
  const char * rotation_axis;
  const char * type;
  CCTK_INT symmetry_xmax;
  CCTK_INT symmetry_xmin;
  CCTK_INT symmetry_ymax;
  CCTK_INT symmetry_ymin;
  CCTK_INT symmetry_zmax;
  CCTK_INT symmetry_zmin;
} RESTRICTED_GRID_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_RESTRICTED_GRID_STRUCT_PARAMS \
  CCTK_REAL const dx = RESTRICTED_GRID_STRUCT.dx; \
  CCTK_REAL const dxyz = RESTRICTED_GRID_STRUCT.dxyz; \
  CCTK_REAL const dy = RESTRICTED_GRID_STRUCT.dy; \
  CCTK_REAL const dz = RESTRICTED_GRID_STRUCT.dz; \
  CCTK_REAL const xmax = RESTRICTED_GRID_STRUCT.xmax; \
  CCTK_REAL const xmin = RESTRICTED_GRID_STRUCT.xmin; \
  CCTK_REAL const xyzmax = RESTRICTED_GRID_STRUCT.xyzmax; \
  CCTK_REAL const xyzmin = RESTRICTED_GRID_STRUCT.xyzmin; \
  CCTK_REAL const ymax = RESTRICTED_GRID_STRUCT.ymax; \
  CCTK_REAL const ymin = RESTRICTED_GRID_STRUCT.ymin; \
  CCTK_REAL const zmax = RESTRICTED_GRID_STRUCT.zmax; \
  CCTK_REAL const zmin = RESTRICTED_GRID_STRUCT.zmin; \
  const char * const bitant_plane = RESTRICTED_GRID_STRUCT.bitant_plane; \
  const char * const domain = RESTRICTED_GRID_STRUCT.domain; \
  const char * const quadrant_direction = RESTRICTED_GRID_STRUCT.quadrant_direction; \
  const char * const rotation_axis = RESTRICTED_GRID_STRUCT.rotation_axis; \
  const char * const type = RESTRICTED_GRID_STRUCT.type; \
  CCTK_INT const symmetry_xmax = RESTRICTED_GRID_STRUCT.symmetry_xmax; \
  CCTK_INT const symmetry_xmin = RESTRICTED_GRID_STRUCT.symmetry_xmin; \
  CCTK_INT const symmetry_ymax = RESTRICTED_GRID_STRUCT.symmetry_ymax; \
  CCTK_INT const symmetry_ymin = RESTRICTED_GRID_STRUCT.symmetry_ymin; \
  CCTK_INT const symmetry_zmax = RESTRICTED_GRID_STRUCT.symmetry_zmax; \
  CCTK_INT const symmetry_zmin = RESTRICTED_GRID_STRUCT.symmetry_zmin; \
  enum { \
      dummy_RESTRICTED_GRID_STRUCT_dx = sizeof( dx ) \
    , dummy_RESTRICTED_GRID_STRUCT_dxyz = sizeof( dxyz ) \
    , dummy_RESTRICTED_GRID_STRUCT_dy = sizeof( dy ) \
    , dummy_RESTRICTED_GRID_STRUCT_dz = sizeof( dz ) \
    , dummy_RESTRICTED_GRID_STRUCT_xmax = sizeof( xmax ) \
    , dummy_RESTRICTED_GRID_STRUCT_xmin = sizeof( xmin ) \
    , dummy_RESTRICTED_GRID_STRUCT_xyzmax = sizeof( xyzmax ) \
    , dummy_RESTRICTED_GRID_STRUCT_xyzmin = sizeof( xyzmin ) \
    , dummy_RESTRICTED_GRID_STRUCT_ymax = sizeof( ymax ) \
    , dummy_RESTRICTED_GRID_STRUCT_ymin = sizeof( ymin ) \
    , dummy_RESTRICTED_GRID_STRUCT_zmax = sizeof( zmax ) \
    , dummy_RESTRICTED_GRID_STRUCT_zmin = sizeof( zmin ) \
    , dummy_RESTRICTED_GRID_STRUCT_bitant_plane = sizeof( bitant_plane ) \
    , dummy_RESTRICTED_GRID_STRUCT_domain = sizeof( domain ) \
    , dummy_RESTRICTED_GRID_STRUCT_quadrant_direction = sizeof( quadrant_direction ) \
    , dummy_RESTRICTED_GRID_STRUCT_rotation_axis = sizeof( rotation_axis ) \
    , dummy_RESTRICTED_GRID_STRUCT_type = sizeof( type ) \
    , dummy_RESTRICTED_GRID_STRUCT_symmetry_xmax = sizeof( symmetry_xmax ) \
    , dummy_RESTRICTED_GRID_STRUCT_symmetry_xmin = sizeof( symmetry_xmin ) \
    , dummy_RESTRICTED_GRID_STRUCT_symmetry_ymax = sizeof( symmetry_ymax ) \
    , dummy_RESTRICTED_GRID_STRUCT_symmetry_ymin = sizeof( symmetry_ymin ) \
    , dummy_RESTRICTED_GRID_STRUCT_symmetry_zmax = sizeof( symmetry_zmax ) \
    , dummy_RESTRICTED_GRID_STRUCT_symmetry_zmin = sizeof( symmetry_zmin ) \
  }; \

