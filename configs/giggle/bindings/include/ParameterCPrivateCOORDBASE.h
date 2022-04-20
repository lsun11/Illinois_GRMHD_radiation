#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL dx;
  CCTK_REAL dy;
  CCTK_REAL dz;
  CCTK_REAL xextent;
  CCTK_REAL xmax;
  CCTK_REAL xmin;
  CCTK_REAL yextent;
  CCTK_REAL ymax;
  CCTK_REAL ymin;
  CCTK_REAL zextent;
  CCTK_REAL zmax;
  CCTK_REAL zmin;
  const char * domainsize;
  const char * spacing;
  CCTK_INT boundary_internal_x_lower;
  CCTK_INT boundary_internal_x_upper;
  CCTK_INT boundary_internal_y_lower;
  CCTK_INT boundary_internal_y_upper;
  CCTK_INT boundary_internal_z_lower;
  CCTK_INT boundary_internal_z_upper;
  CCTK_INT boundary_shiftout_x_lower;
  CCTK_INT boundary_shiftout_x_upper;
  CCTK_INT boundary_shiftout_y_lower;
  CCTK_INT boundary_shiftout_y_upper;
  CCTK_INT boundary_shiftout_z_lower;
  CCTK_INT boundary_shiftout_z_upper;
  CCTK_INT boundary_size_x_lower;
  CCTK_INT boundary_size_x_upper;
  CCTK_INT boundary_size_y_lower;
  CCTK_INT boundary_size_y_upper;
  CCTK_INT boundary_size_z_lower;
  CCTK_INT boundary_size_z_upper;
  CCTK_INT boundary_staggered_x_lower;
  CCTK_INT boundary_staggered_x_upper;
  CCTK_INT boundary_staggered_y_lower;
  CCTK_INT boundary_staggered_y_upper;
  CCTK_INT boundary_staggered_z_lower;
  CCTK_INT boundary_staggered_z_upper;
  CCTK_INT ncells_x;
  CCTK_INT ncells_y;
  CCTK_INT ncells_z;
  CCTK_INT zero_origin_x;
  CCTK_INT zero_origin_y;
  CCTK_INT zero_origin_z;
} PRIVATE_COORDBASE_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_COORDBASE_STRUCT_PARAMS \
  CCTK_REAL const dx = PRIVATE_COORDBASE_STRUCT.dx; \
  CCTK_REAL const dy = PRIVATE_COORDBASE_STRUCT.dy; \
  CCTK_REAL const dz = PRIVATE_COORDBASE_STRUCT.dz; \
  CCTK_REAL const xextent = PRIVATE_COORDBASE_STRUCT.xextent; \
  CCTK_REAL const xmax = PRIVATE_COORDBASE_STRUCT.xmax; \
  CCTK_REAL const xmin = PRIVATE_COORDBASE_STRUCT.xmin; \
  CCTK_REAL const yextent = PRIVATE_COORDBASE_STRUCT.yextent; \
  CCTK_REAL const ymax = PRIVATE_COORDBASE_STRUCT.ymax; \
  CCTK_REAL const ymin = PRIVATE_COORDBASE_STRUCT.ymin; \
  CCTK_REAL const zextent = PRIVATE_COORDBASE_STRUCT.zextent; \
  CCTK_REAL const zmax = PRIVATE_COORDBASE_STRUCT.zmax; \
  CCTK_REAL const zmin = PRIVATE_COORDBASE_STRUCT.zmin; \
  const char * const domainsize = PRIVATE_COORDBASE_STRUCT.domainsize; \
  const char * const spacing = PRIVATE_COORDBASE_STRUCT.spacing; \
  CCTK_INT const boundary_internal_x_lower = PRIVATE_COORDBASE_STRUCT.boundary_internal_x_lower; \
  CCTK_INT const boundary_internal_x_upper = PRIVATE_COORDBASE_STRUCT.boundary_internal_x_upper; \
  CCTK_INT const boundary_internal_y_lower = PRIVATE_COORDBASE_STRUCT.boundary_internal_y_lower; \
  CCTK_INT const boundary_internal_y_upper = PRIVATE_COORDBASE_STRUCT.boundary_internal_y_upper; \
  CCTK_INT const boundary_internal_z_lower = PRIVATE_COORDBASE_STRUCT.boundary_internal_z_lower; \
  CCTK_INT const boundary_internal_z_upper = PRIVATE_COORDBASE_STRUCT.boundary_internal_z_upper; \
  CCTK_INT const boundary_shiftout_x_lower = PRIVATE_COORDBASE_STRUCT.boundary_shiftout_x_lower; \
  CCTK_INT const boundary_shiftout_x_upper = PRIVATE_COORDBASE_STRUCT.boundary_shiftout_x_upper; \
  CCTK_INT const boundary_shiftout_y_lower = PRIVATE_COORDBASE_STRUCT.boundary_shiftout_y_lower; \
  CCTK_INT const boundary_shiftout_y_upper = PRIVATE_COORDBASE_STRUCT.boundary_shiftout_y_upper; \
  CCTK_INT const boundary_shiftout_z_lower = PRIVATE_COORDBASE_STRUCT.boundary_shiftout_z_lower; \
  CCTK_INT const boundary_shiftout_z_upper = PRIVATE_COORDBASE_STRUCT.boundary_shiftout_z_upper; \
  CCTK_INT const boundary_size_x_lower = PRIVATE_COORDBASE_STRUCT.boundary_size_x_lower; \
  CCTK_INT const boundary_size_x_upper = PRIVATE_COORDBASE_STRUCT.boundary_size_x_upper; \
  CCTK_INT const boundary_size_y_lower = PRIVATE_COORDBASE_STRUCT.boundary_size_y_lower; \
  CCTK_INT const boundary_size_y_upper = PRIVATE_COORDBASE_STRUCT.boundary_size_y_upper; \
  CCTK_INT const boundary_size_z_lower = PRIVATE_COORDBASE_STRUCT.boundary_size_z_lower; \
  CCTK_INT const boundary_size_z_upper = PRIVATE_COORDBASE_STRUCT.boundary_size_z_upper; \
  CCTK_INT const boundary_staggered_x_lower = PRIVATE_COORDBASE_STRUCT.boundary_staggered_x_lower; \
  CCTK_INT const boundary_staggered_x_upper = PRIVATE_COORDBASE_STRUCT.boundary_staggered_x_upper; \
  CCTK_INT const boundary_staggered_y_lower = PRIVATE_COORDBASE_STRUCT.boundary_staggered_y_lower; \
  CCTK_INT const boundary_staggered_y_upper = PRIVATE_COORDBASE_STRUCT.boundary_staggered_y_upper; \
  CCTK_INT const boundary_staggered_z_lower = PRIVATE_COORDBASE_STRUCT.boundary_staggered_z_lower; \
  CCTK_INT const boundary_staggered_z_upper = PRIVATE_COORDBASE_STRUCT.boundary_staggered_z_upper; \
  CCTK_INT const ncells_x = PRIVATE_COORDBASE_STRUCT.ncells_x; \
  CCTK_INT const ncells_y = PRIVATE_COORDBASE_STRUCT.ncells_y; \
  CCTK_INT const ncells_z = PRIVATE_COORDBASE_STRUCT.ncells_z; \
  CCTK_INT const zero_origin_x = PRIVATE_COORDBASE_STRUCT.zero_origin_x; \
  CCTK_INT const zero_origin_y = PRIVATE_COORDBASE_STRUCT.zero_origin_y; \
  CCTK_INT const zero_origin_z = PRIVATE_COORDBASE_STRUCT.zero_origin_z; \
  enum { \
      dummy_PRIVATE_COORDBASE_STRUCT_dx = sizeof( dx ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_dy = sizeof( dy ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_dz = sizeof( dz ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_xextent = sizeof( xextent ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_xmax = sizeof( xmax ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_xmin = sizeof( xmin ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_yextent = sizeof( yextent ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_ymax = sizeof( ymax ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_ymin = sizeof( ymin ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_zextent = sizeof( zextent ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_zmax = sizeof( zmax ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_zmin = sizeof( zmin ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_domainsize = sizeof( domainsize ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_spacing = sizeof( spacing ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_boundary_internal_x_lower = sizeof( boundary_internal_x_lower ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_boundary_internal_x_upper = sizeof( boundary_internal_x_upper ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_boundary_internal_y_lower = sizeof( boundary_internal_y_lower ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_boundary_internal_y_upper = sizeof( boundary_internal_y_upper ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_boundary_internal_z_lower = sizeof( boundary_internal_z_lower ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_boundary_internal_z_upper = sizeof( boundary_internal_z_upper ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_boundary_shiftout_x_lower = sizeof( boundary_shiftout_x_lower ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_boundary_shiftout_x_upper = sizeof( boundary_shiftout_x_upper ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_boundary_shiftout_y_lower = sizeof( boundary_shiftout_y_lower ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_boundary_shiftout_y_upper = sizeof( boundary_shiftout_y_upper ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_boundary_shiftout_z_lower = sizeof( boundary_shiftout_z_lower ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_boundary_shiftout_z_upper = sizeof( boundary_shiftout_z_upper ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_boundary_size_x_lower = sizeof( boundary_size_x_lower ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_boundary_size_x_upper = sizeof( boundary_size_x_upper ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_boundary_size_y_lower = sizeof( boundary_size_y_lower ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_boundary_size_y_upper = sizeof( boundary_size_y_upper ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_boundary_size_z_lower = sizeof( boundary_size_z_lower ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_boundary_size_z_upper = sizeof( boundary_size_z_upper ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_boundary_staggered_x_lower = sizeof( boundary_staggered_x_lower ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_boundary_staggered_x_upper = sizeof( boundary_staggered_x_upper ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_boundary_staggered_y_lower = sizeof( boundary_staggered_y_lower ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_boundary_staggered_y_upper = sizeof( boundary_staggered_y_upper ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_boundary_staggered_z_lower = sizeof( boundary_staggered_z_lower ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_boundary_staggered_z_upper = sizeof( boundary_staggered_z_upper ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_ncells_x = sizeof( ncells_x ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_ncells_y = sizeof( ncells_y ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_ncells_z = sizeof( ncells_z ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_zero_origin_x = sizeof( zero_origin_x ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_zero_origin_y = sizeof( zero_origin_y ) \
    , dummy_PRIVATE_COORDBASE_STRUCT_zero_origin_z = sizeof( zero_origin_z ) \
  }; \

