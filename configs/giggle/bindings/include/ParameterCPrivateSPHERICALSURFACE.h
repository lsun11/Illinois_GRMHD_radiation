#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL origin_x[42];
  CCTK_REAL origin_y[42];
  CCTK_REAL origin_z[42];
  CCTK_REAL radius[42];
  CCTK_REAL radius_x[42];
  CCTK_REAL radius_y[42];
  CCTK_REAL radius_z[42];
  CCTK_INT set_elliptic[42];
  CCTK_INT set_spherical[42];
} PRIVATE_SPHERICALSURFACE_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_SPHERICALSURFACE_STRUCT_PARAMS \
  CCTK_REAL const * const origin_x = PRIVATE_SPHERICALSURFACE_STRUCT.origin_x; \
  CCTK_REAL const * const origin_y = PRIVATE_SPHERICALSURFACE_STRUCT.origin_y; \
  CCTK_REAL const * const origin_z = PRIVATE_SPHERICALSURFACE_STRUCT.origin_z; \
  CCTK_REAL const * const radius = PRIVATE_SPHERICALSURFACE_STRUCT.radius; \
  CCTK_REAL const * const radius_x = PRIVATE_SPHERICALSURFACE_STRUCT.radius_x; \
  CCTK_REAL const * const radius_y = PRIVATE_SPHERICALSURFACE_STRUCT.radius_y; \
  CCTK_REAL const * const radius_z = PRIVATE_SPHERICALSURFACE_STRUCT.radius_z; \
  CCTK_INT const * const set_elliptic = PRIVATE_SPHERICALSURFACE_STRUCT.set_elliptic; \
  CCTK_INT const * const set_spherical = PRIVATE_SPHERICALSURFACE_STRUCT.set_spherical; \
  enum { \
      dummy_PRIVATE_SPHERICALSURFACE_STRUCT_origin_x = sizeof( origin_x ) \
    , dummy_PRIVATE_SPHERICALSURFACE_STRUCT_origin_y = sizeof( origin_y ) \
    , dummy_PRIVATE_SPHERICALSURFACE_STRUCT_origin_z = sizeof( origin_z ) \
    , dummy_PRIVATE_SPHERICALSURFACE_STRUCT_radius = sizeof( radius ) \
    , dummy_PRIVATE_SPHERICALSURFACE_STRUCT_radius_x = sizeof( radius_x ) \
    , dummy_PRIVATE_SPHERICALSURFACE_STRUCT_radius_y = sizeof( radius_y ) \
    , dummy_PRIVATE_SPHERICALSURFACE_STRUCT_radius_z = sizeof( radius_z ) \
    , dummy_PRIVATE_SPHERICALSURFACE_STRUCT_set_elliptic = sizeof( set_elliptic ) \
    , dummy_PRIVATE_SPHERICALSURFACE_STRUCT_set_spherical = sizeof( set_spherical ) \
  }; \

