#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL auto_res_ratio[42];
  const char * auto_res_grid[42];
  CCTK_INT auto_res[42];
  CCTK_INT maxnphi;
  CCTK_INT maxntheta;
  CCTK_INT nghostsphi[42];
  CCTK_INT nghoststheta[42];
  CCTK_INT nphi[42];
  CCTK_INT nsurfaces;
  CCTK_INT ntheta[42];
  CCTK_INT symmetric_x[42];
  CCTK_INT symmetric_y[42];
  CCTK_INT symmetric_z[42];
  CCTK_INT verbose;
} RESTRICTED_SPHERICALSURFACE_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_RESTRICTED_SPHERICALSURFACE_STRUCT_PARAMS \
  CCTK_REAL const * const auto_res_ratio = RESTRICTED_SPHERICALSURFACE_STRUCT.auto_res_ratio; \
  const char * const * const auto_res_grid = RESTRICTED_SPHERICALSURFACE_STRUCT.auto_res_grid; \
  CCTK_INT const * const auto_res = RESTRICTED_SPHERICALSURFACE_STRUCT.auto_res; \
  CCTK_INT const maxnphi = RESTRICTED_SPHERICALSURFACE_STRUCT.maxnphi; \
  CCTK_INT const maxntheta = RESTRICTED_SPHERICALSURFACE_STRUCT.maxntheta; \
  CCTK_INT const * const nghostsphi = RESTRICTED_SPHERICALSURFACE_STRUCT.nghostsphi; \
  CCTK_INT const * const nghoststheta = RESTRICTED_SPHERICALSURFACE_STRUCT.nghoststheta; \
  CCTK_INT const * const nphi = RESTRICTED_SPHERICALSURFACE_STRUCT.nphi; \
  CCTK_INT const nsurfaces = RESTRICTED_SPHERICALSURFACE_STRUCT.nsurfaces; \
  CCTK_INT const * const ntheta = RESTRICTED_SPHERICALSURFACE_STRUCT.ntheta; \
  CCTK_INT const * const symmetric_x = RESTRICTED_SPHERICALSURFACE_STRUCT.symmetric_x; \
  CCTK_INT const * const symmetric_y = RESTRICTED_SPHERICALSURFACE_STRUCT.symmetric_y; \
  CCTK_INT const * const symmetric_z = RESTRICTED_SPHERICALSURFACE_STRUCT.symmetric_z; \
  CCTK_INT const verbose = RESTRICTED_SPHERICALSURFACE_STRUCT.verbose; \
  enum { \
      dummy_RESTRICTED_SPHERICALSURFACE_STRUCT_auto_res_ratio = sizeof( auto_res_ratio ) \
    , dummy_RESTRICTED_SPHERICALSURFACE_STRUCT_auto_res_grid = sizeof( auto_res_grid ) \
    , dummy_RESTRICTED_SPHERICALSURFACE_STRUCT_auto_res = sizeof( auto_res ) \
    , dummy_RESTRICTED_SPHERICALSURFACE_STRUCT_maxnphi = sizeof( maxnphi ) \
    , dummy_RESTRICTED_SPHERICALSURFACE_STRUCT_maxntheta = sizeof( maxntheta ) \
    , dummy_RESTRICTED_SPHERICALSURFACE_STRUCT_nghostsphi = sizeof( nghostsphi ) \
    , dummy_RESTRICTED_SPHERICALSURFACE_STRUCT_nghoststheta = sizeof( nghoststheta ) \
    , dummy_RESTRICTED_SPHERICALSURFACE_STRUCT_nphi = sizeof( nphi ) \
    , dummy_RESTRICTED_SPHERICALSURFACE_STRUCT_nsurfaces = sizeof( nsurfaces ) \
    , dummy_RESTRICTED_SPHERICALSURFACE_STRUCT_ntheta = sizeof( ntheta ) \
    , dummy_RESTRICTED_SPHERICALSURFACE_STRUCT_symmetric_x = sizeof( symmetric_x ) \
    , dummy_RESTRICTED_SPHERICALSURFACE_STRUCT_symmetric_y = sizeof( symmetric_y ) \
    , dummy_RESTRICTED_SPHERICALSURFACE_STRUCT_symmetric_z = sizeof( symmetric_z ) \
    , dummy_RESTRICTED_SPHERICALSURFACE_STRUCT_verbose = sizeof( verbose ) \
  }; \

