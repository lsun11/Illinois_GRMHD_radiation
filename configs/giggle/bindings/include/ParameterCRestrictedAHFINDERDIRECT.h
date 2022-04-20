#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL initial_guess__Kerr_KerrSchild__mass[101];
  CCTK_REAL initial_guess__Kerr_KerrSchild__spin[101];
  CCTK_REAL initial_guess__Kerr_KerrSchild__x_posn[101];
  CCTK_REAL initial_guess__Kerr_KerrSchild__y_posn[101];
  CCTK_REAL initial_guess__Kerr_KerrSchild__z_posn[101];
  CCTK_REAL initial_guess__Kerr_Kerr__mass[101];
  CCTK_REAL initial_guess__Kerr_Kerr__spin[101];
  CCTK_REAL initial_guess__Kerr_Kerr__x_posn[101];
  CCTK_REAL initial_guess__Kerr_Kerr__y_posn[101];
  CCTK_REAL initial_guess__Kerr_Kerr__z_posn[101];
  CCTK_REAL initial_guess__coord_ellipsoid__x_center[101];
  CCTK_REAL initial_guess__coord_ellipsoid__x_radius[101];
  CCTK_REAL initial_guess__coord_ellipsoid__y_center[101];
  CCTK_REAL initial_guess__coord_ellipsoid__y_radius[101];
  CCTK_REAL initial_guess__coord_ellipsoid__z_center[101];
  CCTK_REAL initial_guess__coord_ellipsoid__z_radius[101];
  CCTK_REAL initial_guess__coord_sphere__radius[101];
  CCTK_REAL initial_guess__coord_sphere__x_center[101];
  CCTK_REAL initial_guess__coord_sphere__y_center[101];
  CCTK_REAL initial_guess__coord_sphere__z_center[101];
  CCTK_REAL origin_x[101];
  CCTK_REAL origin_y[101];
  CCTK_REAL origin_z[101];
  const char * initial_guess__read_from_named_file__file_name[101];
  const char * initial_guess_method[101];
  CCTK_INT move_origins;
  CCTK_INT predict_origin_movement;
  CCTK_INT reset_horizon_after_not_finding[101];
  CCTK_INT reshape_while_moving;
} RESTRICTED_AHFINDERDIRECT_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_RESTRICTED_AHFINDERDIRECT_STRUCT_PARAMS \
  CCTK_REAL const * const initial_guess__Kerr_KerrSchild__mass = RESTRICTED_AHFINDERDIRECT_STRUCT.initial_guess__Kerr_KerrSchild__mass; \
  CCTK_REAL const * const initial_guess__Kerr_KerrSchild__spin = RESTRICTED_AHFINDERDIRECT_STRUCT.initial_guess__Kerr_KerrSchild__spin; \
  CCTK_REAL const * const initial_guess__Kerr_KerrSchild__x_posn = RESTRICTED_AHFINDERDIRECT_STRUCT.initial_guess__Kerr_KerrSchild__x_posn; \
  CCTK_REAL const * const initial_guess__Kerr_KerrSchild__y_posn = RESTRICTED_AHFINDERDIRECT_STRUCT.initial_guess__Kerr_KerrSchild__y_posn; \
  CCTK_REAL const * const initial_guess__Kerr_KerrSchild__z_posn = RESTRICTED_AHFINDERDIRECT_STRUCT.initial_guess__Kerr_KerrSchild__z_posn; \
  CCTK_REAL const * const initial_guess__Kerr_Kerr__mass = RESTRICTED_AHFINDERDIRECT_STRUCT.initial_guess__Kerr_Kerr__mass; \
  CCTK_REAL const * const initial_guess__Kerr_Kerr__spin = RESTRICTED_AHFINDERDIRECT_STRUCT.initial_guess__Kerr_Kerr__spin; \
  CCTK_REAL const * const initial_guess__Kerr_Kerr__x_posn = RESTRICTED_AHFINDERDIRECT_STRUCT.initial_guess__Kerr_Kerr__x_posn; \
  CCTK_REAL const * const initial_guess__Kerr_Kerr__y_posn = RESTRICTED_AHFINDERDIRECT_STRUCT.initial_guess__Kerr_Kerr__y_posn; \
  CCTK_REAL const * const initial_guess__Kerr_Kerr__z_posn = RESTRICTED_AHFINDERDIRECT_STRUCT.initial_guess__Kerr_Kerr__z_posn; \
  CCTK_REAL const * const initial_guess__coord_ellipsoid__x_center = RESTRICTED_AHFINDERDIRECT_STRUCT.initial_guess__coord_ellipsoid__x_center; \
  CCTK_REAL const * const initial_guess__coord_ellipsoid__x_radius = RESTRICTED_AHFINDERDIRECT_STRUCT.initial_guess__coord_ellipsoid__x_radius; \
  CCTK_REAL const * const initial_guess__coord_ellipsoid__y_center = RESTRICTED_AHFINDERDIRECT_STRUCT.initial_guess__coord_ellipsoid__y_center; \
  CCTK_REAL const * const initial_guess__coord_ellipsoid__y_radius = RESTRICTED_AHFINDERDIRECT_STRUCT.initial_guess__coord_ellipsoid__y_radius; \
  CCTK_REAL const * const initial_guess__coord_ellipsoid__z_center = RESTRICTED_AHFINDERDIRECT_STRUCT.initial_guess__coord_ellipsoid__z_center; \
  CCTK_REAL const * const initial_guess__coord_ellipsoid__z_radius = RESTRICTED_AHFINDERDIRECT_STRUCT.initial_guess__coord_ellipsoid__z_radius; \
  CCTK_REAL const * const initial_guess__coord_sphere__radius = RESTRICTED_AHFINDERDIRECT_STRUCT.initial_guess__coord_sphere__radius; \
  CCTK_REAL const * const initial_guess__coord_sphere__x_center = RESTRICTED_AHFINDERDIRECT_STRUCT.initial_guess__coord_sphere__x_center; \
  CCTK_REAL const * const initial_guess__coord_sphere__y_center = RESTRICTED_AHFINDERDIRECT_STRUCT.initial_guess__coord_sphere__y_center; \
  CCTK_REAL const * const initial_guess__coord_sphere__z_center = RESTRICTED_AHFINDERDIRECT_STRUCT.initial_guess__coord_sphere__z_center; \
  CCTK_REAL const * const origin_x = RESTRICTED_AHFINDERDIRECT_STRUCT.origin_x; \
  CCTK_REAL const * const origin_y = RESTRICTED_AHFINDERDIRECT_STRUCT.origin_y; \
  CCTK_REAL const * const origin_z = RESTRICTED_AHFINDERDIRECT_STRUCT.origin_z; \
  const char * const * const initial_guess__read_from_named_file__file_name = RESTRICTED_AHFINDERDIRECT_STRUCT.initial_guess__read_from_named_file__file_name; \
  const char * const * const initial_guess_method = RESTRICTED_AHFINDERDIRECT_STRUCT.initial_guess_method; \
  CCTK_INT const move_origins = RESTRICTED_AHFINDERDIRECT_STRUCT.move_origins; \
  CCTK_INT const predict_origin_movement = RESTRICTED_AHFINDERDIRECT_STRUCT.predict_origin_movement; \
  CCTK_INT const * const reset_horizon_after_not_finding = RESTRICTED_AHFINDERDIRECT_STRUCT.reset_horizon_after_not_finding; \
  CCTK_INT const reshape_while_moving = RESTRICTED_AHFINDERDIRECT_STRUCT.reshape_while_moving; \
  enum { \
      dummy_RESTRICTED_AHFINDERDIRECT_STRUCT_initial_guess__Kerr_KerrSchild__mass = sizeof( initial_guess__Kerr_KerrSchild__mass ) \
    , dummy_RESTRICTED_AHFINDERDIRECT_STRUCT_initial_guess__Kerr_KerrSchild__spin = sizeof( initial_guess__Kerr_KerrSchild__spin ) \
    , dummy_RESTRICTED_AHFINDERDIRECT_STRUCT_initial_guess__Kerr_KerrSchild__x_posn = sizeof( initial_guess__Kerr_KerrSchild__x_posn ) \
    , dummy_RESTRICTED_AHFINDERDIRECT_STRUCT_initial_guess__Kerr_KerrSchild__y_posn = sizeof( initial_guess__Kerr_KerrSchild__y_posn ) \
    , dummy_RESTRICTED_AHFINDERDIRECT_STRUCT_initial_guess__Kerr_KerrSchild__z_posn = sizeof( initial_guess__Kerr_KerrSchild__z_posn ) \
    , dummy_RESTRICTED_AHFINDERDIRECT_STRUCT_initial_guess__Kerr_Kerr__mass = sizeof( initial_guess__Kerr_Kerr__mass ) \
    , dummy_RESTRICTED_AHFINDERDIRECT_STRUCT_initial_guess__Kerr_Kerr__spin = sizeof( initial_guess__Kerr_Kerr__spin ) \
    , dummy_RESTRICTED_AHFINDERDIRECT_STRUCT_initial_guess__Kerr_Kerr__x_posn = sizeof( initial_guess__Kerr_Kerr__x_posn ) \
    , dummy_RESTRICTED_AHFINDERDIRECT_STRUCT_initial_guess__Kerr_Kerr__y_posn = sizeof( initial_guess__Kerr_Kerr__y_posn ) \
    , dummy_RESTRICTED_AHFINDERDIRECT_STRUCT_initial_guess__Kerr_Kerr__z_posn = sizeof( initial_guess__Kerr_Kerr__z_posn ) \
    , dummy_RESTRICTED_AHFINDERDIRECT_STRUCT_initial_guess__coord_ellipsoid__x_center = sizeof( initial_guess__coord_ellipsoid__x_center ) \
    , dummy_RESTRICTED_AHFINDERDIRECT_STRUCT_initial_guess__coord_ellipsoid__x_radius = sizeof( initial_guess__coord_ellipsoid__x_radius ) \
    , dummy_RESTRICTED_AHFINDERDIRECT_STRUCT_initial_guess__coord_ellipsoid__y_center = sizeof( initial_guess__coord_ellipsoid__y_center ) \
    , dummy_RESTRICTED_AHFINDERDIRECT_STRUCT_initial_guess__coord_ellipsoid__y_radius = sizeof( initial_guess__coord_ellipsoid__y_radius ) \
    , dummy_RESTRICTED_AHFINDERDIRECT_STRUCT_initial_guess__coord_ellipsoid__z_center = sizeof( initial_guess__coord_ellipsoid__z_center ) \
    , dummy_RESTRICTED_AHFINDERDIRECT_STRUCT_initial_guess__coord_ellipsoid__z_radius = sizeof( initial_guess__coord_ellipsoid__z_radius ) \
    , dummy_RESTRICTED_AHFINDERDIRECT_STRUCT_initial_guess__coord_sphere__radius = sizeof( initial_guess__coord_sphere__radius ) \
    , dummy_RESTRICTED_AHFINDERDIRECT_STRUCT_initial_guess__coord_sphere__x_center = sizeof( initial_guess__coord_sphere__x_center ) \
    , dummy_RESTRICTED_AHFINDERDIRECT_STRUCT_initial_guess__coord_sphere__y_center = sizeof( initial_guess__coord_sphere__y_center ) \
    , dummy_RESTRICTED_AHFINDERDIRECT_STRUCT_initial_guess__coord_sphere__z_center = sizeof( initial_guess__coord_sphere__z_center ) \
    , dummy_RESTRICTED_AHFINDERDIRECT_STRUCT_origin_x = sizeof( origin_x ) \
    , dummy_RESTRICTED_AHFINDERDIRECT_STRUCT_origin_y = sizeof( origin_y ) \
    , dummy_RESTRICTED_AHFINDERDIRECT_STRUCT_origin_z = sizeof( origin_z ) \
    , dummy_RESTRICTED_AHFINDERDIRECT_STRUCT_initial_guess__read_from_named_file__file_name = sizeof( initial_guess__read_from_named_file__file_name ) \
    , dummy_RESTRICTED_AHFINDERDIRECT_STRUCT_initial_guess_method = sizeof( initial_guess_method ) \
    , dummy_RESTRICTED_AHFINDERDIRECT_STRUCT_move_origins = sizeof( move_origins ) \
    , dummy_RESTRICTED_AHFINDERDIRECT_STRUCT_predict_origin_movement = sizeof( predict_origin_movement ) \
    , dummy_RESTRICTED_AHFINDERDIRECT_STRUCT_reset_horizon_after_not_finding = sizeof( reset_horizon_after_not_finding ) \
    , dummy_RESTRICTED_AHFINDERDIRECT_STRUCT_reshape_while_moving = sizeof( reshape_while_moving ) \
  }; \

