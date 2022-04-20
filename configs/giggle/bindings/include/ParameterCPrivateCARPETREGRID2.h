#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL min_fraction;
  CCTK_REAL movement_threshold_1;
  CCTK_REAL movement_threshold_10;
  CCTK_REAL movement_threshold_2;
  CCTK_REAL movement_threshold_3;
  CCTK_REAL movement_threshold_4;
  CCTK_REAL movement_threshold_5;
  CCTK_REAL movement_threshold_6;
  CCTK_REAL movement_threshold_7;
  CCTK_REAL movement_threshold_8;
  CCTK_REAL movement_threshold_9;
  CCTK_REAL position_x_1;
  CCTK_REAL position_x_10;
  CCTK_REAL position_x_2;
  CCTK_REAL position_x_3;
  CCTK_REAL position_x_4;
  CCTK_REAL position_x_5;
  CCTK_REAL position_x_6;
  CCTK_REAL position_x_7;
  CCTK_REAL position_x_8;
  CCTK_REAL position_x_9;
  CCTK_REAL position_y_1;
  CCTK_REAL position_y_10;
  CCTK_REAL position_y_2;
  CCTK_REAL position_y_3;
  CCTK_REAL position_y_4;
  CCTK_REAL position_y_5;
  CCTK_REAL position_y_6;
  CCTK_REAL position_y_7;
  CCTK_REAL position_y_8;
  CCTK_REAL position_y_9;
  CCTK_REAL position_z_1;
  CCTK_REAL position_z_10;
  CCTK_REAL position_z_2;
  CCTK_REAL position_z_3;
  CCTK_REAL position_z_4;
  CCTK_REAL position_z_5;
  CCTK_REAL position_z_6;
  CCTK_REAL position_z_7;
  CCTK_REAL position_z_8;
  CCTK_REAL position_z_9;
  CCTK_REAL radius_1[30];
  CCTK_REAL radius_10[30];
  CCTK_REAL radius_2[30];
  CCTK_REAL radius_3[30];
  CCTK_REAL radius_4[30];
  CCTK_REAL radius_5[30];
  CCTK_REAL radius_6[30];
  CCTK_REAL radius_7[30];
  CCTK_REAL radius_8[30];
  CCTK_REAL radius_9[30];
  CCTK_REAL radius_change_threshold_1;
  CCTK_REAL radius_change_threshold_10;
  CCTK_REAL radius_change_threshold_2;
  CCTK_REAL radius_change_threshold_3;
  CCTK_REAL radius_change_threshold_4;
  CCTK_REAL radius_change_threshold_5;
  CCTK_REAL radius_change_threshold_6;
  CCTK_REAL radius_change_threshold_7;
  CCTK_REAL radius_change_threshold_8;
  CCTK_REAL radius_change_threshold_9;
  CCTK_INT active_1;
  CCTK_INT active_10;
  CCTK_INT active_2;
  CCTK_INT active_3;
  CCTK_INT active_4;
  CCTK_INT active_5;
  CCTK_INT active_6;
  CCTK_INT active_7;
  CCTK_INT active_8;
  CCTK_INT active_9;
  CCTK_INT boundary_shiftout;
  CCTK_INT ensure_proper_nesting;
  CCTK_INT freeze_unaligned_levels;
  CCTK_INT freeze_unaligned_parent_levels;
  CCTK_INT min_distance;
  CCTK_INT num_centres;
  CCTK_INT num_levels_1;
  CCTK_INT num_levels_10;
  CCTK_INT num_levels_2;
  CCTK_INT num_levels_3;
  CCTK_INT num_levels_4;
  CCTK_INT num_levels_5;
  CCTK_INT num_levels_6;
  CCTK_INT num_levels_7;
  CCTK_INT num_levels_8;
  CCTK_INT num_levels_9;
  CCTK_INT regrid_every;
  CCTK_INT snap_to_coarse;
  CCTK_INT symmetry_rotating180;
  CCTK_INT symmetry_rotating90;
  CCTK_INT verbose;
} PRIVATE_CARPETREGRID2_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_CARPETREGRID2_STRUCT_PARAMS \
  CCTK_REAL const min_fraction = PRIVATE_CARPETREGRID2_STRUCT.min_fraction; \
  CCTK_REAL const movement_threshold_1 = PRIVATE_CARPETREGRID2_STRUCT.movement_threshold_1; \
  CCTK_REAL const movement_threshold_10 = PRIVATE_CARPETREGRID2_STRUCT.movement_threshold_10; \
  CCTK_REAL const movement_threshold_2 = PRIVATE_CARPETREGRID2_STRUCT.movement_threshold_2; \
  CCTK_REAL const movement_threshold_3 = PRIVATE_CARPETREGRID2_STRUCT.movement_threshold_3; \
  CCTK_REAL const movement_threshold_4 = PRIVATE_CARPETREGRID2_STRUCT.movement_threshold_4; \
  CCTK_REAL const movement_threshold_5 = PRIVATE_CARPETREGRID2_STRUCT.movement_threshold_5; \
  CCTK_REAL const movement_threshold_6 = PRIVATE_CARPETREGRID2_STRUCT.movement_threshold_6; \
  CCTK_REAL const movement_threshold_7 = PRIVATE_CARPETREGRID2_STRUCT.movement_threshold_7; \
  CCTK_REAL const movement_threshold_8 = PRIVATE_CARPETREGRID2_STRUCT.movement_threshold_8; \
  CCTK_REAL const movement_threshold_9 = PRIVATE_CARPETREGRID2_STRUCT.movement_threshold_9; \
  CCTK_REAL const position_x_1 = PRIVATE_CARPETREGRID2_STRUCT.position_x_1; \
  CCTK_REAL const position_x_10 = PRIVATE_CARPETREGRID2_STRUCT.position_x_10; \
  CCTK_REAL const position_x_2 = PRIVATE_CARPETREGRID2_STRUCT.position_x_2; \
  CCTK_REAL const position_x_3 = PRIVATE_CARPETREGRID2_STRUCT.position_x_3; \
  CCTK_REAL const position_x_4 = PRIVATE_CARPETREGRID2_STRUCT.position_x_4; \
  CCTK_REAL const position_x_5 = PRIVATE_CARPETREGRID2_STRUCT.position_x_5; \
  CCTK_REAL const position_x_6 = PRIVATE_CARPETREGRID2_STRUCT.position_x_6; \
  CCTK_REAL const position_x_7 = PRIVATE_CARPETREGRID2_STRUCT.position_x_7; \
  CCTK_REAL const position_x_8 = PRIVATE_CARPETREGRID2_STRUCT.position_x_8; \
  CCTK_REAL const position_x_9 = PRIVATE_CARPETREGRID2_STRUCT.position_x_9; \
  CCTK_REAL const position_y_1 = PRIVATE_CARPETREGRID2_STRUCT.position_y_1; \
  CCTK_REAL const position_y_10 = PRIVATE_CARPETREGRID2_STRUCT.position_y_10; \
  CCTK_REAL const position_y_2 = PRIVATE_CARPETREGRID2_STRUCT.position_y_2; \
  CCTK_REAL const position_y_3 = PRIVATE_CARPETREGRID2_STRUCT.position_y_3; \
  CCTK_REAL const position_y_4 = PRIVATE_CARPETREGRID2_STRUCT.position_y_4; \
  CCTK_REAL const position_y_5 = PRIVATE_CARPETREGRID2_STRUCT.position_y_5; \
  CCTK_REAL const position_y_6 = PRIVATE_CARPETREGRID2_STRUCT.position_y_6; \
  CCTK_REAL const position_y_7 = PRIVATE_CARPETREGRID2_STRUCT.position_y_7; \
  CCTK_REAL const position_y_8 = PRIVATE_CARPETREGRID2_STRUCT.position_y_8; \
  CCTK_REAL const position_y_9 = PRIVATE_CARPETREGRID2_STRUCT.position_y_9; \
  CCTK_REAL const position_z_1 = PRIVATE_CARPETREGRID2_STRUCT.position_z_1; \
  CCTK_REAL const position_z_10 = PRIVATE_CARPETREGRID2_STRUCT.position_z_10; \
  CCTK_REAL const position_z_2 = PRIVATE_CARPETREGRID2_STRUCT.position_z_2; \
  CCTK_REAL const position_z_3 = PRIVATE_CARPETREGRID2_STRUCT.position_z_3; \
  CCTK_REAL const position_z_4 = PRIVATE_CARPETREGRID2_STRUCT.position_z_4; \
  CCTK_REAL const position_z_5 = PRIVATE_CARPETREGRID2_STRUCT.position_z_5; \
  CCTK_REAL const position_z_6 = PRIVATE_CARPETREGRID2_STRUCT.position_z_6; \
  CCTK_REAL const position_z_7 = PRIVATE_CARPETREGRID2_STRUCT.position_z_7; \
  CCTK_REAL const position_z_8 = PRIVATE_CARPETREGRID2_STRUCT.position_z_8; \
  CCTK_REAL const position_z_9 = PRIVATE_CARPETREGRID2_STRUCT.position_z_9; \
  CCTK_REAL const * const radius_1 = PRIVATE_CARPETREGRID2_STRUCT.radius_1; \
  CCTK_REAL const * const radius_10 = PRIVATE_CARPETREGRID2_STRUCT.radius_10; \
  CCTK_REAL const * const radius_2 = PRIVATE_CARPETREGRID2_STRUCT.radius_2; \
  CCTK_REAL const * const radius_3 = PRIVATE_CARPETREGRID2_STRUCT.radius_3; \
  CCTK_REAL const * const radius_4 = PRIVATE_CARPETREGRID2_STRUCT.radius_4; \
  CCTK_REAL const * const radius_5 = PRIVATE_CARPETREGRID2_STRUCT.radius_5; \
  CCTK_REAL const * const radius_6 = PRIVATE_CARPETREGRID2_STRUCT.radius_6; \
  CCTK_REAL const * const radius_7 = PRIVATE_CARPETREGRID2_STRUCT.radius_7; \
  CCTK_REAL const * const radius_8 = PRIVATE_CARPETREGRID2_STRUCT.radius_8; \
  CCTK_REAL const * const radius_9 = PRIVATE_CARPETREGRID2_STRUCT.radius_9; \
  CCTK_REAL const radius_change_threshold_1 = PRIVATE_CARPETREGRID2_STRUCT.radius_change_threshold_1; \
  CCTK_REAL const radius_change_threshold_10 = PRIVATE_CARPETREGRID2_STRUCT.radius_change_threshold_10; \
  CCTK_REAL const radius_change_threshold_2 = PRIVATE_CARPETREGRID2_STRUCT.radius_change_threshold_2; \
  CCTK_REAL const radius_change_threshold_3 = PRIVATE_CARPETREGRID2_STRUCT.radius_change_threshold_3; \
  CCTK_REAL const radius_change_threshold_4 = PRIVATE_CARPETREGRID2_STRUCT.radius_change_threshold_4; \
  CCTK_REAL const radius_change_threshold_5 = PRIVATE_CARPETREGRID2_STRUCT.radius_change_threshold_5; \
  CCTK_REAL const radius_change_threshold_6 = PRIVATE_CARPETREGRID2_STRUCT.radius_change_threshold_6; \
  CCTK_REAL const radius_change_threshold_7 = PRIVATE_CARPETREGRID2_STRUCT.radius_change_threshold_7; \
  CCTK_REAL const radius_change_threshold_8 = PRIVATE_CARPETREGRID2_STRUCT.radius_change_threshold_8; \
  CCTK_REAL const radius_change_threshold_9 = PRIVATE_CARPETREGRID2_STRUCT.radius_change_threshold_9; \
  CCTK_INT const active_1 = PRIVATE_CARPETREGRID2_STRUCT.active_1; \
  CCTK_INT const active_10 = PRIVATE_CARPETREGRID2_STRUCT.active_10; \
  CCTK_INT const active_2 = PRIVATE_CARPETREGRID2_STRUCT.active_2; \
  CCTK_INT const active_3 = PRIVATE_CARPETREGRID2_STRUCT.active_3; \
  CCTK_INT const active_4 = PRIVATE_CARPETREGRID2_STRUCT.active_4; \
  CCTK_INT const active_5 = PRIVATE_CARPETREGRID2_STRUCT.active_5; \
  CCTK_INT const active_6 = PRIVATE_CARPETREGRID2_STRUCT.active_6; \
  CCTK_INT const active_7 = PRIVATE_CARPETREGRID2_STRUCT.active_7; \
  CCTK_INT const active_8 = PRIVATE_CARPETREGRID2_STRUCT.active_8; \
  CCTK_INT const active_9 = PRIVATE_CARPETREGRID2_STRUCT.active_9; \
  CCTK_INT const boundary_shiftout = PRIVATE_CARPETREGRID2_STRUCT.boundary_shiftout; \
  CCTK_INT const ensure_proper_nesting = PRIVATE_CARPETREGRID2_STRUCT.ensure_proper_nesting; \
  CCTK_INT const freeze_unaligned_levels = PRIVATE_CARPETREGRID2_STRUCT.freeze_unaligned_levels; \
  CCTK_INT const freeze_unaligned_parent_levels = PRIVATE_CARPETREGRID2_STRUCT.freeze_unaligned_parent_levels; \
  CCTK_INT const min_distance = PRIVATE_CARPETREGRID2_STRUCT.min_distance; \
  CCTK_INT const num_centres = PRIVATE_CARPETREGRID2_STRUCT.num_centres; \
  CCTK_INT const num_levels_1 = PRIVATE_CARPETREGRID2_STRUCT.num_levels_1; \
  CCTK_INT const num_levels_10 = PRIVATE_CARPETREGRID2_STRUCT.num_levels_10; \
  CCTK_INT const num_levels_2 = PRIVATE_CARPETREGRID2_STRUCT.num_levels_2; \
  CCTK_INT const num_levels_3 = PRIVATE_CARPETREGRID2_STRUCT.num_levels_3; \
  CCTK_INT const num_levels_4 = PRIVATE_CARPETREGRID2_STRUCT.num_levels_4; \
  CCTK_INT const num_levels_5 = PRIVATE_CARPETREGRID2_STRUCT.num_levels_5; \
  CCTK_INT const num_levels_6 = PRIVATE_CARPETREGRID2_STRUCT.num_levels_6; \
  CCTK_INT const num_levels_7 = PRIVATE_CARPETREGRID2_STRUCT.num_levels_7; \
  CCTK_INT const num_levels_8 = PRIVATE_CARPETREGRID2_STRUCT.num_levels_8; \
  CCTK_INT const num_levels_9 = PRIVATE_CARPETREGRID2_STRUCT.num_levels_9; \
  CCTK_INT const regrid_every = PRIVATE_CARPETREGRID2_STRUCT.regrid_every; \
  CCTK_INT const snap_to_coarse = PRIVATE_CARPETREGRID2_STRUCT.snap_to_coarse; \
  CCTK_INT const symmetry_rotating180 = PRIVATE_CARPETREGRID2_STRUCT.symmetry_rotating180; \
  CCTK_INT const symmetry_rotating90 = PRIVATE_CARPETREGRID2_STRUCT.symmetry_rotating90; \
  CCTK_INT const verbose = PRIVATE_CARPETREGRID2_STRUCT.verbose; \
  enum { \
      dummy_PRIVATE_CARPETREGRID2_STRUCT_min_fraction = sizeof( min_fraction ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_movement_threshold_1 = sizeof( movement_threshold_1 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_movement_threshold_10 = sizeof( movement_threshold_10 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_movement_threshold_2 = sizeof( movement_threshold_2 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_movement_threshold_3 = sizeof( movement_threshold_3 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_movement_threshold_4 = sizeof( movement_threshold_4 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_movement_threshold_5 = sizeof( movement_threshold_5 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_movement_threshold_6 = sizeof( movement_threshold_6 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_movement_threshold_7 = sizeof( movement_threshold_7 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_movement_threshold_8 = sizeof( movement_threshold_8 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_movement_threshold_9 = sizeof( movement_threshold_9 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_position_x_1 = sizeof( position_x_1 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_position_x_10 = sizeof( position_x_10 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_position_x_2 = sizeof( position_x_2 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_position_x_3 = sizeof( position_x_3 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_position_x_4 = sizeof( position_x_4 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_position_x_5 = sizeof( position_x_5 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_position_x_6 = sizeof( position_x_6 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_position_x_7 = sizeof( position_x_7 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_position_x_8 = sizeof( position_x_8 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_position_x_9 = sizeof( position_x_9 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_position_y_1 = sizeof( position_y_1 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_position_y_10 = sizeof( position_y_10 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_position_y_2 = sizeof( position_y_2 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_position_y_3 = sizeof( position_y_3 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_position_y_4 = sizeof( position_y_4 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_position_y_5 = sizeof( position_y_5 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_position_y_6 = sizeof( position_y_6 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_position_y_7 = sizeof( position_y_7 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_position_y_8 = sizeof( position_y_8 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_position_y_9 = sizeof( position_y_9 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_position_z_1 = sizeof( position_z_1 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_position_z_10 = sizeof( position_z_10 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_position_z_2 = sizeof( position_z_2 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_position_z_3 = sizeof( position_z_3 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_position_z_4 = sizeof( position_z_4 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_position_z_5 = sizeof( position_z_5 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_position_z_6 = sizeof( position_z_6 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_position_z_7 = sizeof( position_z_7 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_position_z_8 = sizeof( position_z_8 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_position_z_9 = sizeof( position_z_9 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_radius_1 = sizeof( radius_1 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_radius_10 = sizeof( radius_10 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_radius_2 = sizeof( radius_2 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_radius_3 = sizeof( radius_3 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_radius_4 = sizeof( radius_4 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_radius_5 = sizeof( radius_5 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_radius_6 = sizeof( radius_6 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_radius_7 = sizeof( radius_7 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_radius_8 = sizeof( radius_8 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_radius_9 = sizeof( radius_9 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_radius_change_threshold_1 = sizeof( radius_change_threshold_1 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_radius_change_threshold_10 = sizeof( radius_change_threshold_10 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_radius_change_threshold_2 = sizeof( radius_change_threshold_2 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_radius_change_threshold_3 = sizeof( radius_change_threshold_3 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_radius_change_threshold_4 = sizeof( radius_change_threshold_4 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_radius_change_threshold_5 = sizeof( radius_change_threshold_5 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_radius_change_threshold_6 = sizeof( radius_change_threshold_6 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_radius_change_threshold_7 = sizeof( radius_change_threshold_7 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_radius_change_threshold_8 = sizeof( radius_change_threshold_8 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_radius_change_threshold_9 = sizeof( radius_change_threshold_9 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_active_1 = sizeof( active_1 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_active_10 = sizeof( active_10 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_active_2 = sizeof( active_2 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_active_3 = sizeof( active_3 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_active_4 = sizeof( active_4 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_active_5 = sizeof( active_5 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_active_6 = sizeof( active_6 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_active_7 = sizeof( active_7 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_active_8 = sizeof( active_8 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_active_9 = sizeof( active_9 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_boundary_shiftout = sizeof( boundary_shiftout ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_ensure_proper_nesting = sizeof( ensure_proper_nesting ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_freeze_unaligned_levels = sizeof( freeze_unaligned_levels ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_freeze_unaligned_parent_levels = sizeof( freeze_unaligned_parent_levels ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_min_distance = sizeof( min_distance ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_num_centres = sizeof( num_centres ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_num_levels_1 = sizeof( num_levels_1 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_num_levels_10 = sizeof( num_levels_10 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_num_levels_2 = sizeof( num_levels_2 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_num_levels_3 = sizeof( num_levels_3 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_num_levels_4 = sizeof( num_levels_4 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_num_levels_5 = sizeof( num_levels_5 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_num_levels_6 = sizeof( num_levels_6 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_num_levels_7 = sizeof( num_levels_7 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_num_levels_8 = sizeof( num_levels_8 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_num_levels_9 = sizeof( num_levels_9 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_regrid_every = sizeof( regrid_every ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_snap_to_coarse = sizeof( snap_to_coarse ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_symmetry_rotating180 = sizeof( symmetry_rotating180 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_symmetry_rotating90 = sizeof( symmetry_rotating90 ) \
    , dummy_PRIVATE_CARPETREGRID2_STRUCT_verbose = sizeof( verbose ) \
  }; \

