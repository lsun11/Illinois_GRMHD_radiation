#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL out0D_dt;
  CCTK_REAL out0D_point_x;
  CCTK_REAL out0D_point_y;
  CCTK_REAL out0D_point_z;
  CCTK_REAL out1D_dt;
  CCTK_REAL out1D_xline_y;
  CCTK_REAL out1D_xline_z;
  CCTK_REAL out1D_yline_x;
  CCTK_REAL out1D_yline_z;
  CCTK_REAL out1D_zline_x;
  CCTK_REAL out1D_zline_y;
  CCTK_REAL out2D_dt;
  CCTK_REAL out2D_xyplane_z;
  CCTK_REAL out2D_xzplane_y;
  CCTK_REAL out2D_yzplane_x;
  CCTK_REAL out_dt;
  const char * out0D_criterion;
  const char * out0D_dir;
  const char * out0D_vars;
  const char * out1D_criterion;
  const char * out1D_dir;
  const char * out1D_vars;
  const char * out2D_criterion;
  const char * out2D_dir;
  const char * out2D_vars;
  const char * out_criterion;
  const char * out_dir;
  const char * out_extension;
  const char * out_vars;
  CCTK_INT checkpoint;
  CCTK_INT checkpoint_next;
  CCTK_INT compression_level;
  CCTK_INT one_file_per_group;
  CCTK_INT open_one_input_file_at_a_time;
  CCTK_INT out0D_every;
  CCTK_INT out0D_point_xi;
  CCTK_INT out0D_point_yi;
  CCTK_INT out0D_point_zi;
  CCTK_INT out1D_d;
  CCTK_INT out1D_every;
  CCTK_INT out1D_x;
  CCTK_INT out1D_xline_yi;
  CCTK_INT out1D_xline_zi;
  CCTK_INT out1D_y;
  CCTK_INT out1D_yline_xi;
  CCTK_INT out1D_yline_zi;
  CCTK_INT out1D_z;
  CCTK_INT out1D_zline_xi;
  CCTK_INT out1D_zline_yi;
  CCTK_INT out2D_every;
  CCTK_INT out2D_xy;
  CCTK_INT out2D_xyplane_zi;
  CCTK_INT out2D_xz;
  CCTK_INT out2D_xzplane_yi;
  CCTK_INT out2D_yz;
  CCTK_INT out2D_yzplane_xi;
  CCTK_INT out3D_ghosts;
  CCTK_INT out3D_outer_ghosts;
  CCTK_INT out_every;
  CCTK_INT output_all_timelevels;
  CCTK_INT output_index;
  CCTK_INT output_symmetry_points;
  CCTK_INT use_checksums;
  CCTK_INT use_grid_structure_from_checkpoint;
  CCTK_INT use_reflevels_from_checkpoint;
} PRIVATE_CARPETIOHDF5_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_CARPETIOHDF5_STRUCT_PARAMS \
  CCTK_REAL const out0D_dt = PRIVATE_CARPETIOHDF5_STRUCT.out0D_dt; \
  CCTK_REAL const out0D_point_x = PRIVATE_CARPETIOHDF5_STRUCT.out0D_point_x; \
  CCTK_REAL const out0D_point_y = PRIVATE_CARPETIOHDF5_STRUCT.out0D_point_y; \
  CCTK_REAL const out0D_point_z = PRIVATE_CARPETIOHDF5_STRUCT.out0D_point_z; \
  CCTK_REAL const out1D_dt = PRIVATE_CARPETIOHDF5_STRUCT.out1D_dt; \
  CCTK_REAL const out1D_xline_y = PRIVATE_CARPETIOHDF5_STRUCT.out1D_xline_y; \
  CCTK_REAL const out1D_xline_z = PRIVATE_CARPETIOHDF5_STRUCT.out1D_xline_z; \
  CCTK_REAL const out1D_yline_x = PRIVATE_CARPETIOHDF5_STRUCT.out1D_yline_x; \
  CCTK_REAL const out1D_yline_z = PRIVATE_CARPETIOHDF5_STRUCT.out1D_yline_z; \
  CCTK_REAL const out1D_zline_x = PRIVATE_CARPETIOHDF5_STRUCT.out1D_zline_x; \
  CCTK_REAL const out1D_zline_y = PRIVATE_CARPETIOHDF5_STRUCT.out1D_zline_y; \
  CCTK_REAL const out2D_dt = PRIVATE_CARPETIOHDF5_STRUCT.out2D_dt; \
  CCTK_REAL const out2D_xyplane_z = PRIVATE_CARPETIOHDF5_STRUCT.out2D_xyplane_z; \
  CCTK_REAL const out2D_xzplane_y = PRIVATE_CARPETIOHDF5_STRUCT.out2D_xzplane_y; \
  CCTK_REAL const out2D_yzplane_x = PRIVATE_CARPETIOHDF5_STRUCT.out2D_yzplane_x; \
  CCTK_REAL const out_dt = PRIVATE_CARPETIOHDF5_STRUCT.out_dt; \
  const char * const out0D_criterion = PRIVATE_CARPETIOHDF5_STRUCT.out0D_criterion; \
  const char * const out0D_dir = PRIVATE_CARPETIOHDF5_STRUCT.out0D_dir; \
  const char * const out0D_vars = PRIVATE_CARPETIOHDF5_STRUCT.out0D_vars; \
  const char * const out1D_criterion = PRIVATE_CARPETIOHDF5_STRUCT.out1D_criterion; \
  const char * const out1D_dir = PRIVATE_CARPETIOHDF5_STRUCT.out1D_dir; \
  const char * const out1D_vars = PRIVATE_CARPETIOHDF5_STRUCT.out1D_vars; \
  const char * const out2D_criterion = PRIVATE_CARPETIOHDF5_STRUCT.out2D_criterion; \
  const char * const out2D_dir = PRIVATE_CARPETIOHDF5_STRUCT.out2D_dir; \
  const char * const out2D_vars = PRIVATE_CARPETIOHDF5_STRUCT.out2D_vars; \
  const char * const out_criterion = PRIVATE_CARPETIOHDF5_STRUCT.out_criterion; \
  const char * const out_dir = PRIVATE_CARPETIOHDF5_STRUCT.out_dir; \
  const char * const out_extension = PRIVATE_CARPETIOHDF5_STRUCT.out_extension; \
  const char * const out_vars = PRIVATE_CARPETIOHDF5_STRUCT.out_vars; \
  CCTK_INT const checkpoint = PRIVATE_CARPETIOHDF5_STRUCT.checkpoint; \
  CCTK_INT const checkpoint_next = PRIVATE_CARPETIOHDF5_STRUCT.checkpoint_next; \
  CCTK_INT const compression_level = PRIVATE_CARPETIOHDF5_STRUCT.compression_level; \
  CCTK_INT const one_file_per_group = PRIVATE_CARPETIOHDF5_STRUCT.one_file_per_group; \
  CCTK_INT const open_one_input_file_at_a_time = PRIVATE_CARPETIOHDF5_STRUCT.open_one_input_file_at_a_time; \
  CCTK_INT const out0D_every = PRIVATE_CARPETIOHDF5_STRUCT.out0D_every; \
  CCTK_INT const out0D_point_xi = PRIVATE_CARPETIOHDF5_STRUCT.out0D_point_xi; \
  CCTK_INT const out0D_point_yi = PRIVATE_CARPETIOHDF5_STRUCT.out0D_point_yi; \
  CCTK_INT const out0D_point_zi = PRIVATE_CARPETIOHDF5_STRUCT.out0D_point_zi; \
  CCTK_INT const out1D_d = PRIVATE_CARPETIOHDF5_STRUCT.out1D_d; \
  CCTK_INT const out1D_every = PRIVATE_CARPETIOHDF5_STRUCT.out1D_every; \
  CCTK_INT const out1D_x = PRIVATE_CARPETIOHDF5_STRUCT.out1D_x; \
  CCTK_INT const out1D_xline_yi = PRIVATE_CARPETIOHDF5_STRUCT.out1D_xline_yi; \
  CCTK_INT const out1D_xline_zi = PRIVATE_CARPETIOHDF5_STRUCT.out1D_xline_zi; \
  CCTK_INT const out1D_y = PRIVATE_CARPETIOHDF5_STRUCT.out1D_y; \
  CCTK_INT const out1D_yline_xi = PRIVATE_CARPETIOHDF5_STRUCT.out1D_yline_xi; \
  CCTK_INT const out1D_yline_zi = PRIVATE_CARPETIOHDF5_STRUCT.out1D_yline_zi; \
  CCTK_INT const out1D_z = PRIVATE_CARPETIOHDF5_STRUCT.out1D_z; \
  CCTK_INT const out1D_zline_xi = PRIVATE_CARPETIOHDF5_STRUCT.out1D_zline_xi; \
  CCTK_INT const out1D_zline_yi = PRIVATE_CARPETIOHDF5_STRUCT.out1D_zline_yi; \
  CCTK_INT const out2D_every = PRIVATE_CARPETIOHDF5_STRUCT.out2D_every; \
  CCTK_INT const out2D_xy = PRIVATE_CARPETIOHDF5_STRUCT.out2D_xy; \
  CCTK_INT const out2D_xyplane_zi = PRIVATE_CARPETIOHDF5_STRUCT.out2D_xyplane_zi; \
  CCTK_INT const out2D_xz = PRIVATE_CARPETIOHDF5_STRUCT.out2D_xz; \
  CCTK_INT const out2D_xzplane_yi = PRIVATE_CARPETIOHDF5_STRUCT.out2D_xzplane_yi; \
  CCTK_INT const out2D_yz = PRIVATE_CARPETIOHDF5_STRUCT.out2D_yz; \
  CCTK_INT const out2D_yzplane_xi = PRIVATE_CARPETIOHDF5_STRUCT.out2D_yzplane_xi; \
  CCTK_INT const out3D_ghosts = PRIVATE_CARPETIOHDF5_STRUCT.out3D_ghosts; \
  CCTK_INT const out3D_outer_ghosts = PRIVATE_CARPETIOHDF5_STRUCT.out3D_outer_ghosts; \
  CCTK_INT const out_every = PRIVATE_CARPETIOHDF5_STRUCT.out_every; \
  CCTK_INT const output_all_timelevels = PRIVATE_CARPETIOHDF5_STRUCT.output_all_timelevels; \
  CCTK_INT const output_index = PRIVATE_CARPETIOHDF5_STRUCT.output_index; \
  CCTK_INT const output_symmetry_points = PRIVATE_CARPETIOHDF5_STRUCT.output_symmetry_points; \
  CCTK_INT const use_checksums = PRIVATE_CARPETIOHDF5_STRUCT.use_checksums; \
  CCTK_INT const use_grid_structure_from_checkpoint = PRIVATE_CARPETIOHDF5_STRUCT.use_grid_structure_from_checkpoint; \
  CCTK_INT const use_reflevels_from_checkpoint = PRIVATE_CARPETIOHDF5_STRUCT.use_reflevels_from_checkpoint; \
  enum { \
      dummy_PRIVATE_CARPETIOHDF5_STRUCT_out0D_dt = sizeof( out0D_dt ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out0D_point_x = sizeof( out0D_point_x ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out0D_point_y = sizeof( out0D_point_y ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out0D_point_z = sizeof( out0D_point_z ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out1D_dt = sizeof( out1D_dt ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out1D_xline_y = sizeof( out1D_xline_y ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out1D_xline_z = sizeof( out1D_xline_z ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out1D_yline_x = sizeof( out1D_yline_x ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out1D_yline_z = sizeof( out1D_yline_z ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out1D_zline_x = sizeof( out1D_zline_x ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out1D_zline_y = sizeof( out1D_zline_y ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out2D_dt = sizeof( out2D_dt ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out2D_xyplane_z = sizeof( out2D_xyplane_z ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out2D_xzplane_y = sizeof( out2D_xzplane_y ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out2D_yzplane_x = sizeof( out2D_yzplane_x ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out_dt = sizeof( out_dt ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out0D_criterion = sizeof( out0D_criterion ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out0D_dir = sizeof( out0D_dir ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out0D_vars = sizeof( out0D_vars ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out1D_criterion = sizeof( out1D_criterion ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out1D_dir = sizeof( out1D_dir ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out1D_vars = sizeof( out1D_vars ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out2D_criterion = sizeof( out2D_criterion ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out2D_dir = sizeof( out2D_dir ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out2D_vars = sizeof( out2D_vars ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out_criterion = sizeof( out_criterion ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out_dir = sizeof( out_dir ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out_extension = sizeof( out_extension ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out_vars = sizeof( out_vars ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_checkpoint = sizeof( checkpoint ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_checkpoint_next = sizeof( checkpoint_next ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_compression_level = sizeof( compression_level ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_one_file_per_group = sizeof( one_file_per_group ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_open_one_input_file_at_a_time = sizeof( open_one_input_file_at_a_time ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out0D_every = sizeof( out0D_every ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out0D_point_xi = sizeof( out0D_point_xi ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out0D_point_yi = sizeof( out0D_point_yi ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out0D_point_zi = sizeof( out0D_point_zi ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out1D_d = sizeof( out1D_d ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out1D_every = sizeof( out1D_every ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out1D_x = sizeof( out1D_x ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out1D_xline_yi = sizeof( out1D_xline_yi ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out1D_xline_zi = sizeof( out1D_xline_zi ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out1D_y = sizeof( out1D_y ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out1D_yline_xi = sizeof( out1D_yline_xi ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out1D_yline_zi = sizeof( out1D_yline_zi ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out1D_z = sizeof( out1D_z ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out1D_zline_xi = sizeof( out1D_zline_xi ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out1D_zline_yi = sizeof( out1D_zline_yi ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out2D_every = sizeof( out2D_every ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out2D_xy = sizeof( out2D_xy ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out2D_xyplane_zi = sizeof( out2D_xyplane_zi ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out2D_xz = sizeof( out2D_xz ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out2D_xzplane_yi = sizeof( out2D_xzplane_yi ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out2D_yz = sizeof( out2D_yz ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out2D_yzplane_xi = sizeof( out2D_yzplane_xi ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out3D_ghosts = sizeof( out3D_ghosts ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out3D_outer_ghosts = sizeof( out3D_outer_ghosts ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_out_every = sizeof( out_every ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_output_all_timelevels = sizeof( output_all_timelevels ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_output_index = sizeof( output_index ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_output_symmetry_points = sizeof( output_symmetry_points ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_use_checksums = sizeof( use_checksums ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_use_grid_structure_from_checkpoint = sizeof( use_grid_structure_from_checkpoint ) \
    , dummy_PRIVATE_CARPETIOHDF5_STRUCT_use_reflevels_from_checkpoint = sizeof( use_reflevels_from_checkpoint ) \
  }; \

