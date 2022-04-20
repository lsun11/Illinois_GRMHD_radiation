#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL checkpoint_every_walltime_hours;
  CCTK_REAL out_dt;
  CCTK_REAL out_xline_y;
  CCTK_REAL out_xline_z;
  CCTK_REAL out_xyplane_z;
  CCTK_REAL out_xzplane_y;
  CCTK_REAL out_yline_x;
  CCTK_REAL out_yline_z;
  CCTK_REAL out_yzplane_x;
  CCTK_REAL out_zline_x;
  CCTK_REAL out_zline_y;
  const char * checkpoint_ID_file;
  const char * checkpoint_dir;
  const char * checkpoint_file;
  const char * filereader_ID_dir;
  const char * filereader_ID_files;
  const char * filereader_ID_vars;
  const char * out_criterion;
  const char * out_dir;
  const char * out_fileinfo;
  const char * out_mode;
  const char * out_save_parameters;
  const char * parfile_name;
  const char * parfile_write;
  const char * recover;
  const char * recover_dir;
  const char * recover_file;
  const char * verbose;
  CCTK_INT abort_on_io_errors;
  CCTK_INT checkpoint_ID;
  CCTK_INT checkpoint_every;
  CCTK_INT checkpoint_keep;
  CCTK_INT checkpoint_on_terminate;
  CCTK_INT new_filename_scheme;
  CCTK_INT out3D_septimefiles;
  CCTK_INT out_downsample_x;
  CCTK_INT out_downsample_y;
  CCTK_INT out_downsample_z;
  CCTK_INT out_every;
  CCTK_INT out_proc_every;
  CCTK_INT out_single_precision;
  CCTK_INT out_timesteps_per_file;
  CCTK_INT out_unchunked;
  CCTK_INT out_xline_yi;
  CCTK_INT out_xline_zi;
  CCTK_INT out_xyplane_zi;
  CCTK_INT out_xzplane_yi;
  CCTK_INT out_yline_xi;
  CCTK_INT out_yline_zi;
  CCTK_INT out_yzplane_xi;
  CCTK_INT out_zline_xi;
  CCTK_INT out_zline_yi;
  CCTK_INT parfile_update_every;
  CCTK_INT print_timing_info;
  CCTK_INT recover_and_remove;
  CCTK_INT require_empty_output_directory;
  CCTK_INT strict_io_parameter_check;
} RESTRICTED_IO_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_RESTRICTED_IO_STRUCT_PARAMS \
  CCTK_REAL const checkpoint_every_walltime_hours = RESTRICTED_IO_STRUCT.checkpoint_every_walltime_hours; \
  CCTK_REAL const out_dt = RESTRICTED_IO_STRUCT.out_dt; \
  CCTK_REAL const out_xline_y = RESTRICTED_IO_STRUCT.out_xline_y; \
  CCTK_REAL const out_xline_z = RESTRICTED_IO_STRUCT.out_xline_z; \
  CCTK_REAL const out_xyplane_z = RESTRICTED_IO_STRUCT.out_xyplane_z; \
  CCTK_REAL const out_xzplane_y = RESTRICTED_IO_STRUCT.out_xzplane_y; \
  CCTK_REAL const out_yline_x = RESTRICTED_IO_STRUCT.out_yline_x; \
  CCTK_REAL const out_yline_z = RESTRICTED_IO_STRUCT.out_yline_z; \
  CCTK_REAL const out_yzplane_x = RESTRICTED_IO_STRUCT.out_yzplane_x; \
  CCTK_REAL const out_zline_x = RESTRICTED_IO_STRUCT.out_zline_x; \
  CCTK_REAL const out_zline_y = RESTRICTED_IO_STRUCT.out_zline_y; \
  const char * const checkpoint_ID_file = RESTRICTED_IO_STRUCT.checkpoint_ID_file; \
  const char * const checkpoint_dir = RESTRICTED_IO_STRUCT.checkpoint_dir; \
  const char * const checkpoint_file = RESTRICTED_IO_STRUCT.checkpoint_file; \
  const char * const filereader_ID_dir = RESTRICTED_IO_STRUCT.filereader_ID_dir; \
  const char * const filereader_ID_files = RESTRICTED_IO_STRUCT.filereader_ID_files; \
  const char * const filereader_ID_vars = RESTRICTED_IO_STRUCT.filereader_ID_vars; \
  const char * const out_criterion = RESTRICTED_IO_STRUCT.out_criterion; \
  const char * const out_dir = RESTRICTED_IO_STRUCT.out_dir; \
  const char * const out_fileinfo = RESTRICTED_IO_STRUCT.out_fileinfo; \
  const char * const out_mode = RESTRICTED_IO_STRUCT.out_mode; \
  const char * const out_save_parameters = RESTRICTED_IO_STRUCT.out_save_parameters; \
  const char * const parfile_name = RESTRICTED_IO_STRUCT.parfile_name; \
  const char * const parfile_write = RESTRICTED_IO_STRUCT.parfile_write; \
  const char * const recover = RESTRICTED_IO_STRUCT.recover; \
  const char * const recover_dir = RESTRICTED_IO_STRUCT.recover_dir; \
  const char * const recover_file = RESTRICTED_IO_STRUCT.recover_file; \
  const char * const verbose = RESTRICTED_IO_STRUCT.verbose; \
  CCTK_INT const abort_on_io_errors = RESTRICTED_IO_STRUCT.abort_on_io_errors; \
  CCTK_INT const checkpoint_ID = RESTRICTED_IO_STRUCT.checkpoint_ID; \
  CCTK_INT const checkpoint_every = RESTRICTED_IO_STRUCT.checkpoint_every; \
  CCTK_INT const checkpoint_keep = RESTRICTED_IO_STRUCT.checkpoint_keep; \
  CCTK_INT const checkpoint_on_terminate = RESTRICTED_IO_STRUCT.checkpoint_on_terminate; \
  CCTK_INT const new_filename_scheme = RESTRICTED_IO_STRUCT.new_filename_scheme; \
  CCTK_INT const out3D_septimefiles = RESTRICTED_IO_STRUCT.out3D_septimefiles; \
  CCTK_INT const out_downsample_x = RESTRICTED_IO_STRUCT.out_downsample_x; \
  CCTK_INT const out_downsample_y = RESTRICTED_IO_STRUCT.out_downsample_y; \
  CCTK_INT const out_downsample_z = RESTRICTED_IO_STRUCT.out_downsample_z; \
  CCTK_INT const out_every = RESTRICTED_IO_STRUCT.out_every; \
  CCTK_INT const out_proc_every = RESTRICTED_IO_STRUCT.out_proc_every; \
  CCTK_INT const out_single_precision = RESTRICTED_IO_STRUCT.out_single_precision; \
  CCTK_INT const out_timesteps_per_file = RESTRICTED_IO_STRUCT.out_timesteps_per_file; \
  CCTK_INT const out_unchunked = RESTRICTED_IO_STRUCT.out_unchunked; \
  CCTK_INT const out_xline_yi = RESTRICTED_IO_STRUCT.out_xline_yi; \
  CCTK_INT const out_xline_zi = RESTRICTED_IO_STRUCT.out_xline_zi; \
  CCTK_INT const out_xyplane_zi = RESTRICTED_IO_STRUCT.out_xyplane_zi; \
  CCTK_INT const out_xzplane_yi = RESTRICTED_IO_STRUCT.out_xzplane_yi; \
  CCTK_INT const out_yline_xi = RESTRICTED_IO_STRUCT.out_yline_xi; \
  CCTK_INT const out_yline_zi = RESTRICTED_IO_STRUCT.out_yline_zi; \
  CCTK_INT const out_yzplane_xi = RESTRICTED_IO_STRUCT.out_yzplane_xi; \
  CCTK_INT const out_zline_xi = RESTRICTED_IO_STRUCT.out_zline_xi; \
  CCTK_INT const out_zline_yi = RESTRICTED_IO_STRUCT.out_zline_yi; \
  CCTK_INT const parfile_update_every = RESTRICTED_IO_STRUCT.parfile_update_every; \
  CCTK_INT const print_timing_info = RESTRICTED_IO_STRUCT.print_timing_info; \
  CCTK_INT const recover_and_remove = RESTRICTED_IO_STRUCT.recover_and_remove; \
  CCTK_INT const require_empty_output_directory = RESTRICTED_IO_STRUCT.require_empty_output_directory; \
  CCTK_INT const strict_io_parameter_check = RESTRICTED_IO_STRUCT.strict_io_parameter_check; \
  enum { \
      dummy_RESTRICTED_IO_STRUCT_checkpoint_every_walltime_hours = sizeof( checkpoint_every_walltime_hours ) \
    , dummy_RESTRICTED_IO_STRUCT_out_dt = sizeof( out_dt ) \
    , dummy_RESTRICTED_IO_STRUCT_out_xline_y = sizeof( out_xline_y ) \
    , dummy_RESTRICTED_IO_STRUCT_out_xline_z = sizeof( out_xline_z ) \
    , dummy_RESTRICTED_IO_STRUCT_out_xyplane_z = sizeof( out_xyplane_z ) \
    , dummy_RESTRICTED_IO_STRUCT_out_xzplane_y = sizeof( out_xzplane_y ) \
    , dummy_RESTRICTED_IO_STRUCT_out_yline_x = sizeof( out_yline_x ) \
    , dummy_RESTRICTED_IO_STRUCT_out_yline_z = sizeof( out_yline_z ) \
    , dummy_RESTRICTED_IO_STRUCT_out_yzplane_x = sizeof( out_yzplane_x ) \
    , dummy_RESTRICTED_IO_STRUCT_out_zline_x = sizeof( out_zline_x ) \
    , dummy_RESTRICTED_IO_STRUCT_out_zline_y = sizeof( out_zline_y ) \
    , dummy_RESTRICTED_IO_STRUCT_checkpoint_ID_file = sizeof( checkpoint_ID_file ) \
    , dummy_RESTRICTED_IO_STRUCT_checkpoint_dir = sizeof( checkpoint_dir ) \
    , dummy_RESTRICTED_IO_STRUCT_checkpoint_file = sizeof( checkpoint_file ) \
    , dummy_RESTRICTED_IO_STRUCT_filereader_ID_dir = sizeof( filereader_ID_dir ) \
    , dummy_RESTRICTED_IO_STRUCT_filereader_ID_files = sizeof( filereader_ID_files ) \
    , dummy_RESTRICTED_IO_STRUCT_filereader_ID_vars = sizeof( filereader_ID_vars ) \
    , dummy_RESTRICTED_IO_STRUCT_out_criterion = sizeof( out_criterion ) \
    , dummy_RESTRICTED_IO_STRUCT_out_dir = sizeof( out_dir ) \
    , dummy_RESTRICTED_IO_STRUCT_out_fileinfo = sizeof( out_fileinfo ) \
    , dummy_RESTRICTED_IO_STRUCT_out_mode = sizeof( out_mode ) \
    , dummy_RESTRICTED_IO_STRUCT_out_save_parameters = sizeof( out_save_parameters ) \
    , dummy_RESTRICTED_IO_STRUCT_parfile_name = sizeof( parfile_name ) \
    , dummy_RESTRICTED_IO_STRUCT_parfile_write = sizeof( parfile_write ) \
    , dummy_RESTRICTED_IO_STRUCT_recover = sizeof( recover ) \
    , dummy_RESTRICTED_IO_STRUCT_recover_dir = sizeof( recover_dir ) \
    , dummy_RESTRICTED_IO_STRUCT_recover_file = sizeof( recover_file ) \
    , dummy_RESTRICTED_IO_STRUCT_verbose = sizeof( verbose ) \
    , dummy_RESTRICTED_IO_STRUCT_abort_on_io_errors = sizeof( abort_on_io_errors ) \
    , dummy_RESTRICTED_IO_STRUCT_checkpoint_ID = sizeof( checkpoint_ID ) \
    , dummy_RESTRICTED_IO_STRUCT_checkpoint_every = sizeof( checkpoint_every ) \
    , dummy_RESTRICTED_IO_STRUCT_checkpoint_keep = sizeof( checkpoint_keep ) \
    , dummy_RESTRICTED_IO_STRUCT_checkpoint_on_terminate = sizeof( checkpoint_on_terminate ) \
    , dummy_RESTRICTED_IO_STRUCT_new_filename_scheme = sizeof( new_filename_scheme ) \
    , dummy_RESTRICTED_IO_STRUCT_out3D_septimefiles = sizeof( out3D_septimefiles ) \
    , dummy_RESTRICTED_IO_STRUCT_out_downsample_x = sizeof( out_downsample_x ) \
    , dummy_RESTRICTED_IO_STRUCT_out_downsample_y = sizeof( out_downsample_y ) \
    , dummy_RESTRICTED_IO_STRUCT_out_downsample_z = sizeof( out_downsample_z ) \
    , dummy_RESTRICTED_IO_STRUCT_out_every = sizeof( out_every ) \
    , dummy_RESTRICTED_IO_STRUCT_out_proc_every = sizeof( out_proc_every ) \
    , dummy_RESTRICTED_IO_STRUCT_out_single_precision = sizeof( out_single_precision ) \
    , dummy_RESTRICTED_IO_STRUCT_out_timesteps_per_file = sizeof( out_timesteps_per_file ) \
    , dummy_RESTRICTED_IO_STRUCT_out_unchunked = sizeof( out_unchunked ) \
    , dummy_RESTRICTED_IO_STRUCT_out_xline_yi = sizeof( out_xline_yi ) \
    , dummy_RESTRICTED_IO_STRUCT_out_xline_zi = sizeof( out_xline_zi ) \
    , dummy_RESTRICTED_IO_STRUCT_out_xyplane_zi = sizeof( out_xyplane_zi ) \
    , dummy_RESTRICTED_IO_STRUCT_out_xzplane_yi = sizeof( out_xzplane_yi ) \
    , dummy_RESTRICTED_IO_STRUCT_out_yline_xi = sizeof( out_yline_xi ) \
    , dummy_RESTRICTED_IO_STRUCT_out_yline_zi = sizeof( out_yline_zi ) \
    , dummy_RESTRICTED_IO_STRUCT_out_yzplane_xi = sizeof( out_yzplane_xi ) \
    , dummy_RESTRICTED_IO_STRUCT_out_zline_xi = sizeof( out_zline_xi ) \
    , dummy_RESTRICTED_IO_STRUCT_out_zline_yi = sizeof( out_zline_yi ) \
    , dummy_RESTRICTED_IO_STRUCT_parfile_update_every = sizeof( parfile_update_every ) \
    , dummy_RESTRICTED_IO_STRUCT_print_timing_info = sizeof( print_timing_info ) \
    , dummy_RESTRICTED_IO_STRUCT_recover_and_remove = sizeof( recover_and_remove ) \
    , dummy_RESTRICTED_IO_STRUCT_require_empty_output_directory = sizeof( require_empty_output_directory ) \
    , dummy_RESTRICTED_IO_STRUCT_strict_io_parameter_check = sizeof( strict_io_parameter_check ) \
  }; \

