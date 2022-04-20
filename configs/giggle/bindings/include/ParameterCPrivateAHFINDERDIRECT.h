#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL ILUCG__error_tolerance;
  CCTK_REAL Jacobian_perturbation_amplitude;
  CCTK_REAL Theta_norm_for_convergence;
  CCTK_REAL desired_value[101];
  CCTK_REAL desired_value_factor[101];
  CCTK_REAL desired_value_offset[101];
  CCTK_REAL dont_find_after_individual_time[101];
  CCTK_REAL find_after_individual_time[101];
  CCTK_REAL geometry__Schwarzschild_EF__Delta_xyz;
  CCTK_REAL geometry__Schwarzschild_EF__epsilon;
  CCTK_REAL geometry__Schwarzschild_EF__mass;
  CCTK_REAL geometry__Schwarzschild_EF__x_posn;
  CCTK_REAL geometry__Schwarzschild_EF__y_posn;
  CCTK_REAL geometry__Schwarzschild_EF__z_posn;
  CCTK_REAL mask_buffer_thickness;
  CCTK_REAL mask_radius_multiplier;
  CCTK_REAL mask_radius_offset;
  CCTK_REAL max_allowable_Delta_h_over_h;
  CCTK_REAL max_allowable_Theta;
  CCTK_REAL max_allowable_horizon_radius[101];
  CCTK_REAL min_horizon_radius_points_for_mask;
  CCTK_REAL old_style_mask_buffer_value;
  CCTK_REAL old_style_mask_inside_value;
  CCTK_REAL old_style_mask_outside_value;
  CCTK_REAL pretracking_delta[101];
  CCTK_REAL pretracking_maximum_delta[101];
  CCTK_REAL pretracking_maximum_value[101];
  CCTK_REAL pretracking_minimum_delta[101];
  CCTK_REAL pretracking_minimum_value[101];
  CCTK_REAL pretracking_value[101];
  CCTK_REAL shiftout_factor[101];
  CCTK_REAL smoothing_factor[101];
  const char * ASCII_gnuplot_file_name_extension;
  const char * BH_diagnostics_base_file_name;
  const char * BH_diagnostics_directory;
  const char * BH_diagnostics_file_name_extension;
  const char * Delta_h_base_file_name;
  const char * HDF5_file_name_extension;
  const char * Jacobian_base_file_name;
  const char * Jacobian_compute_method;
  const char * Jacobian_store_solve_method;
  const char * OpenDX_control_file_name_extension;
  const char * Theta_base_file_name;
  const char * coordinate_system_name;
  const char * geometry_interpolator_name;
  const char * geometry_interpolator_pars;
  const char * h_base_file_name;
  const char * h_directory;
  const char * integral_method;
  const char * interpatch_interpolator_name;
  const char * interpatch_interpolator_pars;
  const char * mean_curvature_base_file_name;
  const char * method;
  const char * metric_type;
  const char * new_style_mask_bitfield_name;
  const char * new_style_mask_buffer_value;
  const char * new_style_mask_gridfn_name;
  const char * new_style_mask_inside_value;
  const char * new_style_mask_outside_value;
  const char * old_style_mask_gridfn_name;
  const char * patch_system_type[101];
  const char * surface_definition[101];
  const char * surface_interpolator_name;
  const char * surface_interpolator_pars;
  const char * surface_modification[101];
  const char * surface_selection[101];
  const char * verbose_level;
  CCTK_INT ILUCG__limit_CG_iterations;
  CCTK_INT N_horizons;
  CCTK_INT N_zones_per_right_angle[101];
  CCTK_INT UMFPACK__N_II_iterations;
  CCTK_INT check_that_geometry_is_finite;
  CCTK_INT check_that_h_is_finite;
  CCTK_INT debugging_output_at_each_Newton_iteration;
  CCTK_INT depends_on[101];
  CCTK_INT disable_horizon[101];
  CCTK_INT dont_find_after_individual[101];
  CCTK_INT find_after_individual[101];
  CCTK_INT find_every;
  CCTK_INT find_every_individual[101];
  CCTK_INT ghost_zone_width;
  CCTK_INT h_min_digits;
  CCTK_INT hardwire_Schwarzschild_EF_geometry;
  CCTK_INT mask_is_noshrink;
  CCTK_INT max_N_zones_per_right_angle;
  CCTK_INT max_Newton_iterations__initial;
  CCTK_INT max_Newton_iterations__subsequent;
  CCTK_INT max_allowable_Theta_growth_iterations;
  CCTK_INT max_allowable_Theta_nonshrink_iterations;
  CCTK_INT output_ASCII_files;
  CCTK_INT output_BH_diagnostics;
  CCTK_INT output_HDF5_files;
  CCTK_INT output_OpenDX_control_files;
  CCTK_INT output_Theta_every;
  CCTK_INT output_ghost_zones_for_h;
  CCTK_INT output_h_every;
  CCTK_INT output_initial_guess;
  CCTK_INT output_mean_curvature_every;
  CCTK_INT patch_overlap_width;
  CCTK_INT pretracking_max_iterations[101];
  CCTK_INT print_timing_stats;
  CCTK_INT run_at_CCTK_ANALYSIS;
  CCTK_INT run_at_CCTK_POSTINITIAL;
  CCTK_INT run_at_CCTK_POSTSTEP;
  CCTK_INT run_at_CCTK_POST_RECOVER_VARIABLES;
  CCTK_INT set_mask_for_all_horizons;
  CCTK_INT set_mask_for_individual_horizon[101];
  CCTK_INT set_new_style_mask;
  CCTK_INT set_old_style_mask;
  CCTK_INT test_all_Jacobian_compute_methods;
  CCTK_INT use_pretracking[101];
  CCTK_INT want_expansion_gradients;
  CCTK_INT warn_level__gij_not_positive_definite__initial;
  CCTK_INT warn_level__gij_not_positive_definite__subsequent;
  CCTK_INT warn_level__nonfinite_geometry;
  CCTK_INT warn_level__point_outside__initial;
  CCTK_INT warn_level__point_outside__subsequent;
  CCTK_INT warn_level__skipping_finite_check;
  CCTK_INT which_horizon_to_announce_centroid;
  CCTK_INT which_surface_to_store_info[101];
} PRIVATE_AHFINDERDIRECT_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_AHFINDERDIRECT_STRUCT_PARAMS \
  CCTK_REAL const ILUCG__error_tolerance = PRIVATE_AHFINDERDIRECT_STRUCT.ILUCG__error_tolerance; \
  CCTK_REAL const Jacobian_perturbation_amplitude = PRIVATE_AHFINDERDIRECT_STRUCT.Jacobian_perturbation_amplitude; \
  CCTK_REAL const Theta_norm_for_convergence = PRIVATE_AHFINDERDIRECT_STRUCT.Theta_norm_for_convergence; \
  CCTK_REAL const * const desired_value = PRIVATE_AHFINDERDIRECT_STRUCT.desired_value; \
  CCTK_REAL const * const desired_value_factor = PRIVATE_AHFINDERDIRECT_STRUCT.desired_value_factor; \
  CCTK_REAL const * const desired_value_offset = PRIVATE_AHFINDERDIRECT_STRUCT.desired_value_offset; \
  CCTK_REAL const * const dont_find_after_individual_time = PRIVATE_AHFINDERDIRECT_STRUCT.dont_find_after_individual_time; \
  CCTK_REAL const * const find_after_individual_time = PRIVATE_AHFINDERDIRECT_STRUCT.find_after_individual_time; \
  CCTK_REAL const geometry__Schwarzschild_EF__Delta_xyz = PRIVATE_AHFINDERDIRECT_STRUCT.geometry__Schwarzschild_EF__Delta_xyz; \
  CCTK_REAL const geometry__Schwarzschild_EF__epsilon = PRIVATE_AHFINDERDIRECT_STRUCT.geometry__Schwarzschild_EF__epsilon; \
  CCTK_REAL const geometry__Schwarzschild_EF__mass = PRIVATE_AHFINDERDIRECT_STRUCT.geometry__Schwarzschild_EF__mass; \
  CCTK_REAL const geometry__Schwarzschild_EF__x_posn = PRIVATE_AHFINDERDIRECT_STRUCT.geometry__Schwarzschild_EF__x_posn; \
  CCTK_REAL const geometry__Schwarzschild_EF__y_posn = PRIVATE_AHFINDERDIRECT_STRUCT.geometry__Schwarzschild_EF__y_posn; \
  CCTK_REAL const geometry__Schwarzschild_EF__z_posn = PRIVATE_AHFINDERDIRECT_STRUCT.geometry__Schwarzschild_EF__z_posn; \
  CCTK_REAL const mask_buffer_thickness = PRIVATE_AHFINDERDIRECT_STRUCT.mask_buffer_thickness; \
  CCTK_REAL const mask_radius_multiplier = PRIVATE_AHFINDERDIRECT_STRUCT.mask_radius_multiplier; \
  CCTK_REAL const mask_radius_offset = PRIVATE_AHFINDERDIRECT_STRUCT.mask_radius_offset; \
  CCTK_REAL const max_allowable_Delta_h_over_h = PRIVATE_AHFINDERDIRECT_STRUCT.max_allowable_Delta_h_over_h; \
  CCTK_REAL const max_allowable_Theta = PRIVATE_AHFINDERDIRECT_STRUCT.max_allowable_Theta; \
  CCTK_REAL const * const max_allowable_horizon_radius = PRIVATE_AHFINDERDIRECT_STRUCT.max_allowable_horizon_radius; \
  CCTK_REAL const min_horizon_radius_points_for_mask = PRIVATE_AHFINDERDIRECT_STRUCT.min_horizon_radius_points_for_mask; \
  CCTK_REAL const old_style_mask_buffer_value = PRIVATE_AHFINDERDIRECT_STRUCT.old_style_mask_buffer_value; \
  CCTK_REAL const old_style_mask_inside_value = PRIVATE_AHFINDERDIRECT_STRUCT.old_style_mask_inside_value; \
  CCTK_REAL const old_style_mask_outside_value = PRIVATE_AHFINDERDIRECT_STRUCT.old_style_mask_outside_value; \
  CCTK_REAL const * const pretracking_delta = PRIVATE_AHFINDERDIRECT_STRUCT.pretracking_delta; \
  CCTK_REAL const * const pretracking_maximum_delta = PRIVATE_AHFINDERDIRECT_STRUCT.pretracking_maximum_delta; \
  CCTK_REAL const * const pretracking_maximum_value = PRIVATE_AHFINDERDIRECT_STRUCT.pretracking_maximum_value; \
  CCTK_REAL const * const pretracking_minimum_delta = PRIVATE_AHFINDERDIRECT_STRUCT.pretracking_minimum_delta; \
  CCTK_REAL const * const pretracking_minimum_value = PRIVATE_AHFINDERDIRECT_STRUCT.pretracking_minimum_value; \
  CCTK_REAL const * const pretracking_value = PRIVATE_AHFINDERDIRECT_STRUCT.pretracking_value; \
  CCTK_REAL const * const shiftout_factor = PRIVATE_AHFINDERDIRECT_STRUCT.shiftout_factor; \
  CCTK_REAL const * const smoothing_factor = PRIVATE_AHFINDERDIRECT_STRUCT.smoothing_factor; \
  const char * const ASCII_gnuplot_file_name_extension = PRIVATE_AHFINDERDIRECT_STRUCT.ASCII_gnuplot_file_name_extension; \
  const char * const BH_diagnostics_base_file_name = PRIVATE_AHFINDERDIRECT_STRUCT.BH_diagnostics_base_file_name; \
  const char * const BH_diagnostics_directory = PRIVATE_AHFINDERDIRECT_STRUCT.BH_diagnostics_directory; \
  const char * const BH_diagnostics_file_name_extension = PRIVATE_AHFINDERDIRECT_STRUCT.BH_diagnostics_file_name_extension; \
  const char * const Delta_h_base_file_name = PRIVATE_AHFINDERDIRECT_STRUCT.Delta_h_base_file_name; \
  const char * const HDF5_file_name_extension = PRIVATE_AHFINDERDIRECT_STRUCT.HDF5_file_name_extension; \
  const char * const Jacobian_base_file_name = PRIVATE_AHFINDERDIRECT_STRUCT.Jacobian_base_file_name; \
  const char * const Jacobian_compute_method = PRIVATE_AHFINDERDIRECT_STRUCT.Jacobian_compute_method; \
  const char * const Jacobian_store_solve_method = PRIVATE_AHFINDERDIRECT_STRUCT.Jacobian_store_solve_method; \
  const char * const OpenDX_control_file_name_extension = PRIVATE_AHFINDERDIRECT_STRUCT.OpenDX_control_file_name_extension; \
  const char * const Theta_base_file_name = PRIVATE_AHFINDERDIRECT_STRUCT.Theta_base_file_name; \
  const char * const coordinate_system_name = PRIVATE_AHFINDERDIRECT_STRUCT.coordinate_system_name; \
  const char * const geometry_interpolator_name = PRIVATE_AHFINDERDIRECT_STRUCT.geometry_interpolator_name; \
  const char * const geometry_interpolator_pars = PRIVATE_AHFINDERDIRECT_STRUCT.geometry_interpolator_pars; \
  const char * const h_base_file_name = PRIVATE_AHFINDERDIRECT_STRUCT.h_base_file_name; \
  const char * const h_directory = PRIVATE_AHFINDERDIRECT_STRUCT.h_directory; \
  const char * const integral_method = PRIVATE_AHFINDERDIRECT_STRUCT.integral_method; \
  const char * const interpatch_interpolator_name = PRIVATE_AHFINDERDIRECT_STRUCT.interpatch_interpolator_name; \
  const char * const interpatch_interpolator_pars = PRIVATE_AHFINDERDIRECT_STRUCT.interpatch_interpolator_pars; \
  const char * const mean_curvature_base_file_name = PRIVATE_AHFINDERDIRECT_STRUCT.mean_curvature_base_file_name; \
  const char * const method = PRIVATE_AHFINDERDIRECT_STRUCT.method; \
  const char * const metric_type = PRIVATE_AHFINDERDIRECT_STRUCT.metric_type; \
  const char * const new_style_mask_bitfield_name = PRIVATE_AHFINDERDIRECT_STRUCT.new_style_mask_bitfield_name; \
  const char * const new_style_mask_buffer_value = PRIVATE_AHFINDERDIRECT_STRUCT.new_style_mask_buffer_value; \
  const char * const new_style_mask_gridfn_name = PRIVATE_AHFINDERDIRECT_STRUCT.new_style_mask_gridfn_name; \
  const char * const new_style_mask_inside_value = PRIVATE_AHFINDERDIRECT_STRUCT.new_style_mask_inside_value; \
  const char * const new_style_mask_outside_value = PRIVATE_AHFINDERDIRECT_STRUCT.new_style_mask_outside_value; \
  const char * const old_style_mask_gridfn_name = PRIVATE_AHFINDERDIRECT_STRUCT.old_style_mask_gridfn_name; \
  const char * const * const patch_system_type = PRIVATE_AHFINDERDIRECT_STRUCT.patch_system_type; \
  const char * const * const surface_definition = PRIVATE_AHFINDERDIRECT_STRUCT.surface_definition; \
  const char * const surface_interpolator_name = PRIVATE_AHFINDERDIRECT_STRUCT.surface_interpolator_name; \
  const char * const surface_interpolator_pars = PRIVATE_AHFINDERDIRECT_STRUCT.surface_interpolator_pars; \
  const char * const * const surface_modification = PRIVATE_AHFINDERDIRECT_STRUCT.surface_modification; \
  const char * const * const surface_selection = PRIVATE_AHFINDERDIRECT_STRUCT.surface_selection; \
  const char * const verbose_level = PRIVATE_AHFINDERDIRECT_STRUCT.verbose_level; \
  CCTK_INT const ILUCG__limit_CG_iterations = PRIVATE_AHFINDERDIRECT_STRUCT.ILUCG__limit_CG_iterations; \
  CCTK_INT const N_horizons = PRIVATE_AHFINDERDIRECT_STRUCT.N_horizons; \
  CCTK_INT const * const N_zones_per_right_angle = PRIVATE_AHFINDERDIRECT_STRUCT.N_zones_per_right_angle; \
  CCTK_INT const UMFPACK__N_II_iterations = PRIVATE_AHFINDERDIRECT_STRUCT.UMFPACK__N_II_iterations; \
  CCTK_INT const check_that_geometry_is_finite = PRIVATE_AHFINDERDIRECT_STRUCT.check_that_geometry_is_finite; \
  CCTK_INT const check_that_h_is_finite = PRIVATE_AHFINDERDIRECT_STRUCT.check_that_h_is_finite; \
  CCTK_INT const debugging_output_at_each_Newton_iteration = PRIVATE_AHFINDERDIRECT_STRUCT.debugging_output_at_each_Newton_iteration; \
  CCTK_INT const * const depends_on = PRIVATE_AHFINDERDIRECT_STRUCT.depends_on; \
  CCTK_INT const * const disable_horizon = PRIVATE_AHFINDERDIRECT_STRUCT.disable_horizon; \
  CCTK_INT const * const dont_find_after_individual = PRIVATE_AHFINDERDIRECT_STRUCT.dont_find_after_individual; \
  CCTK_INT const * const find_after_individual = PRIVATE_AHFINDERDIRECT_STRUCT.find_after_individual; \
  CCTK_INT const find_every = PRIVATE_AHFINDERDIRECT_STRUCT.find_every; \
  CCTK_INT const * const find_every_individual = PRIVATE_AHFINDERDIRECT_STRUCT.find_every_individual; \
  CCTK_INT const ghost_zone_width = PRIVATE_AHFINDERDIRECT_STRUCT.ghost_zone_width; \
  CCTK_INT const h_min_digits = PRIVATE_AHFINDERDIRECT_STRUCT.h_min_digits; \
  CCTK_INT const hardwire_Schwarzschild_EF_geometry = PRIVATE_AHFINDERDIRECT_STRUCT.hardwire_Schwarzschild_EF_geometry; \
  CCTK_INT const mask_is_noshrink = PRIVATE_AHFINDERDIRECT_STRUCT.mask_is_noshrink; \
  CCTK_INT const max_N_zones_per_right_angle = PRIVATE_AHFINDERDIRECT_STRUCT.max_N_zones_per_right_angle; \
  CCTK_INT const max_Newton_iterations__initial = PRIVATE_AHFINDERDIRECT_STRUCT.max_Newton_iterations__initial; \
  CCTK_INT const max_Newton_iterations__subsequent = PRIVATE_AHFINDERDIRECT_STRUCT.max_Newton_iterations__subsequent; \
  CCTK_INT const max_allowable_Theta_growth_iterations = PRIVATE_AHFINDERDIRECT_STRUCT.max_allowable_Theta_growth_iterations; \
  CCTK_INT const max_allowable_Theta_nonshrink_iterations = PRIVATE_AHFINDERDIRECT_STRUCT.max_allowable_Theta_nonshrink_iterations; \
  CCTK_INT const output_ASCII_files = PRIVATE_AHFINDERDIRECT_STRUCT.output_ASCII_files; \
  CCTK_INT const output_BH_diagnostics = PRIVATE_AHFINDERDIRECT_STRUCT.output_BH_diagnostics; \
  CCTK_INT const output_HDF5_files = PRIVATE_AHFINDERDIRECT_STRUCT.output_HDF5_files; \
  CCTK_INT const output_OpenDX_control_files = PRIVATE_AHFINDERDIRECT_STRUCT.output_OpenDX_control_files; \
  CCTK_INT const output_Theta_every = PRIVATE_AHFINDERDIRECT_STRUCT.output_Theta_every; \
  CCTK_INT const output_ghost_zones_for_h = PRIVATE_AHFINDERDIRECT_STRUCT.output_ghost_zones_for_h; \
  CCTK_INT const output_h_every = PRIVATE_AHFINDERDIRECT_STRUCT.output_h_every; \
  CCTK_INT const output_initial_guess = PRIVATE_AHFINDERDIRECT_STRUCT.output_initial_guess; \
  CCTK_INT const output_mean_curvature_every = PRIVATE_AHFINDERDIRECT_STRUCT.output_mean_curvature_every; \
  CCTK_INT const patch_overlap_width = PRIVATE_AHFINDERDIRECT_STRUCT.patch_overlap_width; \
  CCTK_INT const * const pretracking_max_iterations = PRIVATE_AHFINDERDIRECT_STRUCT.pretracking_max_iterations; \
  CCTK_INT const print_timing_stats = PRIVATE_AHFINDERDIRECT_STRUCT.print_timing_stats; \
  CCTK_INT const run_at_CCTK_ANALYSIS = PRIVATE_AHFINDERDIRECT_STRUCT.run_at_CCTK_ANALYSIS; \
  CCTK_INT const run_at_CCTK_POSTINITIAL = PRIVATE_AHFINDERDIRECT_STRUCT.run_at_CCTK_POSTINITIAL; \
  CCTK_INT const run_at_CCTK_POSTSTEP = PRIVATE_AHFINDERDIRECT_STRUCT.run_at_CCTK_POSTSTEP; \
  CCTK_INT const run_at_CCTK_POST_RECOVER_VARIABLES = PRIVATE_AHFINDERDIRECT_STRUCT.run_at_CCTK_POST_RECOVER_VARIABLES; \
  CCTK_INT const set_mask_for_all_horizons = PRIVATE_AHFINDERDIRECT_STRUCT.set_mask_for_all_horizons; \
  CCTK_INT const * const set_mask_for_individual_horizon = PRIVATE_AHFINDERDIRECT_STRUCT.set_mask_for_individual_horizon; \
  CCTK_INT const set_new_style_mask = PRIVATE_AHFINDERDIRECT_STRUCT.set_new_style_mask; \
  CCTK_INT const set_old_style_mask = PRIVATE_AHFINDERDIRECT_STRUCT.set_old_style_mask; \
  CCTK_INT const test_all_Jacobian_compute_methods = PRIVATE_AHFINDERDIRECT_STRUCT.test_all_Jacobian_compute_methods; \
  CCTK_INT const * const use_pretracking = PRIVATE_AHFINDERDIRECT_STRUCT.use_pretracking; \
  CCTK_INT const want_expansion_gradients = PRIVATE_AHFINDERDIRECT_STRUCT.want_expansion_gradients; \
  CCTK_INT const warn_level__gij_not_positive_definite__initial = PRIVATE_AHFINDERDIRECT_STRUCT.warn_level__gij_not_positive_definite__initial; \
  CCTK_INT const warn_level__gij_not_positive_definite__subsequent = PRIVATE_AHFINDERDIRECT_STRUCT.warn_level__gij_not_positive_definite__subsequent; \
  CCTK_INT const warn_level__nonfinite_geometry = PRIVATE_AHFINDERDIRECT_STRUCT.warn_level__nonfinite_geometry; \
  CCTK_INT const warn_level__point_outside__initial = PRIVATE_AHFINDERDIRECT_STRUCT.warn_level__point_outside__initial; \
  CCTK_INT const warn_level__point_outside__subsequent = PRIVATE_AHFINDERDIRECT_STRUCT.warn_level__point_outside__subsequent; \
  CCTK_INT const warn_level__skipping_finite_check = PRIVATE_AHFINDERDIRECT_STRUCT.warn_level__skipping_finite_check; \
  CCTK_INT const which_horizon_to_announce_centroid = PRIVATE_AHFINDERDIRECT_STRUCT.which_horizon_to_announce_centroid; \
  CCTK_INT const * const which_surface_to_store_info = PRIVATE_AHFINDERDIRECT_STRUCT.which_surface_to_store_info; \
  enum { \
      dummy_PRIVATE_AHFINDERDIRECT_STRUCT_ILUCG__error_tolerance = sizeof( ILUCG__error_tolerance ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_Jacobian_perturbation_amplitude = sizeof( Jacobian_perturbation_amplitude ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_Theta_norm_for_convergence = sizeof( Theta_norm_for_convergence ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_desired_value = sizeof( desired_value ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_desired_value_factor = sizeof( desired_value_factor ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_desired_value_offset = sizeof( desired_value_offset ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_dont_find_after_individual_time = sizeof( dont_find_after_individual_time ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_find_after_individual_time = sizeof( find_after_individual_time ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_geometry__Schwarzschild_EF__Delta_xyz = sizeof( geometry__Schwarzschild_EF__Delta_xyz ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_geometry__Schwarzschild_EF__epsilon = sizeof( geometry__Schwarzschild_EF__epsilon ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_geometry__Schwarzschild_EF__mass = sizeof( geometry__Schwarzschild_EF__mass ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_geometry__Schwarzschild_EF__x_posn = sizeof( geometry__Schwarzschild_EF__x_posn ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_geometry__Schwarzschild_EF__y_posn = sizeof( geometry__Schwarzschild_EF__y_posn ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_geometry__Schwarzschild_EF__z_posn = sizeof( geometry__Schwarzschild_EF__z_posn ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_mask_buffer_thickness = sizeof( mask_buffer_thickness ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_mask_radius_multiplier = sizeof( mask_radius_multiplier ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_mask_radius_offset = sizeof( mask_radius_offset ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_max_allowable_Delta_h_over_h = sizeof( max_allowable_Delta_h_over_h ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_max_allowable_Theta = sizeof( max_allowable_Theta ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_max_allowable_horizon_radius = sizeof( max_allowable_horizon_radius ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_min_horizon_radius_points_for_mask = sizeof( min_horizon_radius_points_for_mask ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_old_style_mask_buffer_value = sizeof( old_style_mask_buffer_value ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_old_style_mask_inside_value = sizeof( old_style_mask_inside_value ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_old_style_mask_outside_value = sizeof( old_style_mask_outside_value ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_pretracking_delta = sizeof( pretracking_delta ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_pretracking_maximum_delta = sizeof( pretracking_maximum_delta ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_pretracking_maximum_value = sizeof( pretracking_maximum_value ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_pretracking_minimum_delta = sizeof( pretracking_minimum_delta ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_pretracking_minimum_value = sizeof( pretracking_minimum_value ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_pretracking_value = sizeof( pretracking_value ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_shiftout_factor = sizeof( shiftout_factor ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_smoothing_factor = sizeof( smoothing_factor ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_ASCII_gnuplot_file_name_extension = sizeof( ASCII_gnuplot_file_name_extension ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_BH_diagnostics_base_file_name = sizeof( BH_diagnostics_base_file_name ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_BH_diagnostics_directory = sizeof( BH_diagnostics_directory ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_BH_diagnostics_file_name_extension = sizeof( BH_diagnostics_file_name_extension ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_Delta_h_base_file_name = sizeof( Delta_h_base_file_name ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_HDF5_file_name_extension = sizeof( HDF5_file_name_extension ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_Jacobian_base_file_name = sizeof( Jacobian_base_file_name ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_Jacobian_compute_method = sizeof( Jacobian_compute_method ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_Jacobian_store_solve_method = sizeof( Jacobian_store_solve_method ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_OpenDX_control_file_name_extension = sizeof( OpenDX_control_file_name_extension ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_Theta_base_file_name = sizeof( Theta_base_file_name ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_coordinate_system_name = sizeof( coordinate_system_name ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_geometry_interpolator_name = sizeof( geometry_interpolator_name ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_geometry_interpolator_pars = sizeof( geometry_interpolator_pars ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_h_base_file_name = sizeof( h_base_file_name ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_h_directory = sizeof( h_directory ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_integral_method = sizeof( integral_method ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_interpatch_interpolator_name = sizeof( interpatch_interpolator_name ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_interpatch_interpolator_pars = sizeof( interpatch_interpolator_pars ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_mean_curvature_base_file_name = sizeof( mean_curvature_base_file_name ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_method = sizeof( method ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_metric_type = sizeof( metric_type ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_new_style_mask_bitfield_name = sizeof( new_style_mask_bitfield_name ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_new_style_mask_buffer_value = sizeof( new_style_mask_buffer_value ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_new_style_mask_gridfn_name = sizeof( new_style_mask_gridfn_name ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_new_style_mask_inside_value = sizeof( new_style_mask_inside_value ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_new_style_mask_outside_value = sizeof( new_style_mask_outside_value ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_old_style_mask_gridfn_name = sizeof( old_style_mask_gridfn_name ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_patch_system_type = sizeof( patch_system_type ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_surface_definition = sizeof( surface_definition ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_surface_interpolator_name = sizeof( surface_interpolator_name ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_surface_interpolator_pars = sizeof( surface_interpolator_pars ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_surface_modification = sizeof( surface_modification ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_surface_selection = sizeof( surface_selection ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_verbose_level = sizeof( verbose_level ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_ILUCG__limit_CG_iterations = sizeof( ILUCG__limit_CG_iterations ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_N_horizons = sizeof( N_horizons ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_N_zones_per_right_angle = sizeof( N_zones_per_right_angle ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_UMFPACK__N_II_iterations = sizeof( UMFPACK__N_II_iterations ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_check_that_geometry_is_finite = sizeof( check_that_geometry_is_finite ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_check_that_h_is_finite = sizeof( check_that_h_is_finite ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_debugging_output_at_each_Newton_iteration = sizeof( debugging_output_at_each_Newton_iteration ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_depends_on = sizeof( depends_on ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_disable_horizon = sizeof( disable_horizon ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_dont_find_after_individual = sizeof( dont_find_after_individual ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_find_after_individual = sizeof( find_after_individual ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_find_every = sizeof( find_every ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_find_every_individual = sizeof( find_every_individual ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_ghost_zone_width = sizeof( ghost_zone_width ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_h_min_digits = sizeof( h_min_digits ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_hardwire_Schwarzschild_EF_geometry = sizeof( hardwire_Schwarzschild_EF_geometry ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_mask_is_noshrink = sizeof( mask_is_noshrink ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_max_N_zones_per_right_angle = sizeof( max_N_zones_per_right_angle ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_max_Newton_iterations__initial = sizeof( max_Newton_iterations__initial ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_max_Newton_iterations__subsequent = sizeof( max_Newton_iterations__subsequent ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_max_allowable_Theta_growth_iterations = sizeof( max_allowable_Theta_growth_iterations ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_max_allowable_Theta_nonshrink_iterations = sizeof( max_allowable_Theta_nonshrink_iterations ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_output_ASCII_files = sizeof( output_ASCII_files ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_output_BH_diagnostics = sizeof( output_BH_diagnostics ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_output_HDF5_files = sizeof( output_HDF5_files ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_output_OpenDX_control_files = sizeof( output_OpenDX_control_files ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_output_Theta_every = sizeof( output_Theta_every ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_output_ghost_zones_for_h = sizeof( output_ghost_zones_for_h ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_output_h_every = sizeof( output_h_every ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_output_initial_guess = sizeof( output_initial_guess ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_output_mean_curvature_every = sizeof( output_mean_curvature_every ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_patch_overlap_width = sizeof( patch_overlap_width ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_pretracking_max_iterations = sizeof( pretracking_max_iterations ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_print_timing_stats = sizeof( print_timing_stats ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_run_at_CCTK_ANALYSIS = sizeof( run_at_CCTK_ANALYSIS ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_run_at_CCTK_POSTINITIAL = sizeof( run_at_CCTK_POSTINITIAL ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_run_at_CCTK_POSTSTEP = sizeof( run_at_CCTK_POSTSTEP ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_run_at_CCTK_POST_RECOVER_VARIABLES = sizeof( run_at_CCTK_POST_RECOVER_VARIABLES ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_set_mask_for_all_horizons = sizeof( set_mask_for_all_horizons ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_set_mask_for_individual_horizon = sizeof( set_mask_for_individual_horizon ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_set_new_style_mask = sizeof( set_new_style_mask ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_set_old_style_mask = sizeof( set_old_style_mask ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_test_all_Jacobian_compute_methods = sizeof( test_all_Jacobian_compute_methods ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_use_pretracking = sizeof( use_pretracking ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_want_expansion_gradients = sizeof( want_expansion_gradients ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_warn_level__gij_not_positive_definite__initial = sizeof( warn_level__gij_not_positive_definite__initial ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_warn_level__gij_not_positive_definite__subsequent = sizeof( warn_level__gij_not_positive_definite__subsequent ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_warn_level__nonfinite_geometry = sizeof( warn_level__nonfinite_geometry ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_warn_level__point_outside__initial = sizeof( warn_level__point_outside__initial ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_warn_level__point_outside__subsequent = sizeof( warn_level__point_outside__subsequent ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_warn_level__skipping_finite_check = sizeof( warn_level__skipping_finite_check ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_which_horizon_to_announce_centroid = sizeof( which_horizon_to_announce_centroid ) \
    , dummy_PRIVATE_AHFINDERDIRECT_STRUCT_which_surface_to_store_info = sizeof( which_surface_to_store_info ) \
  }; \

