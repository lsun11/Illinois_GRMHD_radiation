#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL aspect_ratio_x;
  CCTK_REAL aspect_ratio_y;
  CCTK_REAL aspect_ratio_z;
  const char * base_extents;
  const char * base_outerbounds;
  const char * grid_coordinates_filename;
  const char * grid_structure_filename;
  const char * model;
  const char * processor_topology;
  const char * refinement_centering;
  const char * space_refinement_factors;
  const char * time_refinement_factors;
  const char * timer_file;
  CCTK_INT adaptive_stepsize;
  CCTK_INT additional_buffer_zones;
  CCTK_INT barriers;
  CCTK_INT check_for_poison;
  CCTK_INT checksum_timelevels;
  CCTK_INT constant_load_per_processor;
  CCTK_INT convergence_factor;
  CCTK_INT convergence_level;
  CCTK_INT deadbeef;
  CCTK_INT domain_from_coordbase;
  CCTK_INT domain_from_multipatch;
  CCTK_INT enable_all_storage;
  CCTK_INT ghost_size;
  CCTK_INT ghost_size_x;
  CCTK_INT ghost_size_y;
  CCTK_INT ghost_size_z;
  CCTK_INT global_nsize;
  CCTK_INT global_nx;
  CCTK_INT global_ny;
  CCTK_INT global_nz;
  CCTK_INT init_3_timelevels;
  CCTK_INT init_each_timelevel;
  CCTK_INT init_fill_timelevels;
  CCTK_INT max_core_size_MB;
  CCTK_INT max_memory_size_MB;
  CCTK_INT max_poison_locations;
  CCTK_INT max_refinement_levels;
  CCTK_INT min_points_per_proc;
  CCTK_INT num_convergence_levels;
  CCTK_INT num_integrator_substeps;
  CCTK_INT num_maps;
  CCTK_INT num_threads;
  CCTK_INT output_internal_data;
  CCTK_INT output_timers_every;
  CCTK_INT poison_new_timelevels;
  CCTK_INT poison_value;
  CCTK_INT print_timestats_every;
  CCTK_INT processor_topology_3d_x;
  CCTK_INT processor_topology_3d_y;
  CCTK_INT processor_topology_3d_z;
  CCTK_INT prolongate_initial_data;
  CCTK_INT prolongation_order_space;
  CCTK_INT prolongation_order_time;
  CCTK_INT refine_timestep;
  CCTK_INT refinement_factor;
  CCTK_INT regrid_during_initialisation;
  CCTK_INT regrid_during_recovery;
  CCTK_INT regrid_in_level_mode;
  CCTK_INT schedule_barriers;
  CCTK_INT split_direction;
  CCTK_INT storage_verbose;
  CCTK_INT suppress_restriction;
  CCTK_INT sync_during_time_integration;
  CCTK_INT use_buffer_zones;
  CCTK_INT use_tapered_grids;
  CCTK_INT verbose;
  CCTK_INT veryverbose;
} PRIVATE_CARPET_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_CARPET_STRUCT_PARAMS \
  CCTK_REAL const aspect_ratio_x = PRIVATE_CARPET_STRUCT.aspect_ratio_x; \
  CCTK_REAL const aspect_ratio_y = PRIVATE_CARPET_STRUCT.aspect_ratio_y; \
  CCTK_REAL const aspect_ratio_z = PRIVATE_CARPET_STRUCT.aspect_ratio_z; \
  const char * const base_extents = PRIVATE_CARPET_STRUCT.base_extents; \
  const char * const base_outerbounds = PRIVATE_CARPET_STRUCT.base_outerbounds; \
  const char * const grid_coordinates_filename = PRIVATE_CARPET_STRUCT.grid_coordinates_filename; \
  const char * const grid_structure_filename = PRIVATE_CARPET_STRUCT.grid_structure_filename; \
  const char * const model = PRIVATE_CARPET_STRUCT.model; \
  const char * const processor_topology = PRIVATE_CARPET_STRUCT.processor_topology; \
  const char * const refinement_centering = PRIVATE_CARPET_STRUCT.refinement_centering; \
  const char * const space_refinement_factors = PRIVATE_CARPET_STRUCT.space_refinement_factors; \
  const char * const time_refinement_factors = PRIVATE_CARPET_STRUCT.time_refinement_factors; \
  const char * const timer_file = PRIVATE_CARPET_STRUCT.timer_file; \
  CCTK_INT const adaptive_stepsize = PRIVATE_CARPET_STRUCT.adaptive_stepsize; \
  CCTK_INT const additional_buffer_zones = PRIVATE_CARPET_STRUCT.additional_buffer_zones; \
  CCTK_INT const barriers = PRIVATE_CARPET_STRUCT.barriers; \
  CCTK_INT const check_for_poison = PRIVATE_CARPET_STRUCT.check_for_poison; \
  CCTK_INT const checksum_timelevels = PRIVATE_CARPET_STRUCT.checksum_timelevels; \
  CCTK_INT const constant_load_per_processor = PRIVATE_CARPET_STRUCT.constant_load_per_processor; \
  CCTK_INT const convergence_factor = PRIVATE_CARPET_STRUCT.convergence_factor; \
  CCTK_INT const convergence_level = PRIVATE_CARPET_STRUCT.convergence_level; \
  CCTK_INT const deadbeef = PRIVATE_CARPET_STRUCT.deadbeef; \
  CCTK_INT const domain_from_coordbase = PRIVATE_CARPET_STRUCT.domain_from_coordbase; \
  CCTK_INT const domain_from_multipatch = PRIVATE_CARPET_STRUCT.domain_from_multipatch; \
  CCTK_INT const enable_all_storage = PRIVATE_CARPET_STRUCT.enable_all_storage; \
  CCTK_INT const ghost_size = PRIVATE_CARPET_STRUCT.ghost_size; \
  CCTK_INT const ghost_size_x = PRIVATE_CARPET_STRUCT.ghost_size_x; \
  CCTK_INT const ghost_size_y = PRIVATE_CARPET_STRUCT.ghost_size_y; \
  CCTK_INT const ghost_size_z = PRIVATE_CARPET_STRUCT.ghost_size_z; \
  CCTK_INT const global_nsize = PRIVATE_CARPET_STRUCT.global_nsize; \
  CCTK_INT const global_nx = PRIVATE_CARPET_STRUCT.global_nx; \
  CCTK_INT const global_ny = PRIVATE_CARPET_STRUCT.global_ny; \
  CCTK_INT const global_nz = PRIVATE_CARPET_STRUCT.global_nz; \
  CCTK_INT const init_3_timelevels = PRIVATE_CARPET_STRUCT.init_3_timelevels; \
  CCTK_INT const init_each_timelevel = PRIVATE_CARPET_STRUCT.init_each_timelevel; \
  CCTK_INT const init_fill_timelevels = PRIVATE_CARPET_STRUCT.init_fill_timelevels; \
  CCTK_INT const max_core_size_MB = PRIVATE_CARPET_STRUCT.max_core_size_MB; \
  CCTK_INT const max_memory_size_MB = PRIVATE_CARPET_STRUCT.max_memory_size_MB; \
  CCTK_INT const max_poison_locations = PRIVATE_CARPET_STRUCT.max_poison_locations; \
  CCTK_INT const max_refinement_levels = PRIVATE_CARPET_STRUCT.max_refinement_levels; \
  CCTK_INT const min_points_per_proc = PRIVATE_CARPET_STRUCT.min_points_per_proc; \
  CCTK_INT const num_convergence_levels = PRIVATE_CARPET_STRUCT.num_convergence_levels; \
  CCTK_INT const num_integrator_substeps = PRIVATE_CARPET_STRUCT.num_integrator_substeps; \
  CCTK_INT const num_maps = PRIVATE_CARPET_STRUCT.num_maps; \
  CCTK_INT const num_threads = PRIVATE_CARPET_STRUCT.num_threads; \
  CCTK_INT const output_internal_data = PRIVATE_CARPET_STRUCT.output_internal_data; \
  CCTK_INT const output_timers_every = PRIVATE_CARPET_STRUCT.output_timers_every; \
  CCTK_INT const poison_new_timelevels = PRIVATE_CARPET_STRUCT.poison_new_timelevels; \
  CCTK_INT const poison_value = PRIVATE_CARPET_STRUCT.poison_value; \
  CCTK_INT const print_timestats_every = PRIVATE_CARPET_STRUCT.print_timestats_every; \
  CCTK_INT const processor_topology_3d_x = PRIVATE_CARPET_STRUCT.processor_topology_3d_x; \
  CCTK_INT const processor_topology_3d_y = PRIVATE_CARPET_STRUCT.processor_topology_3d_y; \
  CCTK_INT const processor_topology_3d_z = PRIVATE_CARPET_STRUCT.processor_topology_3d_z; \
  CCTK_INT const prolongate_initial_data = PRIVATE_CARPET_STRUCT.prolongate_initial_data; \
  CCTK_INT const prolongation_order_space = PRIVATE_CARPET_STRUCT.prolongation_order_space; \
  CCTK_INT const prolongation_order_time = PRIVATE_CARPET_STRUCT.prolongation_order_time; \
  CCTK_INT const refine_timestep = PRIVATE_CARPET_STRUCT.refine_timestep; \
  CCTK_INT const refinement_factor = PRIVATE_CARPET_STRUCT.refinement_factor; \
  CCTK_INT const regrid_during_initialisation = PRIVATE_CARPET_STRUCT.regrid_during_initialisation; \
  CCTK_INT const regrid_during_recovery = PRIVATE_CARPET_STRUCT.regrid_during_recovery; \
  CCTK_INT const regrid_in_level_mode = PRIVATE_CARPET_STRUCT.regrid_in_level_mode; \
  CCTK_INT const schedule_barriers = PRIVATE_CARPET_STRUCT.schedule_barriers; \
  CCTK_INT const split_direction = PRIVATE_CARPET_STRUCT.split_direction; \
  CCTK_INT const storage_verbose = PRIVATE_CARPET_STRUCT.storage_verbose; \
  CCTK_INT const suppress_restriction = PRIVATE_CARPET_STRUCT.suppress_restriction; \
  CCTK_INT const sync_during_time_integration = PRIVATE_CARPET_STRUCT.sync_during_time_integration; \
  CCTK_INT const use_buffer_zones = PRIVATE_CARPET_STRUCT.use_buffer_zones; \
  CCTK_INT const use_tapered_grids = PRIVATE_CARPET_STRUCT.use_tapered_grids; \
  CCTK_INT const verbose = PRIVATE_CARPET_STRUCT.verbose; \
  CCTK_INT const veryverbose = PRIVATE_CARPET_STRUCT.veryverbose; \
  enum { \
      dummy_PRIVATE_CARPET_STRUCT_aspect_ratio_x = sizeof( aspect_ratio_x ) \
    , dummy_PRIVATE_CARPET_STRUCT_aspect_ratio_y = sizeof( aspect_ratio_y ) \
    , dummy_PRIVATE_CARPET_STRUCT_aspect_ratio_z = sizeof( aspect_ratio_z ) \
    , dummy_PRIVATE_CARPET_STRUCT_base_extents = sizeof( base_extents ) \
    , dummy_PRIVATE_CARPET_STRUCT_base_outerbounds = sizeof( base_outerbounds ) \
    , dummy_PRIVATE_CARPET_STRUCT_grid_coordinates_filename = sizeof( grid_coordinates_filename ) \
    , dummy_PRIVATE_CARPET_STRUCT_grid_structure_filename = sizeof( grid_structure_filename ) \
    , dummy_PRIVATE_CARPET_STRUCT_model = sizeof( model ) \
    , dummy_PRIVATE_CARPET_STRUCT_processor_topology = sizeof( processor_topology ) \
    , dummy_PRIVATE_CARPET_STRUCT_refinement_centering = sizeof( refinement_centering ) \
    , dummy_PRIVATE_CARPET_STRUCT_space_refinement_factors = sizeof( space_refinement_factors ) \
    , dummy_PRIVATE_CARPET_STRUCT_time_refinement_factors = sizeof( time_refinement_factors ) \
    , dummy_PRIVATE_CARPET_STRUCT_timer_file = sizeof( timer_file ) \
    , dummy_PRIVATE_CARPET_STRUCT_adaptive_stepsize = sizeof( adaptive_stepsize ) \
    , dummy_PRIVATE_CARPET_STRUCT_additional_buffer_zones = sizeof( additional_buffer_zones ) \
    , dummy_PRIVATE_CARPET_STRUCT_barriers = sizeof( barriers ) \
    , dummy_PRIVATE_CARPET_STRUCT_check_for_poison = sizeof( check_for_poison ) \
    , dummy_PRIVATE_CARPET_STRUCT_checksum_timelevels = sizeof( checksum_timelevels ) \
    , dummy_PRIVATE_CARPET_STRUCT_constant_load_per_processor = sizeof( constant_load_per_processor ) \
    , dummy_PRIVATE_CARPET_STRUCT_convergence_factor = sizeof( convergence_factor ) \
    , dummy_PRIVATE_CARPET_STRUCT_convergence_level = sizeof( convergence_level ) \
    , dummy_PRIVATE_CARPET_STRUCT_deadbeef = sizeof( deadbeef ) \
    , dummy_PRIVATE_CARPET_STRUCT_domain_from_coordbase = sizeof( domain_from_coordbase ) \
    , dummy_PRIVATE_CARPET_STRUCT_domain_from_multipatch = sizeof( domain_from_multipatch ) \
    , dummy_PRIVATE_CARPET_STRUCT_enable_all_storage = sizeof( enable_all_storage ) \
    , dummy_PRIVATE_CARPET_STRUCT_ghost_size = sizeof( ghost_size ) \
    , dummy_PRIVATE_CARPET_STRUCT_ghost_size_x = sizeof( ghost_size_x ) \
    , dummy_PRIVATE_CARPET_STRUCT_ghost_size_y = sizeof( ghost_size_y ) \
    , dummy_PRIVATE_CARPET_STRUCT_ghost_size_z = sizeof( ghost_size_z ) \
    , dummy_PRIVATE_CARPET_STRUCT_global_nsize = sizeof( global_nsize ) \
    , dummy_PRIVATE_CARPET_STRUCT_global_nx = sizeof( global_nx ) \
    , dummy_PRIVATE_CARPET_STRUCT_global_ny = sizeof( global_ny ) \
    , dummy_PRIVATE_CARPET_STRUCT_global_nz = sizeof( global_nz ) \
    , dummy_PRIVATE_CARPET_STRUCT_init_3_timelevels = sizeof( init_3_timelevels ) \
    , dummy_PRIVATE_CARPET_STRUCT_init_each_timelevel = sizeof( init_each_timelevel ) \
    , dummy_PRIVATE_CARPET_STRUCT_init_fill_timelevels = sizeof( init_fill_timelevels ) \
    , dummy_PRIVATE_CARPET_STRUCT_max_core_size_MB = sizeof( max_core_size_MB ) \
    , dummy_PRIVATE_CARPET_STRUCT_max_memory_size_MB = sizeof( max_memory_size_MB ) \
    , dummy_PRIVATE_CARPET_STRUCT_max_poison_locations = sizeof( max_poison_locations ) \
    , dummy_PRIVATE_CARPET_STRUCT_max_refinement_levels = sizeof( max_refinement_levels ) \
    , dummy_PRIVATE_CARPET_STRUCT_min_points_per_proc = sizeof( min_points_per_proc ) \
    , dummy_PRIVATE_CARPET_STRUCT_num_convergence_levels = sizeof( num_convergence_levels ) \
    , dummy_PRIVATE_CARPET_STRUCT_num_integrator_substeps = sizeof( num_integrator_substeps ) \
    , dummy_PRIVATE_CARPET_STRUCT_num_maps = sizeof( num_maps ) \
    , dummy_PRIVATE_CARPET_STRUCT_num_threads = sizeof( num_threads ) \
    , dummy_PRIVATE_CARPET_STRUCT_output_internal_data = sizeof( output_internal_data ) \
    , dummy_PRIVATE_CARPET_STRUCT_output_timers_every = sizeof( output_timers_every ) \
    , dummy_PRIVATE_CARPET_STRUCT_poison_new_timelevels = sizeof( poison_new_timelevels ) \
    , dummy_PRIVATE_CARPET_STRUCT_poison_value = sizeof( poison_value ) \
    , dummy_PRIVATE_CARPET_STRUCT_print_timestats_every = sizeof( print_timestats_every ) \
    , dummy_PRIVATE_CARPET_STRUCT_processor_topology_3d_x = sizeof( processor_topology_3d_x ) \
    , dummy_PRIVATE_CARPET_STRUCT_processor_topology_3d_y = sizeof( processor_topology_3d_y ) \
    , dummy_PRIVATE_CARPET_STRUCT_processor_topology_3d_z = sizeof( processor_topology_3d_z ) \
    , dummy_PRIVATE_CARPET_STRUCT_prolongate_initial_data = sizeof( prolongate_initial_data ) \
    , dummy_PRIVATE_CARPET_STRUCT_prolongation_order_space = sizeof( prolongation_order_space ) \
    , dummy_PRIVATE_CARPET_STRUCT_prolongation_order_time = sizeof( prolongation_order_time ) \
    , dummy_PRIVATE_CARPET_STRUCT_refine_timestep = sizeof( refine_timestep ) \
    , dummy_PRIVATE_CARPET_STRUCT_refinement_factor = sizeof( refinement_factor ) \
    , dummy_PRIVATE_CARPET_STRUCT_regrid_during_initialisation = sizeof( regrid_during_initialisation ) \
    , dummy_PRIVATE_CARPET_STRUCT_regrid_during_recovery = sizeof( regrid_during_recovery ) \
    , dummy_PRIVATE_CARPET_STRUCT_regrid_in_level_mode = sizeof( regrid_in_level_mode ) \
    , dummy_PRIVATE_CARPET_STRUCT_schedule_barriers = sizeof( schedule_barriers ) \
    , dummy_PRIVATE_CARPET_STRUCT_split_direction = sizeof( split_direction ) \
    , dummy_PRIVATE_CARPET_STRUCT_storage_verbose = sizeof( storage_verbose ) \
    , dummy_PRIVATE_CARPET_STRUCT_suppress_restriction = sizeof( suppress_restriction ) \
    , dummy_PRIVATE_CARPET_STRUCT_sync_during_time_integration = sizeof( sync_during_time_integration ) \
    , dummy_PRIVATE_CARPET_STRUCT_use_buffer_zones = sizeof( use_buffer_zones ) \
    , dummy_PRIVATE_CARPET_STRUCT_use_tapered_grids = sizeof( use_tapered_grids ) \
    , dummy_PRIVATE_CARPET_STRUCT_verbose = sizeof( verbose ) \
    , dummy_PRIVATE_CARPET_STRUCT_veryverbose = sizeof( veryverbose ) \
  }; \

