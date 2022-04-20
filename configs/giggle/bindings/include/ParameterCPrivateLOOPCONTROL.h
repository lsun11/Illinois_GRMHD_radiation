#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL delayed_overhead_threshold;
  CCTK_REAL immediate_overhead_threshold;
  CCTK_REAL maximum_setup_overhead;
  CCTK_REAL probability_random_jump;
  CCTK_REAL probability_small_jump;
  CCTK_REAL siman_T_initial;
  CCTK_REAL siman_T_min;
  CCTK_REAL siman_k;
  CCTK_REAL siman_mu_T;
  CCTK_REAL siman_probability_change_topology;
  CCTK_REAL siman_step_size;
  CCTK_INT cycle_j_tilings;
  CCTK_INT debug;
  CCTK_INT ignore_initial_overhead;
  CCTK_INT lc_inpoints;
  CCTK_INT lc_inthreads;
  CCTK_INT lc_jnpoints;
  CCTK_INT lc_jnthreads;
  CCTK_INT lc_knpoints;
  CCTK_INT lc_knthreads;
  CCTK_INT legacy_init;
  CCTK_INT max_jump_attempts;
  CCTK_INT nsteps;
  CCTK_INT nx;
  CCTK_INT ny;
  CCTK_INT nz;
  CCTK_INT overhead_threshold_delay;
  CCTK_INT printstats;
  CCTK_INT run_demo;
  CCTK_INT siman_iters_fixed_T;
  CCTK_INT small_jump_distance;
  CCTK_INT use_random_restart_hill_climbing;
  CCTK_INT use_simulated_annealing;
  CCTK_INT verbose;
} PRIVATE_LOOPCONTROL_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_LOOPCONTROL_STRUCT_PARAMS \
  CCTK_REAL const delayed_overhead_threshold = PRIVATE_LOOPCONTROL_STRUCT.delayed_overhead_threshold; \
  CCTK_REAL const immediate_overhead_threshold = PRIVATE_LOOPCONTROL_STRUCT.immediate_overhead_threshold; \
  CCTK_REAL const maximum_setup_overhead = PRIVATE_LOOPCONTROL_STRUCT.maximum_setup_overhead; \
  CCTK_REAL const probability_random_jump = PRIVATE_LOOPCONTROL_STRUCT.probability_random_jump; \
  CCTK_REAL const probability_small_jump = PRIVATE_LOOPCONTROL_STRUCT.probability_small_jump; \
  CCTK_REAL const siman_T_initial = PRIVATE_LOOPCONTROL_STRUCT.siman_T_initial; \
  CCTK_REAL const siman_T_min = PRIVATE_LOOPCONTROL_STRUCT.siman_T_min; \
  CCTK_REAL const siman_k = PRIVATE_LOOPCONTROL_STRUCT.siman_k; \
  CCTK_REAL const siman_mu_T = PRIVATE_LOOPCONTROL_STRUCT.siman_mu_T; \
  CCTK_REAL const siman_probability_change_topology = PRIVATE_LOOPCONTROL_STRUCT.siman_probability_change_topology; \
  CCTK_REAL const siman_step_size = PRIVATE_LOOPCONTROL_STRUCT.siman_step_size; \
  CCTK_INT const cycle_j_tilings = PRIVATE_LOOPCONTROL_STRUCT.cycle_j_tilings; \
  CCTK_INT const debug = PRIVATE_LOOPCONTROL_STRUCT.debug; \
  CCTK_INT const ignore_initial_overhead = PRIVATE_LOOPCONTROL_STRUCT.ignore_initial_overhead; \
  CCTK_INT const lc_inpoints = PRIVATE_LOOPCONTROL_STRUCT.lc_inpoints; \
  CCTK_INT const lc_inthreads = PRIVATE_LOOPCONTROL_STRUCT.lc_inthreads; \
  CCTK_INT const lc_jnpoints = PRIVATE_LOOPCONTROL_STRUCT.lc_jnpoints; \
  CCTK_INT const lc_jnthreads = PRIVATE_LOOPCONTROL_STRUCT.lc_jnthreads; \
  CCTK_INT const lc_knpoints = PRIVATE_LOOPCONTROL_STRUCT.lc_knpoints; \
  CCTK_INT const lc_knthreads = PRIVATE_LOOPCONTROL_STRUCT.lc_knthreads; \
  CCTK_INT const legacy_init = PRIVATE_LOOPCONTROL_STRUCT.legacy_init; \
  CCTK_INT const max_jump_attempts = PRIVATE_LOOPCONTROL_STRUCT.max_jump_attempts; \
  CCTK_INT const nsteps = PRIVATE_LOOPCONTROL_STRUCT.nsteps; \
  CCTK_INT const nx = PRIVATE_LOOPCONTROL_STRUCT.nx; \
  CCTK_INT const ny = PRIVATE_LOOPCONTROL_STRUCT.ny; \
  CCTK_INT const nz = PRIVATE_LOOPCONTROL_STRUCT.nz; \
  CCTK_INT const overhead_threshold_delay = PRIVATE_LOOPCONTROL_STRUCT.overhead_threshold_delay; \
  CCTK_INT const printstats = PRIVATE_LOOPCONTROL_STRUCT.printstats; \
  CCTK_INT const run_demo = PRIVATE_LOOPCONTROL_STRUCT.run_demo; \
  CCTK_INT const siman_iters_fixed_T = PRIVATE_LOOPCONTROL_STRUCT.siman_iters_fixed_T; \
  CCTK_INT const small_jump_distance = PRIVATE_LOOPCONTROL_STRUCT.small_jump_distance; \
  CCTK_INT const use_random_restart_hill_climbing = PRIVATE_LOOPCONTROL_STRUCT.use_random_restart_hill_climbing; \
  CCTK_INT const use_simulated_annealing = PRIVATE_LOOPCONTROL_STRUCT.use_simulated_annealing; \
  CCTK_INT const verbose = PRIVATE_LOOPCONTROL_STRUCT.verbose; \
  enum { \
      dummy_PRIVATE_LOOPCONTROL_STRUCT_delayed_overhead_threshold = sizeof( delayed_overhead_threshold ) \
    , dummy_PRIVATE_LOOPCONTROL_STRUCT_immediate_overhead_threshold = sizeof( immediate_overhead_threshold ) \
    , dummy_PRIVATE_LOOPCONTROL_STRUCT_maximum_setup_overhead = sizeof( maximum_setup_overhead ) \
    , dummy_PRIVATE_LOOPCONTROL_STRUCT_probability_random_jump = sizeof( probability_random_jump ) \
    , dummy_PRIVATE_LOOPCONTROL_STRUCT_probability_small_jump = sizeof( probability_small_jump ) \
    , dummy_PRIVATE_LOOPCONTROL_STRUCT_siman_T_initial = sizeof( siman_T_initial ) \
    , dummy_PRIVATE_LOOPCONTROL_STRUCT_siman_T_min = sizeof( siman_T_min ) \
    , dummy_PRIVATE_LOOPCONTROL_STRUCT_siman_k = sizeof( siman_k ) \
    , dummy_PRIVATE_LOOPCONTROL_STRUCT_siman_mu_T = sizeof( siman_mu_T ) \
    , dummy_PRIVATE_LOOPCONTROL_STRUCT_siman_probability_change_topology = sizeof( siman_probability_change_topology ) \
    , dummy_PRIVATE_LOOPCONTROL_STRUCT_siman_step_size = sizeof( siman_step_size ) \
    , dummy_PRIVATE_LOOPCONTROL_STRUCT_cycle_j_tilings = sizeof( cycle_j_tilings ) \
    , dummy_PRIVATE_LOOPCONTROL_STRUCT_debug = sizeof( debug ) \
    , dummy_PRIVATE_LOOPCONTROL_STRUCT_ignore_initial_overhead = sizeof( ignore_initial_overhead ) \
    , dummy_PRIVATE_LOOPCONTROL_STRUCT_lc_inpoints = sizeof( lc_inpoints ) \
    , dummy_PRIVATE_LOOPCONTROL_STRUCT_lc_inthreads = sizeof( lc_inthreads ) \
    , dummy_PRIVATE_LOOPCONTROL_STRUCT_lc_jnpoints = sizeof( lc_jnpoints ) \
    , dummy_PRIVATE_LOOPCONTROL_STRUCT_lc_jnthreads = sizeof( lc_jnthreads ) \
    , dummy_PRIVATE_LOOPCONTROL_STRUCT_lc_knpoints = sizeof( lc_knpoints ) \
    , dummy_PRIVATE_LOOPCONTROL_STRUCT_lc_knthreads = sizeof( lc_knthreads ) \
    , dummy_PRIVATE_LOOPCONTROL_STRUCT_legacy_init = sizeof( legacy_init ) \
    , dummy_PRIVATE_LOOPCONTROL_STRUCT_max_jump_attempts = sizeof( max_jump_attempts ) \
    , dummy_PRIVATE_LOOPCONTROL_STRUCT_nsteps = sizeof( nsteps ) \
    , dummy_PRIVATE_LOOPCONTROL_STRUCT_nx = sizeof( nx ) \
    , dummy_PRIVATE_LOOPCONTROL_STRUCT_ny = sizeof( ny ) \
    , dummy_PRIVATE_LOOPCONTROL_STRUCT_nz = sizeof( nz ) \
    , dummy_PRIVATE_LOOPCONTROL_STRUCT_overhead_threshold_delay = sizeof( overhead_threshold_delay ) \
    , dummy_PRIVATE_LOOPCONTROL_STRUCT_printstats = sizeof( printstats ) \
    , dummy_PRIVATE_LOOPCONTROL_STRUCT_run_demo = sizeof( run_demo ) \
    , dummy_PRIVATE_LOOPCONTROL_STRUCT_siman_iters_fixed_T = sizeof( siman_iters_fixed_T ) \
    , dummy_PRIVATE_LOOPCONTROL_STRUCT_small_jump_distance = sizeof( small_jump_distance ) \
    , dummy_PRIVATE_LOOPCONTROL_STRUCT_use_random_restart_hill_climbing = sizeof( use_random_restart_hill_climbing ) \
    , dummy_PRIVATE_LOOPCONTROL_STRUCT_use_simulated_annealing = sizeof( use_simulated_annealing ) \
    , dummy_PRIVATE_LOOPCONTROL_STRUCT_verbose = sizeof( verbose ) \
  }; \

