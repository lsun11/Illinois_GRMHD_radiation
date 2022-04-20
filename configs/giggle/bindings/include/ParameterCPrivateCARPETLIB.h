#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  const char * memstat_file;
  const char * timestat_file;
  CCTK_INT barrier_between_stages;
  CCTK_INT barriers;
  CCTK_INT check_bboxes;
  CCTK_INT combine_recompose;
  CCTK_INT combine_sends;
  CCTK_INT commstate_verbose;
  CCTK_INT interleave_communications;
  CCTK_INT max_allowed_memory_MB;
  CCTK_INT output_bboxes;
  CCTK_INT poison_new_memory;
  CCTK_INT poison_value;
  CCTK_INT print_memstats_every;
  CCTK_INT print_timestats_every;
  CCTK_INT reduce_mpi_waitall;
  CCTK_INT use_mpi_send;
  CCTK_INT use_mpi_ssend;
  CCTK_INT vary_tags;
  CCTK_INT verbose;
} PRIVATE_CARPETLIB_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_CARPETLIB_STRUCT_PARAMS \
  const char * const memstat_file = PRIVATE_CARPETLIB_STRUCT.memstat_file; \
  const char * const timestat_file = PRIVATE_CARPETLIB_STRUCT.timestat_file; \
  CCTK_INT const barrier_between_stages = PRIVATE_CARPETLIB_STRUCT.barrier_between_stages; \
  CCTK_INT const barriers = PRIVATE_CARPETLIB_STRUCT.barriers; \
  CCTK_INT const check_bboxes = PRIVATE_CARPETLIB_STRUCT.check_bboxes; \
  CCTK_INT const combine_recompose = PRIVATE_CARPETLIB_STRUCT.combine_recompose; \
  CCTK_INT const combine_sends = PRIVATE_CARPETLIB_STRUCT.combine_sends; \
  CCTK_INT const commstate_verbose = PRIVATE_CARPETLIB_STRUCT.commstate_verbose; \
  CCTK_INT const interleave_communications = PRIVATE_CARPETLIB_STRUCT.interleave_communications; \
  CCTK_INT const max_allowed_memory_MB = PRIVATE_CARPETLIB_STRUCT.max_allowed_memory_MB; \
  CCTK_INT const output_bboxes = PRIVATE_CARPETLIB_STRUCT.output_bboxes; \
  CCTK_INT const poison_new_memory = PRIVATE_CARPETLIB_STRUCT.poison_new_memory; \
  CCTK_INT const poison_value = PRIVATE_CARPETLIB_STRUCT.poison_value; \
  CCTK_INT const print_memstats_every = PRIVATE_CARPETLIB_STRUCT.print_memstats_every; \
  CCTK_INT const print_timestats_every = PRIVATE_CARPETLIB_STRUCT.print_timestats_every; \
  CCTK_INT const reduce_mpi_waitall = PRIVATE_CARPETLIB_STRUCT.reduce_mpi_waitall; \
  CCTK_INT const use_mpi_send = PRIVATE_CARPETLIB_STRUCT.use_mpi_send; \
  CCTK_INT const use_mpi_ssend = PRIVATE_CARPETLIB_STRUCT.use_mpi_ssend; \
  CCTK_INT const vary_tags = PRIVATE_CARPETLIB_STRUCT.vary_tags; \
  CCTK_INT const verbose = PRIVATE_CARPETLIB_STRUCT.verbose; \
  enum { \
      dummy_PRIVATE_CARPETLIB_STRUCT_memstat_file = sizeof( memstat_file ) \
    , dummy_PRIVATE_CARPETLIB_STRUCT_timestat_file = sizeof( timestat_file ) \
    , dummy_PRIVATE_CARPETLIB_STRUCT_barrier_between_stages = sizeof( barrier_between_stages ) \
    , dummy_PRIVATE_CARPETLIB_STRUCT_barriers = sizeof( barriers ) \
    , dummy_PRIVATE_CARPETLIB_STRUCT_check_bboxes = sizeof( check_bboxes ) \
    , dummy_PRIVATE_CARPETLIB_STRUCT_combine_recompose = sizeof( combine_recompose ) \
    , dummy_PRIVATE_CARPETLIB_STRUCT_combine_sends = sizeof( combine_sends ) \
    , dummy_PRIVATE_CARPETLIB_STRUCT_commstate_verbose = sizeof( commstate_verbose ) \
    , dummy_PRIVATE_CARPETLIB_STRUCT_interleave_communications = sizeof( interleave_communications ) \
    , dummy_PRIVATE_CARPETLIB_STRUCT_max_allowed_memory_MB = sizeof( max_allowed_memory_MB ) \
    , dummy_PRIVATE_CARPETLIB_STRUCT_output_bboxes = sizeof( output_bboxes ) \
    , dummy_PRIVATE_CARPETLIB_STRUCT_poison_new_memory = sizeof( poison_new_memory ) \
    , dummy_PRIVATE_CARPETLIB_STRUCT_poison_value = sizeof( poison_value ) \
    , dummy_PRIVATE_CARPETLIB_STRUCT_print_memstats_every = sizeof( print_memstats_every ) \
    , dummy_PRIVATE_CARPETLIB_STRUCT_print_timestats_every = sizeof( print_timestats_every ) \
    , dummy_PRIVATE_CARPETLIB_STRUCT_reduce_mpi_waitall = sizeof( reduce_mpi_waitall ) \
    , dummy_PRIVATE_CARPETLIB_STRUCT_use_mpi_send = sizeof( use_mpi_send ) \
    , dummy_PRIVATE_CARPETLIB_STRUCT_use_mpi_ssend = sizeof( use_mpi_ssend ) \
    , dummy_PRIVATE_CARPETLIB_STRUCT_vary_tags = sizeof( vary_tags ) \
    , dummy_PRIVATE_CARPETLIB_STRUCT_verbose = sizeof( verbose ) \
  }; \

