#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL cctk_final_time;
  CCTK_REAL cctk_initial_time;
  CCTK_REAL max_runtime;
  const char * terminate;
  CCTK_INT cctk_itlast;
  CCTK_INT terminate_next;
} RESTRICTED_CACTUS_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_RESTRICTED_CACTUS_STRUCT_PARAMS \
  CCTK_REAL const cctk_final_time = RESTRICTED_CACTUS_STRUCT.cctk_final_time; \
  CCTK_REAL const cctk_initial_time = RESTRICTED_CACTUS_STRUCT.cctk_initial_time; \
  CCTK_REAL const max_runtime = RESTRICTED_CACTUS_STRUCT.max_runtime; \
  const char * const terminate = RESTRICTED_CACTUS_STRUCT.terminate; \
  CCTK_INT const cctk_itlast = RESTRICTED_CACTUS_STRUCT.cctk_itlast; \
  CCTK_INT const terminate_next = RESTRICTED_CACTUS_STRUCT.terminate_next; \
  enum { \
      dummy_RESTRICTED_CACTUS_STRUCT_cctk_final_time = sizeof( cctk_final_time ) \
    , dummy_RESTRICTED_CACTUS_STRUCT_cctk_initial_time = sizeof( cctk_initial_time ) \
    , dummy_RESTRICTED_CACTUS_STRUCT_max_runtime = sizeof( max_runtime ) \
    , dummy_RESTRICTED_CACTUS_STRUCT_terminate = sizeof( terminate ) \
    , dummy_RESTRICTED_CACTUS_STRUCT_cctk_itlast = sizeof( cctk_itlast ) \
    , dummy_RESTRICTED_CACTUS_STRUCT_terminate_next = sizeof( terminate_next ) \
  }; \

