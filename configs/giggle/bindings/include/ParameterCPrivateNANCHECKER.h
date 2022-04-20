#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  const char * action_if_found;
  const char * check_for;
  const char * check_vars;
  const char * verbose;
  CCTK_INT check_after;
  CCTK_INT check_every;
  CCTK_INT out_NaNmask;
  CCTK_INT report_max;
} PRIVATE_NANCHECKER_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_NANCHECKER_STRUCT_PARAMS \
  const char * const action_if_found = PRIVATE_NANCHECKER_STRUCT.action_if_found; \
  const char * const check_for = PRIVATE_NANCHECKER_STRUCT.check_for; \
  const char * const check_vars = PRIVATE_NANCHECKER_STRUCT.check_vars; \
  const char * const verbose = PRIVATE_NANCHECKER_STRUCT.verbose; \
  CCTK_INT const check_after = PRIVATE_NANCHECKER_STRUCT.check_after; \
  CCTK_INT const check_every = PRIVATE_NANCHECKER_STRUCT.check_every; \
  CCTK_INT const out_NaNmask = PRIVATE_NANCHECKER_STRUCT.out_NaNmask; \
  CCTK_INT const report_max = PRIVATE_NANCHECKER_STRUCT.report_max; \
  enum { \
      dummy_PRIVATE_NANCHECKER_STRUCT_action_if_found = sizeof( action_if_found ) \
    , dummy_PRIVATE_NANCHECKER_STRUCT_check_for = sizeof( check_for ) \
    , dummy_PRIVATE_NANCHECKER_STRUCT_check_vars = sizeof( check_vars ) \
    , dummy_PRIVATE_NANCHECKER_STRUCT_verbose = sizeof( verbose ) \
    , dummy_PRIVATE_NANCHECKER_STRUCT_check_after = sizeof( check_after ) \
    , dummy_PRIVATE_NANCHECKER_STRUCT_check_every = sizeof( check_every ) \
    , dummy_PRIVATE_NANCHECKER_STRUCT_out_NaNmask = sizeof( out_NaNmask ) \
    , dummy_PRIVATE_NANCHECKER_STRUCT_report_max = sizeof( report_max ) \
  }; \

