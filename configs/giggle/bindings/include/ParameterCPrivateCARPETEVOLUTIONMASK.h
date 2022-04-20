#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  const char * enforce_vars;
  CCTK_INT enforce_mask;
  CCTK_INT verbose;
  CCTK_INT writeNaNs;
} PRIVATE_CARPETEVOLUTIONMASK_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_CARPETEVOLUTIONMASK_STRUCT_PARAMS \
  const char * const enforce_vars = PRIVATE_CARPETEVOLUTIONMASK_STRUCT.enforce_vars; \
  CCTK_INT const enforce_mask = PRIVATE_CARPETEVOLUTIONMASK_STRUCT.enforce_mask; \
  CCTK_INT const verbose = PRIVATE_CARPETEVOLUTIONMASK_STRUCT.verbose; \
  CCTK_INT const writeNaNs = PRIVATE_CARPETEVOLUTIONMASK_STRUCT.writeNaNs; \
  enum { \
      dummy_PRIVATE_CARPETEVOLUTIONMASK_STRUCT_enforce_vars = sizeof( enforce_vars ) \
    , dummy_PRIVATE_CARPETEVOLUTIONMASK_STRUCT_enforce_mask = sizeof( enforce_mask ) \
    , dummy_PRIVATE_CARPETEVOLUTIONMASK_STRUCT_verbose = sizeof( verbose ) \
    , dummy_PRIVATE_CARPETEVOLUTIONMASK_STRUCT_writeNaNs = sizeof( writeNaNs ) \
  }; \

