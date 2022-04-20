#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  const char * timestep_method;
  CCTK_INT timestep_outonly;
} RESTRICTED_TIME_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_RESTRICTED_TIME_STRUCT_PARAMS \
  const char * const timestep_method = RESTRICTED_TIME_STRUCT.timestep_method; \
  CCTK_INT const timestep_outonly = RESTRICTED_TIME_STRUCT.timestep_outonly; \
  enum { \
      dummy_RESTRICTED_TIME_STRUCT_timestep_method = sizeof( timestep_method ) \
    , dummy_RESTRICTED_TIME_STRUCT_timestep_outonly = sizeof( timestep_outonly ) \
  }; \

