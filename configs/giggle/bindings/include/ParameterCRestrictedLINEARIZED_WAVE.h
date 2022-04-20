#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  const char * bound;
  const char * initial_data;
} RESTRICTED_LINEARIZED_WAVE_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_RESTRICTED_LINEARIZED_WAVE_STRUCT_PARAMS \
  const char * const bound = RESTRICTED_LINEARIZED_WAVE_STRUCT.bound; \
  const char * const initial_data = RESTRICTED_LINEARIZED_WAVE_STRUCT.initial_data; \
  enum { \
      dummy_RESTRICTED_LINEARIZED_WAVE_STRUCT_bound = sizeof( bound ) \
    , dummy_RESTRICTED_LINEARIZED_WAVE_STRUCT_initial_data = sizeof( initial_data ) \
  }; \

