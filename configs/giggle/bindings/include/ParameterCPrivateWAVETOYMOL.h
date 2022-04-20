#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  const char * bound;
  CCTK_INT estimate_error;
  CCTK_INT order;
} PRIVATE_WAVETOYMOL_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_WAVETOYMOL_STRUCT_PARAMS \
  const char * const bound = PRIVATE_WAVETOYMOL_STRUCT.bound; \
  CCTK_INT const estimate_error = PRIVATE_WAVETOYMOL_STRUCT.estimate_error; \
  CCTK_INT const order = PRIVATE_WAVETOYMOL_STRUCT.order; \
  enum { \
      dummy_PRIVATE_WAVETOYMOL_STRUCT_bound = sizeof( bound ) \
    , dummy_PRIVATE_WAVETOYMOL_STRUCT_estimate_error = sizeof( estimate_error ) \
    , dummy_PRIVATE_WAVETOYMOL_STRUCT_order = sizeof( order ) \
  }; \

