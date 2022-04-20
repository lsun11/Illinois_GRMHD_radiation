#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_INT log_interp_coords;
} PRIVATE_AEILOCALINTERP_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_AEILOCALINTERP_STRUCT_PARAMS \
  CCTK_INT const log_interp_coords = PRIVATE_AEILOCALINTERP_STRUCT.log_interp_coords; \
  enum { \
      dummy_PRIVATE_AEILOCALINTERP_STRUCT_log_interp_coords = sizeof( log_interp_coords ) \
  }; \

