#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_INT iteration_to_compute_temp;
  CCTK_INT iteration_to_insert_rad;
  CCTK_INT num_CO;
} RESTRICTED_BHNS_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_RESTRICTED_BHNS_STRUCT_PARAMS \
  CCTK_INT const iteration_to_compute_temp = RESTRICTED_BHNS_STRUCT.iteration_to_compute_temp; \
  CCTK_INT const iteration_to_insert_rad = RESTRICTED_BHNS_STRUCT.iteration_to_insert_rad; \
  CCTK_INT const num_CO = RESTRICTED_BHNS_STRUCT.num_CO; \
  enum { \
      dummy_RESTRICTED_BHNS_STRUCT_iteration_to_compute_temp = sizeof( iteration_to_compute_temp ) \
    , dummy_RESTRICTED_BHNS_STRUCT_iteration_to_insert_rad = sizeof( iteration_to_insert_rad ) \
    , dummy_RESTRICTED_BHNS_STRUCT_num_CO = sizeof( num_CO ) \
  }; \

