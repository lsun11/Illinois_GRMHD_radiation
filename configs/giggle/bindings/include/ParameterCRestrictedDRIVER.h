#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_INT periodic;
  CCTK_INT periodic_x;
  CCTK_INT periodic_y;
  CCTK_INT periodic_z;
} RESTRICTED_DRIVER_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_RESTRICTED_DRIVER_STRUCT_PARAMS \
  CCTK_INT const periodic = RESTRICTED_DRIVER_STRUCT.periodic; \
  CCTK_INT const periodic_x = RESTRICTED_DRIVER_STRUCT.periodic_x; \
  CCTK_INT const periodic_y = RESTRICTED_DRIVER_STRUCT.periodic_y; \
  CCTK_INT const periodic_z = RESTRICTED_DRIVER_STRUCT.periodic_z; \
  enum { \
      dummy_RESTRICTED_DRIVER_STRUCT_periodic = sizeof( periodic ) \
    , dummy_RESTRICTED_DRIVER_STRUCT_periodic_x = sizeof( periodic_x ) \
    , dummy_RESTRICTED_DRIVER_STRUCT_periodic_y = sizeof( periodic_y ) \
    , dummy_RESTRICTED_DRIVER_STRUCT_periodic_z = sizeof( periodic_z ) \
  }; \

