#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_INT radpower;
  CCTK_INT register_copy;
  CCTK_INT register_flat;
  CCTK_INT register_none;
  CCTK_INT register_radiation;
  CCTK_INT register_robin;
  CCTK_INT register_scalar;
  CCTK_INT register_static;
} RESTRICTED_BOUNDARY_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_RESTRICTED_BOUNDARY_STRUCT_PARAMS \
  CCTK_INT const radpower = RESTRICTED_BOUNDARY_STRUCT.radpower; \
  CCTK_INT const register_copy = RESTRICTED_BOUNDARY_STRUCT.register_copy; \
  CCTK_INT const register_flat = RESTRICTED_BOUNDARY_STRUCT.register_flat; \
  CCTK_INT const register_none = RESTRICTED_BOUNDARY_STRUCT.register_none; \
  CCTK_INT const register_radiation = RESTRICTED_BOUNDARY_STRUCT.register_radiation; \
  CCTK_INT const register_robin = RESTRICTED_BOUNDARY_STRUCT.register_robin; \
  CCTK_INT const register_scalar = RESTRICTED_BOUNDARY_STRUCT.register_scalar; \
  CCTK_INT const register_static = RESTRICTED_BOUNDARY_STRUCT.register_static; \
  enum { \
      dummy_RESTRICTED_BOUNDARY_STRUCT_radpower = sizeof( radpower ) \
    , dummy_RESTRICTED_BOUNDARY_STRUCT_register_copy = sizeof( register_copy ) \
    , dummy_RESTRICTED_BOUNDARY_STRUCT_register_flat = sizeof( register_flat ) \
    , dummy_RESTRICTED_BOUNDARY_STRUCT_register_none = sizeof( register_none ) \
    , dummy_RESTRICTED_BOUNDARY_STRUCT_register_radiation = sizeof( register_radiation ) \
    , dummy_RESTRICTED_BOUNDARY_STRUCT_register_robin = sizeof( register_robin ) \
    , dummy_RESTRICTED_BOUNDARY_STRUCT_register_scalar = sizeof( register_scalar ) \
    , dummy_RESTRICTED_BOUNDARY_STRUCT_register_static = sizeof( register_static ) \
  }; \

