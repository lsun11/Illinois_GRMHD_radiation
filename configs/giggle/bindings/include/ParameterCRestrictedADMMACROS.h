#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_INT spatial_order;
} RESTRICTED_ADMMACROS_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_RESTRICTED_ADMMACROS_STRUCT_PARAMS \
  CCTK_INT const spatial_order = RESTRICTED_ADMMACROS_STRUCT.spatial_order; \
  enum { \
      dummy_RESTRICTED_ADMMACROS_STRUCT_spatial_order = sizeof( spatial_order ) \
  }; \

