#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_INT use_mask;
} RESTRICTED_SPACEMASK_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_RESTRICTED_SPACEMASK_STRUCT_PARAMS \
  CCTK_INT const use_mask = RESTRICTED_SPACEMASK_STRUCT.use_mask; \
  enum { \
      dummy_RESTRICTED_SPACEMASK_STRUCT_use_mask = sizeof( use_mask ) \
  }; \

