#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  const char * conformal_storage;
} RESTRICTED_STATICCONFORMAL_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_RESTRICTED_STATICCONFORMAL_STRUCT_PARAMS \
  const char * const conformal_storage = RESTRICTED_STATICCONFORMAL_STRUCT.conformal_storage; \
  enum { \
      dummy_RESTRICTED_STATICCONFORMAL_STRUCT_conformal_storage = sizeof( conformal_storage ) \
  }; \

