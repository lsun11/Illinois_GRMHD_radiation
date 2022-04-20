#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  const char * bound;
} PRIVATE_WAVEMOL_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_WAVEMOL_STRUCT_PARAMS \
  const char * const bound = PRIVATE_WAVEMOL_STRUCT.bound; \
  enum { \
      dummy_PRIVATE_WAVEMOL_STRUCT_bound = sizeof( bound ) \
  }; \

