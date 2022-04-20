#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL C_ko;
} RESTRICTED_EXCISION_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_RESTRICTED_EXCISION_STRUCT_PARAMS \
  CCTK_REAL const C_ko = RESTRICTED_EXCISION_STRUCT.C_ko; \
  enum { \
      dummy_RESTRICTED_EXCISION_STRUCT_C_ko = sizeof( C_ko ) \
  }; \

