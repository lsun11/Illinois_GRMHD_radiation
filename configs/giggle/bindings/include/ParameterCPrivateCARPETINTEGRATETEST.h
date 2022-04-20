#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL constant;
  CCTK_REAL timefact;
} PRIVATE_CARPETINTEGRATETEST_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_CARPETINTEGRATETEST_STRUCT_PARAMS \
  CCTK_REAL const constant = PRIVATE_CARPETINTEGRATETEST_STRUCT.constant; \
  CCTK_REAL const timefact = PRIVATE_CARPETINTEGRATETEST_STRUCT.timefact; \
  enum { \
      dummy_PRIVATE_CARPETINTEGRATETEST_STRUCT_constant = sizeof( constant ) \
    , dummy_PRIVATE_CARPETINTEGRATETEST_STRUCT_timefact = sizeof( timefact ) \
  }; \

