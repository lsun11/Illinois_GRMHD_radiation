#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL opl_lapse_floor;
} PRIVATE_LAPSE_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_LAPSE_STRUCT_PARAMS \
  CCTK_REAL const opl_lapse_floor = PRIVATE_LAPSE_STRUCT.opl_lapse_floor; \
  enum { \
      dummy_PRIVATE_LAPSE_STRUCT_opl_lapse_floor = sizeof( opl_lapse_floor ) \
  }; \

