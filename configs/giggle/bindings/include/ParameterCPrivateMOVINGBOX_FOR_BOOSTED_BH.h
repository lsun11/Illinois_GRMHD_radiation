#ifdef __cplusplus
extern "C"
{
#endif

extern struct
{
  CCTK_REAL max_dev;
  CCTK_REAL x0_bh;
} PRIVATE_MOVINGBOX_FOR_BOOSTED_BH_STRUCT;

#ifdef __cplusplus
}
#endif

#define DECLARE_PRIVATE_MOVINGBOX_FOR_BOOSTED_BH_STRUCT_PARAMS \
  CCTK_REAL const max_dev = PRIVATE_MOVINGBOX_FOR_BOOSTED_BH_STRUCT.max_dev; \
  CCTK_REAL const x0_bh = PRIVATE_MOVINGBOX_FOR_BOOSTED_BH_STRUCT.x0_bh; \
  enum { \
      dummy_PRIVATE_MOVINGBOX_FOR_BOOSTED_BH_STRUCT_max_dev = sizeof( max_dev ) \
    , dummy_PRIVATE_MOVINGBOX_FOR_BOOSTED_BH_STRUCT_x0_bh = sizeof( x0_bh ) \
  }; \

